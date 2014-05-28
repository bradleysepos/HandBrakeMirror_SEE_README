/* nlmeans.c

   Copyright (c) 2013 Dirk Farin
   Copyright (c) 2003-2014 HandBrake Team
   This file is part of the HandBrake source code
   Homepage: <http://handbrake.fr/>.
   It may be used under the terms of the GNU General Public License v2.
   For full terms see the file COPYING file or visit http://www.gnu.org/licenses/gpl-2.0.html
 */

#include "hb.h"
#include "hbffmpeg.h"

#define NLMEANS_STRENGTH_DEFAULT    8
#define NLMEANS_ORIGIN_TUNE_DEFAULT 1
#define NLMEANS_PATCH_DEFAULT       7
#define NLMEANS_RANGE_DEFAULT       3
#define NLMEANS_FRAMES_DEFAULT      2

#define NLMEANS_FRAMES_MAX 32
#define NLMEANS_EXPSIZE    128

typedef struct
{
    uint8_t *mem;
    uint8_t *image;
    int w;
    int h;
    int border;
} BorderedPlane;

struct PixelSum
{
    float weight_sum;
    float pixel_sum;
};

struct hb_filter_private_s
{
    double strength;    // averaging weight decay, larger produces smoother output
    double origin_tune; // weight tuning for origin patch, 0.00..1.00
    int    patch;       // pixel context region width  (must be odd)
    int    range;       // spatial search window width (must be odd)
    int    frames;      // temporal search depth in frames

    BorderedPlane frame_tmp[3][32];
    int           frame_ready[3][32];
};

static int hb_nlmeans_init(hb_filter_object_t *filter,
                           hb_filter_init_t   *init);

static int hb_nlmeans_work(hb_filter_object_t *filter,
                           hb_buffer_t **buf_in,
                           hb_buffer_t **buf_out);

static void hb_nlmeans_close(hb_filter_object_t *filter);

hb_filter_object_t hb_filter_nlmeans =
{
    .id            = HB_FILTER_NLMEANS,
    .enforce_order = 1,
    .name          = "Denoise (nlmeans)",
    .settings      = NULL,
    .init          = hb_nlmeans_init,
    .work          = hb_nlmeans_work,
    .close         = hb_nlmeans_close,
};

static void nlmeans_copy_bordered(uint8_t *src,
                                  int src_w,
                                  int src_h,
                                  BorderedPlane *dst,
                                  int dst_w,
                                  int dst_h,
                                  int border)
{

    uint8_t *mem   = malloc(dst_w * dst_h * sizeof(uint8_t));
    uint8_t *image = mem + border + dst_w*border;

    // Copy main image
    for (int y = 0; y < src_h; y++)
    {
        memcpy(image + y*dst_w, src + y*src_w, src_w);
    }

    // Copy borders
    for (int k = 0; k < border; k++)
    {
        memcpy(image - (k+1)    *dst_w, src,                   src_w);
        memcpy(image + (src_h+k)*dst_w, src + (src_h-1)*src_w, src_w);
    }
    for (int k = 0; k < border; k++)
    {
        for (int y = -border; y < src_h + border; y++)
        {
            *(image - (k+1)     + y*dst_w) = image[y*dst_w];
            *(image + (k+src_w) + y*dst_w) = image[y*dst_w + (src_w-1)];
        }
    }

    dst->mem    = mem;
    dst->image  = image;
    dst->w      = dst_w;
    dst->h      = dst_h;
    dst->border = border;

}

static void nlmeans_plane(BorderedPlane *plane_tmp,
                          int *plane_ready,
                          uint8_t *dst,
                          int w,
                          int h,
                          double h_param,
                          double origin_tune,
                          int n,
                          int r)
{

    int n_half = (n-1) /2;
    int r_half = (r-1) /2;

    // Source image
    uint8_t *src = plane_tmp[0].image;
    int src_w    = plane_tmp[0].w;

    // Allocate temporary pixel sums
    struct PixelSum *tmp_data = calloc(w * h, sizeof(struct PixelSum));

    // Allocate integral image
    int integral_stride = w + 2*16;
    uint32_t *integral_mem = malloc(integral_stride * (h+1) * sizeof(uint32_t));
    uint32_t *integral     = integral_mem + integral_stride + 16;

    // Precompute exponential table
    float exptable[NLMEANS_EXPSIZE];
    float weight_factor       = 1.0/n/n / (h_param * h_param);
    float min_weight_in_table = 0.0005;
    float stretch             = NLMEANS_EXPSIZE / (-log(min_weight_in_table));
    float weight_fact_table   = weight_factor * stretch;
    int diff_max              = NLMEANS_EXPSIZE / weight_fact_table;
    for (int i = 0; i < NLMEANS_EXPSIZE; i++)
    {
        exptable[i] = exp(-i/stretch);
    }
    exptable[NLMEANS_EXPSIZE-1] = 0;

    // Iterate through available frames
    for (int plane_index = 0; plane_ready[plane_index] == 1; plane_index++)
    {

        // Compare image
        uint8_t *compare = plane_tmp[plane_index].image;
        int compare_w    = plane_tmp[plane_index].w;

        // Iterate through all displacements
        for (int dy = -r_half; dy <= r_half; dy++)
        {
            for (int dx = -r_half; dx <= r_half; dx++)
            {

                // Apply special weight tuning to origin patch
                if (dx == 0 && dy == 0 && plane_index == 0)
                {
                    // TODO: Parallelize this
                    for (int y = n_half; y < h-n + n_half; y++)
                    {
                        for (int x = n_half; x < w-n + n_half; x++)
                        {
                            tmp_data[y*w + x].weight_sum += origin_tune;
                            tmp_data[y*w + x].pixel_sum  += origin_tune * src[y*src_w + x];
                        }
                    }
                    continue;
                }

                // Build integral
                memset(integral-1 - integral_stride, 0, (w+1) * sizeof(uint32_t));
                for (int y = 0; y < h; y++)
                {
                    uint8_t  *p1  = src + y*src_w;
                    uint8_t  *p2  = compare + (y+dy) * compare_w + dx;
                    uint32_t *out = integral + (y*integral_stride) - 1;

                    *out++ = 0;

                    for (int x = 0; x < w; x++)
                    {
                        int diff = *p1++ - *p2++;
                        *out = *(out-1) + diff * diff;
                        out++;
                    }

                    if (y > 0)
                    {
                        out = integral + y*integral_stride;

                        for (int x = 0; x < w; x++)
                        {
                            *out += *(out - integral_stride);
                            out++;
                        }
                    }
                }

                // Average displacements
                // TODO: Parallelize this
                for (int y = 0; y <= h-n; y++)
                {
                    uint32_t *integral_ptr1 = integral + (y  -1)*integral_stride - 1;
                    uint32_t *integral_ptr2 = integral + (y+n-1)*integral_stride - 1;

                    for (int x = 0; x <= w-n; x++)
                    {
                        int xc = x + n_half;
                        int yc = y + n_half;

                        // Difference between patches
                        int diff = (uint32_t)(integral_ptr2[n] - integral_ptr2[0] - integral_ptr1[n] + integral_ptr1[0]);

                        // Sum pixel with weight
                        if (diff < diff_max)
                        {
                            int diffidx = diff * weight_fact_table;

                            //float weight = exp(-diff*weightFact);
                            float weight = exptable[diffidx];

                            tmp_data[yc*w + xc].weight_sum += weight;
                            tmp_data[yc*w + xc].pixel_sum  += weight * compare[(yc+dy)*compare_w + xc + dx];
                        }

                        integral_ptr1++;
                        integral_ptr2++;
                    }
                }
            }
        }
    }

    // Copy border area
    {
        for (int y = 0; y < n_half; y++)
        {
            memcpy(dst + y*w, src + y*src_w, w);
        }
        for (int y = h-n_half; y < h; y++)
        {
            memcpy(dst + y*w, src + y*src_w, w);
        }
        for (int y = n_half; y < h-n_half; y++)
        {
            memcpy(dst + y*w,            src + y*src_w,            n_half);
            memcpy(dst + y*w + w-n_half, src + y*src_w + w-n_half, n_half);
        }
    }

    // Copy main image
    for (int y = n_half; y < h-n_half; y++)
    {
        for (int x = n_half; x < w-n_half; x++)
        {
            *(dst + y*w + x) = tmp_data[y*w + x].pixel_sum / tmp_data[y*w + x].weight_sum;
        }
    }

    free(tmp_data);
    free(integral_mem);

}

static int hb_nlmeans_init(hb_filter_object_t *filter,
                           hb_filter_init_t *init)
{
    filter->private_data = calloc(sizeof(struct hb_filter_private_s), 1);
    hb_filter_private_t *pv = filter->private_data;

    pv->strength    = NLMEANS_STRENGTH_DEFAULT;
    pv->origin_tune = NLMEANS_ORIGIN_TUNE_DEFAULT;
    pv->patch       = NLMEANS_PATCH_DEFAULT;
    pv->range       = NLMEANS_RANGE_DEFAULT;
    pv->frames      = NLMEANS_FRAMES_DEFAULT;

    if (filter->settings)
    {
        sscanf(filter->settings, "%lf:%lf:%d:%d:%d", &pv->strength, &pv->origin_tune, &pv->patch, &pv->range, &pv->frames);
    }

    if (pv->origin_tune < 0.01)
    {
        pv->origin_tune = 0.01; // avoid black artifacts
    }
    if (pv->origin_tune > 1)
    {
        pv->origin_tune = 1;
    }
    if (pv->patch % 2 == 0)
    {
        pv->patch--;
    }
    if (pv->patch < 1)
    {
        pv->patch = 1;
    }
    if (pv->range % 2 == 0)
    {
        pv->range--;
    }
    if (pv->range < 1)
    {
        pv->range = 1;
    }
    if (pv->frames < 1)
    {
        pv->frames = 1;
    }
    if (pv->frames > NLMEANS_FRAMES_MAX)
    {
        pv->frames = NLMEANS_FRAMES_MAX;
    }

    for (int c = 0; c < 3; c++)
    {
        for (int f = 0; f < NLMEANS_FRAMES_MAX; f++)
        {
            pv->frame_ready[c][f] = 0;
        }
    }

    return 0;
}

static void hb_nlmeans_close(hb_filter_object_t *filter)
{
    hb_filter_private_t *pv = filter->private_data;

    if (pv == NULL)
    {
        return;
    }

    for (int c = 0; c < 3; c++)
    {
        for (int f = 0; f < pv->frames; f++)
        {
            if (pv->frame_tmp[c][f].mem != NULL)
            {
                free(pv->frame_tmp[c][f].mem);
                pv->frame_tmp[c][f].mem = NULL;
            }
        }
    }

    free(pv);
    filter->private_data = NULL;
}

static int hb_nlmeans_work(hb_filter_object_t *filter,
                           hb_buffer_t **buf_in,
                           hb_buffer_t **buf_out )
{
    hb_filter_private_t *pv = filter->private_data;
    hb_buffer_t *in = *buf_in, *out;

    if (in->size <= 0)
    {
        *buf_out = in;
        *buf_in  = NULL;
        return HB_FILTER_DONE;
    }

    out = hb_video_buffer_init(in->f.width, in->f.height);

    for (int c = 0; c < 3; c++)
    {

        // Release last frame in buffer
        if (pv->frame_tmp[c][pv->frames-1].mem != NULL)
        {
            free(pv->frame_tmp[c][pv->frames-1].mem);
            pv->frame_tmp[c][pv->frames-1].mem = NULL;
        }
        pv->frame_ready[c][pv->frames-1] = 0;

        // Shift frames in buffer down one level
        for (int f = pv->frames-1; f > 0; f--)
        {
            pv->frame_tmp[c][f]   = pv->frame_tmp[c][f-1];
            pv->frame_ready[c][f] = pv->frame_ready[c][f-1];
        }

        // Extend copy of plane with extra border and place in buffer
        int border = (pv->range/2 + 15) /16*16;
        int tmp_w = in->plane[c].stride + 2*border;
        int tmp_h = in->plane[c].height + 2*border;
        nlmeans_copy_bordered(in->plane[c].data,
                              in->plane[c].stride,
                              in->plane[c].height,
                              &pv->frame_tmp[c][0],
                              tmp_w,
                              tmp_h,
                              border);
        pv->frame_ready[c][0] = 1;

        // Process current plane
        nlmeans_plane(pv->frame_tmp[c],
                      pv->frame_ready[c],
                      out->plane[c].data,
                      in->plane[c].stride,
                      in->plane[c].height,
                      pv->strength,
                      pv->origin_tune,
                      pv->patch,
                      pv->range);

    }

    out->s = in->s;
    hb_buffer_move_subs(out, in);

    *buf_out = out;

    return HB_FILTER_OK;
}
