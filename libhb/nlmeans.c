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

#define NLMEANS_STRENGTH_LUMA_DEFAULT      8
#define NLMEANS_STRENGTH_CHROMA_DEFAULT    8
#define NLMEANS_ORIGIN_TUNE_LUMA_DEFAULT   1
#define NLMEANS_ORIGIN_TUNE_CHROMA_DEFAULT 1
#define NLMEANS_PATCH_SIZE_LUMA_DEFAULT    7
#define NLMEANS_PATCH_SIZE_CHROMA_DEFAULT  7
#define NLMEANS_RANGE_LUMA_DEFAULT         3
#define NLMEANS_RANGE_CHROMA_DEFAULT       3
#define NLMEANS_FRAMES_LUMA_DEFAULT        2
#define NLMEANS_FRAMES_CHROMA_DEFAULT      2
#define NLMEANS_PREFILTER_LUMA_DEFAULT     0
#define NLMEANS_PREFILTER_CHROMA_DEFAULT   0

#define NLMEANS_FRAMES_MAX 32
#define NLMEANS_EXPSIZE    128

typedef struct
{
    uint8_t *mem;
    uint8_t *mem_pre;
    uint8_t *image;
    uint8_t *image_pre;
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
    double strength[3];    // averaging weight decay, larger produces smoother output
    double origin_tune[3]; // weight tuning for origin patch, 0.00..1.00
    int    patch_size[3];  // pixel context region width  (must be odd)
    int    range[3];       // spatial search window width (must be odd)
    int    frames[3];      // temporal search depth in frames
    int    prefilter[3];   // prefilter type, can improve weight analysis

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

static void nlmeans_border(uint8_t *src,
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
            *(image - (k-1)     + y*dst_w) = image[y*dst_w];
            *(image + (k+src_w) + y*dst_w) = image[y*dst_w + (src_w-1)];
        }
    }

    dst->mem       = mem;
    dst->mem_pre   = mem;
    dst->image     = image;
    dst->image_pre = image;
    dst->w         = dst_w;
    dst->h         = dst_h;
    dst->border    = border;

}

static void nlmeans_filter_mean(uint8_t *src,
                                uint8_t *dst,
                                int w,
                                int h,
                                int border,
                                int size)
{

    int w_cropped = w - 2*border;
    int h_cropped = h - 2*border;

    // Mean filter
    uint16_t pixel_sum;
    int offset_min = -((size - 1) /2);
    int offset_max =   (size + 1) /2;
    for (int y = 0; y < h_cropped; y++)
    {
        for (int x = 0; x < w_cropped; x++)
        {
            pixel_sum = 0;
            for (int k = offset_min; k < offset_max; k++)
            {
                for (int j = offset_min; j < offset_max; j++)
                {
                    pixel_sum = pixel_sum + src[w*(y+j) + (x+k)];
                }
            }
            *(dst + w*y + x) = (uint8_t)(pixel_sum / (size * size));
        }
    }

}

static void nlmeans_prefilter(BorderedPlane *src,
                              int filter_type)
{

    // Source image
    uint8_t *mem   = src->mem;
    uint8_t *image = src->image;
    int w          = src->w;
    int h          = src->h;
    int border     = src->border;

    // Duplicate plane
    uint8_t *mem_pre = malloc(w * h * sizeof(uint8_t));
    uint8_t *image_pre = mem_pre + border + w*border;
    for (int y = 0; y < h; y++)
    {
        memcpy(mem_pre + y*w, mem + y*w, w);
    }

    // Filter plane; should already have at least 2px extra border on each side
    switch (filter_type)
    {
        case 1:
            // Mean 3x3
            nlmeans_filter_mean(image, image_pre, w, h, border, 3);
            break;
        case 2:
            // Mean 5x5
            nlmeans_filter_mean(image, image_pre, w, h, border, 5);
            break;
    }

    src->mem_pre   = mem_pre;
    src->image_pre = image_pre;

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
    uint8_t *src     = plane_tmp[0].image;
    uint8_t *src_pre = plane_tmp[0].image_pre;
    int src_w        = plane_tmp[0].w;

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
        uint8_t *compare     = plane_tmp[plane_index].image;
        uint8_t *compare_pre = plane_tmp[plane_index].image_pre;
        int compare_w        = plane_tmp[plane_index].w;

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
                    const uint8_t *p1 = src_pre + y*src_w;
                    const uint8_t *p2 = compare_pre + (y+dy) * compare_w + dx;
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
                    const uint32_t *integral_ptr1 = integral + (y  -1)*integral_stride - 1;
                    const uint32_t *integral_ptr2 = integral + (y+n-1)*integral_stride - 1;

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

    // Mark parameters unset
    for (int c = 0; c < 3; c++)
    {
        pv->strength[c]    = -1;
        pv->origin_tune[c] = -1;
        pv->patch_size[c]  = -1;
        pv->range[c]       = -1;
        pv->frames[c]      = -1;
        pv->prefilter[c]   = -1;
    }

    // Read user parameters
    if (filter->settings != NULL)
    {
        sscanf(filter->settings, "%lf:%lf:%d:%d:%d:%d:%lf:%lf:%d:%d:%d:%d:%lf:%lf:%d:%d:%d:%d",
               &pv->strength[0], &pv->origin_tune[0], &pv->patch_size[0], &pv->range[0], &pv->frames[0], &pv->prefilter[0],
               &pv->strength[1], &pv->origin_tune[1], &pv->patch_size[1], &pv->range[1], &pv->frames[1], &pv->prefilter[1],
               &pv->strength[2], &pv->origin_tune[2], &pv->patch_size[2], &pv->range[2], &pv->frames[2], &pv->prefilter[2]);
    }

    // Cascade values
    // Cr not set; inherit Cb. Cb not set; inherit Y. Y not set; defaults.
    for (int c = 1; c < 3; c++)
    {
        if (pv->strength[c]    == -1) { pv->strength[c]    = pv->strength[c-1]; }
        if (pv->origin_tune[c] == -1) { pv->origin_tune[c] = pv->origin_tune[c-1]; }
        if (pv->patch_size[c]  == -1) { pv->patch_size[c]  = pv->patch_size[c-1]; }
        if (pv->range[c]       == -1) { pv->range[c]       = pv->range[c-1]; }
        if (pv->frames[c]      == -1) { pv->frames[c]      = pv->frames[c-1]; }
        if (pv->prefilter[c]   == -1) { pv->prefilter[c]   = pv->prefilter[c-1]; }
    }

    for (int c = 0; c < 3; c++)
    {
        // Replace NULL values with defaults
        if (pv->strength[c]    == -1) { pv->strength[c]    = c ? NLMEANS_STRENGTH_LUMA_DEFAULT    : NLMEANS_STRENGTH_CHROMA_DEFAULT; }
        if (pv->origin_tune[c] == -1) { pv->origin_tune[c] = c ? NLMEANS_ORIGIN_TUNE_LUMA_DEFAULT : NLMEANS_ORIGIN_TUNE_CHROMA_DEFAULT; }
        if (pv->patch_size[c]  == -1) { pv->patch_size[c]  = c ? NLMEANS_PATCH_SIZE_LUMA_DEFAULT  : NLMEANS_PATCH_SIZE_CHROMA_DEFAULT; }
        if (pv->range[c]       == -1) { pv->range[c]       = c ? NLMEANS_RANGE_LUMA_DEFAULT       : NLMEANS_RANGE_CHROMA_DEFAULT; }
        if (pv->frames[c]      == -1) { pv->frames[c]      = c ? NLMEANS_FRAMES_LUMA_DEFAULT      : NLMEANS_FRAMES_CHROMA_DEFAULT; }
        if (pv->prefilter[c]   == -1) { pv->prefilter[c]   = c ? NLMEANS_PREFILTER_LUMA_DEFAULT   : NLMEANS_PREFILTER_CHROMA_DEFAULT; }

        // Sanitize
        if (pv->origin_tune[c] < 0.01)  { pv->origin_tune[c] = 0.01; } // avoid black artifacts
        if (pv->origin_tune[c] > 1)     { pv->origin_tune[c] = 1; }
        if (pv->patch_size[c] % 2 == 0) { pv->patch_size[c]--; }
        if (pv->patch_size[c] < 1)      { pv->patch_size[c] = 1; }
        if (pv->range[c] % 2 == 0)      { pv->range[c]--; }
        if (pv->range[c] < 1)           { pv->range[c] = 1; }
        if (pv->frames[c] < 1)          { pv->frames[c] = 1; }
        if (pv->frames[c] > NLMEANS_FRAMES_MAX) { pv->frames[c] = NLMEANS_FRAMES_MAX; }

        // Mark buffer empty
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
        for (int f = 0; f < pv->frames[c]; f++)
        {
            if (pv->frame_tmp[c][f].mem_pre != NULL &&
                pv->frame_tmp[c][f].mem_pre != pv->frame_tmp[c][f].mem)
            {
                free(pv->frame_tmp[c][f].mem_pre);
                pv->frame_tmp[c][f].mem_pre = NULL;
            }
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

        if (pv->strength[c] == 0)
        {
            out->plane[c].data = in->plane[c].data;
            continue;
        }

        int frames = pv->frames[c];

        // Release last frame in buffer
        if (pv->frame_tmp[c][frames-1].mem_pre != NULL &&
            pv->frame_tmp[c][frames-1].mem_pre != pv->frame_tmp[c][frames-1].mem)
        {
            free(pv->frame_tmp[c][frames-1].mem_pre);
            pv->frame_tmp[c][frames-1].mem_pre = NULL;
        }
        if (pv->frame_tmp[c][frames-1].mem != NULL)
        {
            free(pv->frame_tmp[c][frames-1].mem);
            pv->frame_tmp[c][frames-1].mem = NULL;
        }
        pv->frame_ready[c][frames-1] = 0;

        // Shift frames in buffer down one level
        for (int f = frames-1; f > 0; f--)
        {
            pv->frame_tmp[c][f]   = pv->frame_tmp[c][f-1];
            pv->frame_ready[c][f] = pv->frame_ready[c][f-1];
        }

        // Extend copy of plane with extra border and place in buffer
        int border = ((pv->range[c] + 2) / 2 + 15) /16*16;
        int tmp_w = in->plane[c].stride + 2*border;
        int tmp_h = in->plane[c].height + 2*border;
        nlmeans_border(in->plane[c].data,
                       in->plane[c].stride,
                       in->plane[c].height,
                       &pv->frame_tmp[c][0],
                       tmp_w,
                       tmp_h,
                       border);
        if (pv->prefilter[c] > 0)
        {
            nlmeans_prefilter(&pv->frame_tmp[c][0],
                              pv->prefilter[c]);
        }
        pv->frame_ready[c][0] = 1;

        // Process current plane
        nlmeans_plane(pv->frame_tmp[c],
                      pv->frame_ready[c],
                      out->plane[c].data,
                      in->plane[c].stride,
                      in->plane[c].height,
                      pv->strength[c],
                      pv->origin_tune[c],
                      pv->patch_size[c],
                      pv->range[c]);

    }

    out->s = in->s;
    hb_buffer_move_subs(out, in);

    *buf_out = out;

    return HB_FILTER_OK;
}
