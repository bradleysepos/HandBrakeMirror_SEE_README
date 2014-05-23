/*
 Copyright (c) 2013 Dirk Farin <dirk.farin@gmail.com>

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along
 with this program; if not, write to the Free Software Foundation, Inc.,
 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
*/

#include "hb.h"
#include "hbffmpeg.h"

#define NLMEANS_PATCH_DEFAULT    7
#define NLMEANS_RANGE_DEFAULT    3
#define NLMEANS_FRAMES_DEFAULT   2
#define NLMEANS_STRENGTH_DEFAULT 8.00
#define NLMEANS_ORIGIN_DEFAULT   1.00

#define NLMEANS_FRAMES_MAX 32
#define EXP_TABLE_SIZE     128

#define ABS(A) ( (A) > 0 ? (A) : -(A) )
#define MIN( a, b ) ( (a) > (b) ? (b) : (a) )
#define CEIL_RSHIFT(a,b) (-((-(a)) >> (b)))

typedef struct
{
    unsigned char * mem;
    unsigned char * data;
    int w;
    int h;
    int border;
} BorderedPlane;

struct hb_filter_private_s
{
    int    patch;    // pixel context region width  (must be odd, even--)
    int    range;    // spatial search window width (must be odd, even--)
    int    frames;   // temporal search depth in frames
    double strength; // averaging weight decay, larger produces smoother output
    double origin;   // weight tuning for origin patch, 0.00..1.00

    BorderedPlane frame_tmp[3][32];
    int           frame_ready[3][32];
};

static int hb_nlmeans_init( hb_filter_object_t * filter,
                            hb_filter_init_t * init );

static int hb_nlmeans_work( hb_filter_object_t * filter,
                            hb_buffer_t ** buf_in,
                            hb_buffer_t ** buf_out );

static void hb_nlmeans_close( hb_filter_object_t * filter );

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

struct PixelSum
{
    float weight_sum;
    float pixel_sum;
};

static void nlmeans_copy_bordered( unsigned char * src,
                                   int src_w,
                                   int src_h,
                                   BorderedPlane * dst,
                                   int dst_w,
                                   int dst_h,
                                   int border )
{

//    hb_log("copy_bordered");

    unsigned char * mem = malloc(dst_w * dst_h * sizeof(unsigned char));
    unsigned char * data = mem + border + dst_w*border;
//    hb_log(" malloc %p", mem);

    // Copy main image
    int y;
    for (y = 0; y < src_h; y++)
    {
        memcpy(data + y*dst_w, src + y*src_w, src_w);
    }

    // Copy borders
    int k;
    for (k = 0; k < border; k++)
    {
        memcpy(data - (k+1)*dst_w, src, src_w);
        memcpy(data + (src_h+k)*dst_w, src + (src_h-1)*src_w, src_w);
    }
    for (k = 0; k < border; k++)
    {
        for (y = -border; y < src_h + border; y++)
        {
            *(data - (k+1) + y*dst_w) = data[y*dst_w];
            *(data + (k+src_w) + y*dst_w) = data[y*dst_w + (src_w-1)];
        }
    }

    dst->mem = mem;
//    hb_log(" mem %p", dst->mem);
    dst->data = data;
    dst->w = dst_w;
    dst->h = dst_h;
    dst->border = border;

}

static void nlmeans_plane( BorderedPlane * plane_tmp,
                           int * plane_ready,
                           unsigned char * dst,
                           int w,
                           int h,
                           int n,
                           int r,
                           int f,
                           double h_param,
                           double origin_tune )
{

    // Source image
    //BorderedPlane * plane_tmp0 = &plane_tmp[0];
    //unsigned char * src = plane_tmp0->data;
    //int src_w = plane_tmp0->w;
    unsigned char * src = plane_tmp[0].data;
    int src_w = plane_tmp[0].w;
//    hb_log("process %p", plane_tmp[0].mem);

    int n_half = (n-1)/2;
    int r_half = (r-1)/2;

    // Allocate temporary pixel sums
    struct PixelSum* tmp_data = calloc(w*h,sizeof(struct PixelSum));

    // Allocate integral image
    int integral_stride = w + 2*16;
    uint32_t* integral_mem = malloc( integral_stride * ( h + 1 ) * sizeof(uint32_t) );
    uint32_t* integral = integral_mem + integral_stride + 16;

    // Precompute exponential table
    float exptable[EXP_TABLE_SIZE];
    float weight_factor = 1.0/n/n / (h_param * h_param);
    float min_weight_in_table = 0.0005;
    float stretch = EXP_TABLE_SIZE / (-log(min_weight_in_table));
    float weight_fact_table = weight_factor*stretch;
    int diff_max = EXP_TABLE_SIZE/weight_fact_table;
    // TODO: Parallelize this?
    for (int i=0;i<EXP_TABLE_SIZE;i++)
    {
        exptable[i] = exp(-i/stretch);
    }
    exptable[EXP_TABLE_SIZE-1] = 0;

    // Iterate through available frames
    for (int plane_index = 0; plane_ready[plane_index]; plane_index++)
    {

        // Compare image
        //BorderedPlane * plane_tmpN = &plane_tmp[plane_index];
        //unsigned char * compare = plane_tmpN->data;
        //int compare_w = plane_tmpN->w;
        unsigned char * compare = plane_tmp[plane_index].data;
        int compare_w = plane_tmp[plane_index].w;

        // Iterate through all displacements
        for (int dy=-r_half ; dy<=r_half ; dy++)
        {
            for (int dx=-r_half ; dx<=r_half ; dx++)
            {

                // Apply special weight tuning to origin patch
                if (dx==0 && dy==0 && plane_index==0)
                {
                    // TODO: Parallelize this
                    for (int y=n_half;y<h-n+n_half;y++)
                    {
                        for (int x=n_half;x<w-n+n_half;x++)
                        {
                            tmp_data[y*w+x].weight_sum += origin_tune;
                            tmp_data[y*w+x].pixel_sum  += origin_tune * src[y*src_w+x];
                        }
                    }
                    continue;
                }

                // Build integral
                memset(integral -1 -integral_stride, 0, (w+1)*sizeof(uint32_t));
                for (int y=0;y<h;y++)
                {
                    unsigned char* p1 = src +  y    *src_w;
                    unsigned char* p2 = compare + (y+dy)*compare_w + dx;
                    uint32_t* out = integral + y*integral_stride -1;

                    *out++ = 0;

                    for (int x=0;x<w;x++)
                    {
                        int diff = *p1++ - *p2++;
                        *out = *(out-1) + diff * diff;
                        out++;
                    }

                    if (y>0)
                    {
                        out = integral + y*integral_stride;

                        for (int x=0;x<w;x++)
                        {
                            *out += *(out - integral_stride);
                            out++;
                        }
                    }
                }
                // Average displacements
                // TODO: Parallelize this
                for (int y=0;y<=h-n;y++)
                {
                    uint32_t* integral_ptr1 = integral+(y  -1)*integral_stride-1;
                    uint32_t* integral_ptr2 = integral+(y+n-1)*integral_stride-1;

                    for (int x=0;x<=w-n;x++)
                    {
                        int xc = x+n_half;
                        int yc = y+n_half;

                        // Difference between patches
                        int diff = (uint32_t)(integral_ptr2[n] - integral_ptr2[0] - integral_ptr1[n] + integral_ptr1[0]);

                        // Sum pixel with weight
                        if (diff<diff_max)
                        {
                            int diffidx = diff*weight_fact_table;

                            //float weight = exp(-diff*weightFact);
                            float weight = exptable[diffidx];

                            tmp_data[yc*w+xc].weight_sum += weight;
                            tmp_data[yc*w+xc].pixel_sum  += weight * compare[(yc+dy)*compare_w+xc+dx];
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
            memcpy(dst+y*w, src+y*src_w, w);
        }
        for (int y = h-n_half; y < h; y++)
        {
            memcpy(dst+y*w, src+y*src_w, w);
        }
        for (int y = n_half; y < h-n_half; y++)
        {
            memcpy(dst+y*w,          src+y*src_w,          n_half);
            memcpy(dst+y*w+w-n_half, src+y*src_w+w-n_half, n_half);
        }
    }

    // Copy main image
    for (int y=n_half;y<h-n_half;y++)
    {
        for (int x=n_half;x<w-n_half;x++)
        {
            *(dst+y*w+x) = tmp_data[y*w+x].pixel_sum / tmp_data[y*w+x].weight_sum;
        }
    }

    free(tmp_data);
    free(integral_mem);

}

static int hb_nlmeans_init( hb_filter_object_t * filter,
                            hb_filter_init_t * init )
{
    filter->private_data = calloc( sizeof(struct hb_filter_private_s), 1 );
    hb_filter_private_t * pv = filter->private_data;

    pv->patch    = NLMEANS_PATCH_DEFAULT;
    pv->range    = NLMEANS_RANGE_DEFAULT;
    pv->frames   = NLMEANS_FRAMES_DEFAULT;
    pv->strength = NLMEANS_STRENGTH_DEFAULT;
    pv->origin   = NLMEANS_ORIGIN_DEFAULT;

    if( filter->settings )
    {
        sscanf( filter->settings, "%d:%d:%d:%lf:%lf", &pv->patch, &pv->range, &pv->frames, &pv->strength, &pv->origin );
    }

    if ( pv->patch < 1 )
    {
        pv->patch = 1;
    }
    if ( pv->range < 1 )
    {
        pv->range = 1;
    }
    if ( pv->frames < 1 )
    {
        pv->frames = 1;
    }
    if ( pv->frames > NLMEANS_FRAMES_MAX )
    {
        pv->frames = NLMEANS_FRAMES_MAX;
    }
    if ( pv->origin < 0.01 )
    {
        pv->origin = 0.01;
    }
    if ( pv->origin > 1.00 )
    {
        pv->origin = 1.00;
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

static void hb_nlmeans_close( hb_filter_object_t * filter )
{
    hb_filter_private_t * pv = filter->private_data;

    if (!pv)
    {
        return;
    }

    for (int c = 0; c < 3; c++)
    {
        for (int f = 0; f < pv->frames; f++)
        {
            if (pv->frame_tmp[c][f].mem != NULL)
            {
                hb_log("free %p", pv->frame_tmp[c][f].mem);
                free( pv->frame_tmp[c][f].mem );
                pv->frame_tmp[c][f].mem = NULL;
            }
        }
    }

    free( pv );
    filter->private_data = NULL;
}

static int hb_nlmeans_work( hb_filter_object_t * filter,
                            hb_buffer_t ** buf_in,
                            hb_buffer_t ** buf_out )
{
    hb_filter_private_t * pv = filter->private_data;
    hb_buffer_t * in = *buf_in, * out;

    if ( in->size <= 0 )
    {
        *buf_out = in;
        *buf_in = NULL;
        return HB_FILTER_DONE;
    }

//    hb_log("frame");

    out = hb_video_buffer_init( in->f.width, in->f.height );

    for (int c = 0; c < 3; c++)
    {
//        hb_log(" plane %d", c);
//        hb_log("  mem %p", pv->frame_tmp[c][0].mem);
        // Release last frame in buffer
        int frames_buffered = 0;
        for (int f = 0; f < pv->frames; f++)
        {
            if (pv->frame_ready[c][f])
            {
                frames_buffered++;
            }
        }
//        hb_log("  frames_buffered %d", frames_buffered);
        if (frames_buffered == pv->frames && pv->frame_tmp[c][frames_buffered].mem != NULL)
        {
//            hb_log("   free %p", pv->frame_tmp[c][frames_buffered].mem);
            free( pv->frame_tmp[c][frames_buffered].mem );
            pv->frame_tmp[c][frames_buffered].mem = NULL;
            pv->frame_ready[c][frames_buffered] = 0;
        }
        // Shift frames in buffer
        for (int f = pv->frames-1; f > 0; f--)
        {
//            hb_log("  shift %d", f);
//            hb_log("   cur  %p", pv->frame_tmp[c][f].mem);
//            hb_log("   prev %p", pv->frame_tmp[c][f-1].mem);
            pv->frame_tmp[c][f] = pv->frame_tmp[c][f-1];
//            hb_log("   cur  %p", pv->frame_tmp[c][f].mem);
            pv->frame_ready[c][f] = pv->frame_ready[c][f-1];
        }

        // Extend copy of plane with extra border and place in buffer
        int border = (pv->range/2 + 15) /16*16;
        int tmp_w = in->plane[c].stride + 2*border;
        int tmp_h = in->plane[c].height + 2*border;
        nlmeans_copy_bordered( in->plane[c].data,
                               in->plane[c].stride,
                               in->plane[c].height,
                               &pv->frame_tmp[c][0],
                               tmp_w,
                               tmp_h,
                               border );
        pv->frame_ready[c][0] = 1;

        // Process plane
        nlmeans_plane( pv->frame_tmp[c],
                       pv->frame_ready[c],
                       out->plane[c].data,
                       in->plane[c].stride,
                       in->plane[c].height,
                       pv->patch,
                       pv->range,
                       pv->frames,
                       pv->strength,
                       pv->origin );

    }

    out->s = in->s;
    hb_buffer_move_subs( out, in );

    *buf_out = out;

    return HB_FILTER_OK;
}
