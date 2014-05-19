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

#define NLMEANS_H_DEFAULT          8
#define NLMEANS_RANGE_DEFAULT      3
#define NLMEANS_FRAMES_DEFAULT     1
#define NLMEANS_PATCH_DEFAULT 7

#define NLMEANS_MAX_IMAGES 1 //32
#define EXP_TABLE_SIZE     128

#define ABS(A) ( (A) > 0 ? (A) : -(A) )
#define MIN( a, b ) ( (a) > (b) ? (b) : (a) )
#define CEIL_RSHIFT(a,b) (-((-(a)) >> (b)))

struct hb_filter_private_s
{
    unsigned char * out_tmp[3];
    double h_param;  // averaging weight decay (larger == smoother)
    int    range;    // spatial search window width (must be odd)
                     //   e.g. 15 = 15x15 search window
    int    frames;   // temporal search depth (1 == static only)
    int    patch;    // pixel context region width (must be odd)
                     //   e.g. 3 = 3x3 context window
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

static void nlmeans_plane(unsigned char * src,
                          int sw,
                          int sh,
                          unsigned char * dst,
                          int w,
                          int h,
                          int c,
                          int h_param,
                          int range,
                          int frames,
                          int patch )
{

    int n = (patch|1);
    int r = (range     |1);

    int n_half = (n-1)/2;
    int r_half = (r-1)/2;

    // Allocate temporary pixel sums
    struct PixelSum* tmp_data = calloc(w*h,sizeof(struct PixelSum));

    // Allocate integral image
    int integral_stride = w + 2*16;
    uint32_t* integral_mem = malloc( integral_stride * ( h + 1 ) * sizeof(uint32_t) );
    uint32_t* integral = integral_mem + integral_stride + 16;

    // Precompute exponential table
    float weight_factor = 1.0/n/n / (h_param * h_param);

    int table_size=EXP_TABLE_SIZE;
    float min_weight_in_table = 0.0005;

    float exptable[EXP_TABLE_SIZE];

    float stretch = table_size/ (-log(min_weight_in_table));
    float weight_fact_table = weight_factor*stretch;
    int diff_max = table_size/weight_fact_table;

    for (int i=0;i<table_size;i++) {
        exptable[i] = exp(-i/stretch);
    }
    exptable[table_size-1]=0;

    // Iterate through available frames
    int image_idx=0;
    //for (int image_idx=0; image_idx < n_refs; image_idx++)
    //{

        // Source image
        unsigned char* current_image = src;
        int current_image_stride = sw;

        // Compare image
        //unsigned char* compare_image = ref[image_idx].plane[c].data;
        //int compare_image_stride = ref[image_idx].plane[c].stride;
        unsigned char* compare_image = src;
        int compare_image_stride = sw;

        // Iterate through all displacements
        for (int dy=-r_half ; dy<=r_half ; dy++)
        {
            for (int dx=-r_half ; dx<=r_half ; dx++)
            {

                // Special case: origin, weight = 1
                if (dx==0 && dy==0 && image_idx==0) {
                    // TODO: Parallelize this
                    for (int y=n_half;y<h-n+n_half;y++) {
                        for (int x=n_half;x<w-n+n_half;x++) {
                            tmp_data[y*w+x].weight_sum += 1;
                            tmp_data[y*w+x].pixel_sum  += current_image[y*current_image_stride+x];
                        }
                    }
                    continue;
                }

                // Normal case
                // Build integral
                memset(integral -1 -integral_stride, 0, (w+1)*sizeof(uint32_t));
                for (int y=0;y<h;y++) {
                    unsigned char* p1 = current_image +  y    *current_image_stride;
                    unsigned char* p2 = compare_image + (y+dy)*compare_image_stride + dx;
                    uint32_t* out = integral + y*integral_stride -1;

                    *out++ = 0;

                    for (int x=0;x<w;x++)
                    {
                        int diff = *p1++ - *p2++;
                        *out = *(out-1) + diff * diff;
                        out++;
                    }

                    if (y>0) {
                        out = integral + y*integral_stride;

                        for (int x=0;x<w;x++) {
                            *out += *(out - integral_stride);
                            out++;
                        }
                    }
                }
                // Average displacements
                // TODO: Parallelize this
                for (int y=0;y<=h-n;y++) {
                    uint32_t* integral_ptr1 = integral+(y  -1)*integral_stride-1;
                    uint32_t* integral_ptr2 = integral+(y+n-1)*integral_stride-1;

                    for (int x=0;x<=w-n;x++) {
                        int xc = x+n_half;
                        int yc = y+n_half;

                        // Difference between patches
                        int diff = (uint32_t)(integral_ptr2[n] - integral_ptr2[0] - integral_ptr1[n] + integral_ptr1[0]);

                        // Sum pixel with weight
                        if (diff<diff_max) {
                            int diffidx = diff*weight_fact_table;

                            //float weight = exp(-diff*weightFact);
                            float weight = exptable[diffidx];

                            tmp_data[yc*w+xc].weight_sum += weight;
                            tmp_data[yc*w+xc].pixel_sum  += weight * compare_image[(yc+dy)*compare_image_stride+xc+dx];
                        }

                        integral_ptr1++;
                        integral_ptr2++;
                    }
                }
            }
        }
    //}

    // Output border area
    {

        for (int y=0;       y<n_half  ;y++)
        {
            memcpy(dst+y*w, src+y*sw, w);
        }
        for (int y=h-n_half;y<h       ;y++)
        {
            memcpy(dst+y*w, src+y*sw, w);
        }
        for (int y=n_half  ;y<h-n_half;y++)
        {
            memcpy(dst+y*w,          src+y*sw,          n_half);
            memcpy(dst+y*w+w-n_half, src+y*sw+w-n_half, n_half);
        }
    }

    // Output main image
    for (int y=n_half;y<h-n_half;y++) \
    {
        for (int x=n_half;x<w-n_half;x++) \
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

    pv->h_param = NLMEANS_H_DEFAULT;
    pv->range   = NLMEANS_RANGE_DEFAULT;
    pv->frames  = NLMEANS_FRAMES_DEFAULT;
    pv->patch   = NLMEANS_PATCH_DEFAULT;

    if( filter->settings )
    {
        sscanf( filter->settings, "%lf:%d:%d:%d", &pv->h_param, &pv->range, &pv->frames, &pv->patch );
    }

    if ( pv->frames < 1 )
    {
        pv->frames = 1;
    }
    if ( pv->frames > NLMEANS_MAX_IMAGES )
    {
        pv->frames = NLMEANS_MAX_IMAGES;
    }

    return 0;
}

static void hb_nlmeans_close( hb_filter_object_t * filter )
{
    hb_filter_private_t * pv = filter->private_data;

    if( !pv )
    {
        return;
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

    out = hb_video_buffer_init( in->f.width, in->f.height );

    // Iterate through planes
    int c;
    for (c = 0; c < 3; c++)
    {
        // Allocate temporary image
        if ( !pv->out_tmp[c] )
        {
            pv->out_tmp[c] = malloc( in->plane[c].stride * in->plane[c].height * sizeof(unsigned char) );
        }

        // Source data, width, height
        unsigned char * src = in->plane[c].data;
        int w = in->plane[c].stride;
        int h = in->plane[c].height;

        // Extra border, new stride, new height
        int border = (pv->range/2 + 15) /16*16;
        int out_w = w + 2*border;
        int out_h = h + 2*border;

        // Allocate bordered image
        unsigned char * mem_start = malloc(out_w * out_h);
        unsigned char * bordered = mem_start + border + out_w*border;

        // Copy main image
        int y;
        for (y = 0; y < h; y++)
        {
            memcpy(bordered + y*out_w, src + y*w, w);
        }

        // Copy borders
        int k;
        // Top, bottom
        for (k = 0; k < border; k++)
        {
            memcpy(bordered - (k+1)*out_w, src, w);
            memcpy(bordered + (h+k)*out_w, src + (h-1)*w, w);
        }
        // Left, right
        for (k = 0; k < border; k++) {
            for (y = -border; y < h + border; y++)
            {
                *(bordered - (k+1) + y*out_w) = bordered[y*out_w];
                *(bordered + (k+w) + y*out_w) = bordered[y*out_w + (w-1)];
            }
        }

        // Process plane
        nlmeans_plane( bordered,
                       out_w,
                       out_h,
                       pv->out_tmp[c],
                       w,
                       h,
                       c,
                       pv->h_param,
                       pv->range,
                       pv->frames,
                       pv->patch );

        out->plane[c].data = pv->out_tmp[c];
        pv->out_tmp[c] = NULL;
    }

    out->s = in->s;
    hb_buffer_move_subs( out, in );

    *buf_out = out;

    return HB_FILTER_OK;
}
