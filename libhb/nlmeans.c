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
#include "nlmeans.h"

#define NLMEANS_H_DEFAULT          8
#define NLMEANS_RANGE_DEFAULT      3
#define NLMEANS_FRAMES_DEFAULT     1
#define NLMEANS_PATCH_SIZE_DEFAULT 7

#define NLMEANS_MAX_IMAGES 1 //32
#define EXP_TABLE_SIZE     128

#define ABS(A) ( (A) > 0 ? (A) : -(A) )
#define MIN( a, b ) ( (a) > (b) ? (b) : (a) )
#define CEIL_RSHIFT(a,b) (-((-(a)) >> (b)))

struct hb_filter_private_s
{
    unsigned char * out_tmp[3];
    double h_param;
    int    range;      // search range (must be odd number)
    int    n_frames;   // temporal search depth
    int    patch_size;

    //hb_buffer_t * ref[NLMEANS_MAX_IMAGES]; // frame cache
    //int           n_refs;                  // frame cache count

    int hsub;
    int vsub;
    //NLMeansFunctions func;
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

typedef struct
{
    unsigned char* img; // point to logical origin (0;0)
    int stride;

    int w,h;
    int border; // extra border on all four sides

    unsigned char* mem_start; // start of allocated memory
} MonoImage;

/*
static void store_ref(hb_filter_private_t * pv, hb_buffer_t * b)
{
    for (int i = NLMEANS_MAX_IMAGES-1; i > 0; i--) {
        hb_buffer_close(&pv->ref[i]);
        memmove(&pv->ref[i], &pv->ref[i-1], sizeof(hb_buffer_t *) * 2 );
    }
    pv->ref[0] = b;
}
*/

struct PixelSum
{
    float weight_sum;
    float pixel_sum;
};

/*
static void rawplane2bordered( MonoImage* out,
                               unsigned char * src,
                               unsigned char * dst,
                               int w,
                               int h,
                               hb_filter_private_t * pv )
{

    //hb_filter_private_t * pv = filter->private_data;

    hb_log("rawplane2bordered");
    int border = (pv->range/2 + 15) /16*16;
    int out_w = w + 2*border;
    int out_h = h + 2*border;

    unsigned char * mem_start = malloc(out_w * out_h);
    unsigned char * image = mem_start + border + out_w*border;

    int y;
    int k;

    // main image copy
    for (y = 0; y < h; y++)
    {
        memcpy(image + y*out_w, src + y*w, w);
    }

    // top/bottom borders
    for (k = 0; k < border; k++)
    {
        memcpy(image - (k+1)*out_w, src, w);
    }

    // left/right borders
    for (k = 0; k < border; k++) {
        for (y = -border; y < h + border; y++)
        {
            *(image - (k+1) + y*out_w) = image[y*out_w];
            *(image + (k+w) + y*out_w) = image[y*out_w + (w-1)];
        }
    }

    out->img = image;
    out->mem_start = mem_start;
    hb_log(" out->mem_start %p", out->mem_start);
    out->stride = out_w;
    out->w = w;
    hb_log(" out->w %d", out->w);
    out->h = h;
    hb_log(" out->h %d", out->h);
    out->border = border;

}
*/

/*
static void buildIntegralImage_scalar(uint32_t* integral_image,   int integral_stride,
				      unsigned char* current_image, int current_image_stride,
				      unsigned char* compare_image, int compare_image_stride,
				      int  w,int  h,
				      int dx,int dy)
{
    hb_log("buildintegral");
    memset(integral_image -1 -integral_stride, 0, (w+1)*sizeof(uint32_t));

    for (int y=0;y<h;y++) {
        unsigned char* p1 = current_image +  y    *current_image_stride;
        unsigned char* p2 = compare_image + (y+dy)*compare_image_stride + dx;
        uint32_t* out = integral_image + y*integral_stride -1;

        *out++ = 0;

        for (int x=0;x<w;x++)
        {
            int diff = *p1++ - *p2++;
            *out = *(out-1) + diff * diff;
            out++;
        }

        if (y>0) {
            out = integral_image + y*integral_stride;

            for (int x=0;x<w;x++) {
                *out += *(out - integral_stride);
                out++;
            }
        }
    }
}
*/

static void nlmeans_plane(unsigned char * src,
                          int sw,
                          int sh,
                          unsigned char * dst,
                          int w,
                          int h,
                          int c,
                          //hb_buffer_t * ref,
                          //int n_refs,
                          int h_param,
                          int range,
                          int n_frames,
                          int patch_size )
{

    hb_log("plane");

    int n = (patch_size|1);
    int r = (range     |1);
    hb_log(" patch_size %d", patch_size);
    hb_log(" range %d", range);
    hb_log(" n %d", n);
    hb_log(" r %d", r);

    int n_half = (n-1)/2;
    int r_half = (r-1)/2;


    // alloc memory for temporary pixel sums

    struct PixelSum* tmp_data = calloc(w*h,sizeof(struct PixelSum));
    hb_log(" pixelsum tmp_data %p", tmp_data);

    // allocate integral image

    int integral_stride = w + 2*16;
    uint32_t* integral_mem = malloc( integral_stride * ( h + 1 ) * sizeof(uint32_t) );
    uint32_t* integral = integral_mem + integral_stride + 16;
    hb_log(" integral_stride %d", integral_stride);
    hb_log(" integral_mem %d", integral_mem);
    hb_log(" sizeof(uint32_t) %d", sizeof(uint32_t));
    hb_log(" integral %d", integral);


    // precompute exponential table

    float weight_factor = 1.0/n/n / (h_param * h_param);
    hb_log(" weight_factor %f", weight_factor);

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

    //return;

    int image_idx=0;
    //for (int image_idx=0; image_idx < n_refs; image_idx++)
    //{
        // copy input image

        unsigned char* current_image = src;
        int current_image_stride = sw;

        //unsigned char* compare_image = ref[image_idx].plane[c].data;
        //int compare_image_stride = ref[image_idx].plane[c].stride;
        unsigned char* compare_image = src;
        int compare_image_stride = sw;

        // --- iterate through all displacements ---

        for (int dy=-r_half ; dy<=r_half ; dy++)
            for (int dx=-r_half ; dx<=r_half ; dx++)
            {
                // special, simple implementation for no shift (no difference -> weight 1)

                if (dx==0 && dy==0 && image_idx==0) {
//#pragma omp parallel for
                    for (int y=n_half;y<h-n+n_half;y++) {
                        for (int x=n_half;x<w-n+n_half;x++) {
                            tmp_data[y*w+x].weight_sum += 1;
                            tmp_data[y*w+x].pixel_sum  += current_image[y*current_image_stride+x];
                        }
                    }

                    continue;
                }


                // --- regular case ---

                /*
                pv->func.buildIntegralImage(integral,integral_stride,
                                             current_image, current_image_stride,
                                             compare_image, compare_image_stride,
                                             w,h, dx,dy);
                hb_log("buildintegral->mono");
                */

                hb_log("buildintegral-inline");
                memset(integral -1 -integral_stride, 0, (w+1)*sizeof(uint32_t));

                hb_log(" loop 1");
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

                hb_log(" loop 2");
//#pragma omp parallel for
                for (int y=0;y<=h-n;y++) {
                    uint32_t* integral_ptr1 = integral+(y  -1)*integral_stride-1;
                    uint32_t* integral_ptr2 = integral+(y+n-1)*integral_stride-1;

                    for (int x=0;x<=w-n;x++) {
                        int xc = x+n_half;
                        int yc = y+n_half;

                        // patch difference

                        int diff = (uint32_t)(integral_ptr2[n] - integral_ptr2[0] - integral_ptr1[n] + integral_ptr1[0]);


                        // sum pixel with weight

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
    //}



    // --- fill output image ---

    // copy border area

    {

        for (int y=0;       y<n_half  ;y++)
        {
            hb_log(" out %p", dst+y*w);
            hb_log(" in  %p", src+y*sw);
            hb_log(" bytes %d", sw);
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

    // output main image

    for (int y=n_half;y<h-n_half;y++) \
    {
        for (int x=n_half;x<w-n_half;x++) \
        {
            *(dst+y*w+x) = tmp_data[y*w+x].pixel_sum / tmp_data[y*w+x].weight_sum;
        }
    }

    free(tmp_data);
    free(integral_mem);

    hb_log("return");

}

static int hb_nlmeans_init( hb_filter_object_t * filter,
                            hb_filter_init_t * init )
{
    filter->private_data = calloc( sizeof(struct hb_filter_private_s), 1 );
    hb_filter_private_t * pv = filter->private_data;

    pv->h_param    = NLMEANS_H_DEFAULT;
    pv->range      = NLMEANS_RANGE_DEFAULT;
    pv->n_frames   = NLMEANS_FRAMES_DEFAULT;
    pv->patch_size = NLMEANS_PATCH_SIZE_DEFAULT;

    if( filter->settings )
    {
        sscanf( filter->settings, "%lf:%d:%d:%d", &pv->h_param, &pv->range, &pv->n_frames, &pv->patch_size );
    }

    if ( pv->n_frames < 1 )
    {
        pv->n_frames = 1;
    }
    if ( pv->n_frames > NLMEANS_MAX_IMAGES )
    {
        pv->n_frames = NLMEANS_MAX_IMAGES;
    }

    //pv->func.buildIntegralImage = buildIntegralImage_scalar;

    // amount to shift width/height right to find the corresponding chroma dimension
    // since handbrake only operates in one colorspace, for now hard code these values
    pv->hsub = 1;
    pv->vsub = 1;

    //pv->n_refs = 0;

    return 0;
}

static void hb_nlmeans_close( hb_filter_object_t * filter )
{
    hb_log("close");
    hb_filter_private_t * pv = filter->private_data;

    if( !pv )
    {
        return;
    }
/*
    // Close reference buffers
    int ii;
    for (ii = 0; ii < NLMEANS_MAX_IMAGES; ii++)
    {
        hb_buffer_close(&pv->ref[ii]);
    }
*/
    free( pv );
    filter->private_data = NULL;
}

static int hb_nlmeans_work( hb_filter_object_t * filter,
                            hb_buffer_t ** buf_in,
                            hb_buffer_t ** buf_out )
{
    hb_log("work");
    hb_filter_private_t * pv = filter->private_data;
    hb_buffer_t * in = *buf_in, * out;

    hb_log(" pv->h_param %f", pv->h_param);
    hb_log(" pv->range %d", pv->range);
    hb_log(" pv->n_frames %d", pv->n_frames);
    hb_log(" pv->patch_size %d", pv->patch_size);

    if ( in->size <= 0 )
    {
        *buf_out = in;
        *buf_in = NULL;
        return HB_FILTER_DONE;
    }

    /* Store current frame in nlmeans cache */
    //*buf_in = NULL;
    //store_ref(pv, in);
    //if (pv->n_refs < NLMEANS_MAX_IMAGES)
    //{
    //    pv->n_refs = pv->n_refs + 1;
    //}

    out = hb_video_buffer_init( in->f.width, in->f.height );

    int c;
    for (c = 0; c < 3; c++) {
        if ( !pv->out_tmp[c] )
        {
            pv->out_tmp[c] = malloc( in->plane[c].stride * in->plane[c].height * sizeof(unsigned char) );
        }

        hb_log("border");

        unsigned char * src = in->plane[c].data;
        int w = in->plane[c].stride;
        int h = in->plane[c].height;

        int border = (pv->range/2 + 15) /16*16;
        int out_w = w + 2*border;
        int out_h = h + 2*border;
        unsigned char * mem_start = malloc(out_w * out_h);
        unsigned char * bordered = mem_start + border + out_w*border;

        int y;
        int k;

        // main image copy
        for (y = 0; y < h; y++)
        {
            memcpy(bordered + y*out_w, src + y*w, w);
        }

        // top/bottom borders
        for (k = 0; k < border; k++)
        {
            memcpy(bordered - (k+1)*out_w, src, w);
        }

        // left/right borders
        for (k = 0; k < border; k++) {
            for (y = -border; y < h + border; y++)
            {
                *(bordered - (k+1) + y*out_w) = bordered[y*out_w];
                *(bordered + (k+w) + y*out_w) = bordered[y*out_w + (w-1)];
            }
        }

        hb_log(" mem_start %p", mem_start);
        hb_log(" w %d", w);
        hb_log(" h %d", h);
        hb_log(" border %d", border);

        nlmeans_plane( bordered,
                       out_w,
                       out_h,
                       pv->out_tmp[c],
                       w,
                       h,
                       c,
                       pv->h_param,
                       pv->range,
                       pv->n_frames,
                       pv->patch_size );

/*
        nlmeans_plane( in->plane[c].data,
                       pv->out_tmp[c],
                       in->plane[c].stride,
                       in->plane[c].height,
                       c,
                       //pv->ref,
                       //pv->n_refs,
                       pv->h_param,
                       pv->range,
                       pv->n_frames,
                       pv->patch_size );
*/

        //MonoImage* images[NLMEANS_MAX_IMAGES];
        //int images[NLMEANS_MAX_IMAGES];
        //int i;
        //for (i=0; i < pv->n_refs; i++) {
            //rawplane2bordered(&images[i], pv->ref[i]->plane[c].data, pv->out_tmp[c], in->plane[c].stride, in->plane[c].height, &pv);
        //    rawplane2bordered(&images[i], pv->ref[i]->plane[c].data, pv->out_tmp[c], in->plane[c].stride, in->plane[c].height, &pv);
        //}
        //hb_log("images[0]->w %d", images[0]->w);

        //NLMeans_mono_multi(pv->out_tmp[c], in->plane[c].stride,
        //                   images, i, &pv);

        out->plane[c].data = pv->out_tmp[c];
        //out->plane[c].data = in->plane[c].data; // DEBUG
        pv->out_tmp[c] = NULL;
    }

    //hb_buffer_move_subs( out, pv->ref[0] );
    //out->s = pv->ref[0]->s;
    out->s = in->s;
    hb_buffer_move_subs( out, in );

    *buf_out = out;

    return HB_FILTER_OK;
}
