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
#include "mpeg2dec/mpeg2.h"
#include "nlmeans.h"

#define NLMEANS_H_DEFAULT          8
#define NLMEANS_RANGE_DEFAULT      3
#define NLMEANS_FRAMES_DEFAULT     2
#define NLMEANS_PATCH_SIZE_DEFAULT 7

#define NLMEANS_MAX_IMAGES 32
#define EXP_TABLE_SIZE     128

#define ABS(A) ( (A) > 0 ? (A) : -(A) )
#define MIN( a, b ) ( (a) > (b) ? (b) : (a) )

typedef struct
{
    double h_param;
    int    range;      // search range (must be odd number)
    int    n_frames;   // temporal search depth
    int    patch_size;
} NLMeansParams;

struct hb_filter_private_s
{
    //short            nlmeans_coef[4][512*16];
    unsigned short * nlmeans_line;
    unsigned short * nlmeans_frame[3];
    //double h_param;
    //int    range;      // search range (must be odd number)
    //int    n_frames;   // temporal search depth
    //int    patch_size;
    NLMeansParams param;
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
    uint8_t* img; // point to logical origin (0;0)
    int stride;

    int w,h;
    int border; // extra border on all four sides

    uint8_t* mem_start; // start of allocated memory
} MonoImage;

typedef struct
{
    MonoImage   plane[3];
} ColorImage;

static void free_mono_image(MonoImage* img)
{
    if (img->mem_start) {
        free(img->mem_start);
        img->mem_start=NULL;
        img->img=NULL;
    }
}

static void free_color_image(ColorImage* img)
{
    for (int c=0;c<3;c++) {
        free_mono_image(&(img->plane[c]));
    }
}

typedef struct {
    const AVClass *class;

    int hsub,vsub;

    NLMeansParams param;

    ColorImage images[NLMEANS_MAX_IMAGES];
    int        image_available[NLMEANS_MAX_IMAGES];

    NLMeansFunctions func;

} NLMContext;

struct PixelSum
{
    float weight_sum;
    float pixel_sum;
};

static void alloc_and_copy_image_with_border(MonoImage* out_img,
                                             const uint8_t* input_image, int input_stride,
                                             int w,int h, int req_border)
{
    const int border = (req_border+15)/16*16;
    const int out_stride = (w+2*border);
    const int out_total_height = (h+2*border);

    uint8_t* const memory_ptr   = (uint8_t*)malloc(out_stride*out_total_height);
    uint8_t* const output_image = memory_ptr + border + border*out_stride; // logical output image origin (0,0)


    // copy main image content

    for (int y=0;y<h;y++) {
        memcpy(output_image + y*out_stride, input_image+y*input_stride, w);
    }

    // top/bottom borders

    for (int k=0;k<border;k++) {
        memcpy(output_image-(k+1)*out_stride, input_image, w);
        memcpy(output_image+(h+k)*out_stride, input_image+(h-1)*input_stride, w);
    }

    // left/right borders

    for (int k=0;k<border;k++) {
        for (int y=-border;y<h+border;y++)
        {
            *(output_image  -k-1+y*out_stride) = output_image[y*out_stride];
            *(output_image+w+k  +y*out_stride) = output_image[y*out_stride+w-1];
        }
    }

    out_img->img = output_image;
    out_img->mem_start = memory_ptr;
    out_img->stride = out_stride;
    out_img->w = w;
    out_img->h = h;
    out_img->border = border;
}

static void buildIntegralImage_scalar(uint32_t* integral_image,   int integral_stride,
				      const uint8_t* current_image, int current_image_stride,
				      const uint8_t* compare_image, int compare_image_stride,
				      int  w,int  h,
				      int dx,int dy)
{
    memset(integral_image -1 -integral_stride, 0, (w+1)*sizeof(uint32_t));

    for (int y=0;y<h;y++) {
        const uint8_t* p1 = current_image +  y    *current_image_stride;
        const uint8_t* p2 = compare_image + (y+dy)*compare_image_stride + dx;
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

static void NLMeans_mono_multi(uint8_t* out, int out_stride,
			       const MonoImage*const* images, int n_images,
			       const NLMContext* ctx)
{
    const int w = images[0]->w;
    const int h = images[0]->h;

    const int n = (ctx->param.patch_size|1);
    const int r = (ctx->param.range     |1);

    const int n_half = (n-1)/2;
    const int r_half = (r-1)/2;


    // alloc memory for temporary pixel sums

    struct PixelSum* const tmp_data = (struct PixelSum*)calloc(w*h,sizeof(struct PixelSum));


    // allocate integral image

    const int integral_stride = w+2*16;
    uint32_t* const integral_mem = (uint32_t*)malloc( integral_stride*(h+1)*sizeof(uint32_t) );
    uint32_t* const integral = integral_mem + integral_stride + 16;


    // precompute exponential table

    const float weight_factor = 1.0/n/n / (ctx->param.h_param * ctx->param.h_param);

    const int table_size=EXP_TABLE_SIZE;
    const float min_weight_in_table = 0.0005;

    float exptable[EXP_TABLE_SIZE];

    const float stretch = table_size/ (-log(min_weight_in_table));
    const float weight_fact_table = weight_factor*stretch;
    const int diff_max = table_size/weight_fact_table;

    for (int i=0;i<table_size;i++) {
        exptable[i] = exp(-i/stretch);
    }
    exptable[table_size-1]=0;



    for (int image_idx=0; image_idx<n_images; image_idx++)
    {
        // copy input image

        const uint8_t* current_image = images[0]->img;
        int current_image_stride = images[0]->stride;

        const uint8_t* compare_image = images[image_idx]->img;
        int compare_image_stride = images[image_idx]->stride;


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

                ctx->func.buildIntegralImage(integral,integral_stride,
                                             current_image, current_image_stride,
                                             compare_image, compare_image_stride,
                                             w,h, dx,dy);

//#pragma omp parallel for
                for (int y=0;y<=h-n;y++) {
                    const uint32_t* integral_ptr1 = integral+(y  -1)*integral_stride-1;
                    const uint32_t* integral_ptr2 = integral+(y+n-1)*integral_stride-1;

                    for (int x=0;x<=w-n;x++) {
                        const int xc = x+n_half;
                        const int yc = y+n_half;

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
    }



    // --- fill output image ---

    // copy border area

    {
        const uint8_t* in  = images[0]->img;
        int orig_in_stride = images[0]->stride;

        for (int y=0;       y<n_half  ;y++) { memcpy(out+y*out_stride, in+y*orig_in_stride, w); }
        for (int y=h-n_half;y<h       ;y++) { memcpy(out+y*out_stride, in+y*orig_in_stride, w); }
        for (int y=n_half  ;y<h-n_half;y++) {
            memcpy(out+y*out_stride,          in+y*orig_in_stride,          n_half);
            memcpy(out+y*out_stride+w-n_half, in+y*orig_in_stride+w-n_half, n_half);
        }
    }

    // output main image

    for (int y=n_half;y<h-n_half;y++) {
        for (int x=n_half;x<w-n_half;x++) {
            *(out+y*out_stride+x) = tmp_data[y*w+x].pixel_sum / tmp_data[y*w+x].weight_sum;
        }
    }

    free(tmp_data);
    free(integral_mem);
}


static void NLMeans_color_auto(uint8_t** out, int* out_stride,
			       const ColorImage* img, // function takes ownership
			       NLMContext* ctx)
{
    //assert(ctx->param.n_frames >= 1);
    //assert(ctx->param.n_frames <= NLMEANS_MAX_IMAGES);

    // free oldest image
/*
    free_color_image(&ctx->images[ctx->param.n_frames-1]);

    // shift old images one down and put new image into entry [0]

    for (int i=ctx->param.n_frames-1; i>0; i--) {
        ctx->images[i] = ctx->images[i-1];
        ctx->image_available[i] = ctx->image_available[i-1];
    }

    ctx->images[0] = *img;
    ctx->image_available[0] = 1;


    // process color planes separately

    for (int c=0;c<3;c++)
        if (ctx->images[0].plane[c].img != NULL)
        {
            const MonoImage* images[NLMEANS_MAX_IMAGES];
            int i;
            for (i=0; ctx->image_available[i]; i++) {
                images[i] = &ctx->images[i].plane[c];
            }

            NLMeans_mono_multi(out[c], out_stride[c],
                               images, i, ctx);
        }
*/
}

static void nlmeans_precalc_coef( short * ct,
                                 double dist25 )
{
    int i;
    double gamma, simil, c;

    gamma = log( 0.25 ) / log( 1.0 - MIN(dist25,252.0)/255.0 - 0.00001 );

    for( i = -255*16; i <= 255*16; i++ )
    {
        /* nlmeans_lowpass_mul() truncates (not rounds) the diff, use +15/32 as midpoint */
        double f = (i + 15.0/32.0) / 16.0;
        simil = 1.0 - ABS(f) / 255.0;
        c = pow(simil, gamma) * 256.0 * f;
        ct[16*256+i] = (c<0) ? (c-0.5) : (c+0.5);
    }

    ct[0] = (dist25 != 0);
}

static inline unsigned int nlmeans_lowpass_mul( int prev_mul,
                                               int curr_mul,
                                               short * coef )
{
    int d = (prev_mul - curr_mul)>>4;
    return curr_mul + coef[d];
}

static void nlmeans_denoise_temporal( unsigned char * frame_src,
                                     unsigned char * frame_dst,
                                     unsigned short * frame_ant,
                                     int w, int h,
                                     short * temporal)
{
    int x, y;
    unsigned int tmp;

    temporal += 0x1000;

    for( y = 0; y < h; y++ )
    {
        for( x = 0; x < w; x++ )
        {
            frame_ant[x] = tmp = nlmeans_lowpass_mul( frame_ant[x],
                                                     frame_src[x]<<8,
                                                     temporal );
            frame_dst[x] = (tmp+0x7F)>>8;
        }

        frame_src += w;
        frame_dst += w;
        frame_ant += w;
    }
}

static void nlmeans_denoise_spatial( unsigned char * frame_src,
                                    unsigned char * frame_dst,
                                    unsigned short * line_ant,
                                    unsigned short * frame_ant,
                                    int w, int h,
                                    short * spatial,
                                    short * temporal )
{
    int x, y;
    unsigned int pixel_ant;
    unsigned int tmp;

    spatial  += 0x1000;
    temporal += 0x1000;

    /* First line has no top neighbor. Only left one for each tmp and last frame */
    pixel_ant = frame_src[0]<<8;
    for ( x = 0; x < w; x++)
    {
        line_ant[x] = tmp = pixel_ant = nlmeans_lowpass_mul( pixel_ant,
                                                            frame_src[x]<<8,
                                                            spatial );
        frame_ant[x] = tmp = nlmeans_lowpass_mul( frame_ant[x],
                                                 tmp,
                                                 temporal );
        frame_dst[x] = (tmp+0x7F)>>8;
    }

    for( y = 1; y < h; y++ )
    {
        frame_src += w;
        frame_dst += w;
        frame_ant += w;
        pixel_ant = frame_src[0]<<8;
        for ( x = 0; x < w-1; x++ )
        {
            line_ant[x] = tmp =  nlmeans_lowpass_mul( line_ant[x],
                                                     pixel_ant,
                                                     spatial );
            pixel_ant =          nlmeans_lowpass_mul( pixel_ant,
                                                     frame_src[x+1]<<8,
                                                     spatial );
            frame_ant[x] = tmp = nlmeans_lowpass_mul( frame_ant[x],
                                                     tmp,
                                                     temporal );
            frame_dst[x] = (tmp+0x7F)>>8;
        }
        line_ant[x] = tmp =  nlmeans_lowpass_mul( line_ant[x],
                                                 pixel_ant,
                                                 spatial );
        frame_ant[x] = tmp = nlmeans_lowpass_mul( frame_ant[x],
                                                 tmp,
                                                 temporal );
        frame_dst[x] = (tmp+0x7F)>>8;
    }
}

static void nlmeans_denoise( unsigned char * frame_src,
                            unsigned char * frame_dst,
                            unsigned short * line_ant,
                            unsigned short ** frame_ant_ptr,
                            int w,
                            int h,
                            short * spatial,
                            short * temporal )
{
    int x, y;
    unsigned short* frame_ant = (*frame_ant_ptr);

    if( !frame_ant)
    {
        unsigned char * src = frame_src;
        (*frame_ant_ptr) = frame_ant = malloc( w*h*sizeof(unsigned short) );
        for ( y = 0; y < h; y++, frame_src += w, frame_ant += w )
        {
            for( x = 0; x < w; x++ )
            {
                frame_ant[x] = frame_src[x]<<8;
            }
        }
        frame_src = src;
        frame_ant = *frame_ant_ptr;
    }

    /* If no spatial coefficients, do temporal denoise only */
    if( spatial[0] )
    {
        nlmeans_denoise_spatial( frame_src,
                                frame_dst,
                                line_ant,
                                frame_ant,
                                w, h,
                                spatial,
                                temporal );
    }
    else
    {
        nlmeans_denoise_temporal( frame_src,
                                 frame_dst,
                                 frame_ant,
                                 w, h,
                                 temporal);
    }
}

static int hb_nlmeans_init( hb_filter_object_t * filter,
                            hb_filter_init_t * init )
{
    filter->private_data = calloc( sizeof(struct hb_filter_private_s), 1 );
    hb_filter_private_t * pv = filter->private_data;

    pv->param.h_param    = NLMEANS_H_DEFAULT;
    pv->param.range      = NLMEANS_RANGE_DEFAULT;
    pv->param.n_frames   = NLMEANS_FRAMES_DEFAULT;
    pv->param.patch_size = NLMEANS_PATCH_SIZE_DEFAULT;

    if( filter->settings )
    {
        sscanf( filter->settings, "%lf:%d:%d:%d", &pv->param.h_param, &pv->param.range, &pv->param.n_frames, &pv->param.patch_size );
    }

    if ( pv->param.n_frames < 1 )
    {
        pv->param.n_frames = 1;
    }
    if ( pv->param.n_frames > NLMEANS_MAX_IMAGES )
    {
        pv->param.n_frames = NLMEANS_MAX_IMAGES;
    }

    NLMContext * nlm = filter->private_data;
    
    int i;
    for ( i = 0; i < NLMEANS_MAX_IMAGES; i++ )
    {
        nlm->image_available[i] = 0;
    }
    nlm->func.buildIntegralImage = buildIntegralImage_scalar;

    return 0;
}

static void hb_nlmeans_close( hb_filter_object_t * filter )
{
    hb_filter_private_t * pv = filter->private_data;

    if( !pv )
    {
        return;
    }

	if( pv->nlmeans_frame[0] )
    {
        free( pv->nlmeans_frame[0] );
        pv->nlmeans_frame[0] = NULL;
    }
	if( pv->nlmeans_frame[1] )
    {
        free( pv->nlmeans_frame[1] );
        pv->nlmeans_frame[1] = NULL;
    }
	if( pv->nlmeans_frame[2] )
    {
        free( pv->nlmeans_frame[2] );
        pv->nlmeans_frame[2] = NULL;
    }

    NLMContext * nlm = filter->private_data;

    for (int i=0;i<NLMEANS_MAX_IMAGES;i++) {
        if (nlm->image_available[i]) {
            free_color_image(&(nlm->images[i]));

            nlm->image_available[i] = 0;
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

    out = hb_video_buffer_init( in->f.width, in->f.height );

    if( !pv->nlmeans_line )
    {
        pv->nlmeans_line = malloc( in->plane[0].stride * sizeof(unsigned short) );
    }

    NLMContext * nlm = filter->private_data;

    ColorImage bordered_image;

    // extend image with border
    int c;
    for (c = 0; c < 3; c++) {
        //int w = FF_CEIL_RSHIFT(in->width,  (!!c * nlm->hsub));
        //int h = FF_CEIL_RSHIFT(in->height, (!!c * nlm->vsub));
        int border = nlm->param.range/2;

        /*
        alloc_and_copy_image_with_border(&bordered_image.plane[c],
                                         in->plane[c].data, pv->nlmeans_line,
                                         in->plane[c].stride,in->plane[c].height,border);
        */

        
    }

    const ColorImage * img = &bordered_image;

    // free oldest image
    free_color_image(&nlm->images[nlm->param.n_frames-1]);

    // shift old images one down and put new image into entry [0]
    for (int i=nlm->param.n_frames-1; i>0; i--) {
        nlm->images[i] = nlm->images[i-1];
        nlm->image_available[i] = nlm->image_available[i-1];
    }

    nlm->images[0] = *img;
    nlm->image_available[0] = 1;


    // process color planes separately
    for (int c=0;c<3;c++)
        if (nlm->images[0].plane[c].img != NULL)
        {
            const MonoImage* images[NLMEANS_MAX_IMAGES];
            int i;
            for (i=0; nlm->image_available[i]; i++) {
                images[i] = &nlm->images[i].plane[c];
            }

            NLMeans_mono_multi(out->plane[c].data, in->plane[c].stride,
                               images, i, nlm);
        }

/*
    NLMeans_color_auto(out->data, out->linesize,
		       &bordered_image,
		       nlm);
*/

/*
    int c;

    for ( c = 0; c < 3; c++ )
    {
        nlmeans_denoise( in->plane[c].data,
                        out->plane[c].data,
                        pv->nlmeans_line,
                        &pv->nlmeans_frame[c],
                        in->plane[c].stride,
                        in->plane[c].height,
                        pv->nlmeans_coef[c?2:0],
                        pv->nlmeans_coef[c?3:1] );
    }
*/
    out->s = in->s;
    hb_buffer_move_subs( out, in );

    *buf_out = out;

    return HB_FILTER_OK;
}
