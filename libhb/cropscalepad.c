/* cropscalepad.c

   Copyright (c) 2003-2015 HandBrake Team
   This file is part of the HandBrake source code
   Homepage: <http://handbrake.fr/>.
   It may be used under the terms of the GNU General Public License v2.
   For full terms see the file COPYING file or visit http://www.gnu.org/licenses/gpl-2.0.html
 */

#include "hb.h"
#include "hbffmpeg.h"
#include "common.h"
#include "opencl.h"

struct hb_filter_private_s
{
    hb_job_t            *job;
    int                 width_in;
    int                 height_in;
    int                 pix_fmt;
    int                 pix_fmt_out;
    int                 crop[4];
    int                 width_scaled;
    int                 height_scaled;
    int                 pad[4];
    int                 width_out;
    int                 height_out;

    /* OpenCL/DXVA2 */
    int                 use_dxva;
    int                 use_decomb;
    int                 use_detelecine;
    hb_oclscale_t      *os; //ocl scaler handler

    struct SwsContext * context;
};

static int hb_crop_scale_pad_init(hb_filter_object_t * filter,
                                  hb_filter_init_t * init);

static int hb_crop_scale_pad_work(hb_filter_object_t * filter,
                                  hb_buffer_t ** buf_in,
                                  hb_buffer_t ** buf_out);

static int hb_crop_scale_pad_info(hb_filter_object_t * filter,
                                  hb_filter_info_t * info);

static void hb_crop_scale_pad_close(hb_filter_object_t * filter);

hb_filter_object_t hb_filter_crop_scale_pad =
{
    .id            = HB_FILTER_CROP_SCALE_PAD,
    .enforce_order = 1,
    .name          = "Crop, Scale, Pad",
    .settings      = NULL,
    .init          = hb_crop_scale_pad_init,
    .work          = hb_crop_scale_pad_work,
    .close         = hb_crop_scale_pad_close,
    .info          = hb_crop_scale_pad_info,
};

static int hb_crop_scale_pad_init(hb_filter_object_t * filter,
                                  hb_filter_init_t * init)
{
    filter->private_data = calloc(1, sizeof(struct hb_filter_private_s));
    hb_filter_private_t * pv = filter->private_data;

    // TODO: add pix format option to settings
    pv->job = init->job;
    pv->pix_fmt_out = init->pix_fmt;
    pv->width_in      = init->geometry.width;
    pv->height_in     = init->geometry.height;
    pv->width_scaled  = init->geometry.width  - (init->crop[2] + init->crop[3]);
    pv->height_scaled = init->geometry.height - (init->crop[0] + init->crop[1]);
    pv->width_out     = init->geometry.width  - (init->crop[2] + init->crop[3]) + (init->pad[2] + init->pad[3]);
    pv->height_out    = init->geometry.height - (init->crop[0] + init->crop[1]) + (init->pad[0] + init->pad[1]);

    /* OpenCL/DXVA2 */
    pv->use_dxva       = hb_hwd_enabled(init->job->h);
    pv->use_decomb     = init->job->use_decomb;
    pv->use_detelecine = init->job->use_detelecine;

    if (pv->job->use_opencl && pv->job->title->opencl_support)
    {
        pv->os = (hb_oclscale_t *)(malloc(sizeof(hb_oclscale_t)));
        memset(pv->os, 0, sizeof(hb_oclscale_t));
    }

    memcpy(pv->crop, init->crop, sizeof(int[4]));
    memcpy(pv->pad,  init->pad,  sizeof(int[4]));
    if (filter->settings != NULL)
    {
        sscanf(filter->settings, "%d:%d:%d:%d:%d:%d:%d:%d",
               &pv->width_scaled, &pv->height_scaled,
               &pv->crop[0], &pv->crop[1], &pv->crop[2], &pv->crop[3],
               &pv->width_out, &pv->height_out);
    }

    if (pv->width_out > 0 && pv->width_out > pv->width_scaled)
    {
        pv->pad[2] = (pv->width_out  - pv->width_scaled)  / 2;
        pv->pad[3] = pv->pad[2];
    }
    else
    {
        pv->pad[2] = pv->pad[3] = 0;
        pv->width_out = pv->width_scaled;
    }
    if (pv->height_out > 0 && pv->height_out > pv->height_scaled)
    {
        pv->pad[0] = (pv->height_out - pv->height_scaled) / 2;
        pv->pad[1] = pv->pad[0];
    }
    else
    {
        pv->pad[0] = pv->pad[1] = 0;
        pv->height_out = pv->height_scaled;
    }

    // Set init values so the next stage in the pipline
    // knows what it will be getting
    init->pix_fmt = pv->pix_fmt;
    init->geometry.width = pv->width_out;
    init->geometry.height = pv->height_out;
    memcpy(init->crop, pv->crop, sizeof(int[4]));
    memcpy(init->pad,  pv->pad,  sizeof(int[4]));

    return 0;
}

static int hb_crop_scale_pad_info(hb_filter_object_t * filter,
                                  hb_filter_info_t * info)
{
    hb_filter_private_t * pv = filter->private_data;

    if (pv == NULL)
    {
        return 0;
    }

    // Set info values so the next stage in the pipline
    // knows what it will be getting
    memset(info, 0, sizeof(hb_filter_info_t));
    info->out.pix_fmt = pv->pix_fmt;
    info->out.geometry.width = pv->width_out;
    info->out.geometry.height = pv->height_out;
    memcpy(info->out.crop, pv->crop, sizeof(int[4]));
    memcpy(info->out.pad,  pv->pad,  sizeof(int[4]));

    int width_cropped  = pv->width_in   - (pv->crop[2] + pv->crop[3]);
    int height_cropped = pv->height_in  - (pv->crop[0] + pv->crop[1]);

    sprintf(info->human_readable_desc,
            "source: %d * %d, crop (%d/%d/%d/%d): %d * %d, scale: %d * %d, pad (%d/%d/%d/%d): %d * %d",
            pv->width_in, pv->height_in,
            pv->crop[0], pv->crop[1], pv->crop[2], pv->crop[3],
            width_cropped, height_cropped, pv->width_scaled, pv->height_scaled,
            pv->pad[0], pv->pad[1], pv->pad[2], pv->pad[3],
            pv->width_out, pv->height_out);

    return 0;
}

static void hb_crop_scale_pad_close( hb_filter_object_t * filter )
{
    hb_filter_private_t * pv = filter->private_data;

    if (pv == NULL)
    {
        return;
    }

    /* OpenCL */
    if (pv->job->use_opencl && pv->job->title->opencl_support && pv->os)
    {
        if (hb_ocl != NULL)
        {
            HB_OCL_BUF_FREE(hb_ocl, pv->os->bicubic_x_weights);
            HB_OCL_BUF_FREE(hb_ocl, pv->os->bicubic_y_weights);
        }
        free(pv->os);
    }

    if (pv->context != NULL)
    {
        sws_freeContext(pv->context);
    }

    free(pv);
    filter->private_data = NULL;
}

/* OpenCL */
static hb_buffer_t* crop_scale_pad(hb_filter_private_t * pv, hb_buffer_t * in)
{
    AVPicture           pic_in;
    AVPicture           pic_cropped;
    AVPicture           pic_scaled;
    AVPicture           pic_out;
    hb_buffer_t * scaled;
    scaled = hb_video_buffer_init(pv->width_scaled, pv->height_scaled);
    hb_buffer_t * out;
    out = hb_video_buffer_init(pv->width_out, pv->height_out);

    int pad_color[3];
    pad_color[0] = 0;    // 0 for black, 255 for white
    pad_color[1] = -128; // -128 for no color
    pad_color[2] = -128; // -128 for no color

    hb_avpicture_fill(&pic_in, in);
    hb_avpicture_fill(&pic_scaled, scaled);
    hb_avpicture_fill(&pic_out, out);

    // Crop; this alters the pointer to the data to point to the
    // correct place for cropped frame
    av_picture_crop(&pic_cropped, &pic_in, in->f.fmt,
                    pv->crop[0], pv->crop[2]);

    // Use bicubic OpenCL scaling when selected and when downsampling < 4:1;
    if ((pv->job->use_opencl && pv->job->title->opencl_support) &&
        (pv->width_scaled * 4 > pv->width_in) &&
        (in->cl.buffer != NULL) && (scaled->cl.buffer != NULL))
    {
        /* OpenCL */
        hb_ocl_scale(in, scaled, pv->crop, pv->os);
    }
    else
    {
        if (pv->context   == NULL         ||
            pv->width_in  != in->f.width  ||
            pv->height_in != in->f.height ||
            pv->pix_fmt   != in->f.fmt)
        {
            // Something changed, need a new scaling context.
            if (pv->context != NULL)
            {
                sws_freeContext(pv->context);
            }

            pv->context = hb_sws_get_context(in->f.width  - (pv->crop[2] + pv->crop[3]),
                                             in->f.height - (pv->crop[0] + pv->crop[1]),
                                             in->f.fmt,
                                             scaled->f.width, scaled->f.height, scaled->f.fmt,
                                             SWS_LANCZOS|SWS_ACCURATE_RND);
            pv->width_in  = in->f.width;
            pv->height_in = in->f.height;
            pv->pix_fmt   = in->f.fmt;
        }

        // Scale pic_cropped into pic_render according to the
        // context set up above
        sws_scale(pv->context,
                  (const uint8_t* const*)(pic_cropped.data), pic_cropped.linesize,
                  0, in->f.height - (pv->crop[0] + pv->crop[1]),
                  pic_scaled.data, pic_scaled.linesize);
    }

    // Pad
    av_picture_pad(&pic_out, &pic_scaled,
                   pv->height_out, pv->width_out, scaled->f.fmt,
                   pv->pad[0], pv->pad[1], pv->pad[2], pv->pad[3],
                   pad_color);
    avpicture_free(&pic_scaled);
    free(scaled);

    out->s = in->s;
    hb_buffer_move_subs(out, in);
    return out;
}

static int hb_crop_scale_pad_work(hb_filter_object_t * filter,
                                  hb_buffer_t ** buf_in,
                                  hb_buffer_t ** buf_out)
{
    hb_filter_private_t * pv = filter->private_data;
    hb_buffer_t * in = *buf_in;

    if (in->s.flags & HB_BUF_FLAG_EOF)
    {
        *buf_out = in;
        *buf_in = NULL;
        return HB_FILTER_DONE;
    }

    if (pv == NULL)
    {
        *buf_out = in;
        *buf_in = NULL;
        return HB_FILTER_OK;
    }

    // If scaled width or scaled height were not set, set them now based on the
    // input width & height
    if (pv->width_scaled <= 0)
    {
        pv->width_scaled  = in->f.width  - (pv->crop[2] + pv->crop[3]);
    }
    if (pv->height_scaled <= 0)
    {
        pv->height_scaled = in->f.height - (pv->crop[0] + pv->crop[1]);
    }

    // If padded width or padded height were not set, set them now based on the
    // input width & height
    if (pv->width_out <= 0)
    {
        pv->width_out  = pv->width_scaled;
    }
    if (pv->height_out <= 0)
    {
        pv->height_out = pv->height_scaled;
    }

    /* OpenCL/DXVA2 */
    if ((!pv->use_dxva &&
         !pv->crop[0] && !pv->crop[1] && !pv->crop[2] && !pv->crop[3] &&
         !pv->pad[0]  && !pv->pad[1]  && !pv->pad[2]  && !pv->pad[3]  &&
         in->f.fmt == pv->pix_fmt_out && in->f.width == pv->width_out &&
         in->f.height == pv->height_out) ||
        (pv->use_dxva && !pv->use_decomb && !pv->use_detelecine &&
         !pv->pad[0]  && !pv->pad[1]  && !pv->pad[2]  && !pv->pad[3] &&
         in->f.width  == pv->width_out && in->f.height == pv->height_out))
    {
        *buf_out = in;
        *buf_in  = NULL;
        return HB_FILTER_OK;
    }

    *buf_out = crop_scale_pad(pv, in);

    return HB_FILTER_OK;
}
