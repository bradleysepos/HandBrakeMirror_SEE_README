﻿// --------------------------------------------------------------------------------------------------------------------
// <copyright file="AnamorphicFactory.cs" company="HandBrake Project (http://handbrake.fr)">
//   This file is part of the HandBrake source code - It may be used under the terms of the GNU General Public License.
// </copyright>
// <summary>
//   The Anamorphic factory.
// </summary>
// --------------------------------------------------------------------------------------------------------------------

namespace HandBrake.Interop.Json.Factories
{
    using System;
    using System.Collections.Generic;
    using System.Runtime.InteropServices;

    using HandBrake.Interop.HbLib;
    using HandBrake.Interop.Json.Anamorphic;
    using HandBrake.Interop.Model;
    using HandBrake.Interop.Model.Encoding;
    using HandBrake.Interop.Model.Scan;

    using Newtonsoft.Json;

    /// <summary>
    /// The anamorphic factory.
    /// </summary>
    public class AnamorphicFactory
    {
        /// <summary>
        /// The keep setting.
        /// </summary>
        public enum KeepSetting
        {
            HB_KEEP_WIDTH = 0x01,
            HB_KEEP_HEIGHT = 0x02,
            HB_KEEP_DISPLAY_ASPECT = 0x04
        }

        /// <summary>
        /// The create geometry.
        /// </summary>
        /// <param name="job">
        /// The job.
        /// </param>
        /// <param name="title">
        /// The current title.
        /// </param>
        /// <param name="keepWidthOrHeight">
        /// Keep Width or Height. (Not Display Aspect)
        /// </param>
        /// <returns>
        /// The <see cref="Scan.Geometry"/>.
        /// </returns>
        public static Geometry CreateGeometry(EncodeJob job, Title title, KeepSetting keepWidthOrHeight) // Todo remove the need for these objects. Should use simpler objects.
        {
            int settingMode = (int)keepWidthOrHeight + (job.KeepDisplayAspect ? 0x04 : 0);

            // Sanatise the Geometry First.
            AnamorphicGeometry anamorphicGeometry = new AnamorphicGeometry
            {
                SourceGeometry = new SourceGeometry
                                 {
                                    Width = title.Resolution.Width,
                                    Height = title.Resolution.Height,
                                    PAR = new PAR { Num = title.ParVal.Width, Den = title.ParVal.Height }
                                 },
                DestSettings = new DestSettings
                               {
                                    AnamorphicMode = (int)job.Anamorphic,
                                    Geometry = { 
                                                Width = job.Width, Height = job.Height,
                                                PAR = new PAR
                                                      {
                                                          Num = job.Anamorphic != Anamorphic.Custom ? title.ParVal.Width : job.PixelAspectX,
                                                          Den = job.Anamorphic != Anamorphic.Custom ? title.ParVal.Height : job.PixelAspectY,
                                                      } 
                                               },
                                    Keep = settingMode,
                                    Crop = new List<int> { job.Cropping.Top, job.Cropping.Bottom, job.Cropping.Left, job.Cropping.Right },
                                    Modulus = job.Modulus,
                                    MaxWidth = job.MaxWidth,
                                    MaxHeight = job.MaxHeight,
                                    ItuPAR = false
                               }
            };

            if (job.Anamorphic == Anamorphic.Custom)
            {
                anamorphicGeometry.DestSettings.Geometry.PAR = new PAR { Num = job.PixelAspectX, Den = job.PixelAspectY };
            }
            else
            {
                anamorphicGeometry.DestSettings.Geometry.PAR = new PAR { Num = title.ParVal.Width, Den = title.ParVal.Height };
            }

            string encode = JsonConvert.SerializeObject(anamorphicGeometry, Formatting.Indented, new JsonSerializerSettings { NullValueHandling = NullValueHandling.Ignore });
            IntPtr json = HBFunctions.hb_set_anamorphic_size_json(Marshal.StringToHGlobalAnsi(encode));
            string result = Marshal.PtrToStringAnsi(json);
            AnamorphicResult resultGeometry = JsonConvert.DeserializeObject<AnamorphicResult>(result);

            // Setup the Destination Gemotry.
            Geometry geometry = new Geometry
                {
                    Width = resultGeometry.Width,
                    Height = resultGeometry.Height,
                    PAR = resultGeometry.PAR
                };
            return geometry;
        }
    }
}
