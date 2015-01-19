﻿// --------------------------------------------------------------------------------------------------------------------
// <copyright file="AudioList.cs" company="HandBrake Project (http://handbrake.fr)">
//   This file is part of the HandBrake source code - It may be used under the terms of the GNU General Public License.
// </copyright>
// <summary>
//   The audio list.
// </summary>
// --------------------------------------------------------------------------------------------------------------------

namespace HandBrake.Interop.Json.Encode
{
    /// <summary>
    /// The audio list.
    /// </summary>
    public class AudioList
    {
        /// <summary>
        /// Gets or sets the bitrate.
        /// </summary>
        public int? Bitrate { get; set; }

        /// <summary>
        /// Gets or sets the compression level.
        /// </summary>
        public double? CompressionLevel { get; set; }

        /// <summary>
        /// Gets or sets the drc.
        /// </summary>
        public double DRC { get; set; }

        /// <summary>
        /// Gets or sets the encoder.
        /// </summary>
        public int Encoder { get; set; }

        /// <summary>
        /// Gets or sets the gain.
        /// </summary>
        public double Gain { get; set; }

        /// <summary>
        /// Gets or sets the mixdown.
        /// </summary>
        public int Mixdown { get; set; }

        /// <summary>
        /// Gets or sets a value indicating whether normalize mix level.
        /// </summary>
        public bool NormalizeMixLevel { get; set; }

        /// <summary>
        /// Gets or sets the quality.
        /// </summary>
        public double? Quality { get; set; }

        /// <summary>
        /// Gets or sets the samplerate.
        /// </summary>
        public int Samplerate { get; set; }

        /// <summary>
        /// Gets or sets the Name of the audio track.
        /// </summary>
        public string Name { get; set; }

        /// <summary>
        /// Gets or sets the track.
        /// </summary>
        public int Track { get; set; }
    }
}