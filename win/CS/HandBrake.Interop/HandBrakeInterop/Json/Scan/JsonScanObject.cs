﻿// --------------------------------------------------------------------------------------------------------------------
// <copyright file="JsonScanObject.cs" company="HandBrake Project (http://handbrake.fr)">
//   This file is part of the HandBrake source code - It may be used under the terms of the GNU General Public License.
// </copyright>
// <summary>
//   The root object.
// </summary>
// --------------------------------------------------------------------------------------------------------------------

namespace HandBrake.Interop.Json.Scan
{
    using System.Collections.Generic;

    /// <summary>
    /// The root object.
    /// </summary>
    internal class JsonScanObject
    {
        /// <summary>
        /// Gets or sets the main feature.
        /// </summary>
        public int MainFeature { get; set; }

        /// <summary>
        /// Gets or sets the title list.
        /// </summary>
        public List<TitleList> TitleList { get; set; }
    }
}