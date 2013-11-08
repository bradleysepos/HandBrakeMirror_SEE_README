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

#include <stdint.h>

typedef struct
{
    void (*buildIntegralImage)(uint32_t* integral,   int integral_stride32,
                               const uint8_t* currimage, int currstride,
                               const uint8_t* image, int stride,
                               int  w,int  h,
                               int dx,int dy);
} NLMeansFunctions;