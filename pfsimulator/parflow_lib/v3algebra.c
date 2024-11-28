/*BHEADER**********************************************************************
*
*  Copyright (c) 1995-2024, Lawrence Livermore National Security,
*  LLC. Produced at the Lawrence Livermore National Laboratory. Written
*  by the Parflow Team (see the CONTRIBUTORS file)
*  <parflow@lists.llnl.gov> CODE-OCEC-08-103. All rights reserved.
*
*  This file is part of Parflow. For details, see
*  http://www.llnl.gov/casc/parflow
*
*  Please read the COPYRIGHT file or Our Notice and the LICENSE file
*  for the GNU Lesser General Public License.
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License (as published
*  by the Free Software Foundation) version 2.1 dated February 1999.
*
*  This program is distributed in the hope that it will be useful, but
*  WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms
*  and conditions of the GNU General Public License for more details.
*
*  You should have received a copy of the GNU Lesser General Public
*  License along with this program; if not, write to the Free Software
*  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
*  USA
**********************************************************************EHEADER*/

#include "v3algebra.h"

v3 v3_scale(const double c, const v3 A)
{
  v3 res = {c * A.u, c * A.v, c * A.w};
  return res;
}

v3 v3_add(const v3 A, const v3 B)
{
  v3 res = {A.u + B.u, A.v + B.v, A.w + B.w};
  return res;
}

v3 v3_subtract(const v3 A, const v3 B)
{
  v3 res = {A.u - B.u, A.v - B.v, A.w - B.w};
  return res;
}

double v3_dot(const v3 A, const v3 B)
{
  return (A.u * B.u + A.v * B.v + A.w * B.w);
}

v3 v3_cross(const v3 A, const v3 B)
{
  v3 res = {
    A.v * B.w - A.w * B.v,
    A.w * B.u - A.u * B.w,
    A.u * B.v - A.v * B.u
  };
  return res;
}


m3 m3_add(const m3 A, const m3 B)
{
  m3 res = {
    A.uu + B.uu, A.uv + B.uv, A.uw + B.uw,
    A.vu + B.vu, A.vv + B.vv, A.vw + B.vw,
    A.wu + B.wu, A.wv + B.wv, A.ww + B.ww
  };
  return res;
}

m3 m3_subtract(const m3 A, const m3 B)
{
  m3 res = {
    A.uu - B.uu, A.uv - B.uv, A.uw - B.uw,
    A.vu - B.vu, A.vv - B.vv, A.vw - B.vw,
    A.wu - B.wu, A.wv - B.wv, A.ww - B.ww
  };
  return res;
}

m3 m3_multiply(const m3 A, const m3 B)
{
  m3 res = {
    A.uu * B.uu + A.uv * B.vu + A.uw * B.wu, // uu
    A.uu * B.uv + A.uv * B.vv + A.uw * B.wv, // uv
    A.uu * B.uw + A.uv * B.vw + A.uw * B.ww, // uw
    A.vu * B.uu + A.vv * B.vu + A.vw * B.wu, // vu
    A.vu * B.uv + A.vv * B.vv + A.vw * B.wv, // vv
    A.vu * B.uw + A.vv * B.vw + A.vw * B.ww, // vw
    A.wu * B.uu + A.wv * B.vu + A.ww * B.wu, // wu
    A.wu * B.uv + A.wv * B.vv + A.ww * B.wv, // wv
    A.wu * B.uw + A.wv * B.vw + A.ww * B.ww  // ww
  };
  return res;
}

m3 m3_inverse(const m3 A)
{
  double idet = 1.0 / m3_determinant(A);
  m3 res = {
    idet * (A.vv * A.ww - A.vw * A.wv), // uu
    idet * (A.uw * A.wv - A.uv * A.ww), // uv
    idet * (A.uv * A.vw - A.uw * A.vv), // uw
    idet * (A.vw * A.wu - A.vu * A.ww), // vu
    idet * (A.uu * A.ww - A.uw * A.wu), // vv
    idet * (A.uw * A.vu - A.uu * A.vw), // vw
    idet * (A.vu * A.wv - A.vv * A.wu), // wu
    idet * (A.uv * A.wu - A.uu * A.wv), // wv
    idet * (A.uu * A.vv - A.uv * A.vu)  // ww
  };
  return res;
}

m3 m3_transpose(const m3 A)
{
  m3 res = {
    A.uu, // uu
    A.vu, // uv
    A.wu, // uw
    A.uv, // vu
    A.vv, // vv
    A.wv, // vw
    A.uw, // wu
    A.vw, // wv
    A.ww  // ww
  };
  return res;
}

double m3_determinant(const m3 A) 
{
  return (A.uu * (A.vv * A.ww - A.vw * A.wv)
    + A.uv * (A.vw * A.wu - A.vu * A.ww)
    + A.uw * (A.vu * A.wv - A.vv * A.wu));
}

v3 m3v3_contraction(const m3 A, const v3 B)
{
  v3 res = {
    A.uu * B.u + A.uv * B.v + A.uw * B.w, // u
    A.vu * B.u + A.vv * B.v + A.vw * B.w, // v
    A.wu * B.u + A.wv * B.v + A.ww * B.w  // w
  };
  return res;
}