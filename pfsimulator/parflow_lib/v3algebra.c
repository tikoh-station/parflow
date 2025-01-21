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
#include "math.h"

v3 v3_init(double u, double v, double w)
{
  v3 res = {u, v, w};
  return res;
}

v3 v3_scale(double c, v3 A)
{
  return v3_init(c * A.u, c * A.v, c * A.w);
}

v3 v3_add(v3 A, v3 B)
{
  return v3_init(A.u + B.u, A.v + B.v, A.w + B.w);
}

v3 v3_subtract(v3 A, v3 B)
{
  return v3_init(A.u - B.u, A.v - B.v, A.w - B.w);
}

v3 v3_linear_combo(double a, v3 A, double b, v3 B)
{
  return v3_add(v3_scale(a, A), v3_scale(b, B));
}

double v3_dot(v3 A, v3 B)
{
  return (A.u * B.u + A.v * B.v + A.w * B.w);
}

v3 v3_cross(v3 A, v3 B)
{
  return v3_init(
    A.v * B.w - A.w * B.v,
    A.w * B.u - A.u * B.w,
    A.u * B.v - A.v * B.u
  );
}

double v3_norm(v3 A)
{
  return sqrt(v3_dot(A, A));
}


m3 m3_diagonal(double uu, double vv, double ww)
{
  m3 M = {uu, 0.0, 0.0, 0.0, vv, 0.0, 0.0, 0.0, ww};
  return M;
}

m3 m3_factor(double c, m3 M)
{
  m3 cM = {
    c * M.uu, // uu
    c * M.uv, // uv
    c * M.uw, // uw
    c * M.vu, // vu
    c * M.vv, // vv
    c * M.vw, // vw
    c * M.wu, // wu
    c * M.wv, // wv
    c * M.ww  // ww
  };
  return cM;
}

m3 m3_add(m3 M, m3 N)
{
  m3 M_plus_N = {
    M.uu + N.uu, M.uv + N.uv, M.uw + N.uw,
    M.vu + N.vu, M.vv + N.vv, M.vw + N.vw,
    M.wu + N.wu, M.wv + N.wv, M.ww + N.ww
  };
  return M_plus_N;
}

m3 m3_subtract(m3 M, m3 N)
{
  m3 M_minus_N = {
    M.uu - N.uu, M.uv - N.uv, M.uw - N.uw,
    M.vu - N.vu, M.vv - N.vv, M.vw - N.vw,
    M.wu - N.wu, M.wv - N.wv, M.ww - N.ww
  };
  return M_minus_N;
}

m3 m3_multiply(m3 M, m3 N)
{
  m3 M_times_N = {
    M.uu * N.uu + M.uv * N.vu + M.uw * N.wu, // uu
    M.uu * N.uv + M.uv * N.vv + M.uw * N.wv, // uv
    M.uu * N.uw + M.uv * N.vw + M.uw * N.ww, // uw
    M.vu * N.uu + M.vv * N.vu + M.vw * N.wu, // vu
    M.vu * N.uv + M.vv * N.vv + M.vw * N.wv, // vv
    M.vu * N.uw + M.vv * N.vw + M.vw * N.ww, // vw
    M.wu * N.uu + M.wv * N.vu + M.ww * N.wu, // wu
    M.wu * N.uv + M.wv * N.vv + M.ww * N.wv, // wv
    M.wu * N.uw + M.wv * N.vw + M.ww * N.ww  // ww
  };
  return M_times_N;
}

m3 m3_inverse(m3 M)
{
  double idet = 1.0 / m3_determinant(M);
  m3 iM = {
    idet * (M.vv * M.ww - M.vw * M.wv), // uu
    idet * (M.uw * M.wv - M.uv * M.ww), // uv
    idet * (M.uv * M.vw - M.uw * M.vv), // uw
    idet * (M.vw * M.wu - M.vu * M.ww), // vu
    idet * (M.uu * M.ww - M.uw * M.wu), // vv
    idet * (M.uw * M.vu - M.uu * M.vw), // vw
    idet * (M.vu * M.wv - M.vv * M.wu), // wu
    idet * (M.uv * M.wu - M.uu * M.wv), // wv
    idet * (M.uu * M.vv - M.uv * M.vu)  // ww
  };
  return iM;
}

m3 m3_transpose(m3 M)
{
  m3 res = {
    M.uu, // uu
    M.vu, // uv
    M.wu, // uw
    M.uv, // vu
    M.vv, // vv
    M.wv, // vw
    M.uw, // wu
    M.vw, // wv
    M.ww  // ww
  };
  return res;
}

double m3_determinant(m3 M) 
{
  return (M.uu * (M.vv * M.ww - M.vw * M.wv)
    + M.uv * (M.vw * M.wu - M.vu * M.ww)
    + M.uw * (M.vu * M.wv - M.vv * M.wu));
}

v3 m3v3_left_dot(v3 A, m3 M)
{
  return v3_init(
    M.uu * A.u + M.vu * A.v + M.wu * A.w, // u
    M.uv * A.u + M.vv * A.v + M.wv * A.w, // v
    M.uw * A.u + M.vw * A.v + M.ww * A.w  // w
  );
}

v3 m3v3_right_dot(m3 M, v3 A)
{
  return v3_init(
    M.uu * A.u + M.uv * A.v + M.uw * A.w, // u
    M.vu * A.u + M.vv * A.v + M.vw * A.w, // v
    M.wu * A.u + M.wv * A.v + M.ww * A.w  // w
  );
}