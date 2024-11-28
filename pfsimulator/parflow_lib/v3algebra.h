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

#ifndef _V3ALGEBRA_HEADER
#define _V3ALGEBRA_HEADER

typedef struct {
  double u, v, w;
} v3;

typedef struct {
  double uu, uv, uw, vv, vw, ww;
} sm3;

typedef struct {
  double uu, uv, uw, vu, vv, vw, wu, wv, ww;
} m3;

v3 v3_scale(const double c, const v3 A);
v3 v3_add(const v3 A, const v3 B);
v3 v3_subtract(const v3 A, const v3 B);
double v3_dot(const v3 A, const v3 B);
v3 v3_cross(const v3 A, const v3 B);

m3 m3_add(const m3 A, const m3 B);
m3 m3_subtract(const m3 A, const m3 B);
m3 m3_multiply(const m3 A, const m3 B);
m3 m3_inverse(const m3 A);
m3 m3_transpose(const m3 A);
double m3_determinant(const m3 A);

v3 m3v3_contraction(const m3 A, const v3 B);

#endif // _V3ALGEBRA_HEADER