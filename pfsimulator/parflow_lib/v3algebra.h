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
  double uu, uv, uw, vu, vv, vw, wu, wv, ww;
} m3;

v3 v3_init(double u, double v, double w);
v3 v3_scale(double c, v3 A);
v3 v3_add(v3 A, v3 B);
v3 v3_subtract(v3 A, v3 B);
v3 v3_linear_combo(double a, v3 A, double b, v3 B);
double v3_dot(v3 A, v3 B);
v3 v3_cross(v3 A, v3 B);
double v3_norm(v3 A);

m3 m3_diagonal(double uu, double vv, double ww);
m3 m3_factor(double c, m3 M);
m3 m3_add(m3 M, m3 N);
m3 m3_subtract(m3 M, m3 N);
m3 m3_multiply(m3 M, m3 N);
m3 m3_inverse(m3 M);
m3 m3_transpose(m3 M);
double m3_determinant(m3 M);

v3 m3v3_left_dot(v3 A, m3 M);
v3 m3v3_right_dot(m3 M, v3 A);

#endif // _V3ALGEBRA_HEADER