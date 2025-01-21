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

#include "v3basis.h"

v3 v3basis_projection_onto(v3basis bas, v3 vec)
{
  return v3_init(v3_dot(bas.u, vec), v3_dot(bas.v, vec), v3_dot(bas.w, vec));
}

v3 v3basis_projection_from(v3basis bas, v3 vec)
{
  return v3_add(v3_scale(vec.u, bas.u), 
    v3_add(v3_scale(vec.v, bas.v), v3_scale(vec.w, bas.w)));
}

v3basis v3basis_reciprocal(v3basis bas)
{
  double ivol = 1.0 / v3basis_volume(bas);
  v3basis rec = {
    v3_scale(ivol, v3_cross(bas.v, bas.w)),
    v3_scale(ivol, v3_cross(bas.w, bas.u)),
    v3_scale(ivol, v3_cross(bas.u, bas.v))
  };
  return rec;
}

double v3basis_volume(v3basis bas)
{
  return v3_dot(bas.u, v3_cross(bas.v, bas.w));
}

m3 m3v3basis_left_dot(v3basis bas, m3 M)
{
  v3 ru = m3v3_left_dot(bas.u, M);
  v3 rv = m3v3_left_dot(bas.v, M);
  v3 rw = m3v3_left_dot(bas.w, M);
  m3 res = {ru.u, ru.v, ru.w, rv.u, rv.v, rv.w, rw.u, rw.v, rw.w};
  return res;
}

m3 m3v3basis_right_dot(m3 M, v3basis bas)
{
  v3 ru = m3v3_right_dot(M, bas.u);
  v3 rv = m3v3_right_dot(M, bas.v);
  v3 rw = m3v3_right_dot(M, bas.w);
  m3 res = {ru.u, rv.u, rw.u, ru.v, rv.v, rw.v, ru.w, rv.w, rw.w};
  return res;
}