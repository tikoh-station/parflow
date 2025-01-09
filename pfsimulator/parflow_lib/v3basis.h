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

#ifndef _V3BASIS_HEADER
#define _V3BASIS_HEADER

#include "v3algebra.h"

typedef struct {
  v3 u, v, w;
} v3basis;

v3 v3basis_projection_onto(v3basis bas, v3 vec);

v3 v3basis_projection_from(v3basis bas, v3 vec);

v3basis v3basis_reciprocal(v3basis bas);

double v3basis_volume(v3basis bas);

m3 m3v3basis_left_dot(v3basis bas, m3 M);

m3 m3v3basis_right_dot(m3 M, v3basis bas);

#endif // _V3BASIS_HEADER
