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

#ifndef _MORPHISM_TERRAIN_FOLLOWING_HEADER
#define _MORPHISM_TERRAIN_FOLLOWING_HEADER

#include "morphism.h"

typedef void (*TerrainFollowingInvoke)(Morphism **my_morphism);
void TerrainFollowing(Morphism **my_morphism);
PFModule* TerrainFollowingInitInstanceXtra();
void TerrainFollowingFreeInstanceXtra();
PFModule* TerrainFollowingNewPublicXtra();
void  TerrainFollowingFreePublicXtra();
int  TerrainFollowingSizeOfTempData();

v3 TerrainFollowingEval(PFModule *morphism_core, v3 zeta);
v3 TerrainFollowingInverse(PFModule *morphism_core, v3 x);
v3basis TerrainFollowingDel(PFModule *morphism_core, v3 zeta);
v3basis TerrainFollowingDelInverse(PFModule *morphism_core, v3 zeta);
double TerrainFollowingJacobian(PFModule *morphism_core, v3 zeta);

v3 TerrainFollowingToContravariant(v3 vec, v3basis basis_contravariant);
v3 TerrainFollowingToCovariant(v3 vec, v3basis basis_covariant);
v3 TerrainFollowingFromContravariant(v3 vec, v3basis basis_contravariant);
v3 TerrainFollowingFromCovariant(v3 vec, v3basis basis_covariant);

#endif // _MORPHISM_TERRAIN_FOLLOWING_HEADER