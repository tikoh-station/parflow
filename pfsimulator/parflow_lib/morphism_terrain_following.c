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

#include "morphism_terrain_following.h"

/*--------------------------------------------------------------------------
 * Structures
 *--------------------------------------------------------------------------*/

typedef struct {
  Morphism *my_morphism;
} PublicXtra;

typedef void InstanceXtra;

/*-------------------------------------------------------------------------
 * TerrainFollowing
 *-------------------------------------------------------------------------*/
void TerrainFollowing(Morphism *my_morphism)
{
  PFModule *this_module = ThisPFModule;
  PublicXtra *public_xtra = (PublicXtra*)PFModulePublicXtra(this_module);
  my_morphism = public_xtra->my_morphism;
  return;
}

/*--------------------------------------------------------------------------
 * TerrainFollowingInitInstanceXtra
 *--------------------------------------------------------------------------*/
PFModule* TerrainFollowingInitInstanceXtra()
{
  PFModule      *this_module = ThisPFModule;
  InstanceXtra  *instance_xtra = ctalloc(InstanceXtra, 1);

  PFModuleInstanceXtra(this_module) = instance_xtra;
  return this_module;
}

/*--------------------------------------------------------------------------
 * TerrainFollowingFreeInstanceXtra
 *--------------------------------------------------------------------------*/
void TerrainFollowingFreeInstanceXtra()
{
  PFModule      *this_module = ThisPFModule;
  InstanceXtra  *instance_xtra = 
    (InstanceXtra*)PFModuleInstanceXtra(this_module);

  if (instance_xtra)
  {
    tfree(instance_xtra);
  }
}

/*--------------------------------------------------------------------------
 * TerrainFollowingNewPublicXtra
 *--------------------------------------------------------------------------*/
PFModule* TerrainFollowingNewPublicXtra()
{
  PFModule *this_module = ThisPFModule;
  PublicXtra *public_xtra = ctalloc(PublicXtra, 1);
  Morphism *my_morphism = ctalloc(Morphism, 1);

  my_morphism->morphism_core = ThisPFModule;
  my_morphism->eval = TerrainFollowingEval;
  my_morphism->inverse = TerrainFollowingInverse;
  my_morphism->del = TerrainFollowingDel;
  my_morphism->del_inverse = TerrainFollowingDelInverse;
  my_morphism->jacobian = TerrainFollowingJacobian;

  public_xtra->my_morphism = my_morphism;

  PFModulePublicXtra(this_module) = public_xtra;
  return this_module;
}

/*-------------------------------------------------------------------------
 * TerrainFollowingFreePublicXtra
 *-------------------------------------------------------------------------*/
void  TerrainFollowingFreePublicXtra()
{
  PFModule *this_module = ThisPFModule;
  PublicXtra *public_xtra = (PublicXtra*)PFModulePublicXtra(this_module);

  if (public_xtra)
  {
    tfree(public_xtra->my_morphism);
    tfree(public_xtra);
  }
}

/*--------------------------------------------------------------------------
 * TerrainFollowingSizeOfTempData
 *--------------------------------------------------------------------------*/
int  TerrainFollowingSizeOfTempData()
{
  return 0;
}



/*--------------------------------------------------------------------------
 * Morphism Functions
 *--------------------------------------------------------------------------*/
v3 TerrainFollowingEval(const v3 zeta)
{
  double elevation = 0;
  v3 x = {zeta.u, zeta.v, zeta.w + elevation};
  return x;
}

v3 TerrainFollowingInverse(const v3 x)
{
  double elevation = 0;
  v3 zeta = {x.u, x.v, x.w - elevation};
  return zeta;
}
