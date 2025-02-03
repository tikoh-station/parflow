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
  // Vector *slope_x;
  // Vector *slope_y;
  double slope_x;
  double slope_y;
} PublicXtra;

typedef void InstanceXtra;

/*-------------------------------------------------------------------------
 * TerrainFollowing
 *-------------------------------------------------------------------------*/
void TerrainFollowing(Morphism **my_morphism)
{
  PFModule *this_module = ThisPFModule;
  PublicXtra *public_xtra = (PublicXtra*)PFModulePublicXtra(this_module);
  *my_morphism = public_xtra->my_morphism;
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

  my_morphism->to_contravariant = MorphismToContravariant;
  my_morphism->to_covariant = MorphismToCovariant;
  my_morphism->from_contravariant = MorphismFromContravariant;
  my_morphism->from_covariant = MorphismFromCovariant;

  public_xtra->my_morphism = my_morphism;

  // read constant slopes from script
  char keyname[IDB_MAX_KEY_LEN];
  NameArray na_types = NA_NewNameArray("Constant");

  strcpy(keyname, "Solver.TerrainFollowingGrid.SlopesX.Type");
  char *xswitch_name = GetString(keyname);
  int xtype = NA_NameToIndexExitOnError(na_types, xswitch_name, keyname);
  switch(xtype)
  {
    case 0:
    {
      strcpy(keyname, "Solver.TerrainFollowingGrid.SlopesX.Value");
      public_xtra->slope_x = GetDouble(keyname);
      break;
    }
    default:
    {
      public_xtra->slope_x = 0.0;
      break;
    }
  }

  strcpy(keyname, "Solver.TerrainFollowingGrid.SlopesY.Type");
  char *yswitch_name = GetString(keyname);
  int ytype = NA_NameToIndexExitOnError(na_types, yswitch_name, keyname);
  switch(ytype)
  {
    case 0:
    {
      strcpy(keyname, "Solver.TerrainFollowingGrid.SlopesY.Value");
      public_xtra->slope_y = GetDouble(keyname);
      break;
    }
    default:
    {
      public_xtra->slope_y = 0.0;
      break;
    }
  }
  
  // public_xtra->slope_x = ProblemDataTSlopeX(problem_data);
  // public_xtra->slope_y = ProblemDataTSlopeY(problem_data);

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
v3 TerrainFollowingEval(PFModule *morphism_core, const v3 zeta)
{
  // for now, this methods aren't really required,
  // which is why this isn't yet implemented.
  double elevation = 0; 
  v3 x = {zeta.u, zeta.v, zeta.w + elevation};
  return x;
}

v3 TerrainFollowingInverse(PFModule *morphism_core, const v3 x)
{
  // for now, this methods aren't really required,
  // which is why this isn't yet implemented.
  double elevation = 0; 
  v3 zeta = {x.u, x.v, x.w - elevation};
  return zeta;
}

v3basis TerrainFollowingDel(PFModule *morphism_core, v3 zeta)
{
  PublicXtra *public_xtra = PFModulePublicXtra(morphism_core);
  // int i = IndexSpaceX(zeta.u);
  // int j = IndexSpaceY(zeta.v);
  // double lx = zeta.u-ceil(zeta.u);
  // double ly = zeta.v-ceil(zeta.v);
  double sx = public_xtra->slope_x;
  double sy = public_xtra->slope_y;

  v3basis del = {{1, 0, sx}, {0, 1, sy}, {0, 0, 1}};
  return del;
}

v3basis TerrainFollowingDelInverse(PFModule *morphism_core, v3 zeta)
{
  PublicXtra *public_xtra = PFModulePublicXtra(morphism_core);
  double sx = public_xtra->slope_x;
  double sy = public_xtra->slope_y;

  v3basis del_inverse = {{1, 0, 0}, {0, 1, 0}, {-sx, -sy, 1}};
  return del_inverse;
}

double TerrainFollowingJacobian(PFModule *morphism_core, v3 zeta)
{
  return 1.0;
}



/*--------------------------------------------------------------------------
 * Morphism Functions
 *--------------------------------------------------------------------------*/
// v3 TerrainFollowingToContravariant(v3 vec, v3basis basis_contravariant)
// {
//   return vec;
// }

// v3 TerrainFollowingToCovariant(v3 vec, v3basis basis_covariant)
// {
//   return vec;
// }

// v3 TerrainFollowingFromContravariant(v3 vec, v3basis basis_contravariant)
// {
//   return vec;
// }

// v3 TerrainFollowingFromCovariant(v3 vec, v3basis basis_covariant)
// {
//   return vec;
// }
