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
/*****************************************************************************
*
*  This module computes the contributions for the spatial discretization of the
*  kinematic equation for the overland flow boundary condition:KE,KW,KN,KS.
*
*  It also computes the derivatives of these terms for inclusion in the Jacobian.
*
* Could add a switch statement to handle the diffusion wave also.
* -DOK
*****************************************************************************/

#include "morphism_variable_dz.h"
#include "math.h"
// #include "globals.h"

/*--------------------------------------------------------------------------
 * Structures
 *--------------------------------------------------------------------------*/

typedef struct {
  Morphism *my_morphism;
} PublicXtra;

typedef void InstanceXtra;

/*-------------------------------------------------------------------------
 * MorphismVariableDZ
 *-------------------------------------------------------------------------*/
void MorphismVariableDZ(Morphism **my_morphism)
{
  PFModule *this_module = ThisPFModule;
  PublicXtra *public_xtra = (PublicXtra*)PFModulePublicXtra(this_module);
  *my_morphism = public_xtra->my_morphism;
  return;
}

/*--------------------------------------------------------------------------
 * MorphismVariableDZInitInstanceXtra
 *--------------------------------------------------------------------------*/
PFModule* MorphismVariableDZInitInstanceXtra()
{
  PFModule      *this_module = ThisPFModule;
  InstanceXtra  *instance_xtra = ctalloc(InstanceXtra, 1);

  PFModuleInstanceXtra(this_module) = instance_xtra;
  return this_module;
}

/*--------------------------------------------------------------------------
 * MorphismVariableDZFreeInstanceXtra
 *--------------------------------------------------------------------------*/
void MorphismVariableDZFreeInstanceXtra()
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
 * MorphismVariableDZNewPublicXtra
 *--------------------------------------------------------------------------*/
PFModule* MorphismVariableDZNewPublicXtra()
{
  PFModule *this_module = ThisPFModule;
  PublicXtra *public_xtra = ctalloc(PublicXtra, 1);
  Morphism *my_morphism = ctalloc(Morphism, 1);

  my_morphism->morphism_core = ThisPFModule;
  my_morphism->eval = MorphismVariableDZEval;
  my_morphism->inverse = MorphismVariableDZInverse;
  my_morphism->del = MorphismVariableDZDel;
  my_morphism->del_inverse = MorphismVariableDZDelInverse;
  my_morphism->jacobian = MorphismVariableDZJacobian;

  my_morphism->to_contravariant = MorphismToContravariant;
  my_morphism->to_covariant = MorphismToCovariant;
  my_morphism->from_contravariant = MorphismFromContravariant;
  my_morphism->from_covariant = MorphismFromCovariant;

  public_xtra->my_morphism = my_morphism;

  PFModulePublicXtra(this_module) = public_xtra;
  return this_module;
}

/*-------------------------------------------------------------------------
 * MorphismVariableDZFreePublicXtra
 *-------------------------------------------------------------------------*/
void  MorphismVariableDZFreePublicXtra()
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
 * MorphismVariableDZSizeOfTempData
 *--------------------------------------------------------------------------*/
int  MorphismVariableDZSizeOfTempData()
{
  return 0;
}



/*--------------------------------------------------------------------------
 * Morphism Functions
 *--------------------------------------------------------------------------*/
v3 MorphismVariableDZEval(PFModule *morphism_core, v3 zeta)
{
  return v3_init(zeta.u, zeta.v, ZFunction(zeta.w));
}

v3 MorphismVariableDZInverse(PFModule *morphism_core, v3 x)
{
  return v3_init(x.u, x.v, ZFunctionInverse(x.w));
}

v3basis MorphismVariableDZDel(PFModule *morphism_core, v3 zeta)
{
  v3basis del = {{1, 0, 0}, {0, 1, 0}, {0, 0, ZFunctionDel(zeta.w)}};
  return del;
}

v3basis MorphismVariableDZDelInverse(PFModule *morphism_core, v3 zeta)
{
  v3basis del_inverse = {{1, 0, 0}, {0, 1, 0}, 
    {0, 0, 1.0 / ZFunctionDel(zeta.w)}};
  return del_inverse;
}

double MorphismVariableDZJacobian(PFModule *morphism_core, v3 zeta)
{
  return ZFunctionDel(zeta.w);
}


// defined for a zlabel varying between 0 and 1, and gives z between 0 and 2
double ZFunction(double zlabel)
{
  double c = 4.0 * zlabel - 2.0;
  double cubed = c * c * c;
  return (0.25 * cubed + 2.0 * zlabel + 2.0) / 3.0;
}

double ZFunctionInverse(double z)
{
  double p = 0.125;
  double q = 0.1875 * (1 - z);
  double sqrt_Delta = sqrt(0.25 * q * q + p * p * p / 27.0);
  double half_q = 0.5 * q;
  double zlabel = 0.5 + cbrt(-half_q + sqrt_Delta) + cbrt(-half_q - sqrt_Delta);
  return zlabel;
}

double ZFunctionDel(double zlabel)
{
  double c = 4.0 * zlabel - 2.0;
  return (c * c + 2.0 / 3.0);
}