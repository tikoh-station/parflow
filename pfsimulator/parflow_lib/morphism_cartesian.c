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

#include "morphism_cartesian.h"
// #include "globals.h"

/*--------------------------------------------------------------------------
 * Structures
 *--------------------------------------------------------------------------*/

typedef struct {
  Morphism *my_morphism;
} PublicXtra;

typedef void InstanceXtra;

/*-------------------------------------------------------------------------
 * Cartesian
 *-------------------------------------------------------------------------*/
void Cartesian(Morphism **my_morphism)
{
  PFModule *this_module = ThisPFModule;
  PublicXtra *public_xtra = (PublicXtra*)PFModulePublicXtra(this_module);
  *my_morphism = public_xtra->my_morphism;
  return;
}

/*--------------------------------------------------------------------------
 * CartesianInitInstanceXtra
 *--------------------------------------------------------------------------*/
PFModule* CartesianInitInstanceXtra()
{
  PFModule      *this_module = ThisPFModule;
  InstanceXtra  *instance_xtra = ctalloc(InstanceXtra, 1);

  PFModuleInstanceXtra(this_module) = instance_xtra;
  return this_module;
}

/*--------------------------------------------------------------------------
 * CartesianFreeInstanceXtra
 *--------------------------------------------------------------------------*/
void CartesianFreeInstanceXtra()
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
 * CartesianNewPublicXtra
 *--------------------------------------------------------------------------*/
PFModule* CartesianNewPublicXtra()
{
  PFModule *this_module = ThisPFModule;
  PublicXtra *public_xtra = ctalloc(PublicXtra, 1);
  Morphism *my_morphism = ctalloc(Morphism, 1);

  my_morphism->morphism_core = ThisPFModule;
  my_morphism->eval = CartesianEval;
  my_morphism->inverse = CartesianInverse;
  my_morphism->del = CartesianDel;
  my_morphism->del_inverse = CartesianDelInverse;
  my_morphism->jacobian = CartesianJacobian;

  my_morphism->to_contravariant = CartesianToContravariant;
  my_morphism->to_covariant = CartesianToCovariant;
  my_morphism->from_contravariant = CartesianFromContravariant;
  my_morphism->from_covariant = CartesianFromCovariant;

  public_xtra->my_morphism = my_morphism;

  PFModulePublicXtra(this_module) = public_xtra;
  return this_module;
}

/*-------------------------------------------------------------------------
 * CartesianFreePublicXtra
 *-------------------------------------------------------------------------*/
void  CartesianFreePublicXtra()
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
 * CartesianSizeOfTempData
 *--------------------------------------------------------------------------*/
int  CartesianSizeOfTempData()
{
  return 0;
}



/*--------------------------------------------------------------------------
 * Morphism Functions
 *--------------------------------------------------------------------------*/
v3 CartesianEval(PFModule *morphism_core, v3 zeta)
{
  return zeta;
}

v3 CartesianInverse(PFModule *morphism_core, v3 x)
{
  return x;
}

v3basis CartesianDel(PFModule *morphism_core, v3 zeta)
{
  v3basis del = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  return del;
}

v3basis CartesianDelInverse(PFModule *morphism_core, v3 zeta)
{
  v3basis del_inverse = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  return del_inverse;
}

double CartesianJacobian(PFModule *morphism_core, v3 zeta)
{
  return 1.0;
}



/*--------------------------------------------------------------------------
 * Morphism Functions
 *--------------------------------------------------------------------------*/
v3 CartesianToContravariant(v3 vec, v3basis basis_contravariant)
{
  return vec;
}

v3 CartesianToCovariant(v3 vec, v3basis basis_covariant)
{
  return vec;
}

v3 CartesianFromContravariant(v3 vec, v3basis basis_contravariant)
{
  return vec;
}

v3 CartesianFromCovariant(v3 vec, v3basis basis_covariant)
{
  return vec;
}
