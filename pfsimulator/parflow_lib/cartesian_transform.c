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

#include "parflow.h"

typedef void PublicXtra;

typedef void InstanceXtra;

/*--------------------------------------------------------------------------
 * CartesianTransform
 *--------------------------------------------------------------------------*/

void CartesianTransform()
{
  return;
}

/*--------------------------------------------------------------------------
 * CartesianTransformEval
 *--------------------------------------------------------------------------*/

void CartesianTransformEval(int pos)
{
  amps_Printf("[CartesianTransform::Eval] pos: %d\n", pos);
  return;
}

/*--------------------------------------------------------------------------
 * CartesianTransformZCoordinate
 *--------------------------------------------------------------------------*/

void CartesianTransformZCoordinate(Vector *z_coords)
{
  PFModule     *this_module = ThisPFModule;
  InstanceXtra *instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);

  Grid *grid = VectorGrid(z_coords);

  int is = 0;
  ForSubgridI(is, GridSubgrids(grid))
  {
    Subgrid *subgrid = GridSubgrid(grid, is);
    Subvector *z_sub = VectorSubvector(z_coords, is);
    double *z_dat = SubvectorData(z_sub);

    /* RDF: assumes resolutions are the same in all 3 directions */
    int r = SubgridRX(subgrid);

    int iu = SubgridIX(subgrid);
    int iv = SubgridIY(subgrid);
    int iw = SubgridIZ(subgrid);

    int nu = SubgridNX(subgrid);
    int nv = SubgridNY(subgrid);
    int nw = SubgridNZ(subgrid);

    int nu_v = SubvectorNX(z_sub);
    int nv_v = SubvectorNY(z_sub);
    int nw_v = SubvectorNZ(z_sub);

    int i = 0, j = 0, k = 0;
    int idx = SubvectorEltIndex(z_sub, i, j, k);
    BoxLoopI1(i, j, k, iu, iv, iw, nu, nv, nw,
              idx, nu_v, nv_v, nw_v, 1, 1, 1,
    {
      z_dat[idx] = RealSpaceZ(k, r);
    });
  }

  VectorUpdateCommHandle *handle = InitVectorUpdate(z_coords, VectorUpdateAll);
  FinalizeVectorUpdate(handle);

  return;
}

/*--------------------------------------------------------------------------
 * CartesianTransformCurveLength
 *--------------------------------------------------------------------------*/

void CartesianTransformCurveLength(Vector *length, Coordinate coordinate,
                                   PositionInCell position)
{
  InitVectorAll(length, 1.0);
  return;
}

/*--------------------------------------------------------------------------
 * CartesianTransformJacobian
 *--------------------------------------------------------------------------*/

void CartesianTransformJacobian(Vector *jacobian)
{
  InitVectorAll(jacobian, 1.0);
  return;
}

/*--------------------------------------------------------------------------
 * CartesianTransformJacobian
 *--------------------------------------------------------------------------*/

void CartesianTransformTranslationFactors(Vector *T_uu,
                                          Vector *T_uv,
                                          Vector *T_uw,
                                          Coordinate coordinate,
                                          PositionInCell position)
{
  InitVectorAll(T_uu, 1.0);
  InitVectorAll(T_uv, 0.0);
  InitVectorAll(T_uw, 0.0);
  return;
}

/*--------------------------------------------------------------------------
 * CartesianTransformMetricContravariant
 *--------------------------------------------------------------------------*/

void CartesianTransformMetricContravariant(Vector *g_uu, Vector *g_uv,
                                           Vector *g_uw, Vector *g_vv,
                                           Vector *g_vw, Vector *g_ww)
{
  InitVectorAll(g_uu, 1.0);
  InitVectorAll(g_uv, 0.0);
  InitVectorAll(g_uw, 0.0);
  InitVectorAll(g_vv, 1.0);
  InitVectorAll(g_vw, 0.0);
  InitVectorAll(g_ww, 1.0);

  return;
}

/*--------------------------------------------------------------------------
 * CartesianTransformInitInstanceXtra
 *--------------------------------------------------------------------------*/

PFModule  *CartesianTransformInitInstanceXtra()
{
  PFModule      *this_module = ThisPFModule;
  InstanceXtra  *instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);

  if (instance_xtra == NULL)
  {
    instance_xtra = ctalloc(InstanceXtra, 1);
  }

  PFModuleInstanceXtra(this_module) = instance_xtra;
  return this_module;
}


/*--------------------------------------------------------------------------
 * CartesianTransformFreeInstanceXtra
 *--------------------------------------------------------------------------*/

void  CartesianTransformFreeInstanceXtra()
{
  PFModule      *this_module = ThisPFModule;
  InstanceXtra  *instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);

  if (instance_xtra)
  {
    tfree(instance_xtra);
  }
}

/*--------------------------------------------------------------------------
 * CartesianTransformNewPublicXtra
 *--------------------------------------------------------------------------*/

PFModule  *CartesianTransformNewPublicXtra(CoordinateTransformMethods *coordinate_transform_methods)
{
  PFModule      *this_module = ThisPFModule;
  PublicXtra    *public_xtra = NULL;

  // coordinate_transform_methods->eval = CartesianTransformEval;
  coordinate_transform_methods->z_coordinate = CartesianTransformZCoordinate;
  coordinate_transform_methods->curve_length = CartesianTransformCurveLength;
  coordinate_transform_methods->jacobian = CartesianTransformJacobian;
  coordinate_transform_methods->metric_contravariant = CartesianTransformMetricContravariant;
  coordinate_transform_methods->translation_factors = CartesianTransformTranslationFactors;
  
  PFModulePublicXtra(this_module) = public_xtra;
  return this_module;
}

/*-------------------------------------------------------------------------
 * CartesianTransformFreePublicXtra
 *-------------------------------------------------------------------------*/

void  CartesianTransformFreePublicXtra()
{
  PFModule    *this_module = ThisPFModule;
  PublicXtra  *public_xtra = (PublicXtra*)PFModulePublicXtra(this_module);

  if (public_xtra)
  {
    tfree(public_xtra);
  }
}

/*--------------------------------------------------------------------------
 * CartesianTransformSizeOfTempData
 *--------------------------------------------------------------------------*/

int  CartesianTransformSizeOfTempData()
{
  return 0;
}