/*BHEADER**********************************************************************
*
*  Copyright (c) 1995-2025, Lawrence Livermore National Security,
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

typedef struct {

  PFModule *coordinate_transform_module;

  CoordinateTransformMethods *coordinate_transform_methods;

  int transform_type;

} PublicXtra;

typedef struct {

  PFModule *coordinate_transform_module;

} InstanceXtra;

enum Transform {
  Cartesian = 0,
  TerrainFollowing = 1
};

/*--------------------------------------------------------------------------
 * CoordinateTransform
 *--------------------------------------------------------------------------*/

void CoordinateTransform(ProblemData *problem_data)
{
  PFModule *this_module = ThisPFModule;

  InitVectorAll(ProblemDataLengthUA(problem_data), 1.0);
  InitVectorAll(ProblemDataLengthUB(problem_data), 1.0);
  InitVectorAll(ProblemDataLengthVA(problem_data), 1.0);
  InitVectorAll(ProblemDataLengthVB(problem_data), 1.0);
  InitVectorAll(ProblemDataLengthWA(problem_data), 1.0);
  InitVectorAll(ProblemDataLengthWB(problem_data), 1.0);

  CoordinateTransformCurveLength(this_module, ProblemDataLengthUA(problem_data), U, OneQuarter);
  CoordinateTransformCurveLength(this_module, ProblemDataLengthUB(problem_data), U, ThreeQuarter);
  CoordinateTransformCurveLength(this_module, ProblemDataLengthVA(problem_data), V, OneQuarter);
  CoordinateTransformCurveLength(this_module, ProblemDataLengthVB(problem_data), V, ThreeQuarter);
  CoordinateTransformCurveLength(this_module, ProblemDataLengthWA(problem_data), W, OneQuarter);
  CoordinateTransformCurveLength(this_module, ProblemDataLengthWB(problem_data), W, ThreeQuarter);

  CoordinateTransformJacobian(this_module, ProblemDataJacobian(problem_data));

  CoordinateTransformZCoordinate(this_module, ProblemDataZCoordinate(problem_data));

  return;
}

// /*--------------------------------------------------------------------------
//  * CoordinateTransformEval
//  *--------------------------------------------------------------------------*/

// void CoordinateTransformEval(PFModule* this_module, int pos)
// {
//   PublicXtra *public_xtra = (PublicXtra*)PFModulePublicXtra(this_module);
//   CoordinateTransformMethods *methods = public_xtra->coordinate_transform_methods;

//   // Set ThisPFModule ahead of member function call
//   ThisPFModule = (public_xtra->coordinate_transform_module);
  
//   // Call Eval Method
//   (*(methods->eval))(pos);

//   return;
// }

/*--------------------------------------------------------------------------
 * CoordinateTransformZCoordinate
 *--------------------------------------------------------------------------*/

void CoordinateTransformZCoordinate(PFModule* this_module, Vector *z_coords)
{
  PublicXtra *public_xtra = (PublicXtra*)PFModulePublicXtra(this_module);
  InstanceXtra *instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);
  CoordinateTransformMethods *methods = public_xtra->coordinate_transform_methods;

  // Set ThisPFModule ahead of member function call
  ThisPFModule = (instance_xtra->coordinate_transform_module);
  
  // Call ZCoordinate Method
  (*(methods->z_coordinate))(z_coords);

  return;
}

/*--------------------------------------------------------------------------
 * CoordinateTransformCurveLength
 *--------------------------------------------------------------------------*/

void CoordinateTransformCurveLength(PFModule* this_module, Vector *length,
                                    Coordinate coordinate,
                                    PositionInCell position)
{
  PublicXtra *public_xtra = (PublicXtra*)PFModulePublicXtra(this_module);
  InstanceXtra *instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);
  CoordinateTransformMethods *methods = public_xtra->coordinate_transform_methods;

  // Set ThisPFModule ahead of member function call
  ThisPFModule = (instance_xtra->coordinate_transform_module);
  
  // Call CurveLength Method
  (*(methods->curve_length))(length, coordinate, position);

  return;
}

/*--------------------------------------------------------------------------
 * CoordinateTransformJacobian
 *--------------------------------------------------------------------------*/

void CoordinateTransformJacobian(PFModule* this_module, Vector *jacobian)
{
  PublicXtra *public_xtra = (PublicXtra*)PFModulePublicXtra(this_module);
  InstanceXtra *instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);
  CoordinateTransformMethods *methods = public_xtra->coordinate_transform_methods;

  // Set ThisPFModule ahead of member function call
  ThisPFModule = (instance_xtra->coordinate_transform_module);
  
  // Call Jacobian Method
  (*(methods->jacobian))(jacobian);

  return;
}

/*--------------------------------------------------------------------------
 * CoordinateTransformTranslationFactors
 *--------------------------------------------------------------------------*/

void CoordinateTransformTranslationFactors(PFModule *this_module,
                                           Vector *T_uu,
                                           Vector *T_uv,
                                           Vector *T_uw,
                                           Coordinate coordinate,
                                           PositionInCell position)
{
  PublicXtra *public_xtra = (PublicXtra*)PFModulePublicXtra(this_module);
  InstanceXtra *instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);
  CoordinateTransformMethods *methods = public_xtra->coordinate_transform_methods;

  // Set ThisPFModule ahead of member function call
  ThisPFModule = (instance_xtra->coordinate_transform_module);
  
  // Call CurveLength Method
  (*(methods->translation_factors))(T_uu, T_uv, T_uw, coordinate, position);

  return;
}

/*--------------------------------------------------------------------------
 * CoordinateTransformMetricContravariant
 *--------------------------------------------------------------------------*/

void CoordinateTransformMetricContravariant(PFModule *this_module,
                                            Vector *g_uu, Vector *g_uv,
                                            Vector *g_uw, Vector *g_vv,
                                            Vector *g_vw, Vector *g_ww)
{
  PublicXtra *public_xtra = (PublicXtra*)PFModulePublicXtra(this_module);
  InstanceXtra *instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);
  CoordinateTransformMethods *methods = public_xtra->coordinate_transform_methods;

  // Set ThisPFModule ahead of member function call
  ThisPFModule = (instance_xtra->coordinate_transform_module);
  
  // Call MetricContravariant Method
  (*(methods->metric_contravariant))(g_uu, g_uv, g_uw, g_vv, g_vw, g_ww);

  return;
}

/*--------------------------------------------------------------------------
 * CoordinateTransformInitInstanceXtra
 *--------------------------------------------------------------------------*/

PFModule  *CoordinateTransformInitInstanceXtra(ProblemData *problem_data)
{
  PFModule      *this_module = ThisPFModule;
  InstanceXtra  *instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);
  PublicXtra    *public_xtra = (PublicXtra*)PFModulePublicXtra(this_module);

  if (instance_xtra == NULL)
  {
    instance_xtra = ctalloc(InstanceXtra, 1);

    switch(public_xtra->transform_type)
    {
      case Cartesian:
      {
        instance_xtra->coordinate_transform_module = 
          PFModuleNewInstanceType(CartesianTransformInitInstanceXtraInvoke,
                              public_xtra->coordinate_transform_module,
                              (problem_data));
        break;
      }

      case TerrainFollowing:
      {
        instance_xtra->coordinate_transform_module = 
          PFModuleNewInstanceType(TerrainFollowingTransformInitInstanceXtraInvoke,
                              public_xtra->coordinate_transform_module,
                              (problem_data));
        break;
      }

      default:
      {
        InputError("Invalid Transform Type.", "", "");
      }
    }
  }
  else
  {
    switch(public_xtra->transform_type)
    {
      case Cartesian:
      {
        PFModuleReNewInstanceType(CartesianTransformInitInstanceXtraInvoke,
                              instance_xtra->coordinate_transform_module, ());
        break;
      }

      case TerrainFollowing:
      {
        PFModuleReNewInstanceType(TerrainFollowingTransformInitInstanceXtraInvoke,
                              instance_xtra->coordinate_transform_module,
                              (problem_data));
        break;
      }

      default:
      {
        InputError("Invalid Transform Type.", "", "");
      }
    }
  }

  PFModuleInstanceXtra(this_module) = instance_xtra;
  return this_module;
}


/*--------------------------------------------------------------------------
 * CoordinateTransformFreeInstanceXtra
 *--------------------------------------------------------------------------*/

void  CoordinateTransformFreeInstanceXtra()
{
  PFModule      *this_module = ThisPFModule;
  InstanceXtra  *instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);

  if (instance_xtra)
  {
    PFModuleFreeInstance(instance_xtra->coordinate_transform_module);

    tfree(instance_xtra);
  }
}

/*--------------------------------------------------------------------------
 * CoordinateTransformNewPublicXtra
 *--------------------------------------------------------------------------*/

PFModule  *CoordinateTransformNewPublicXtra()
{
  PFModule      *this_module = ThisPFModule;
  PublicXtra    *public_xtra = ctalloc(PublicXtra, 1);

  (public_xtra->coordinate_transform_methods) = ctalloc(CoordinateTransformMethods, 1);

  NameArray transform_name_array = NA_NewNameArray("Cartesian TerrainFollowing");

  char key[IDB_MAX_KEY_LEN];
  sprintf(key, "ComputationalGrid.Transform.Type");
  char *transform_name = GetStringDefault(key, "Cartesian");
  int transform_type = NA_NameToIndexExitOnError(transform_name_array, transform_name, key);

  public_xtra->transform_type = transform_type;
  switch (transform_type)
  {
    case Cartesian:
    {
      (public_xtra->coordinate_transform_module) = PFModuleNewModuleType(CartesianTransformNewPublicXtraInvoke, CartesianTransform, (public_xtra->coordinate_transform_methods));
      break;
    }

    case TerrainFollowing:
    {
      (public_xtra->coordinate_transform_module) = PFModuleNewModuleType(TerrainFollowingTransformNewPublicXtraInvoke, TerrainFollowingTransform, (public_xtra->coordinate_transform_methods));
      break;
    }

    default:
    {
      InputError("Invalid value <%s> for key <%s>", transform_name, key);
    }
  }

  NA_FreeNameArray(transform_name_array);

  PFModulePublicXtra(this_module) = public_xtra;
  return this_module;
}

/*-------------------------------------------------------------------------
 * CoordinateTransformFreePublicXtra
 *-------------------------------------------------------------------------*/

void  CoordinateTransformFreePublicXtra()
{
  PFModule    *this_module = ThisPFModule;
  PublicXtra  *public_xtra = (PublicXtra*)PFModulePublicXtra(this_module);

  if (public_xtra)
  {
    tfree(public_xtra->coordinate_transform_methods);
    PFModuleFreeModule(public_xtra->coordinate_transform_module);
    tfree(public_xtra);
  }
}

/*--------------------------------------------------------------------------
 * CoordinateTransformSizeOfTempData
 *--------------------------------------------------------------------------*/

int  CoordinateTransformSizeOfTempData()
{
  return 0;
}