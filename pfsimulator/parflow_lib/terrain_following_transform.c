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

#define PredefinedFunction 0
#define PFBFile 1

typedef struct {

  int input_type;
  char *filename;

} PublicXtra;

typedef struct {

  ProblemData *problem_data;
  
  Vector *elevation;
  Vector *delevation_du;
  Vector *delevation_dv;

} InstanceXtra;

/*--------------------------------------------------------------------------
 * TerrainFollowingTransform
 *--------------------------------------------------------------------------*/

void TerrainFollowingTransform()
{
  return;
}

/*--------------------------------------------------------------------------
 * TerrainFollowingTransformEval
 *--------------------------------------------------------------------------*/

void TerrainFollowingTransformEval(int pos)
{
  amps_Printf("[TerrainFollowingTransform::Eval] pos: %d\n", pos);
  return;
}

/*--------------------------------------------------------------------------
 * TerrainFollowingTransformZCoordinate
 *--------------------------------------------------------------------------*/

void TerrainFollowingTransformZCoordinate(Vector *z_coords)
{
  PFModule     *this_module = ThisPFModule;
  InstanceXtra *instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);

  Vector *elevation = instance_xtra->elevation;

  Grid *grid = VectorGrid(z_coords);
  Grid *grid2d = VectorGrid(elevation);

  int is = 0;
  ForSubgridI(is, GridSubgrids(grid))
  {
    Subgrid *subgrid = GridSubgrid(grid, is);

    Subvector *z_sub = VectorSubvector(z_coords, is);
    double *z_dat = SubvectorData(z_sub);
    Subvector *elevation_sub = VectorSubvector(elevation, is);
    double *elevation_dat = SubvectorData(elevation_sub);

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
      int idx2d = SubvectorEltIndex(elevation_sub, i, j, 0);

      z_dat[idx] = RealSpaceZ(k, r) + elevation_dat[idx2d];
    });
  }

  VectorUpdateCommHandle *handle = InitVectorUpdate(z_coords, VectorUpdateAll);
  FinalizeVectorUpdate(handle);

  return;
}

/*--------------------------------------------------------------------------
 * TerrainFollowingTransformCurveLength
 *--------------------------------------------------------------------------*/

void TerrainFollowingTransformCurveLength(Vector *length, Coordinate coordinate,
                                          PositionInCell position)
{
  if(coordinate == W) {
    InitVectorAll(length, 1.0);
    return;
  }

  PFModule     *this_module = ThisPFModule;
  InstanceXtra *instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);

  Vector *del_elevation = NULL;
  if (coordinate == U)
    del_elevation = instance_xtra->delevation_du;
  else if (coordinate == V)
    del_elevation = instance_xtra->delevation_dv;

  Grid *grid = VectorGrid(length);

  int is = 0;
  ForSubgridI(is, GridSubgrids(grid))
  {
    Subgrid *subgrid = GridSubgrid(grid, is);

    Subvector *length_sub = VectorSubvector(length, is);
    double *length_dat = SubvectorData(length_sub);
    Subvector *del_elevation_sub = VectorSubvector(del_elevation, is);
    double *del_elevation_dat = SubvectorData(del_elevation_sub);

    int iu = SubgridIX(subgrid);
    int iv = SubgridIY(subgrid);
    int iw = SubgridIZ(subgrid);

    int nu = SubgridNX(subgrid);
    int nv = SubgridNY(subgrid);
    int nw = SubgridNZ(subgrid);

    int nu_v = SubvectorNX(length_sub);
    int nv_v = SubvectorNY(length_sub);
    int nw_v = SubvectorNZ(length_sub);

    int stride_u = 1;
    int stride_v = nu_v;

    int stride = 0;
    switch (coordinate)
    {
      case U:
        stride = stride_u;
        break;
      case V:
        stride = stride_v;
        break;
    }

    int i = -1, j = -1, k = -1;
    int idx = SubvectorEltIndex(length_sub, i, j, k);
    BoxLoopI1(i, j, k, (iu-1), (iv-1), (iw-1), (nu+2), (nv+2), (nw+2),
              idx, nu_v, nv_v, nw_v, 1, 1, 1,
    {
      int idx2d = SubvectorEltIndex(del_elevation_sub, i, j, 0);
      
      double del_elev = 0.0;
      switch(position) {
        case CellCenter:
          del_elev = del_elevation_dat[idx2d];
          break;
        case Face:
          del_elev = 0.5 * (del_elevation_dat[idx2d] + del_elevation_dat[idx2d + stride]);
          break;
        case OneQuarter:
          del_elev = 0.75 * del_elevation_dat[idx2d] + 0.25 * del_elevation_dat[idx2d + stride];
          break;
        case ThreeQuarter:
          del_elev = 0.25 * del_elevation_dat[idx2d] + 0.75 * del_elevation_dat[idx2d + stride];
          break;
      }

      length_dat[idx] = sqrt(1.0 + del_elev * del_elev);
    });
  }

  VectorUpdateCommHandle *handle = InitVectorUpdate(length, VectorUpdateAll);
  FinalizeVectorUpdate(handle);

  return;
}

/*--------------------------------------------------------------------------
 * TerrainFollowingTransformJacobian
 *--------------------------------------------------------------------------*/

void TerrainFollowingTransformJacobian(Vector *jacobian)
{
  InitVectorAll(jacobian, 1.0);
  return;
}

/*--------------------------------------------------------------------------
 * TerrainFollowingTransformJacobian
 *--------------------------------------------------------------------------*/

void TerrainFollowingTransformTranslationFactors(Vector *T_uu,
                                                 Vector *T_uv,
                                                 Vector *T_uw,
                                                 Coordinate coordinate,
                                                 PositionInCell position)
{
  if(coordinate == W) {
    InitVectorAll(T_uu, 1.0);
    InitVectorAll(T_uv, 0.0);
    InitVectorAll(T_uw, 0.0);
    return;
  }

  PFModule     *this_module = ThisPFModule;
  InstanceXtra *instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);

  Vector *del_elevation = NULL;
  if (coordinate == U)
    del_elevation = instance_xtra->delevation_du;
  else if (coordinate == V)
    del_elevation = instance_xtra->delevation_dv;

  Grid *grid = VectorGrid(T_uu);

  int is = 0;
  ForSubgridI(is, GridSubgrids(grid))
  {
    Subgrid *subgrid = GridSubgrid(grid, is);

    Subvector *T_uu_sub = VectorSubvector(T_uu, is);
    double *T_uu_dat = SubvectorData(T_uu_sub);
    Subvector *T_uv_sub = VectorSubvector(T_uv, is);
    double *T_uv_dat = SubvectorData(T_uv_sub);
    Subvector *T_uw_sub = VectorSubvector(T_uw, is);
    double *T_uw_dat = SubvectorData(T_uw_sub);

    Subvector *del_elevation_sub = VectorSubvector(del_elevation, is);
    double *del_elevation_dat = SubvectorData(del_elevation_sub);

    int iu = SubgridIX(subgrid);
    int iv = SubgridIY(subgrid);
    int iw = SubgridIZ(subgrid);

    int nu = SubgridNX(subgrid);
    int nv = SubgridNY(subgrid);
    int nw = SubgridNZ(subgrid);

    int nu_v = SubvectorNX(T_uu_sub);
    int nv_v = SubvectorNY(T_uu_sub);
    int nw_v = SubvectorNZ(T_uu_sub);

    int stride_u = 1;
    int stride_v = nu_v;

    int stride = 0;
    switch (coordinate)
    {
      case U:
        stride = stride_u;
        break;
      case V:
        stride = stride_v;
        break;
    }

    int i = -1, j = -1, k = -1;
    int idx = SubvectorEltIndex(T_uu_sub, i, j, k);
    BoxLoopI1(i, j, k, (iu-1), (iv-1), (iw-1), (nu+1), (nv+1), (nw+1),
              idx, nu_v, nv_v, nw_v, 1, 1, 1,
    {
      int idx2d = SubvectorEltIndex(del_elevation_sub, i, j, 0);

      double diff = del_elevation_dat[idx2d] - del_elevation_dat[idx2d + stride];

      double t = 0.0;
      switch(position) {
        case CellCenter:
          t = 0.5 * diff;
          break;
        case Face:
          t = 0.0;
          break;
        case OneQuarter:
          t = 0.25 * diff;
          break;
        case ThreeQuarter:
          t = -0.25 * diff;
          break;
      }

      T_uu_dat[idx] = 1.0;
      if (coordinate == U) {
        T_uv_dat[idx] = 0.0;
        T_uw_dat[idx] = t;
      }
      else if (coordinate == V) {
        T_uv_dat[idx] = t;
        T_uw_dat[idx] = 0.0;
      }
    });
  }

  VectorUpdateCommHandle *handle = NULL; 
  handle = InitVectorUpdate(T_uu, VectorUpdateAll);
  FinalizeVectorUpdate(handle);
  handle = InitVectorUpdate(T_uv, VectorUpdateAll);
  FinalizeVectorUpdate(handle);
  handle = InitVectorUpdate(T_uw, VectorUpdateAll);
  FinalizeVectorUpdate(handle);

  return;
}

/*--------------------------------------------------------------------------
 * TerrainFollowingTransformMetricContravariant
 *--------------------------------------------------------------------------*/

void TerrainFollowingTransformMetricContravariant(Vector *g_uu, Vector *g_uv,
                                                  Vector *g_uw, Vector *g_vv,
                                                  Vector *g_vw, Vector *g_ww)
{
  PFModule     *this_module = ThisPFModule;
  InstanceXtra *instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);

  Vector *delevation_du = instance_xtra->delevation_du;
  Vector *delevation_dv = instance_xtra->delevation_dv;

  Grid *grid_3d = VectorGrid(g_uu);
  Grid *grid_2d = VectorGrid(delevation_du);

  int is = 0;
  ForSubgridI(is, GridSubgrids(grid_3d))
  {
    Subgrid *subgrid_3d = GridSubgrid(grid_3d, is);
    Subgrid *subgrid_2d = GridSubgrid(grid_2d, is);

    Subvector *g_uu_sub = VectorSubvector(g_uu, is);
    double *g_uu_dat = SubvectorData(g_uu_sub);
    Subvector *g_uv_sub = VectorSubvector(g_uv, is);
    double *g_uv_dat = SubvectorData(g_uv_sub);
    Subvector *g_uw_sub = VectorSubvector(g_uw, is);
    double *g_uw_dat = SubvectorData(g_uw_sub);
    Subvector *g_vv_sub = VectorSubvector(g_vv, is);
    double *g_vv_dat = SubvectorData(g_vv_sub);
    Subvector *g_vw_sub = VectorSubvector(g_vw, is);
    double *g_vw_dat = SubvectorData(g_vw_sub);
    Subvector *g_ww_sub = VectorSubvector(g_ww, is);
    double *g_ww_dat = SubvectorData(g_ww_sub);

    Subvector *delevation_du_sub = VectorSubvector(delevation_du, is);
    double *delevation_du_dat = SubvectorData(delevation_du_sub);
    Subvector *delevation_dv_sub = VectorSubvector(delevation_dv, is);
    double *delevation_dv_dat = SubvectorData(delevation_dv_sub);

    int iu = SubgridIX(subgrid_3d);
    int iv = SubgridIY(subgrid_3d);
    int iw = SubgridIZ(subgrid_3d);
    int ww = SubgridIZ(subgrid_2d);

    int nu = SubgridNX(subgrid_3d);
    int nv = SubgridNY(subgrid_3d);
    int nw = SubgridNZ(subgrid_3d);

    int nu_v = SubvectorNX(g_uu_sub);
    int nv_v = SubvectorNY(g_uu_sub);
    int nw_v = SubvectorNZ(g_uu_sub);

    int i = -1, j = -1, k = -1;
    int idx = SubvectorEltIndex(g_uu_sub, i, j, k);
    BoxLoopI1(i, j, k, (iu-1), (iv-1), (iw-1), (nu+2), (nv+2), (nw+2),
              idx, nu_v, nv_v, nw_v, 1, 1, 1,
    {
      int idx2d = SubvectorEltIndex(delevation_du_sub, i, j, ww);

      double delev_du = delevation_du_dat[idx2d];
      double delev_dv = delevation_dv_dat[idx2d];

      g_uu_dat[idx] = 1.0;
      g_uv_dat[idx] = 0.0;
      g_uw_dat[idx] = - delev_du;
      g_vv_dat[idx] = 1.0;
      g_vw_dat[idx] = - delev_dv;
      g_ww_dat[idx] = 1.0 + delev_du * delev_du + delev_dv * delev_dv;
    });
  }

  VectorUpdateCommHandle *handle = NULL;
  handle = InitVectorUpdate(g_uu, VectorUpdateAll);
  FinalizeVectorUpdate(handle);
  handle = InitVectorUpdate(g_uv, VectorUpdateAll);
  FinalizeVectorUpdate(handle);
  handle = InitVectorUpdate(g_uw, VectorUpdateAll);
  FinalizeVectorUpdate(handle);
  handle = InitVectorUpdate(g_vv, VectorUpdateAll);
  FinalizeVectorUpdate(handle);
  handle = InitVectorUpdate(g_vw, VectorUpdateAll);
  FinalizeVectorUpdate(handle);
  handle = InitVectorUpdate(g_ww, VectorUpdateAll);
  FinalizeVectorUpdate(handle);

  return;
}

/*--------------------------------------------------------------------------
 * TerrainFollowingTransformInitInstanceXtra
 *--------------------------------------------------------------------------*/

PFModule  *TerrainFollowingTransformInitInstanceXtra(ProblemData *problem_data)
{
  PFModule      *this_module = ThisPFModule;
  InstanceXtra  *instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);
  PublicXtra  *public_xtra = (PublicXtra*)PFModulePublicXtra(this_module);

  if (instance_xtra == NULL)
  {
    instance_xtra = ctalloc(InstanceXtra, 1);
    instance_xtra->problem_data = problem_data;

    Grid *grid2d = VectorGrid(ProblemDataIndexOfDomainTop(problem_data));
    instance_xtra->elevation = NewVectorType(grid2d, 1, 1, vector_cell_centered_2D);
    instance_xtra->delevation_du = NewVectorType(grid2d, 1, 1, vector_cell_centered_2D);
    instance_xtra->delevation_dv = NewVectorType(grid2d, 1, 1, vector_cell_centered_2D);

    switch (public_xtra->input_type)
    {
      case PFBFile:
      {
        InitVectorAll(instance_xtra->elevation, 0.0);
        ReadPFBinary((public_xtra->filename), instance_xtra->elevation);
        VectorUpdateCommHandle *handle = InitVectorUpdate(instance_xtra->elevation, VectorUpdateAll);
        FinalizeVectorUpdate(handle);
        break;
      }

      default:
      {
        InputError("Something went wrong.%s%s\n", "", "");
      }
    }

    ComputeElevationDerivatives(instance_xtra->delevation_du,
                                instance_xtra->delevation_dv,
                                instance_xtra->elevation, problem_data);
  }

  PFModuleInstanceXtra(this_module) = instance_xtra;
  return this_module;
}


/*--------------------------------------------------------------------------
 * TerrainFollowingTransformFreeInstanceXtra
 *--------------------------------------------------------------------------*/

void  TerrainFollowingTransformFreeInstanceXtra()
{
  PFModule      *this_module = ThisPFModule;
  InstanceXtra  *instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);

  if (instance_xtra)
  {
    FreeVector(instance_xtra->elevation);
    FreeVector(instance_xtra->delevation_du);
    FreeVector(instance_xtra->delevation_dv);
    tfree(instance_xtra);
  }
}

/*--------------------------------------------------------------------------
 * TerrainFollowingTransformNewPublicXtra
 *--------------------------------------------------------------------------*/

PFModule  *TerrainFollowingTransformNewPublicXtra(CoordinateTransformMethods *coordinate_transform_methods)
{
  PFModule      *this_module = ThisPFModule;
  PublicXtra    *public_xtra = ctalloc(PublicXtra, 1);

  // coordinate_transform_methods->eval = TerrainFollowingTransformEval;
  coordinate_transform_methods->z_coordinate = TerrainFollowingTransformZCoordinate;
  coordinate_transform_methods->curve_length = TerrainFollowingTransformCurveLength;
  coordinate_transform_methods->jacobian = TerrainFollowingTransformJacobian;
  coordinate_transform_methods->metric_contravariant = TerrainFollowingTransformMetricContravariant;
  coordinate_transform_methods->translation_factors = TerrainFollowingTransformTranslationFactors;

  NameArray input_name_array = NA_NewNameArray("PredefinedFunction PFBFile");

  char key[IDB_MAX_KEY_LEN];
  sprintf(key, "Transform.TerrainFollowing.Elevation.Type");
  char *input_name = GetString(key);
  int input_type = NA_NameToIndexExitOnError(input_name_array, input_name, key);

  switch (input_type)
  {
    case PredefinedFunction:
    {
      InputError("Value <%s> for key <%s> not yet implemented.\n", input_name, key);
      break;
    }

    case PFBFile:
    {
      public_xtra->input_type = input_type;
      sprintf(key, "Transform.TerrainFollowing.Elevation.FileName");
      public_xtra->filename = GetString(key);
      break;
    }

    default:
    {
      InputError("Invalid value <%s> for key <%s>", input_name, key);
      break;
    }
  }

  NA_FreeNameArray(input_name_array);


  PFModulePublicXtra(this_module) = public_xtra;
  return this_module;
}

/*-------------------------------------------------------------------------
 * TerrainFollowingTransformFreePublicXtra
 *-------------------------------------------------------------------------*/

void  TerrainFollowingTransformFreePublicXtra()
{
  PFModule    *this_module = ThisPFModule;
  PublicXtra  *public_xtra = (PublicXtra*)PFModulePublicXtra(this_module);

  if (public_xtra)
  {
    tfree(public_xtra);
  }
}

/*--------------------------------------------------------------------------
 * TerrainFollowingTransformSizeOfTempData
 *--------------------------------------------------------------------------*/

int  TerrainFollowingTransformSizeOfTempData()
{
  return 0;
}

/*--------------------------------------------------------------------------
 * ComputeElevationDerivatives
 *--------------------------------------------------------------------------*/

void ComputeElevationDerivatives(Vector *delevation_du,
                                 Vector *delevation_dv,
                                 Vector *elevation,
                                 ProblemData *problem_data)
{
  Vector *top = ProblemDataIndexOfDomainTop(problem_data);
  Grid *grid2d = VectorGrid(top);

  int isubgrid = 0;
  ForSubgridI(isubgrid, GridSubgrids(grid2d))
  {
    Subgrid *subgrid = GridSubgrid(grid2d, isubgrid);
    Subvector *top_sub = VectorSubvector(top, isubgrid);
    Subvector *elevation_sub = VectorSubvector(elevation, isubgrid);

    int iu = SubgridIX(subgrid);
    int iv = SubgridIY(subgrid);
    int iw = SubgridIZ(subgrid);

    int nu = SubgridNX(subgrid);
    int nv = SubgridNY(subgrid);
    int nw = SubgridNZ(subgrid);

    double du = SubgridDX(subgrid);
    double dv = SubgridDY(subgrid);
    double dw = SubgridDZ(subgrid);

    int nu_v = SubvectorNX(top_sub);
    int nv_v = SubvectorNY(top_sub);
    int nw_v = SubvectorNZ(top_sub);

    int stride_u = 1;
    int stride_v = nu_v;

    double *top_dat = SubvectorData(top_sub);
    double *elevation_dat = SubvectorData(elevation_sub);
    double *delevation_du_dat = SubvectorData(VectorSubvector(delevation_du, isubgrid));
    double *delevation_dv_dat = SubvectorData(VectorSubvector(delevation_dv, isubgrid));

    int i = 0, j = 0, k = 0;
    int idx = SubvectorEltIndex(elevation_sub, i, j, k);
    BoxLoopI1(i, j, k, iu, iv, iw, nu, nv, nw,
              idx, nu_v, nv_v, nw_v, 1, 1, 1,
    {
      int idx_negu = idx - stride_u;
      int idx_posu = idx + stride_u;
      int idx_negv = idx - stride_v;
      int idx_posv = idx + stride_v;

      // int is_face_negu = top_dat[idx_negu] < 0;
      // int is_face_posu = top_dat[idx_posu] < 0;
      // int is_face_negv = top_dat[idx_negv] < 0;
      // int is_face_posv = top_dat[idx_posv] < 0;
      int is_face_negu = i - 1 < iu;
      int is_face_posu = i + 1 >= iu + nu;
      int is_face_negv = j - 1 < iv;
      int is_face_posv = j + 1 >= iv + nv;

      double elevation_negu = is_face_negu ? elevation_dat[idx] : elevation_dat[idx_negu];
      double elevation_posu = is_face_posu ? elevation_dat[idx] : elevation_dat[idx_posu];
      double elevation_negv = is_face_negv ? elevation_dat[idx] : elevation_dat[idx_negv];
      double elevation_posv = is_face_posv ? elevation_dat[idx] : elevation_dat[idx_posv];

      double idel_u = (is_face_negu || is_face_posu) ? 1.0 / du : 0.5 / du;
      double idel_v = (is_face_negv || is_face_posv) ? 1.0 / dv : 0.5 / dv;

      delevation_du_dat[idx] = idel_u * (elevation_posu - elevation_negu);
      delevation_dv_dat[idx] = idel_v * (elevation_posv - elevation_negv);

      if(is_face_negu) {
        delevation_du_dat[idx_negu] = delevation_du_dat[idx];
        delevation_dv_dat[idx_negu] = delevation_dv_dat[idx];
      }
      if(is_face_posu) {
        delevation_du_dat[idx_posu] = delevation_du_dat[idx];
        delevation_dv_dat[idx_posu] = delevation_dv_dat[idx];
      }
      if(is_face_negv) {
        delevation_du_dat[idx_negv] = delevation_du_dat[idx];
        delevation_dv_dat[idx_negv] = delevation_dv_dat[idx];
      }
      if(is_face_posv) {
        delevation_du_dat[idx_posv] = delevation_du_dat[idx];
        delevation_dv_dat[idx_posv] = delevation_dv_dat[idx];
      }
      if(is_face_negu && is_face_negv) {
        delevation_du_dat[idx_negv - stride_u] = delevation_du_dat[idx];
        delevation_dv_dat[idx_negv - stride_u] = delevation_dv_dat[idx];
      }
      if(is_face_posu && is_face_negv) {
        delevation_du_dat[idx_negv + stride_u] = delevation_du_dat[idx];
        delevation_dv_dat[idx_negv + stride_u] = delevation_dv_dat[idx];
      }
      if(is_face_negu && is_face_posv) {
        delevation_du_dat[idx_posv - stride_u] = delevation_du_dat[idx];
        delevation_dv_dat[idx_posv - stride_u] = delevation_dv_dat[idx];
      }
      if(is_face_posu && is_face_posv) {
        delevation_du_dat[idx_posv + stride_u] = delevation_du_dat[idx];
        delevation_dv_dat[idx_posv + stride_u] = delevation_dv_dat[idx];
      }
    });
  }

  VectorUpdateCommHandle *handle = NULL;
  handle = InitVectorUpdate(delevation_du, VectorUpdateAll);
  FinalizeVectorUpdate(handle);
  handle = InitVectorUpdate(delevation_dv, VectorUpdateAll);
  FinalizeVectorUpdate(handle);

  return;
}