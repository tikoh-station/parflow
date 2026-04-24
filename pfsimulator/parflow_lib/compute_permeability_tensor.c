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

#define Isotropic 0
#define Anisotropic 1

#define Constant 0
#define PFBFile 1


PermeabilityTensor* NewPermeabilityTensor(Grid *grid,
                                          int nc,
                                          int num_ghost,
                                          enum vector_type type)
{
  PermeabilityTensor *tensor = ctalloc(PermeabilityTensor, 1);

  PermeabilityTensorUU(tensor) = NewVectorType(grid, nc, num_ghost, type);
  PermeabilityTensorUV(tensor) = NewVectorType(grid, nc, num_ghost, type);
  PermeabilityTensorUW(tensor) = NewVectorType(grid, nc, num_ghost, type);
  PermeabilityTensorVU(tensor) = NewVectorType(grid, nc, num_ghost, type);
  PermeabilityTensorVV(tensor) = NewVectorType(grid, nc, num_ghost, type);
  PermeabilityTensorVW(tensor) = NewVectorType(grid, nc, num_ghost, type);
  PermeabilityTensorWU(tensor) = NewVectorType(grid, nc, num_ghost, type);
  PermeabilityTensorWV(tensor) = NewVectorType(grid, nc, num_ghost, type);
  PermeabilityTensorWW(tensor) = NewVectorType(grid, nc, num_ghost, type);

  return tensor;
}

void FreePermeabilityTensor(PermeabilityTensor *tensor)
{
  FreeVector(PermeabilityTensorUU(tensor));
  FreeVector(PermeabilityTensorUV(tensor));
  FreeVector(PermeabilityTensorUW(tensor));
  FreeVector(PermeabilityTensorVU(tensor));
  FreeVector(PermeabilityTensorVV(tensor));
  FreeVector(PermeabilityTensorVW(tensor));
  FreeVector(PermeabilityTensorWU(tensor));
  FreeVector(PermeabilityTensorWV(tensor));
  FreeVector(PermeabilityTensorWW(tensor));

  return;
}


typedef void PublicXtra;

typedef struct {

  PFModule *coordinate_transform;

} InstanceXtra;

/*--------------------------------------------------------------------------
 * ComputePermeabilityTensor
 *--------------------------------------------------------------------------*/

void ComputePermeabilityTensor(ProblemData *problem_data)
{
  PFModule      *this_module = ThisPFModule;
  InstanceXtra  *instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);
  PFModule      *coordinate_transform = (instance_xtra->coordinate_transform);

  char key[IDB_MAX_KEY_LEN];
  
  NameArray medium_type_na = NA_NewNameArray("Isotropic Anisotropic");

  /* For now, geometries are not implemented. */

  sprintf(key, "Geom.Perm.Medium");
  char *medium_type_name = GetString(key);
  int medium_type = NA_NameToIndexExitOnError(medium_type_na, medium_type_name, key);
  NA_FreeNameArray(medium_type_na);

  Grid *grid = VectorGrid(ProblemDataPorosity(problem_data));
  PermeabilityTensor* tensor = NewPermeabilityTensor(grid, 1, 1, vector_cell_centered);
  switch (medium_type)
  {
    case Isotropic:
    {
      InitIsotropicPermeabilityTensor(tensor, coordinate_transform, problem_data);
      break;
    }

    case Anisotropic:
    {
      amps_Printf("Anisotropic permeability tensor not yet implemented.\n");
      break;
    }
  }

  Vector *TA_u = NewVectorType(grid, 1, 1, vector_cell_centered);
  Vector *TA_v = NewVectorType(grid, 1, 1, vector_cell_centered);
  Vector *TA_w = NewVectorType(grid, 1, 1, vector_cell_centered);
  Vector *TB_u = NewVectorType(grid, 1, 1, vector_cell_centered);
  Vector *TB_v = NewVectorType(grid, 1, 1, vector_cell_centered);
  Vector *TB_w = NewVectorType(grid, 1, 1, vector_cell_centered);

  Vector *KC_uu = PermeabilityTensorUU(tensor);
  Vector *KC_uv = PermeabilityTensorUV(tensor);
  Vector *KC_uw = PermeabilityTensorUW(tensor);
  Vector *KC_vu = PermeabilityTensorVU(tensor);
  Vector *KC_vv = PermeabilityTensorVV(tensor);
  Vector *KC_vw = PermeabilityTensorVW(tensor);
  Vector *KC_wu = PermeabilityTensorWU(tensor);
  Vector *KC_wv = PermeabilityTensorWV(tensor);
  Vector *KC_ww = PermeabilityTensorWW(tensor);

  PermeabilityTensor *K = ProblemDataPermeabilityTensor(problem_data);
  Vector *K_uu = PermeabilityTensorUU(K);
  Vector *K_uv = PermeabilityTensorUV(K);
  Vector *K_uw = PermeabilityTensorUW(K);
  Vector *K_vu = PermeabilityTensorVU(K);
  Vector *K_vv = PermeabilityTensorVV(K);
  Vector *K_vw = PermeabilityTensorVW(K);
  Vector *K_wu = PermeabilityTensorWU(K);
  Vector *K_wv = PermeabilityTensorWV(K);
  Vector *K_ww = PermeabilityTensorWW(K);

  // U Face Permutation: (u', v', w') -> (u, v, w)
  CoordinateTransformTranslationFactors(coordinate_transform, 
                                        TA_u, TA_v, TA_w, U, OneQuarter);

  CoordinateTransformTranslationFactors(coordinate_transform, 
                                        TB_u, TB_v, TB_w, U, ThreeQuarter);

  ComputeFacePermeabilityTensor(K_uu, K_uv, K_uw, KC_uu, KC_uv, KC_uw,
                                TA_u, TA_v, TA_w, TB_u, TB_v, TB_w,
                                U, problem_data);

  // V Face Permutation: (u', v', w') -> (v, w, u)
  CoordinateTransformTranslationFactors(coordinate_transform, 
                                        TA_v, TA_w, TA_u, V, OneQuarter);

  CoordinateTransformTranslationFactors(coordinate_transform, 
                                        TB_v, TB_w, TB_u, V, ThreeQuarter);

  ComputeFacePermeabilityTensor(K_vv, K_vw, K_vu, KC_vv, KC_vw, KC_vu,
                                TA_v, TA_w, TA_u, TB_v, TB_w, TB_u,
                                V, problem_data);

  // W Face Permutation: (u', v', w') -> (w, u, v)
  CoordinateTransformTranslationFactors(coordinate_transform, 
                                        TA_w, TA_u, TA_v, W, OneQuarter);

  CoordinateTransformTranslationFactors(coordinate_transform, 
                                        TB_w, TB_u, TB_v, W, ThreeQuarter);

  ComputeFacePermeabilityTensor(K_ww, K_wu, K_wv, KC_ww, KC_wu, KC_wv,
                                TA_w, TA_u, TA_v, TB_w, TB_u, TB_v,
                                W, problem_data);

  FreeVector(TA_u);
  FreeVector(TA_v);
  FreeVector(TA_w);
  FreeVector(TB_u);
  FreeVector(TB_v);
  FreeVector(TB_w);

  FreePermeabilityTensor(tensor);

  return;
}


/*--------------------------------------------------------------------------
 * ComputePermeabilityTensorInitInstanceXtra
 *--------------------------------------------------------------------------*/

PFModule  *ComputePermeabilityTensorInitInstanceXtra(PFModule *coordinate_transform)
{
  PFModule      *this_module = ThisPFModule;
  InstanceXtra  *instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);

  if (instance_xtra == NULL)
  {
    instance_xtra = ctalloc(InstanceXtra, 1);

    instance_xtra->coordinate_transform = coordinate_transform;
  }

  PFModuleInstanceXtra(this_module) = instance_xtra;
  return this_module;
}


/*--------------------------------------------------------------------------
 * ComputePermeabilityTensorFreeInstanceXtra
 *--------------------------------------------------------------------------*/

void  ComputePermeabilityTensorFreeInstanceXtra()
{
  PFModule      *this_module = ThisPFModule;
  InstanceXtra  *instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);

  if (instance_xtra)
  {
    tfree(instance_xtra);
  }
}

/*--------------------------------------------------------------------------
 * ComputePermeabilityTensorNewPublicXtra
 *--------------------------------------------------------------------------*/

PFModule  *ComputePermeabilityTensorNewPublicXtra()
{
  PFModule      *this_module = ThisPFModule;
  PublicXtra    *public_xtra = ctalloc(PublicXtra, 1);

  PFModulePublicXtra(this_module) = public_xtra;
  return this_module;
}

/*-------------------------------------------------------------------------
 * ComputePermeabilityTensorFreePublicXtra
 *-------------------------------------------------------------------------*/

void  ComputePermeabilityTensorFreePublicXtra()
{
  PFModule    *this_module = ThisPFModule;
  PublicXtra  *public_xtra = (PublicXtra*)PFModulePublicXtra(this_module);

  if (public_xtra)
  {
    tfree(public_xtra);
  }
}

/*--------------------------------------------------------------------------
 * ComputePermeabilityTensorSizeOfTempData
 *--------------------------------------------------------------------------*/

int  ComputePermeabilityTensorSizeOfTempData()
{
  return 0;
}


/*--------------------------------------------------------------------------
 * InitIsotropicPermeabilityTensor
 *--------------------------------------------------------------------------*/

void InitIsotropicPermeabilityTensor(PermeabilityTensor *perm_tensor,
                                     PFModule *coordinate_transform,
                                     ProblemData *problem_data)
{
  char key[IDB_MAX_KEY_LEN];

  NameArray input_type_na = NA_NewNameArray("Constant PFBFile");

  sprintf(key, "Geom.Perm.Type");
  char *input_type_name = GetString(key);
  int input_type = NA_NameToIndexExitOnError(input_type_na, input_type_name, key);
  NA_FreeNameArray(input_type_na);

  Grid *grid = VectorGrid(ProblemDataPorosity(problem_data));
  Vector *permeability = NewVectorType(grid, 1, 1, vector_cell_centered);

  switch (input_type)
  {
    case Constant:
    {
      double value = GetDouble("Geom.Perm.Value");
      InitVectorAll(permeability, value);

      VectorUpdateCommHandle *handle = InitVectorUpdate(permeability, VectorUpdateAll);
      FinalizeVectorUpdate(handle);
      break;
    }

    case PFBFile:
    {
      char *filename = GetString("Geom.Perm.FileName");
      ReadPFBinary(filename, permeability);

      VectorUpdateCommHandle *handle = InitVectorUpdate(permeability, VectorUpdateAll);
      FinalizeVectorUpdate(handle);
      break;
    }
  }

  /* Set Permeability Tensor */

  Vector *K_uu = PermeabilityTensorUU(perm_tensor);
  Vector *K_uv = PermeabilityTensorUV(perm_tensor);
  Vector *K_uw = PermeabilityTensorUW(perm_tensor);
  Vector *K_vu = PermeabilityTensorVU(perm_tensor);
  Vector *K_vv = PermeabilityTensorVV(perm_tensor);
  Vector *K_vw = PermeabilityTensorVW(perm_tensor);
  Vector *K_wu = PermeabilityTensorWU(perm_tensor);
  Vector *K_wv = PermeabilityTensorWV(perm_tensor);
  Vector *K_ww = PermeabilityTensorWW(perm_tensor);

  CoordinateTransformMetricContravariant(coordinate_transform, 
                                         K_uu, K_uv, K_uw, K_vv, K_vw, K_ww);

  int is = 0;
  ForSubgridI(is, GridSubgrids(grid))
  {
    Subvector *perm_sub = VectorSubvector(permeability, is);
    double *perm_dat = SubvectorData(perm_sub);

    int iu_v = SubvectorIX(perm_sub);
    int iv_v = SubvectorIY(perm_sub);
    int iw_v = SubvectorIZ(perm_sub);

    int nu_v = SubvectorNX(perm_sub);
    int nv_v = SubvectorNY(perm_sub);
    int nw_v = SubvectorNZ(perm_sub);

    double *K_uu_dat = SubvectorData(VectorSubvector(K_uu, is));
    double *K_uv_dat = SubvectorData(VectorSubvector(K_uv, is));
    double *K_uw_dat = SubvectorData(VectorSubvector(K_uw, is));
    double *K_vu_dat = SubvectorData(VectorSubvector(K_vu, is));
    double *K_vv_dat = SubvectorData(VectorSubvector(K_vv, is));
    double *K_vw_dat = SubvectorData(VectorSubvector(K_vw, is));
    double *K_wu_dat = SubvectorData(VectorSubvector(K_wu, is));
    double *K_wv_dat = SubvectorData(VectorSubvector(K_wv, is));
    double *K_ww_dat = SubvectorData(VectorSubvector(K_ww, is));

    int i = 0, j = 0, k = 0, idx = 0;
    BoxLoopI1(i, j, k, iu_v, iv_v, iw_v, nu_v, nv_v, nw_v,
              idx, nu_v, nv_v, nw_v, 1, 1, 1,
    {
      K_uu_dat[idx] *= perm_dat[idx];
      K_uv_dat[idx] *= perm_dat[idx];
      K_uw_dat[idx] *= perm_dat[idx];
      K_vu_dat[idx] = K_uv_dat[idx];
      K_vv_dat[idx] *= perm_dat[idx];
      K_vw_dat[idx] *= perm_dat[idx];
      K_wu_dat[idx] = K_uw_dat[idx];
      K_wv_dat[idx] = K_vw_dat[idx];
      K_ww_dat[idx] *= perm_dat[idx];
    });
  }

  VectorUpdateCommHandle *handle = NULL;
  handle = InitVectorUpdate(K_uu, VectorUpdateAll);
  FinalizeVectorUpdate(handle);
  handle = InitVectorUpdate(K_uv, VectorUpdateAll);
  FinalizeVectorUpdate(handle);
  handle = InitVectorUpdate(K_uw, VectorUpdateAll);
  FinalizeVectorUpdate(handle);
  handle = InitVectorUpdate(K_vu, VectorUpdateAll);
  FinalizeVectorUpdate(handle);
  handle = InitVectorUpdate(K_vv, VectorUpdateAll);
  FinalizeVectorUpdate(handle);
  handle = InitVectorUpdate(K_vw, VectorUpdateAll);
  FinalizeVectorUpdate(handle);
  handle = InitVectorUpdate(K_wu, VectorUpdateAll);
  FinalizeVectorUpdate(handle);
  handle = InitVectorUpdate(K_wv, VectorUpdateAll);
  FinalizeVectorUpdate(handle);
  handle = InitVectorUpdate(K_ww, VectorUpdateAll);
  FinalizeVectorUpdate(handle);

  FreeVector(permeability);

  return;
}


/*--------------------------------------------------------------------------
 * ComputeFacePermeabilityTensor
 *--------------------------------------------------------------------------*/

void ComputeFacePermeabilityTensor(Vector *K_uu, Vector *K_uv, Vector *K_uw,
                                   Vector *KC_uu, Vector *KC_uv, Vector *KC_uw,
                                   Vector *TA_uu, Vector *TA_uv, Vector *TA_uw,
                                   Vector *TB_uu, Vector *TB_uv, Vector *TB_uw,
                                   Coordinate coord, ProblemData *problem_data)
{
  Grid *grid = VectorGrid(K_uu);

  int isubgrid = 0;
  ForSubgridI(isubgrid, GridSubgrids(grid))
  {
    Subgrid *subgrid = GridSubgrid(grid, isubgrid);
    Subvector *K_uu_sub = VectorSubvector(K_uu, isubgrid);

    int iu = SubgridIX(subgrid);
    int iv = SubgridIY(subgrid);
    int iw = SubgridIZ(subgrid);

    int nu = SubgridNX(subgrid);
    int nv = SubgridNY(subgrid);
    int nw = SubgridNZ(subgrid);

    int nu_v = SubvectorNX(K_uu_sub);
    int nv_v = SubvectorNY(K_uu_sub);
    int nw_v = SubvectorNZ(K_uu_sub);


    int stride_u = 1;
    int stride_v = nu_v;
    int stride_w = nu_v * stride_v;
    int stride = 0;
    switch (coord)
    {
      case U:
        stride = stride_u;
        break;
      case V:
        stride = stride_v;
        break;
      case W:
        stride = stride_w;
        break;
    }

    // Hydraulic Conductivity in face
    double *K_uu_dat = SubvectorData(VectorSubvector(K_uu, isubgrid));
    double *K_uv_dat = SubvectorData(VectorSubvector(K_uv, isubgrid));
    double *K_uw_dat = SubvectorData(VectorSubvector(K_uw, isubgrid));

    // Hydraulic Conductivity in cell centers
    double *KC_uu_dat = SubvectorData(VectorSubvector(KC_uu, isubgrid));
    double *KC_uv_dat = SubvectorData(VectorSubvector(KC_uv, isubgrid));
    double *KC_uw_dat = SubvectorData(VectorSubvector(KC_uw, isubgrid));

    // Translation Coefficients in cell A
    double *TA_uu_dat = SubvectorData(VectorSubvector(TA_uu, isubgrid));
    double *TA_uv_dat = SubvectorData(VectorSubvector(TA_uv, isubgrid));
    double *TA_uw_dat = SubvectorData(VectorSubvector(TA_uw, isubgrid));

    // Translation Coefficients in cell B
    double *TB_uu_dat = SubvectorData(VectorSubvector(TB_uu, isubgrid));
    double *TB_uv_dat = SubvectorData(VectorSubvector(TB_uv, isubgrid));
    double *TB_uw_dat = SubvectorData(VectorSubvector(TB_uw, isubgrid));

    int i = -1, j = -1, k = -1;
    int idx = SubvectorEltIndex(K_uu_sub, i, j, k);
    BoxLoopI1(i, j, k, (iu-1), (iv-1), (iw-1), (nu+1), (nv+1), (nw+1),
              idx, nu_v, nv_v, nw_v, 1, 1, 1,
    {
      // Hydraulic Conductivity in cell A
      double KA_uu_ = KC_uu_dat[idx];
      double KA_uv_ = KC_uv_dat[idx];
      double KA_uw_ = KC_uw_dat[idx];
      // Hydraulic Conductivity in cell B
      double KB_uu_ = KC_uu_dat[idx + stride];
      double KB_uv_ = KC_uv_dat[idx + stride];
      double KB_uw_ = KC_uw_dat[idx + stride];

      // Translation Coefficients in cell A
      double TA_uu_ = TA_uu_dat[idx];
      double TA_uv_ = TA_uv_dat[idx];
      double TA_uw_ = TA_uw_dat[idx];
      // Translation Coefficients in cell B
      double TB_uu_ = TB_uu_dat[idx];
      double TB_uv_ = TB_uv_dat[idx];
      double TB_uw_ = TB_uw_dat[idx];

      double idenominator = 1.0 / (KA_uu_ * TB_uu_ + KB_uu_ * TA_uu_);

      K_uu_dat[idx] = 2 * KA_uu_ * KB_uu_ * idenominator;

      K_uv_dat[idx] = (KA_uv_ * TA_uu_ * KB_uu_ + KB_uv_ * TB_uu_ * KA_uu_ - KA_uu_ * KB_uu_ * (TA_uv_ + TB_uv_)) * idenominator;

      K_uw_dat[idx] = (KA_uw_ * TA_uu_ * KB_uu_ + KB_uw_ * TB_uu_ * KA_uu_ - KA_uu_ * KB_uu_ * (TA_uw_ + TB_uw_)) * idenominator;
    });
  }

  // Update Vectors
  VectorUpdateCommHandle  *handle = NULL;

  handle = InitVectorUpdate(K_uu, VectorUpdateAll);
  FinalizeVectorUpdate(handle);
  handle = InitVectorUpdate(K_uv, VectorUpdateAll);
  FinalizeVectorUpdate(handle);
  handle = InitVectorUpdate(K_uw, VectorUpdateAll);
  FinalizeVectorUpdate(handle);

  return;
}