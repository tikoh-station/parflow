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

#ifndef _HYDRAULIC_CONDUCTIVITY_HEADER
#define _HYDRAULIC_CONDUCTIVITY_HEADER

typedef struct {

  Vector *K_uu;
  Vector *K_uv;
  Vector *K_uw;
  Vector *K_vu;
  Vector *K_vv;
  Vector *K_vw;
  Vector *K_wu;
  Vector *K_wv;
  Vector *K_ww;

} PermeabilityTensor;

PermeabilityTensor* NewPermeabilityTensor(Grid *grid,
                                          int nc,
                                          int num_ghost,
                                          enum vector_type type);

void FreePermeabilityTensor(PermeabilityTensor *tensor);

#define PermeabilityTensorUU(tensor) ((tensor)->K_uu)
#define PermeabilityTensorUV(tensor) ((tensor)->K_uv)
#define PermeabilityTensorUW(tensor) ((tensor)->K_uw)
#define PermeabilityTensorVU(tensor) ((tensor)->K_vu)
#define PermeabilityTensorVV(tensor) ((tensor)->K_vv)
#define PermeabilityTensorVW(tensor) ((tensor)->K_vw)
#define PermeabilityTensorWU(tensor) ((tensor)->K_wu)
#define PermeabilityTensorWV(tensor) ((tensor)->K_wv)
#define PermeabilityTensorWW(tensor) ((tensor)->K_ww)


// void (*ComputePermeabilityTensorInvoke)(ProblemData *problem_data);
// void ComputePermeabilityTensor(ProblemData *problem_data);
// PFModule  *(*ComputePermeabilityTensorInitInstanceXtraInvoke)();
// PFModule  *ComputePermeabilityTensorInitInstanceXtra();
// void  ComputePermeabilityTensorFreeInstanceXtra();
// PFModule  *ComputePermeabilityTensorNewPublicXtra();
// void  ComputePermeabilityTensorFreePublicXtra();
// int  ComputePermeabilityTensorSizeOfTempData();

// void ComputeFacePermeabilityTensor(Vector *K_uu, Vector *K_uv, Vector *K_uw,
//                                    Vector *KC_uu, Vector *KC_uv, Vector *KC_uw,
//                                    Vector *TA_uu, Vector *TA_uv, Vector *TA_uw,
//                                    Vector *TB_uu, Vector *TB_uv, Vector *TB_uw,
//                                    Coordinate coord, ProblemData *problem_data)

#endif // _HYDRAULIC_CONDUCTIVITY_HEADER