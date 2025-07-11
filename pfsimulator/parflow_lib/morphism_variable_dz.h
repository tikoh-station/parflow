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

#ifndef _MORPHISM_VARIABLE_DZ_HEADER
#define _MORPHISM_VARIABLE_DZ_HEADER

#include "morphism.h"

typedef void (*MorphismVariableDZInvoke)(Morphism **my_morphism);
void MorphismVariableDZ(Morphism **my_morphism);
PFModule* MorphismVariableDZInitInstanceXtra();
void MorphismVariableDZFreeInstanceXtra();
PFModule* MorphismVariableDZNewPublicXtra();
void MorphismVariableDZFreePublicXtra();
int MorphismVariableDZSizeOfTempData();

v3 MorphismVariableDZEval(PFModule *morphism_core, v3 zeta);
v3 MorphismVariableDZInverse(PFModule *morphism_core, v3 x);
v3basis MorphismVariableDZDel(PFModule *morphism_core, v3 zeta);
v3basis MorphismVariableDZDelInverse(PFModule *morphism_core, v3 zeta);
double MorphismVariableDZJacobian(PFModule *morphism_core, v3 zeta);

double ZFunction(double zlabel);
double ZFunctionInverse(double z);
double ZFunctionDel(double zlabel);

#endif // _MORPHISM_VARIABLE_DZ_HEADER