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

#ifndef _TERRAIN_FOLLOWING_TRANSFORM_HEADER
#define _TERRAIN_FOLLOWING_TRANSFORM_HEADER

void TerrainFollowingTransform();

void TerrainFollowingTransformEval(int pos);

void TerrainFollowingTransformZCoordinate(Vector *z_coords);

void TerrainFollowingTransformCurveLength(Vector *length, Coordinate coordinate,
                                   PositionInCell position);

void TerrainFollowingTransformJacobian(Vector *jacobian);

void TerrainFollowingTransformTranslationFactors(Vector *T_uu,
                                          Vector *T_uv,
                                          Vector *T_uw,
                                          Coordinate coordinate,
                                          PositionInCell position);

void TerrainFollowingTransformMetricContravariant(Vector *g_uu, Vector *g_uv,
                                           Vector *g_uw, Vector *g_vv,
                                           Vector *g_vw, Vector *g_ww);

typedef PFModule *(*TerrainFollowingTransformInitInstanceXtraInvoke)(ProblemData *problem_data);

PFModule *TerrainFollowingTransformInitInstanceXtra(ProblemData *problem_data);

void TerrainFollowingTransformFreeInstanceXtra();

typedef PFModule *(*TerrainFollowingTransformNewPublicXtraInvoke)(CoordinateTransformMethods *coordinate_transform_methods);

PFModule *TerrainFollowingTransformNewPublicXtra(CoordinateTransformMethods *coordinate_transform_methods);

void TerrainFollowingTransformFreePublicXtra();

int TerrainFollowingTransformSizeOfTempData();

void ComputeElevationDerivatives(Vector *delevation_du,
                                 Vector *delevation_dv,
                                 Vector *elevation,
                                 ProblemData *problem_data);

#endif // _TERRAIN_FOLLOWING_TRANSFORM_HEADER