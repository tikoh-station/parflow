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

#ifndef _COORDINATE_TRANSFORM_HEADER
#define _COORDINATE_TRANSFORM_HEADER

/* Fill in the ghost points at the edges with values from the interior */

typedef enum { U = 0, V = 1, W = 2 } Coordinate;

typedef enum { // position relative to cell center
  CellCenter = 0,  // at \alpha
  Face = 1,        // at \alpha + (1/2) \Delta/alpha
  OneQuarter = 2,  // at \alpha + (1/4) \Delta/alpha
  ThreeQuarter = 3 // at \alpha + (3/4) \Delta/alpha
} PositionInCell;

// typedef void (*CoordinateTransformEvalInvoke)(int pos);

// void CoordinateTransformEval(PFModule *coordinate_transform, int pos);

typedef void (*CoordinateTransformZCoordinateInvoke)(Vector *z_coords);

void CoordinateTransformZCoordinate(PFModule *coordinate_transform,
                                    Vector *z_coords);

typedef void (*CoordinateTransformCurveLengthInvoke)(Vector *length,
                                                     Coordinate coordinate,
                                                     PositionInCell position);

void CoordinateTransformCurveLength(PFModule *coordinate_transform,
                                    Vector *length,
                                    Coordinate coordinate,
                                    PositionInCell position);

typedef void (*CoordinateTransformJacobianInvoke)(Vector *jacobian);

void CoordinateTransformJacobian(PFModule *coordinate_transform,
                                 Vector *jacobian);

typedef void (*CoordinateTransformMetricContravariantInvoke)(Vector *g_uu,
                                                             Vector *g_uv,
                                                             Vector *g_uw,
                                                             Vector *g_vv,
                                                             Vector *g_vw,
                                                             Vector *g_ww);

void CoordinateTransformMetricContravariant(PFModule *coordinate_transform,
                                            Vector *g_uu, Vector *g_uv,
                                            Vector *g_uw, Vector *g_vv,
                                            Vector *g_vw, Vector *g_ww);


typedef void (*CoordinateTransformTranslationFactorsInvoke)(Vector *T_uu,
                                                            Vector *T_uv,
                                                            Vector *T_uw,
                                                            Coordinate coordinate,
                                                            PositionInCell position);

void CoordinateTransformTranslationFactors(PFModule *coordinate_transform,
                                           Vector *T_uu,
                                           Vector *T_uv,
                                           Vector *T_uw,
                                           Coordinate coordinate,
                                           PositionInCell position);

typedef struct {

  // CoordinateTransformEvalInvoke eval;

  CoordinateTransformZCoordinateInvoke z_coordinate;

  CoordinateTransformCurveLengthInvoke curve_length;

  CoordinateTransformJacobianInvoke jacobian;

  CoordinateTransformMetricContravariantInvoke metric_contravariant;

  CoordinateTransformTranslationFactorsInvoke translation_factors;

} CoordinateTransformMethods;


typedef void (*CoordinateTransformInvoke)(ProblemData *problem_data);

void CoordinateTransform(ProblemData *problem_data);

typedef PFModule  *(*CoordinateTransformInitInstanceXtraInvoke)(ProblemData *problem_data);

PFModule  *CoordinateTransformInitInstanceXtra(ProblemData *problem_data);

void  CoordinateTransformFreeInstanceXtra();

PFModule  *CoordinateTransformNewPublicXtra();

void  CoordinateTransformFreePublicXtra();

int  CoordinateTransformSizeOfTempData();

#endif // _COORDINATE_TRANSFORM_HEADER