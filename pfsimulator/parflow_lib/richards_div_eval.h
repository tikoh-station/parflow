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

#ifndef _RICHARDS_DIV_EVAL
#define _RICHARDS_DIV_EVAL

#include "parflow.h"
#include "v3algebra.h"
#include "morphism.h"

// structure to store element indexes of a subvector
typedef struct {

  v3 pos; // position of the grid cell

  int mid; // cell in the middle
  int rgt; // cell to the right
  int lft; // cell to the left
  int top; // cell to the top
  int bot; // cell to the bottom
  int frt; // cell to the front
  int bck; // cell to the back

  int rgt_top; // cell to the top right diagonal
  int rgt_bot; // cell to the bottom right diagonal
  int lft_top; // cell to the top left diagonal
  int lft_bot; // cell to the bottom left diagonal

  int rgt_frt; // cell to the front right diagonal
  int rgt_bck; // cell to the back right diagonal
  int lft_frt; // cell to the front left diagonal
  int lft_bck; // cell to the back left diagonal

  int top_frt; // cell to the top front diagonal
  int top_bck; // cell to the top back diagonal
  int bot_frt; // cell to the bottom front diagonal
  int bot_bck; // cell to the bottom back diagonal

} StencilIndx;

void SubvectorStencilIndx(StencilIndx *indx, Subvector *sub,
                          Subgrid *subgrid, int i, int j, int k);
void SubvectorStencilIndxRgt(StencilIndx *indx, Subvector *sub,
                             Subgrid *subgrid, int i, int j, int k);
void SubvectorStencilIndxTop(StencilIndx *indx, Subvector *sub,
                             Subgrid *subgrid, int i, int j, int k);
void SubvectorStencilIndxFrt(StencilIndx *indx, Subvector *sub,
                             Subgrid *subgrid, int i, int j, int k);



// richards divergence term from the right cell
double RichardsDivergenceRgt(StencilIndx *S, Morphism *my_morphism, double *pp, double *dp, double *rpp, double *permxp, double *permyp, double *permzp, double gravity, double viscosity, double du, double dv, double dw);

// richards divergence term from the top cell
double RichardsDivergenceTop(StencilIndx *S, Morphism *my_morphism, double *pp, double *dp, double *rpp, double *permxp, double *permyp, double *permzp, double gravity, double viscosity, double du, double dv, double dw);

// richards divergence term from the front cell
double RichardsDivergenceFrt(StencilIndx *S, Morphism *my_morphism, double *pp, double *dp, double *rpp, double *permxp, double *permyp, double *permzp, double gravity, double viscosity, double du, double dv, double dw);

#endif // _RICHARDS_DIV_EVAL