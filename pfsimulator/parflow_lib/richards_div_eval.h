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

#ifndef _RICHARDS_DIV_EVAL_HEADER
#define _RICHARDS_DIV_EVAL_HEADER

#include "parflow.h"

/* This function computes the density at a cell face. It makes a weighted
 * average with the values in each cell. */
LoopFriendlyFunction
double DensityMean(int idx, int stride,
                   double *lA_dat, double *lB_dat,
                   double *density_dat)
{
  return (density_dat[idx] * lA_dat[idx] +
          density_dat[idx + stride] * lB_dat[idx]) / 
          (lA_dat[idx] + lB_dat[idx]);
}

/* This function computes the derivatives of the hydraulic head. */
LoopFriendlyFunction
void HydraulicHeadDel(int idx, int stride_u, int stride_v, int stride_w,
                      double du, double dv, double dw, double *dh,
                      double *lengthA_dat, double *lengthB_dat, 
                      double * pressure_dat, double *z_dat,
                      double *density_dat, double gravity)
{
  double density = DensityMean(idx, stride_u, lengthA_dat, lengthB_dat, density_dat);

  double *p_ = pressure_dat;
  double *z_ = z_dat;
  int su = stride_u, sv = stride_v, sw = stride_w;

  double dh_du = ((p_[idx + su] - p_[idx]) - density * gravity * (z_[idx + su] - z_[idx])) / du;
  double dh_dv = (0.25 / dv) * ((p_[idx + su + sv] - p_[idx + su - sv] + p_[idx + sv] - p_[idx - sv]) + density * gravity * (z_[idx + su + sv] - z_[idx + su - sv] + z_[idx + sv] - z_[idx - sv]));
  double dh_dw = (0.25 / dw) * ((p_[idx + su + sw] - p_[idx + su - sw] + p_[idx + sw] - p_[idx - sw]) + density * gravity * (z_[idx + su + sw] - z_[idx + su - sw] + z_[idx + sw] - z_[idx - sw]));

  dh[U] = dh_du;
  dh[V] = dh_dv;
  dh[W] = dh_dw;

  return;
}


LoopFriendlyFunction
double RelativePermeabilityMean(int idx, int stride, double dh,
                                double *rel_perm_dat, double *density_dat)
{
  return dh < 0 ? rel_perm_dat[idx] * density_dat[idx] :
                  rel_perm_dat[idx + stride] * density_dat[idx + stride];
}

/* This function computes the Darcy flux at a cell face. It takes al quantities 
 * evaluated at the respective cell face. */
LoopFriendlyFunction
double DarcyFlux(double dArea, double rel_perm_times_density, double viscosity,
                 double dh_du, double dh_dv, double dh_dw,
                 double K_uu, double K_uv, double K_uw)
{
  return -(dArea * rel_perm_times_density / viscosity) * (K_uu * dh_du + K_uv * dh_dv + K_uw * dh_dw);
}

/* This function computes the flux at a cell face.
 * It is important to note that the flux computation in each face:
 * posx, posy and posz is related by a cyclic permutation of the
 * expression. So we can reuse the function for each face, as
 * long as we do a cyclic permutation of the cell dimensions du, dv, dw
 * and the same cyclic permutation of the strides and hydraulic conductivity.
 * */
LoopFriendlyFunction
double FluxAtFace(int idx, int stride_u, double dArea, double *dh,
                  double *density_dat, double *rel_perm_dat, double viscosity,
                  double *K_uu, double *K_uv, double *K_uw)
{
  double rel_perm_times_density = RelativePermeabilityMean(idx, stride_u,
                                                           dh[U],
                                                           rel_perm_dat,
                                                           density_dat);

  return DarcyFlux(dArea, rel_perm_times_density, viscosity,
                   dh[U], dh[V], dh[W], K_uu[idx], K_uv[idx], K_uw[idx]);
}

#endif // _RICHARDS_DIV_EVAL_HEADER