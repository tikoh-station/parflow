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

#include "richards_div_eval.h"


#define PMean(a, b, c, d)    HarmonicMean(c, d)
#define PMeanDZ(a, b, c, d)     HarmonicMeanDZ(a, b, c, d)
#define RPMean(a, b, c, d)   UpstreamMean(a, b, c, d)
#define Mean(a, b)            ArithmeticMean(a, b)


void SubvectorStencilIndx(StencilIndx *indx, Subvector *sub,
                          Subgrid *subgrid, int i, int j, int k)
{
  double u = RealSpaceX(i, SubgridRX(subgrid));
  double v = RealSpaceY(j, SubgridRY(subgrid));
  double w = RealSpaceZ(k, SubgridRZ(subgrid));
  indx->pos = v3_init(u, v, w);

  indx->mid = SubvectorEltIndex(sub, i, j, k);
  indx->rgt = SubvectorEltIndex(sub, i+1, j, k);
  indx->lft = SubvectorEltIndex(sub, i-1, j, k);
  indx->top = SubvectorEltIndex(sub, i, j+1, k);
  indx->bot = SubvectorEltIndex(sub, i, j-1, k);
  indx->frt = SubvectorEltIndex(sub, i, j, k+1);
  indx->bck = SubvectorEltIndex(sub, i, j, k-1);

  indx->rgt_top = SubvectorEltIndex(sub, i+1, j+1, k);
  indx->rgt_bot = SubvectorEltIndex(sub, i+1, j-1, k);
  indx->lft_top = SubvectorEltIndex(sub, i-1, j+1, k);
  indx->lft_bot = SubvectorEltIndex(sub, i-1, j-1, k);

  indx->rgt_frt = SubvectorEltIndex(sub, i+1, j, k+1);
  indx->rgt_bck = SubvectorEltIndex(sub, i+1, j, k-1);
  indx->lft_frt = SubvectorEltIndex(sub, i-1, j, k+1);
  indx->lft_bck = SubvectorEltIndex(sub, i-1, j, k-1);

  indx->top_frt = SubvectorEltIndex(sub, i, j+1, k+1);
  indx->top_bck = SubvectorEltIndex(sub, i, j+1, k-1);
  indx->bot_frt = SubvectorEltIndex(sub, i, j-1, k+1);
  indx->bot_bck = SubvectorEltIndex(sub, i, j-1, k-1);

  return;
}

void SubvectorStencilIndxRgt(StencilIndx *indx, Subvector *sub,
                             Subgrid *subgrid, int i, int j, int k)
{
  double u = RealSpaceX(i, SubgridRX(subgrid));
  double v = RealSpaceY(j, SubgridRY(subgrid));
  double w = RealSpaceZ(k, SubgridRZ(subgrid));
  indx->pos = v3_init(u, v, w);

  indx->mid = SubvectorEltIndex(sub, i, j, k);
  indx->rgt = SubvectorEltIndex(sub, i+1, j, k);
  indx->top = SubvectorEltIndex(sub, i, j+1, k);
  indx->bot = SubvectorEltIndex(sub, i, j-1, k);
  indx->frt = SubvectorEltIndex(sub, i, j, k+1);
  indx->bck = SubvectorEltIndex(sub, i, j, k-1);

  indx->rgt_top = SubvectorEltIndex(sub, i+1, j+1, k);
  indx->rgt_bot = SubvectorEltIndex(sub, i+1, j-1, k);
  indx->rgt_frt = SubvectorEltIndex(sub, i+1, j, k+1);
  indx->rgt_bck = SubvectorEltIndex(sub, i+1, j, k-1);

  return;
}

void SubvectorStencilIndxTop(StencilIndx *indx, Subvector *sub,
                             Subgrid *subgrid, int i, int j, int k)
{
  double u = RealSpaceX(i, SubgridRX(subgrid));
  double v = RealSpaceY(j, SubgridRY(subgrid));
  double w = RealSpaceZ(k, SubgridRZ(subgrid));
  indx->pos = v3_init(u, v, w);

  indx->mid = SubvectorEltIndex(sub, i, j, k);
  indx->rgt = SubvectorEltIndex(sub, i+1, j, k);
  indx->lft = SubvectorEltIndex(sub, i-1, j, k);
  indx->top = SubvectorEltIndex(sub, i, j+1, k);
  indx->frt = SubvectorEltIndex(sub, i, j, k+1);
  indx->bck = SubvectorEltIndex(sub, i, j, k-1);

  indx->rgt_top = SubvectorEltIndex(sub, i+1, j+1, k);
  indx->lft_top = SubvectorEltIndex(sub, i-1, j+1, k);
  indx->top_frt = SubvectorEltIndex(sub, i, j+1, k+1);
  indx->top_bck = SubvectorEltIndex(sub, i, j+1, k-1);

  return;
}

void SubvectorStencilIndxFrt(StencilIndx *indx, Subvector *sub, 
                             Subgrid *subgrid, int i, int j, int k)
{
  double u = RealSpaceX(i, SubgridRX(subgrid));
  double v = RealSpaceY(j, SubgridRY(subgrid));
  double w = RealSpaceZ(k, SubgridRZ(subgrid));
  indx->pos = v3_init(u, v, w);

  indx->mid = SubvectorEltIndex(sub, i, j, k);
  indx->rgt = SubvectorEltIndex(sub, i+1, j, k);
  indx->lft = SubvectorEltIndex(sub, i-1, j, k);
  indx->top = SubvectorEltIndex(sub, i, j+1, k);
  indx->bot = SubvectorEltIndex(sub, i, j-1, k);
  indx->frt = SubvectorEltIndex(sub, i, j, k+1);

  indx->rgt_frt = SubvectorEltIndex(sub, i+1, j, k+1);
  indx->lft_frt = SubvectorEltIndex(sub, i-1, j, k+1);
  indx->top_frt = SubvectorEltIndex(sub, i, j+1, k+1);
  indx->bot_frt = SubvectorEltIndex(sub, i, j-1, k+1);

  return;
}



double RichardsDivergenceRgt(StencilIndx *S, Morphism *my_morphism, double *pp, double *dp, double *rpp, double *permxp, double *permyp, double *permzp, double gravity, double viscosity, double du, double dv, double dw)
{
  double u = S->pos.u;
  double v = S->pos.v;
  double w = S->pos.w;

  double idu = 1.0 / du;
  double idv = 1.0 / dv;
  double idw = 1.0 / dw;

  v3 zeta_rgt = v3_init(u + 0.5 * du, v, w);
  v3basis basis_cov_rgt = MorphismDel(my_morphism, zeta_rgt);

  v3 grad_z = v3_init(0, 0, 1);
  v3 grad_z_cov_rgt = MorphismToCovariant(grad_z, basis_cov_rgt);

  v3 grad_p_cov_rgt = 
    v3_init(
      idu * (pp[S->rgt] - pp[S->mid]),
      0.25 * idv * (pp[S->rgt_top] - pp[S->rgt_bot] + pp[S->top] - pp[S->bot]),
      0.25 * idw * (pp[S->rgt_frt] - pp[S->rgt_bck] + pp[S->frt] - pp[S->bck])
    );

  double upstreams_rgt = gravity * RPMean(pp[S->mid], pp[S->rgt], 
                                          dp[S->mid] * rpp[S->mid], 
                                          dp[S->rgt] * rpp[S->rgt]) 
                          / viscosity;
  double weight_density_rgt = gravity * (dp[S->rgt] - dp[S->mid]) * idu;
  v3 grad_h_cov_rgt = v3_subtract(v3_scale(upstreams_rgt, grad_p_cov_rgt),
                                  v3_scale(weight_density_rgt, grad_z_cov_rgt));

  v3basis basis_cov_rgt_onequarter = 
      MorphismDel(my_morphism, v3_init(u+0.25*du, v, w));
  v3basis basis_cov_rgt_threequarter = 
      MorphismDel(my_morphism, v3_init(u+0.75*du, v, w));

  m3 iK_mid = m3_diagonal(1.0 / permxp[S->mid],
                          1.0 / permyp[S->mid],
                          1.0 / permzp[S->mid]);
  m3 iK_rgt = m3_diagonal(1.0 / permxp[S->rgt],
                          1.0 / permyp[S->rgt],
                          1.0 / permzp[S->rgt]);

  m3 K_dot_basis_con_rgt = m3_factor(2, m3_inverse(m3_add(
    m3v3basis_left_dot(basis_cov_rgt_onequarter, iK_mid),
    m3v3basis_left_dot(basis_cov_rgt_threequarter, iK_mid)
  )));

  v3basis basis_con_rgt = MorphismDelInverse(my_morphism, zeta_rgt);

  v3 K_con_rgt = m3v3_left_dot(basis_con_rgt.u, K_dot_basis_con_rgt);

  double q_con_rgt = -1 * v3_dot(K_con_rgt, grad_h_cov_rgt);

  v3 zeta_mid = S->pos;
  double jac_mid = MorphismJacobian(my_morphism, zeta_mid);
  double jac_rgt = MorphismJacobian(my_morphism, zeta_rgt);

  double u_rgt = dv * dw * jac_rgt * q_con_rgt / jac_mid;
  return u_rgt;
}


double RichardsDivergenceTop(StencilIndx *S, Morphism *my_morphism, double *pp, double *dp, double *rpp, double *permxp, double *permyp, double *permzp, double gravity, double viscosity, double du, double dv, double dw)
{
  double u = S->pos.u;
  double v = S->pos.v;
  double w = S->pos.w;

  double idu = 1.0 / du;
  double idv = 1.0 / dv;
  double idw = 1.0 / dw;

  v3 zeta_top = v3_init(u, v + 0.5 * dv, w);
  v3basis basis_cov_top = MorphismDel(my_morphism, zeta_top);

  v3 grad_z = v3_init(0, 0, 1);
  v3 grad_z_cov_top = MorphismToCovariant(grad_z, basis_cov_top);

  v3 grad_p_cov_top = 
    v3_init(
      0.25 * idu * (pp[S->rgt_top] - pp[S->lft_top] + pp[S->rgt] - pp[S->lft]),
      idv * (pp[S->top] - pp[S->mid]),
      0.25 * idw * (pp[S->top_frt] - pp[S->top_bck] + pp[S->frt] - pp[S->bck])
    );

  double upstreams_top = gravity * RPMean(pp[S->mid], pp[S->top], 
                                          dp[S->mid] * rpp[S->mid], 
                                          dp[S->top] * rpp[S->top]) 
                          / viscosity;
  double weight_density_top = gravity * (dp[S->top] - dp[S->mid]) * idv;
  v3 grad_h_cov_top = v3_subtract(v3_scale(upstreams_top, grad_p_cov_top),
                                  v3_scale(weight_density_top, grad_z_cov_top));

  v3basis basis_cov_top_onequarter = 
      MorphismDel(my_morphism, v3_init(u, v+0.25*dv, w));
  v3basis basis_cov_top_threequarter = 
      MorphismDel(my_morphism, v3_init(u, v+0.75*dv, w));

  m3 iK_mid = m3_diagonal(1.0 / permxp[S->mid],
                          1.0 / permyp[S->mid],
                          1.0 / permzp[S->mid]);
  m3 iK_top = m3_diagonal(1.0 / permxp[S->top],
                          1.0 / permyp[S->top],
                          1.0 / permzp[S->top]);

  m3 K_dot_basis_con_top = m3_factor(2, m3_inverse(m3_add(
    m3v3basis_left_dot(basis_cov_top_onequarter, iK_mid),
    m3v3basis_left_dot(basis_cov_top_threequarter, iK_mid)
  )));

  v3basis basis_con_top = MorphismDelInverse(my_morphism, zeta_top);

  v3 K_con_top = m3v3_left_dot(basis_con_top.v, K_dot_basis_con_top);

  double q_con_top = -1 * v3_dot(K_con_top, grad_h_cov_top);

  v3 zeta_mid = S->pos;
  double jac_mid = MorphismJacobian(my_morphism, zeta_mid);
  double jac_top = MorphismJacobian(my_morphism, zeta_top);

  double u_top = dw * du * jac_top * q_con_top / jac_mid;
  return u_top;
}


double RichardsDivergenceFrt(StencilIndx *S, Morphism *my_morphism, double *pp, double *dp, double *rpp, double *permxp, double *permyp, double *permzp, double gravity, double viscosity, double du, double dv, double dw)
{
  double u = S->pos.u;
  double v = S->pos.v;
  double w = S->pos.w;

  double idu = 1.0 / du;
  double idv = 1.0 / dv;
  double idw = 1.0 / dw;

  v3 zeta_frt = v3_init(u, v, w + 0.5 * dw);
  v3basis basis_cov_frt = MorphismDel(my_morphism, zeta_frt);

  v3 grad_z = v3_init(0, 0, 1);
  v3 grad_z_cov_frt = MorphismToCovariant(grad_z, basis_cov_frt);

  v3 grad_p_cov_frt = 
    v3_init(
      0.25 * idu * (pp[S->rgt_frt] - pp[S->lft_frt] + pp[S->rgt] - pp[S->lft]),
      0.25 * idv * (pp[S->top_frt] - pp[S->bot_frt] + pp[S->top] - pp[S->bot]),
      idw * (pp[S->frt] - pp[S->mid])
    );

  double upstreams_frt = gravity * RPMean(pp[S->mid], pp[S->frt], 
                                          dp[S->mid] * rpp[S->mid], 
                                          dp[S->frt] * rpp[S->frt]) 
                          / viscosity;
  double weight_density_frt = gravity * (dp[S->frt] - dp[S->mid]) * idw;
  v3 grad_h_cov_frt = v3_subtract(v3_scale(upstreams_frt, grad_p_cov_frt),
                                  v3_scale(weight_density_frt, grad_z_cov_frt));

  v3basis basis_cov_frt_onequarter = 
      MorphismDel(my_morphism, v3_init(u, v, w+0.25*dw));
  v3basis basis_cov_frt_threequarter = 
      MorphismDel(my_morphism, v3_init(u, v, w+0.75*dw));

  m3 iK_mid = m3_diagonal(1.0 / permxp[S->mid],
                          1.0 / permyp[S->mid],
                          1.0 / permzp[S->mid]);
  m3 iK_frt = m3_diagonal(1.0 / permxp[S->frt],
                          1.0 / permyp[S->frt],
                          1.0 / permzp[S->frt]);

  m3 K_dot_basis_con_frt = m3_factor(2, m3_inverse(m3_add(
    m3v3basis_left_dot(basis_cov_frt_onequarter, iK_mid),
    m3v3basis_left_dot(basis_cov_frt_threequarter, iK_mid)
  )));

  v3basis basis_con_frt = MorphismDelInverse(my_morphism, zeta_frt);

  v3 K_con_frt = m3v3_left_dot(basis_con_frt.w, K_dot_basis_con_frt);

  double q_con_frt = -1 * v3_dot(K_con_frt, grad_h_cov_frt);

  v3 zeta_mid = S->pos;
  double jac_mid = MorphismJacobian(my_morphism, zeta_mid);
  double jac_frt = MorphismJacobian(my_morphism, zeta_frt);

  double u_frt = du * dv * jac_frt * q_con_frt / jac_mid;
  return u_frt;
}
