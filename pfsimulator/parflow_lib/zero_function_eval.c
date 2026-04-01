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
#include "llnlmath.h"
#include "llnltyps.h"
//#include "math.h"
#include "float.h"
#include "richards_div_eval.h"
// #include "coordinate_transform.h"
// #include "compute_permeability_tensor.h"

#define NegUFace LeftFace
#define PosUFace RightFace
#define NegVFace DownFace
#define PosVFace UpFace
#define NegWFace BackFace
#define PosWFace FrontFace

/*---------------------------------------------------------------------
 * Define module structures
 *---------------------------------------------------------------------*/

typedef struct {
  int time_index;
  double SpinupDampP1;      // NBE
  double SpinupDampP2;      // NBE
  int tfgupwind;           //@RMM added for TFG formulation switch
} PublicXtra;

typedef struct {
  Problem      *problem;

  PFModule     *density_module;
  PFModule     *saturation_module;
  PFModule     *rel_perm_module;
  PFModule     *phase_source;
  PFModule     *bc_pressure;
  PFModule     *bc_internal;
  PFModule     *overlandflow_module;  //DOK
  PFModule     *overlandflow_module_diff;  //@RMM
  PFModule     *overlandflow_module_kin;
} InstanceXtra;

/*---------------------------------------------------------------------
 * Define macros for function evaluation
 *---------------------------------------------------------------------*/

#define PMean(a, b, c, d)    HarmonicMean(c, d)
#define PMeanDZ(a, b, c, d)     HarmonicMeanDZ(a, b, c, d)
#define RPMean(a, b, c, d)   UpstreamMean(a, b, c, d)
#define Mean(a, b)            ArithmeticMean(a, b)

/*  This routine provides the interface between KINSOL and ParFlow
 *  for function evaluations.  */

void     KINSolFunctionEval(
                            int      size,
                            N_Vector pressure,
                            N_Vector fval,
                            void *   current_state)
{
  PFModule  *nl_function_eval = StateFunc(((State*)current_state));
  ProblemData *problem_data = StateProblemData(((State*)current_state));
  Vector      *old_pressure = StateOldPressure(((State*)current_state));
  Vector      *saturation = StateSaturation(((State*)current_state));
  Vector      *old_saturation = StateOldSaturation(((State*)current_state));
  Vector      *density = StateDensity(((State*)current_state));
  Vector      *old_density = StateOldDensity(((State*)current_state));
  double dt = StateDt(((State*)current_state));
  double time = StateTime(((State*)current_state));
  Vector       *evap_trans = StateEvapTrans(((State*)current_state));
  Vector       *ovrl_bc_flx = StateOvrlBcFlx(((State*)current_state));

  /* velocity vectors jjb */
  Vector       *x_velocity = StateXvel(((State*)current_state));
  Vector       *y_velocity = StateYvel(((State*)current_state));
  Vector       *z_velocity = StateZvel(((State*)current_state));

  (void)size;

  PFModuleInvokeType(NlFunctionEvalInvoke, nl_function_eval,
                     (pressure, fval, problem_data, saturation, old_saturation,
                      density, old_density, dt, time, old_pressure, evap_trans,
                      ovrl_bc_flx, x_velocity, y_velocity, z_velocity));

  return;
}


/*  This routine evaluates the nonlinear function based on the current
 *  pressure values.  This evaluation is basically an application
 *  of the stencil to the pressure array. */

void NlFunctionEval(Vector *     pressure, /* Current pressure values */
                    Vector *     fval, /* Return values of the nonlinear function */
                    ProblemData *problem_data,  /* Geometry data for problem */
                    Vector *     saturation, /* Saturation / work vector */
                    Vector *     old_saturation, /* Saturation values at previous time step */
                    Vector *     density, /* Density vector */
                    Vector *     old_density, /* Density values at previous time step */
                    double       dt, /* Time step size */
                    double       time, /* New time value */
                    Vector *     old_pressure,
                    Vector *     evap_trans, /*sk sink term from land surface model*/
                    Vector *     ovrl_bc_flx, /*sk overland flow boundary fluxes*/
                    Vector *     u_velocity, /* velocity vectors jjb */
                    Vector *     v_velocity,
                    Vector *     w_velocity)
{
  PUSH_NVTX("NlFunctionEval", 0)

  PFModule      *this_module = ThisPFModule;
  InstanceXtra  *instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);
  PublicXtra    *public_xtra = (PublicXtra*)PFModulePublicXtra(this_module);

  Problem     *problem = (instance_xtra->problem);

  PFModule    *density_module = (instance_xtra->density_module);
  PFModule    *saturation_module = (instance_xtra->saturation_module);
  PFModule    *rel_perm_module = (instance_xtra->rel_perm_module);
  PFModule    *phase_source = (instance_xtra->phase_source);
  PFModule    *bc_pressure = (instance_xtra->bc_pressure);
  PFModule    *bc_internal = (instance_xtra->bc_internal);

  /* Reuse saturation vector to save memory */
  Vector      *rel_perm = saturation;
  Vector      *source = saturation;

  Vector      *porosity = ProblemDataPorosity(problem_data);
  Vector      *permeability_x = ProblemDataPermeabilityX(problem_data);
  Vector      *permeability_y = ProblemDataPermeabilityY(problem_data);
  Vector      *permeability_z = ProblemDataPermeabilityZ(problem_data);
  Vector      *sstorage = ProblemDataSpecificStorage(problem_data);
  Vector      *top = ProblemDataIndexOfDomainTop(problem_data);
  Vector      *bottom = ProblemDataIndexOfDomainBottom(problem_data);
  
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

  Vector *lengthA_u = ProblemDataLengthUA(problem_data);
  Vector *lengthB_u = ProblemDataLengthUB(problem_data);
  Vector *lengthA_v = ProblemDataLengthVA(problem_data);
  Vector *lengthB_v = ProblemDataLengthVB(problem_data);
  Vector *lengthA_w = ProblemDataLengthWA(problem_data);
  Vector *lengthB_w = ProblemDataLengthWB(problem_data);

  Vector *z_coordinate = ProblemDataZCoordinate(problem_data);

  double gravity = ProblemGravity(problem);
  double viscosity = ProblemPhaseViscosity(problem, 0);

  Grid        *grid = VectorGrid(pressure);
  Grid        *grid2d = VectorGrid(top);

  // Velocity data jjb
  int nx_vy, sy_v;
  int nx_vz, ny_vz, sz_v;

  GrGeomSolid *gr_domain = ProblemDataGrDomain(problem_data);

  /* Get Data Vectors */
  Vector      *jacobian = ProblemDataJacobian(problem_data);

  BeginTiming(public_xtra->time_index);
  /* Initialize function values to zero. */
  InitVectorAll(fval, 0.0);

  /* Pass pressure values to neighbors.  */
  VectorUpdateCommHandle *handle = NULL;
  handle = InitVectorUpdate(pressure, VectorUpdateAll);
  FinalizeVectorUpdate(handle);

  /* Calculate pressure dependent properties: density and saturation */

  double dtmp = 0;
  PFModuleInvokeType(PhaseDensityInvoke, density_module, (0, pressure, density, &dtmp, &dtmp,
                                                          CALCFCN));

  PFModuleInvokeType(SaturationInvoke, saturation_module, (saturation, pressure, density,
                                                           gravity, problem_data, CALCFCN));


  /* Calculate accumulation terms for the function values */

  int is = 0;
  ForSubgridI(is, GridSubgrids(grid))
  {
    Subgrid *subgrid = GridSubgrid(grid, is);

    /* RDF: assumes resolutions are the same in all 3 directions */
    int r = SubgridRX(subgrid);

    int iu = SubgridIX(subgrid);
    int iv = SubgridIY(subgrid);
    int iw = SubgridIZ(subgrid);

    int nu = SubgridNX(subgrid);
    int nv = SubgridNY(subgrid);
    int nw = SubgridNZ(subgrid);

    double du = SubgridDX(subgrid);
    double dv = SubgridDY(subgrid);
    double dw = SubgridDZ(subgrid);

    double dudvdw = du * dv * dw;

    Subvector *fval_sub = VectorSubvector(fval, is);
    double *fval_dat = SubvectorData(fval_sub);

    /* Parameters from coordinate transform */
    double *jacobian_dat = SubvectorData(VectorSubvector(jacobian, is));

    /* Parameters to compute Accumulation Terms */
    double *density_dat = SubvectorData(VectorSubvector(density, is));
    double *old_density_dat = SubvectorData(VectorSubvector(old_density, is));
    double *saturation_dat = SubvectorData(VectorSubvector(saturation, is));
    double *old_saturation_dat = SubvectorData(VectorSubvector(old_saturation, is));
    double *porosity_dat = SubvectorData(VectorSubvector(porosity, is));

    /* Parameters to compute Compressive Storage Terms */
    double *sstor_dat = SubvectorData(VectorSubvector(sstorage, is));
    double *pressure_dat = SubvectorData(VectorSubvector(pressure, is));
    double *old_pressure_dat = SubvectorData(VectorSubvector(old_pressure, is));

    int i = 0, j = 0, k = 0;
    GrGeomInLoop(i, j, k, gr_domain, r, iu, iv, iw, nu, nv, nw,
    {
      int idx = SubvectorEltIndex(fval_sub, i, j, k);

      double jac_ = jacobian_dat[idx];
      double satur_ = saturation_dat[idx];
      double osatur_ = old_saturation_dat[idx];
      double dens_ = density_dat[idx];
      double odens_ = old_density_dat[idx];
      double poros_ = porosity_dat[idx];
      
      double sstor_ = sstor_dat[idx];
      double press_ = pressure_dat[idx];
      double opress_ = old_pressure_dat[idx];

      /* Accuumulation Terms */
      fval_dat[idx] = dudvdw * jac_ * (satur_ * dens_ - osatur_ * odens_) * poros_;
      
      /* Compressive Storage Terms */
      fval_dat[idx] += dudvdw * jac_ * sstor_ * (press_ * satur_ * dens_ - opress_ * osatur_ * odens_);
    });
  }

  /* Add in contributions from source terms - user specified sources and
   * flux wells.  Calculate phase source values overwriting current
   * saturation vector */
  PFModuleInvokeType(PhaseSourceInvoke, phase_source, (source, 0, problem, problem_data,
                                                       time));

  is = 0;
  ForSubgridI(is, GridSubgrids(grid))
  {
    Subgrid *subgrid = GridSubgrid(grid, is);

    /* RDF: assumes resolutions are the same in all 3 directions */
    int r = SubgridRX(subgrid);

    int iu = SubgridIX(subgrid);
    int iv = SubgridIY(subgrid);
    int iw = SubgridIZ(subgrid);

    int nu = SubgridNX(subgrid);
    int nv = SubgridNY(subgrid);
    int nw = SubgridNZ(subgrid);

    double du = SubgridDX(subgrid);
    double dv = SubgridDY(subgrid);
    double dw = SubgridDZ(subgrid);

    double dudvdw = du * dv * dw;

    Subvector *fval_sub = VectorSubvector(fval, is);
    double *fval_dat = SubvectorData(fval_sub);

    double *jacobian_dat = SubvectorData(VectorSubvector(jacobian, is));

    double *source_dat = SubvectorData(VectorSubvector(source, is));
    double *evap_trans_dat = SubvectorData(VectorSubvector(evap_trans, is));

    int i = 0, j = 0, k = 0;
    GrGeomInLoop(i, j, k, gr_domain, r, iu, iv, iw, nu, nv, nw,
    {
      int idx = SubvectorEltIndex(fval_sub, i, j, k);

      double jac_ = jacobian_dat[idx];
      double src_ = source_dat[idx];
      double evap_trans_ = evap_trans_dat[idx];

      fval_dat[idx] -= dudvdw * jac_ * dt * (src_ + evap_trans_);
    });
  }

  BCStruct *bc_struct = PFModuleInvokeType(BCPressureInvoke, bc_pressure,
                                 (problem_data, grid, gr_domain, time));

  /*
   * Temporarily insert boundary pressure values for Dirichlet
   * boundaries into cells that are in the inactive region but next
   * to a Dirichlet boundary condition.  These values are required
   * for use in the rel_perm_module to compute rel_perm values for
   * these cells. They needed for upstream weighting in mobilities.
   *
   * NOTES:
   *
   * These values must be later removed from the pressure field and
   * fval needs to be adjusted for these cells to make the inactive
   * region problem decoupled from the active region cells for the
   * solver.
   *
   * Densities are currently defined everywhere so should be valid for
   * these boundary cells.
   *
   * SGS not sure if this will work or not so left it here for later
   * exploration.  This is a little hacky in the sense that we are
   * inserting values and then need to overwrite them again.  It
   * might be more clean to rewrite the Dirichlet boundary condition
   * code to not require the values be in the pressure field for
   * these cells but instead grab the values out of the
   * BCStructPatchValues as was done here.  In other words use
   * bc_patch_values[ival] in rel_perm_module code and remove this
   * loop.
   */

  is = 0;
  ForSubgridI(is, GridSubgrids(grid))
  {
    Subgrid *subgrid = GridSubgrid(grid, is);

    Subvector *pressure_sub = VectorSubvector(pressure, is);
    double *pressure_dat = SubvectorData(pressure_sub);

    int stride_u = 1;
    int stride_v = SubvectorNX(pressure_sub);
    int stride_w = SubvectorNY(pressure_sub) * stride_v;

    int ipatch = 0;
    ForBCStructNumPatches(ipatch, bc_struct)
    {
      double *bc_patch_values = BCStructPatchValues(bc_struct, ipatch, is);

      int icounter = 0, i = 0, j = 0, k = 0;
      ForPatchCellsPerFace(DirichletBC,
                           BeforeAllCells(DoNothing),
                           LoopVars(i, j, k, icounter, bc_struct, ipatch, is),
                           Locals(int press_idx, idx; double value; ),
                           CellSetup({
        press_idx = 0;
        idx = SubvectorEltIndex(pressure_sub, i, j, k);
        value = bc_patch_values[icounter];
      }),
                           FACE(NegUFace, { press_idx = idx - stride_u; }),
                           FACE(PosUFace, { press_idx = idx + stride_u; }),
                           FACE(NegVFace, { press_idx = idx - stride_v; }),
                           FACE(PosVFace, { press_idx = idx + stride_v; }),
                           FACE(NegWFace, { press_idx = idx - stride_w; }),
                           FACE(PosWFace, { press_idx = idx + stride_w; }),
                           CellFinalize({ pressure_dat[press_idx] = value; }),
                           AfterAllCells(DoNothing)
                           ); /* End DirichletBC */
    }          /* End ipatch loop */
  }            /* End subgrid loop */

  /* Calculate relative permeability values overwriting current
   * phase source values */

  PFModuleInvokeType(PhaseRelPermInvoke, rel_perm_module,
                     (rel_perm, pressure, density, gravity, problem_data,
                      CALCFCN));

  /* Calculate contributions from second order derivatives and gravity */
  is = 0;
  ForSubgridI(is, GridSubgrids(grid))
  {
    Subgrid *subgrid = GridSubgrid(grid, is);

    /* RDF: assumes resolutions are the same in all 3 directions */
    int r = SubgridRX(subgrid);

    int iu = SubgridIX(subgrid);
    int iv = SubgridIY(subgrid);
    int iw = SubgridIZ(subgrid);

    int nu = SubgridNX(subgrid);
    int nv = SubgridNY(subgrid);
    int nw = SubgridNZ(subgrid);

    double du = SubgridDX(subgrid);
    double dv = SubgridDY(subgrid);
    double dw = SubgridDZ(subgrid);

    double dvdw = dv * dw;
    double dwdu = du * dw;
    double dudv = du * dv;

    Subvector *fval_sub = VectorSubvector(fval, is);
    double *fval_dat = SubvectorData(fval_sub);

    int nu_v = SubvectorNX(fval_sub);
    int nv_v = SubvectorNY(fval_sub);
    int nw_v = SubvectorNZ(fval_sub);

    int stride_u = 1;
    int stride_v = SubvectorNX(fval_sub);
    int stride_w = SubvectorNY(fval_sub) * stride_v;

    double *pressure_dat = SubvectorData(VectorSubvector(pressure, is));
    double *z_dat = SubvectorData(VectorSubvector(z_coordinate, is));
    double *density_dat = SubvectorData(VectorSubvector(density, is));
    double *rel_perm_dat = SubvectorData(VectorSubvector(rel_perm, is));
    double *jacobian_dat = SubvectorData(VectorSubvector(jacobian, is));

    double *lengthA_u_dat = SubvectorData(VectorSubvector(lengthA_u, is));
    double *lengthB_u_dat = SubvectorData(VectorSubvector(lengthB_u, is));
    double *lengthA_v_dat = SubvectorData(VectorSubvector(lengthA_v, is));
    double *lengthB_v_dat = SubvectorData(VectorSubvector(lengthB_v, is));
    double *lengthA_w_dat = SubvectorData(VectorSubvector(lengthA_w, is));
    double *lengthB_w_dat = SubvectorData(VectorSubvector(lengthB_w, is));

    double *K_uu_dat = SubvectorData(VectorSubvector(K_uu, is));
    double *K_uv_dat = SubvectorData(VectorSubvector(K_uv, is));
    double *K_uw_dat = SubvectorData(VectorSubvector(K_uw, is));
    double *K_vv_dat = SubvectorData(VectorSubvector(K_vv, is));
    double *K_vw_dat = SubvectorData(VectorSubvector(K_vw, is));
    double *K_vu_dat = SubvectorData(VectorSubvector(K_vu, is));
    double *K_ww_dat = SubvectorData(VectorSubvector(K_ww, is));
    double *K_wu_dat = SubvectorData(VectorSubvector(K_wu, is));
    double *K_wv_dat = SubvectorData(VectorSubvector(K_wv, is));

    /* velocity accessors jjb */
    Subvector *vel_u_sub = VectorSubvector(u_velocity, is);
    Subvector *vel_v_sub = VectorSubvector(v_velocity, is);
    Subvector *vel_w_sub = VectorSubvector(w_velocity, is);
    double *vel_u_dat = SubvectorData(vel_u_sub);
    double *vel_v_dat = SubvectorData(vel_v_sub);
    double *vel_w_dat = SubvectorData(vel_w_sub);

    int i = 0, j = 0, k = 0;
    int idx = SubvectorEltIndex(fval_sub, i, j, k);
    BoxLoopI1(i, j, k, iu, iv, iw, nu, nv, nw,
              idx, nu_v, nv_v, nw_v, 1, 1, 1,
    {
      double dh[3]; dh[U] = 0.0; dh[V] = 0.0; dh[W] = 0.0;

      /* velocity subvector indices jjb */
      int vu_idx = SubvectorEltIndex(vel_u_sub, i + 1, j, k);
      int vv_idx = SubvectorEltIndex(vel_v_sub, i, j + 1, k);
      int vw_idx = SubvectorEltIndex(vel_w_sub, i, j, k + 1);

      // double neg_dh_dz = (pp[ip] - pp[ip + stride_z]) / dz - 0.5 * (dp[ip] + dp[ip + stride_z]) * gravity * (z[ip + stride_z] - z[ip]) / dz;

      // double rel_perm = neg_dh_dz > 0 ? rpp[idx] * dp[idx] : rpp[idx + stride_z] * dp[idx + stride_z];

      // double u_upper = dxdy * HarmonicMeanDZ(permzp[ip], permzp[ip + sz_p])
      //                  * neg_dh_dz * rel_perm / viscosity;

      // U Face Permutation: (u', v', w') -> (u, v, w)
      double jacobian_u = 0.5 * (jacobian_dat[idx] + jacobian_dat[idx + stride_u]);
      double dArea_u = jacobian_u * dvdw;
      HydraulicHeadDel(idx, stride_u, stride_v, stride_w, du, dv, dw, dh,
                       lengthA_u_dat, lengthB_u_dat, pressure_dat, z_dat,
                       density_dat, gravity);
      double u_posu = FluxAtFace(idx, stride_u, dArea_u, dh,
                                 density_dat, rel_perm_dat, viscosity,
                                 K_uu_dat, K_uv_dat, K_uw_dat);

      // V Face Permutation: (u', v', w') -> (v, w, u)
      double jacobian_v = 0.5 * (jacobian_dat[idx] + jacobian_dat[idx + stride_v]);
      double dArea_v = jacobian_v * dwdu;
      HydraulicHeadDel(idx, stride_v, stride_w, stride_u, dv, dw, du, dh,
                       lengthA_v_dat, lengthB_v_dat, pressure_dat, z_dat,
                       density_dat, gravity);
      double u_posv = FluxAtFace(idx, stride_v, dArea_v, dh,
                                 density_dat, rel_perm_dat, viscosity,
                                 K_vv_dat, K_vw_dat, K_vu_dat);

      // W Face Permutation: (u', v', w') -> (w, u, v)
      double jacobian_w = 0.5 * (jacobian_dat[idx] + jacobian_dat[idx + stride_w]);
      double dArea_w = jacobian_w * dudv;
      HydraulicHeadDel(idx, stride_w, stride_u, stride_v, dw, du, dv, dh,
                       lengthA_w_dat, lengthB_w_dat, pressure_dat, z_dat,
                       density_dat, gravity);
      double u_posw = FluxAtFace(idx, stride_w, dArea_w, dh,
                                 density_dat, rel_perm_dat, viscosity,
                                 K_ww_dat, K_wu_dat, K_wv_dat);

      /* velocity data jjb */
      vel_u_dat[vu_idx] = u_posu / dArea_u;
      vel_v_dat[vv_idx] = u_posv / dArea_v;
      vel_w_dat[vw_idx] = u_posw / dArea_w;

      PlusEquals(fval_dat[idx], dt * (u_posu + u_posv + u_posw));
      PlusEquals(fval_dat[idx + stride_u], -dt * u_posu);
      PlusEquals(fval_dat[idx + stride_v], -dt * u_posv);
      PlusEquals(fval_dat[idx + stride_w], -dt * u_posw);
    });

    /* We need to add the flux values coming from the ghost nodes
     * on the negative sides. This has to be done in a separate loop
     * because of the change in stencil, which would access memory
     * out-of-bounds in the ghost nodes.
     *
     * For example, in the nodes in the NegU side, the HydraulicHeadDel
     * for faces V and W will access values outside the grid. */

    /* Ghost nodes on NegU side: */
    i = -1, j = 0, k = 0;
    idx = SubvectorEltIndex(fval_sub, i, j, k);
    BoxLoopI1(i, j, k, (iu-1), iv, iw, 1, nv, nw,
              idx, nu_v, nv_v, nw_v, 1, 1, 1,
    {
      double dh[3]; dh[U] = 0.0; dh[V] = 0.0; dh[W] = 0.0;

      /* velocity subvector indices jjb */
      int vu_idx = SubvectorEltIndex(vel_u_sub, i + 1, j, k);

      // U Face Permutation: (u', v', w') -> (u, v, w)
      double jacobian_u = 0.5 * (jacobian_dat[idx] + jacobian_dat[idx + stride_u]);
      double dArea_u = jacobian_u * dvdw;
      HydraulicHeadDel(idx, stride_u, stride_v, stride_w, du, dv, dw, dh,
                       lengthA_u_dat, lengthB_u_dat, pressure_dat, z_dat,
                       density_dat, gravity);
      double u_posu = FluxAtFace(idx, stride_u, dArea_u, dh,
                                 density_dat, rel_perm_dat, viscosity,
                                 K_uu_dat, K_uv_dat, K_uw_dat);

      /* velocity data jjb */
      vel_u_dat[vu_idx] = u_posu / dArea_u;

      PlusEquals(fval_dat[idx + stride_u], -dt * u_posu);
    });

    /* Ghost nodes on NegV side: */
    i = 0, j = -1, k = 0;
    idx = SubvectorEltIndex(fval_sub, i, j, k);
    BoxLoopI1(i, j, k, iu, (iv-1), iw, nu, 1, nw,
              idx, nu_v, nv_v, nw_v, 1, 1, 1,
    {
      double dh[3]; dh[U] = 0.0; dh[V] = 0.0; dh[W] = 0.0;

      /* velocity subvector indices jjb */
      int vv_idx = SubvectorEltIndex(vel_v_sub, i, j + 1, k);

      // V Face Permutation: (u', v', w') -> (v, w, u)
      double jacobian_v = 0.5 * (jacobian_dat[idx] + jacobian_dat[idx + stride_v]);
      double dArea_v = jacobian_v * dwdu;
      HydraulicHeadDel(idx, stride_v, stride_w, stride_u, dv, dw, du, dh,
                       lengthA_v_dat, lengthB_v_dat, pressure_dat, z_dat,
                       density_dat, gravity);
      double u_posv = FluxAtFace(idx, stride_v, dArea_v, dh,
                                 density_dat, rel_perm_dat, viscosity,
                                 K_vv_dat, K_vw_dat, K_vu_dat);

      /* velocity data jjb */
      vel_v_dat[vv_idx] = u_posv / dArea_v;

      PlusEquals(fval_dat[idx + stride_v], -dt * u_posv);
    });

    /* Ghost nodes on NegW side: */
    i = 0, j = 0, k = -1;
    idx = SubvectorEltIndex(fval_sub, i, j, k);
    BoxLoopI1(i, j, k, iu, iv, (iw-1), nu, nv, 1,
              idx, nu_v, nv_v, nw_v, 1, 1, 1,
    {
      double dh[3]; dh[U] = 0.0; dh[V] = 0.0; dh[W] = 0.0;

      /* velocity subvector indices jjb */
      int vw_idx = SubvectorEltIndex(vel_w_sub, i, j, k + 1);

      // W Face Permutation: (u', v', w') -> (w, u, v)
      double jacobian_w = 0.5 * (jacobian_dat[idx] + jacobian_dat[idx + stride_w]);
      double dArea_w = jacobian_w * dudv;
      HydraulicHeadDel(idx, stride_w, stride_u, stride_v, dw, du, dv, dh,
                       lengthA_w_dat, lengthB_w_dat, pressure_dat, z_dat,
                       density_dat, gravity);
      double u_posw = FluxAtFace(idx, stride_w, dArea_w, dh,
                                 density_dat, rel_perm_dat, viscosity,
                                 K_ww_dat, K_wu_dat, K_wv_dat);

      /* velocity data jjb */
      vel_w_dat[vw_idx] = u_posw / dArea_w;

      PlusEquals(fval_dat[idx + stride_w], -dt * u_posw);
    });
  }

  /* Calculate correction for boundary conditions */
  is = 0;
  ForSubgridI(is, GridSubgrids(grid))
  {
    Subgrid *subgrid = GridSubgrid(grid, is);

    double du = SubgridDX(subgrid);
    double dv = SubgridDY(subgrid);
    double dw = SubgridDZ(subgrid);

    double dvdw = dv * dw;
    double dwdu = du * dw;
    double dudv = du * dv;

    Subvector *fval_sub = VectorSubvector(fval, is);
    double *fval_dat = SubvectorData(fval_sub);

    int stride_u = 1;
    int stride_v = SubvectorNX(fval_sub);
    int stride_w = SubvectorNY(fval_sub) * stride_v;

    double *pressure_dat = SubvectorData(VectorSubvector(pressure, is));
    double *z_dat = SubvectorData(VectorSubvector(z_coordinate, is));
    double *density_dat = SubvectorData(VectorSubvector(density, is));
    double *rel_perm_dat = SubvectorData(VectorSubvector(rel_perm, is));
    double *jacobian_dat = SubvectorData(VectorSubvector(jacobian, is));

    Subvector *top_sub = VectorSubvector(top, is);
    double *top_dat = SubvectorData(top_sub);
    double *bot_dat = SubvectorData(VectorSubvector(bottom, is));

    double *lengthA_u_dat = SubvectorData(VectorSubvector(lengthA_u, is));
    double *lengthB_u_dat = SubvectorData(VectorSubvector(lengthB_u, is));
    double *lengthA_v_dat = SubvectorData(VectorSubvector(lengthA_v, is));
    double *lengthB_v_dat = SubvectorData(VectorSubvector(lengthB_v, is));
    double *lengthA_w_dat = SubvectorData(VectorSubvector(lengthA_w, is));
    double *lengthB_w_dat = SubvectorData(VectorSubvector(lengthB_w, is));

    double *K_uu_dat = SubvectorData(VectorSubvector(K_uu, is));
    double *K_uv_dat = SubvectorData(VectorSubvector(K_uv, is));
    double *K_uw_dat = SubvectorData(VectorSubvector(K_uw, is));
    double *K_vv_dat = SubvectorData(VectorSubvector(K_vv, is));
    double *K_vw_dat = SubvectorData(VectorSubvector(K_vw, is));
    double *K_vu_dat = SubvectorData(VectorSubvector(K_vu, is));
    double *K_ww_dat = SubvectorData(VectorSubvector(K_ww, is));
    double *K_wu_dat = SubvectorData(VectorSubvector(K_wu, is));
    double *K_wv_dat = SubvectorData(VectorSubvector(K_wv, is));

    /* velocity accessors jjb */
    Subvector *vel_u_sub = VectorSubvector(u_velocity, is);
    Subvector *vel_v_sub = VectorSubvector(v_velocity, is);
    Subvector *vel_w_sub = VectorSubvector(w_velocity, is);
    double *vel_u_dat = SubvectorData(vel_u_sub);
    double *vel_v_dat = SubvectorData(vel_v_sub);
    double *vel_w_dat = SubvectorData(vel_w_sub);


    int ipatch = 0;
    ForBCStructNumPatches(ipatch, bc_struct)
    {
      double *bc_patch_values = BCStructPatchValues(bc_struct, ipatch, is);
      int icounter = 0, i = 0, j = 0, k = 0;

      // Remove incorrect flux values from boundary cells
      /* I need to remove almsot all fluxes in each face because most
       * take values from the boundary, except the flux that is aligned
       * with the face.
       */
      ForPatchCellsPerFace(BC_ALL,
                           BeforeAllCells(DoNothing),
                           LoopVars(i, j, k, icounter, bc_struct, ipatch, is),
                           Locals(int idx, idx_2d, idx_2d_negu, idx_2d_posu;
                                  int idx_2d_negv, idx_2d_posv;
                                  int vu_negu_idx, vv_negv_idx, vw_negw_idx;
                                  int vu_posu_idx, vv_posv_idx, vw_posw_idx;
                                  int make_correction_except_negu;
                                  int make_correction_except_posu;
                                  int make_correction_except_negv;
                                  int make_correction_except_posv;
                                  int make_correction_except_negw;
                                  int make_correction_except_posw;
                                  int make_full_correction, is_edge;
                                  int is_face_negu, is_face_posu;
                                  int is_face_negv, is_face_posv;
                                  int is_face_negw, is_face_posw;
                                  int is_face_u, is_face_v, is_face_w;),
                           CellSetup({
        idx = SubvectorEltIndex(fval_sub, i, j, k);
        idx_2d = SubvectorEltIndex(top_sub, i, j, 0);
        idx_2d_negu = SubvectorEltIndex(top_sub, i - 1, j, 0);
        idx_2d_posu = SubvectorEltIndex(top_sub, i + 1, j, 0);
        idx_2d_negv = SubvectorEltIndex(top_sub, i, j - 1, 0);
        idx_2d_posv = SubvectorEltIndex(top_sub, i, j + 1, 0);

        /* velocity subvector indices jjb */
        vu_negu_idx = SubvectorEltIndex(vel_u_sub, i, j, k);
        vv_negv_idx = SubvectorEltIndex(vel_v_sub, i, j, k);
        vw_negw_idx = SubvectorEltIndex(vel_w_sub, i, j, k);
        vu_posu_idx = SubvectorEltIndex(vel_u_sub, i + 1, j, k);
        vv_posv_idx = SubvectorEltIndex(vel_v_sub, i, j + 1, k);
        vw_posw_idx = SubvectorEltIndex(vel_w_sub, i, j, k + 1);

        // since we are using coordinate transformations, we assume that the
        // grid is regular and that there are no inactive cells.
        make_correction_except_negu = FALSE;
        make_correction_except_posu = FALSE;
        make_correction_except_negv = FALSE;
        make_correction_except_posv = FALSE;
        make_correction_except_negw = FALSE;
        make_correction_except_posw = FALSE;
        make_full_correction = FALSE;
        // special case in top and bottom layers to ensure
        // that the flux is not removed twice
        is_face_negu = top_dat[idx_2d_negu] < 0;
        is_face_posu = top_dat[idx_2d_posu] < 0;
        is_face_negv = top_dat[idx_2d_negv] < 0;
        is_face_posv = top_dat[idx_2d_posv] < 0;
        is_face_negw = k == lrint(bot_dat[idx_2d]);
        is_face_posw = k == lrint(top_dat[idx_2d]);

        is_face_u = (is_face_negu || is_face_posu);
        is_face_v = (is_face_negv || is_face_posv);
        is_face_w = (is_face_negw || is_face_posw);

        is_edge = (is_face_u && is_face_v) || (is_face_u && is_face_w) || (is_face_v && is_face_w);
      }),
                           FACE(NegUFace, {
        make_correction_except_posu = !(is_face_v || is_face_w);
        make_full_correction = FALSE;
      }),
                           FACE(PosUFace, {
        make_correction_except_negu = !(is_face_v || is_face_w);
        make_full_correction = FALSE;
      }),
                           FACE(NegVFace, {
        make_correction_except_posv = !(is_face_u || is_face_w);
        make_full_correction = (is_face_u && !is_face_w);
      }),
                           FACE(PosVFace, {
        make_correction_except_negv = !(is_face_u || is_face_w);
        make_full_correction = (is_face_u && !is_face_w);
      }),
                           FACE(NegWFace, {
        make_correction_except_posw = !(is_face_u || is_face_v);
        make_full_correction = (is_face_u || is_face_v);
      }),
                           FACE(PosWFace, {
        make_correction_except_negw = !(is_face_u || is_face_v);
        make_full_correction = (is_face_u || is_face_v);
      }),
                           CellFinalize(
      {
        double dh[3]; dh[U] = 0.0; dh[V] = 0.0; dh[W] = 0.0;

        // U Face Permutation: (u', v', w') -> (u, v, w)
        double jacobian_negu = 0.5 * (jacobian_dat[idx] + jacobian_dat[idx - stride_u]);
        double dArea_negu = jacobian_negu * dvdw;
        HydraulicHeadDel((idx-stride_u), stride_u, stride_v, stride_w,
                        du, dv, dw, dh, lengthA_u_dat, lengthB_u_dat,
                        pressure_dat, z_dat, density_dat, gravity);
        double u_negu = FluxAtFace((idx-stride_u), stride_u, dArea_negu, dh,
                                  density_dat, rel_perm_dat, viscosity,
                                  K_uu_dat, K_uv_dat, K_uw_dat);

        double jacobian_posu = 0.5 * (jacobian_dat[idx] + jacobian_dat[idx + stride_u]);
        double dArea_posu = jacobian_posu * dvdw;
        HydraulicHeadDel(idx, stride_u, stride_v, stride_w, du, dv, dw, dh,
                        lengthA_u_dat, lengthB_u_dat, pressure_dat, z_dat,
                        density_dat, gravity);
        double u_posu = FluxAtFace(idx, stride_u, dArea_posu, dh,
                                  density_dat, rel_perm_dat, viscosity,
                                  K_uu_dat, K_uv_dat, K_uw_dat);

        // V Face Permutation: (u', v', w') -> (v, w, u)
        double jacobian_negv = 0.5 * (jacobian_dat[idx] + jacobian_dat[idx - stride_v]);
        double dArea_negv = jacobian_negv * dwdu;
        HydraulicHeadDel((idx-stride_v), stride_v, stride_w, stride_u,
                        dv, dw, du, dh, lengthA_v_dat, lengthB_v_dat,
                        pressure_dat, z_dat, density_dat, gravity);
        double u_negv = FluxAtFace((idx-stride_v), stride_v, dArea_negv, dh,
                                  density_dat, rel_perm_dat, viscosity,
                                  K_vv_dat, K_vw_dat, K_vu_dat);

        double jacobian_posv = 0.5 * (jacobian_dat[idx] + jacobian_dat[idx + stride_v]);
        double dArea_posv = jacobian_posv * dwdu;
        HydraulicHeadDel(idx, stride_v, stride_w, stride_u, dv, dw, du, dh,
                        lengthA_v_dat, lengthB_v_dat, pressure_dat, z_dat,
                        density_dat, gravity);
        double u_posv = FluxAtFace(idx, stride_v, dArea_posv, dh,
                                  density_dat, rel_perm_dat, viscosity,
                                  K_vv_dat, K_vw_dat, K_vu_dat);

        // W Face Permutation: (u', v', w') -> (w, u, v)
        double jacobian_negw = 0.5 * (jacobian_dat[idx] + jacobian_dat[idx - stride_w]);
        double dArea_negw = jacobian_negw * dudv;
        HydraulicHeadDel((idx-stride_w), stride_w, stride_u, stride_v,
                        dw, du, dv, dh, lengthA_w_dat, lengthB_w_dat,
                        pressure_dat, z_dat, density_dat, gravity);
        double u_negw = FluxAtFace((idx-stride_w), stride_w, dArea_negw, dh,
                                  density_dat, rel_perm_dat, viscosity,
                                  K_ww_dat, K_wu_dat, K_wv_dat);

        double jacobian_posw = 0.5 * (jacobian_dat[idx] + jacobian_dat[idx + stride_w]);
        double dArea_posw = jacobian_posw * dudv;
        HydraulicHeadDel(idx, stride_w, stride_u, stride_v, dw, du, dv, dh,
                        lengthA_w_dat, lengthB_w_dat, pressure_dat, z_dat,
                        density_dat, gravity);
        double u_posw = FluxAtFace(idx, stride_w, dArea_posw, dh,
                                  density_dat, rel_perm_dat, viscosity,
                                  K_ww_dat, K_wu_dat, K_wv_dat);

        vel_u_dat[vu_negu_idx] = (make_correction_except_negu) * u_negu / dArea_negu;
        vel_u_dat[vu_posu_idx] = (make_correction_except_posu) * u_posu / dArea_posu;
        vel_v_dat[vv_negv_idx] = (make_correction_except_negv) * u_negv / dArea_negv;
        vel_v_dat[vv_posv_idx] = (make_correction_except_posv) * u_posv / dArea_posv;
        vel_w_dat[vw_negw_idx] = (make_correction_except_negw) * u_negw / dArea_negw;
        vel_w_dat[vw_posw_idx] = (make_correction_except_posw) * u_posw / dArea_posw;

        double u_correction = - u_negu + u_posu - u_negv
                              + u_posv - u_negw + u_posw;

        u_correction += (make_correction_except_negu) * u_negu;
        u_correction -= (make_correction_except_posu) * u_posu;
        u_correction += (make_correction_except_negv) * u_negv;
        u_correction -= (make_correction_except_posv) * u_posv;
        u_correction += (make_correction_except_negw) * u_negw;
        u_correction -= (make_correction_except_posw) * u_posw;

        PlusEquals(fval_dat[idx], (!is_edge) * (- dt * u_correction));

        // if it is edge, we need to ensure that there are no redundant
        // corrections at the cells that share faces from multiple
        // boundary conditions.
        u_correction *= make_full_correction;
        PlusEquals(fval_dat[idx], - dt * u_correction);
      }),
                           AfterAllCells(DoNothing)
                           ); /* End BC_ALL */

      // Fix flux values from boundary conditions
      ForPatchCellsPerFace(FluxBC,
                           BeforeAllCells(DoNothing),
                           LoopVars(i, j, k, icounter, bc_struct, ipatch, is),
                           Locals(int idx, idx_2d, idx_2d_negu, idx_2d_posu;
                                  int idx_2d_negv, idx_2d_posv;
                                  int su, sv, sw;
                                  double *p_, *z_;
                                  double flow_magnitude;
                                  int vu_negu_idx, vv_negv_idx, vw_negw_idx;
                                  int vu_posu_idx, vv_posv_idx, vw_posw_idx;
                                  int make_correction_except_negu;
                                  int make_correction_except_posu;
                                  int make_correction_except_negv;
                                  int make_correction_except_posv;
                                  int make_correction_except_negw;
                                  int make_correction_except_posw;
                                  int is_patch_negu;
                                  int is_patch_posu;
                                  int is_patch_negv;
                                  int is_patch_posv;
                                  int is_patch_negw;
                                  int is_patch_posw;
                                  int make_full_correction, is_edge;
                                  int is_face_negu, is_face_posu;
                                  int is_face_negv, is_face_posv;
                                  int is_face_negw, is_face_posw;
                                  int is_face_u, is_face_v, is_face_w;
                                  double dh_face_negu[3], dh_face_posu[3];
                                  double dh_face_negv[3], dh_face_posv[3];
                                  double dh_face_negw[3], dh_face_posw[3];
                                  double u_fluxbc, u_correction;
                                  double jacobian_negu, jacobian_posu;
                                  double jacobian_negv, jacobian_posv;
                                  double jacobian_negw, jacobian_posw;
                                  double dArea_negu, dArea_posu;
                                  double dArea_negv, dArea_posv;
                                  double dArea_negw, dArea_posw;
                                  double density_negu, density_posu;
                                  double density_negv, density_posv;
                                  double density_negw, density_posw;
                                ),
                           CellSetup({
        idx = SubvectorEltIndex(fval_sub, i, j, k);
        idx_2d = SubvectorEltIndex(top_sub, i, j, 0);
        idx_2d_negu = SubvectorEltIndex(top_sub, i - 1, j, 0);
        idx_2d_posu = SubvectorEltIndex(top_sub, i + 1, j, 0);
        idx_2d_negv = SubvectorEltIndex(top_sub, i, j - 1, 0);
        idx_2d_posv = SubvectorEltIndex(top_sub, i, j + 1, 0);

        su = stride_u; sv = stride_v; sw = stride_w;
        p_ = pressure_dat; z_ = z_dat;

        flow_magnitude = bc_patch_values[icounter];

        /* velocity subvector indices jjb */
        vu_negu_idx = SubvectorEltIndex(vel_u_sub, i, j, k);
        vv_negv_idx = SubvectorEltIndex(vel_v_sub, i, j, k);
        vw_negw_idx = SubvectorEltIndex(vel_w_sub, i, j, k);
        vu_posu_idx = SubvectorEltIndex(vel_u_sub, i + 1, j, k);
        vv_posv_idx = SubvectorEltIndex(vel_v_sub, i, j + 1, k);
        vw_posw_idx = SubvectorEltIndex(vel_w_sub, i, j, k + 1);

        // since we are using coordinate transformations, we assume that the
        // grid is regular and that there are no inactive cells.
        make_correction_except_negu = FALSE;
        make_correction_except_posu = FALSE;
        make_correction_except_negv = FALSE;
        make_correction_except_posv = FALSE;
        make_correction_except_negw = FALSE;
        make_correction_except_posw = FALSE;
        is_patch_negu = FALSE;
        is_patch_posu = FALSE;
        is_patch_negv = FALSE;
        is_patch_posv = FALSE;
        is_patch_negw = FALSE;
        is_patch_posw = FALSE;
        make_full_correction = FALSE;
        // special case in top and bottom layers to ensure
        // that the flux is not removed twice
        is_face_negu = top_dat[idx_2d_negu] < 0;
        is_face_posu = top_dat[idx_2d_posu] < 0;
        is_face_negv = top_dat[idx_2d_negv] < 0;
        is_face_posv = top_dat[idx_2d_posv] < 0;
        is_face_negw = k == lrint(bot_dat[idx_2d]);
        is_face_posw = k == lrint(top_dat[idx_2d]);

        is_face_u = (is_face_negu || is_face_posu);
        is_face_v = (is_face_negv || is_face_posv);
        is_face_w = (is_face_negw || is_face_posw);

        u_fluxbc = 0.0;
        u_correction = 0.0;

        dh_face_negu[U] = 0.0; dh_face_negu[V] = 0.0; dh_face_negu[W] = 0.0;
        dh_face_posu[U] = 0.0; dh_face_posu[V] = 0.0; dh_face_posu[W] = 0.0;
        dh_face_negv[U] = 0.0; dh_face_negv[V] = 0.0; dh_face_negv[W] = 0.0;
        dh_face_posv[U] = 0.0; dh_face_posv[V] = 0.0; dh_face_posv[W] = 0.0;
        dh_face_negw[U] = 0.0; dh_face_negw[V] = 0.0; dh_face_negw[W] = 0.0;
        dh_face_posw[U] = 0.0; dh_face_posw[V] = 0.0; dh_face_posw[W] = 0.0;

        // U Face Permutation: (u', v', w') -> (u, v, w)
        HydraulicHeadDel((idx-stride_u), stride_u, stride_v, stride_w,
                        du, dv, dw, dh_face_negu, lengthA_u_dat, lengthB_u_dat,
                        pressure_dat, z_dat, density_dat, gravity);

        HydraulicHeadDel(idx, stride_u, stride_v, stride_w, du, dv, dw,
                        dh_face_posu, lengthA_u_dat, lengthB_u_dat,
                        pressure_dat, z_dat, density_dat, gravity);

        // V Face Permutation: (u', v', w') -> (v, w, u)
        HydraulicHeadDel((idx-stride_v), stride_v, stride_w, stride_u,
                        dv, dw, du, dh_face_negv, lengthA_v_dat, lengthB_v_dat,
                        pressure_dat, z_dat, density_dat, gravity);

        HydraulicHeadDel(idx, stride_v, stride_w, stride_u, dv, dw, du,
                        dh_face_posv, lengthA_v_dat, lengthB_v_dat,
                        pressure_dat, z_dat, density_dat, gravity);

        // W Face Permutation: (u', v', w') -> (w, u, v)
        HydraulicHeadDel((idx-stride_w), stride_w, stride_u, stride_v,
                        dw, du, dv, dh_face_negw, lengthA_w_dat, lengthB_w_dat,
                        pressure_dat, z_dat, density_dat, gravity);

        HydraulicHeadDel(idx, stride_w, stride_u, stride_v, dw, du, dv,
                        dh_face_posw, lengthA_w_dat, lengthB_w_dat,
                        pressure_dat, z_dat, density_dat, gravity);

        jacobian_negu = 0.5 * (jacobian_dat[idx] + jacobian_dat[idx - stride_u]);
        dArea_negu = jacobian_negu * dvdw;

        jacobian_posu = 0.5 * (jacobian_dat[idx] + jacobian_dat[idx + stride_u]);
        dArea_posu = jacobian_posu * dvdw;

        jacobian_negv = 0.5 * (jacobian_dat[idx] + jacobian_dat[idx - stride_v]);
        dArea_negv = jacobian_negv * dwdu;

        jacobian_posv = 0.5 * (jacobian_dat[idx] + jacobian_dat[idx + stride_v]);
        dArea_posv = jacobian_posv * dwdu;

        jacobian_negw = 0.5 * (jacobian_dat[idx] + jacobian_dat[idx - stride_w]);
        dArea_negw = jacobian_negw * dudv;

        jacobian_posw = 0.5 * (jacobian_dat[idx] + jacobian_dat[idx + stride_w]);
        dArea_posw = jacobian_posw * dudv;

        density_negu = DensityMean((idx - su), su, lengthA_u_dat, lengthB_u_dat, density_dat);
        density_posu = DensityMean(idx, su, lengthA_u_dat, lengthB_u_dat, density_dat);
        density_negv = DensityMean((idx - sv), sv, lengthA_v_dat, lengthB_v_dat, density_dat);
        density_posv = DensityMean(idx, sv, lengthA_v_dat, lengthB_v_dat, density_dat);
        density_negw = DensityMean((idx - sw), sw, lengthA_w_dat, lengthB_w_dat, density_dat);
        density_posw = DensityMean(idx, sw, lengthA_w_dat, lengthB_w_dat, density_dat);
      }),
                           FACE(NegUFace, {
        double norm_negu = 0.5 * (lengthA_u_dat[idx - stride_u] + lengthB_u_dat[idx - stride_u]);
        u_fluxbc = dArea_negu * flow_magnitude / norm_negu;
        
        is_patch_negu = TRUE;
        make_correction_except_posu = !(is_face_v || is_face_w);
        make_full_correction = FALSE;

        /* Override derivatives in u */
        // U Face Permutation: (u', v', w') -> (u, v, w)
        // dh_face_negu[U] = 0.0; // does not matter
        // dh_face_posu[U] = 0.0; // does not matter
        // V Face Permutation: (u', v', w') -> (v, w, u)
        dh_face_negv[W] = (0.5 / du) * ((p_[idx + su] - p_[idx] + p_[idx - sv + su] - p_[idx - sv]) + density_negv * gravity * (z_[idx + su] - z_[idx] + z_[idx - sv + su] - z_[idx - sv]));
        dh_face_posv[W] = (0.5 / du) * ((p_[idx + sv + su] - p_[idx + sv] + p_[idx + su] - p_[idx]) + density_posv * gravity * (z_[idx + sv + su] - z_[idx + sv] + z_[idx + su] - z_[idx]));
        // W Face Permutation: (u', v', w') -> (w, u, v)
        dh_face_negw[V] = (0.5 / du) * ((p_[idx + su] - p_[idx] + p_[idx - sw + su] - p_[idx - sw]) + density_negw * gravity * (z_[idx + su] - z_[idx] + z_[idx - sw + su] - z_[idx - sw]));
        dh_face_posw[V] = (0.5 / du) * ((p_[idx + sw + su] - p_[idx + sw] + p_[idx + su] - p_[idx]) + density_posw * gravity * (z_[idx + sw + su] - z_[idx + sw] + z_[idx + su] - z_[idx]));
      }),
                           FACE(PosUFace, {
        double norm_posu = 0.5 * (lengthA_u_dat[idx] + lengthB_u_dat[idx]);
        u_fluxbc = dArea_posu * flow_magnitude / norm_posu;

        is_patch_posu = TRUE;
        make_correction_except_negu = !(is_face_v || is_face_w);
        make_full_correction = FALSE;

        /* Override derivatives in u */
        // U Face Permutation: (u', v', w') -> (u, v, w)
        // dh_face_negu[U] = 0.0; // does not matter
        // dh_face_posu[U] = 0.0; // does not matter
        // V Face Permutation: (u', v', w') -> (v, w, u)
        dh_face_negv[W] = (0.5 / du) * ((p_[idx] - p_[idx - su] + p_[idx - sv] - p_[idx - sv - su]) + density_negv * gravity * (z_[idx] - z_[idx - su] + z_[idx - sv] - z_[idx - sv - su]));
        dh_face_posv[W] = (0.5 / du) * ((p_[idx + sv] - p_[idx + sv - su] + p_[idx] - p_[idx - su]) + density_posv * gravity * (z_[idx + sv] - z_[idx + sv - su] + z_[idx] - z_[idx - su]));
        // W Face Permutation: (u', v', w') -> (w, u, v)
        dh_face_negw[V] = (0.5 / du) * ((p_[idx] - p_[idx - su] + p_[idx - sw] - p_[idx - sw - su]) + density_negw * gravity * (z_[idx] - z_[idx - su] + z_[idx - sw] - z_[idx - sw - su]));
        dh_face_posw[V] = (0.5 / du) * ((p_[idx + sw] - p_[idx + sw - su] + p_[idx] - p_[idx - su]) + density_posw * gravity * (z_[idx + sw] - z_[idx + sw - su] + z_[idx] - z_[idx - su]));
      }),
                           FACE(NegVFace, {
        double norm_negv = 0.5 * (lengthA_v_dat[idx - stride_v] + lengthB_v_dat[idx - stride_v]);
        u_fluxbc = dArea_negv * flow_magnitude / norm_negv;

        is_patch_negv = TRUE;
        make_correction_except_posv = !(is_face_u || is_face_w);
        make_full_correction = (is_face_u && !is_face_w);

        /* Override derivatives in v */
        // U Face Permutation: (u', v', w') -> (u, v, w)
        dh_face_negu[V] = (0.5 / dv) * ((p_[idx + sv] - p_[idx] + p_[idx - su + sv] - p_[idx - su]) + density_negu * gravity * (z_[idx + sv] - z_[idx] + z_[idx - su + sv] - z_[idx - su]));
        dh_face_posu[V] = (0.5 / dv) * ((p_[idx + su + sv] - p_[idx + su] + p_[idx + sv] - p_[idx]) + density_posu * gravity * (z_[idx + su + sv] - z_[idx + su] + z_[idx + sv] - z_[idx]));
        // V Face Permutation: (u', v', w') -> (v, w, u)
        // dh_face_negv[U] = 0.0; // does not matter
        // dh_face_posv[U] = 0.0; // does not matter
        // W Face Permutation: (u', v', w') -> (w, u, v)
        dh_face_negw[W] = (0.5 / dv) * ((p_[idx + sv] - p_[idx] + p_[idx - sw + sv] - p_[idx - sw]) + density_negw * gravity * (z_[idx + sv] - z_[idx] + z_[idx - sw + sv] - z_[idx - sw]));
        dh_face_posw[W] = (0.5 / dv) * ((p_[idx + sw + sv] - p_[idx + sw] + p_[idx + sv] - p_[idx]) + density_posw * gravity * (z_[idx + sw + sv] - z_[idx + sw] + z_[idx + sv] - z_[idx]));

        // corrections for the edges
        if(is_face_negu)
        {
          // dh_face_negu is given by bc
          // dh_face_posu nothing changes
          // dh_face_negw[V] change
          dh_face_negw[V] = (0.5 / du) * ((p_[idx + su] - p_[idx] + p_[idx - sw + su] - p_[idx - sw]) + density_negw * gravity * (z_[idx + su] - z_[idx] + z_[idx - sw + su] - z_[idx - sw]));
          // dh_face_posw[V] change
          dh_face_posw[V] = (0.5 / du) * ((p_[idx + sw + su] - p_[idx + sw] + p_[idx + su] - p_[idx]) + density_posw * gravity * (z_[idx + sw + su] - z_[idx + sw] + z_[idx + su] - z_[idx]));
        }
        if(is_face_posu)
        {
          // dh_face_negu nothing changes
          // dh_face_posu is given by bc
          // dh_face_negw[V] change
          dh_face_negw[V] = (0.5 / du) * ((p_[idx] - p_[idx - su] + p_[idx - sw] - p_[idx - sw - su]) + density_negw * gravity * (z_[idx] - z_[idx - su] + z_[idx - sw] - z_[idx - sw - su]));
          // dh_face_posw[V] change
          dh_face_posw[V] = (0.5 / du) * ((p_[idx + sw] - p_[idx + sw - su] + p_[idx] - p_[idx - su]) + density_posw * gravity * (z_[idx + sw] - z_[idx + sw - su] + z_[idx] - z_[idx - su]));
        }
      }),
                           FACE(PosVFace, {
        double norm_posv = 0.5 * (lengthA_v_dat[idx] + lengthB_v_dat[idx]);
        u_fluxbc = dArea_posv * flow_magnitude / norm_posv;

        is_patch_posv = TRUE;
        make_correction_except_negv = !(is_face_u || is_face_w);
        make_full_correction = (is_face_u && !is_face_w);

        /* Override derivatives in v */
        // U Face Permutation: (u', v', w') -> (u, v, w)
        dh_face_negu[V] = (0.5 / dv) * ((p_[idx] - p_[idx - sv] + p_[idx - su] - p_[idx - su - sv]) + density_negu * gravity * (z_[idx] - z_[idx - sv] + z_[idx - su] - z_[idx - su - sv]));
        dh_face_posu[V] = (0.5 / dv) * ((p_[idx + su] - p_[idx + su - sv] + p_[idx] - p_[idx - sv]) + density_posu * gravity * (z_[idx + su] - z_[idx + su - sv] + z_[idx] - z_[idx - sv]));
        // V Face Permutation: (u', v', w') -> (v, w, u)
        // dh_face_negv[U] = 0.0; // does not matter
        // dh_face_posv[U] = 0.0; // does not matter
        // W Face Permutation: (u', v', w') -> (w, u, v)
        dh_face_negw[W] = (0.5 / dv) * ((p_[idx] - p_[idx - sv] + p_[idx - sw] - p_[idx - sw - sv]) + density_negw * gravity * (z_[idx] - z_[idx - sv] + z_[idx - sw] - z_[idx - sw - sv]));
        dh_face_posw[W] = (0.5 / dv) * ((p_[idx + sw] - p_[idx + sw - sv] + p_[idx] - p_[idx - sv]) + density_posw * gravity * (z_[idx + sw] - z_[idx + sw - sv] + z_[idx] - z_[idx - sv]));

        // corrections for the edges
        if(is_face_negu)
        {
          // dh_face_negu is given by bc
          // dh_face_posu nothing changes
          // dh_face_negw[V] change
          dh_face_negw[V] = (0.5 / du) * ((p_[idx + su] - p_[idx] + p_[idx - sw + su] - p_[idx - sw]) + density_negw * gravity * (z_[idx + su] - z_[idx] + z_[idx - sw + su] - z_[idx - sw]));
          // dh_face_posw[V] change
          dh_face_posw[V] = (0.5 / du) * ((p_[idx + sw + su] - p_[idx + sw] + p_[idx + su] - p_[idx]) + density_posw * gravity * (z_[idx + sw + su] - z_[idx + sw] + z_[idx + su] - z_[idx]));
        }
        if(is_face_posu)
        {
          // dh_face_negu nothing changes
          // dh_face_posu is given by bc
          // dh_face_negw[V] change
          dh_face_negw[V] = (0.5 / du) * ((p_[idx] - p_[idx - su] + p_[idx - sw] - p_[idx - sw - su]) + density_negw * gravity * (z_[idx] - z_[idx - su] + z_[idx - sw] - z_[idx - sw - su]));
          // dh_face_posw[V] change
          dh_face_posw[V] = (0.5 / du) * ((p_[idx + sw] - p_[idx + sw - su] + p_[idx] - p_[idx - su]) + density_posw * gravity * (z_[idx + sw] - z_[idx + sw - su] + z_[idx] - z_[idx - su]));
        }
      }),
                           FACE(NegWFace, {
        double norm_negw = 0.5 * (lengthA_w_dat[idx - stride_w] + lengthB_w_dat[idx - stride_w]);
        u_fluxbc = dArea_negw * flow_magnitude / norm_negw;

        is_patch_negw = TRUE;
        make_correction_except_posw = !(is_face_u || is_face_v);
        make_full_correction = (is_face_u || is_face_v);

        /* Override derivatives in w */
        // U Face Permutation: (u', v', w') -> (u, v, w)
        dh_face_negu[W] = (0.5 / dw) * ((p_[idx + sw] - p_[idx] + p_[idx - su + sw] - p_[idx - su]) + density_negu * gravity * (z_[idx + sw] - z_[idx] + z_[idx - su + sw] - z_[idx - su]));
        dh_face_posu[W] = (0.5 / dw) * ((p_[idx + su + sw] - p_[idx + su] + p_[idx + sw] - p_[idx]) + density_posu * gravity * (z_[idx + su + sw] - z_[idx + su] + z_[idx + sw] - z_[idx]));
        // V Face Permutation: (u', v', w') -> (v, w, u)
        dh_face_negv[V] = (0.5 / dw) * ((p_[idx + sw] - p_[idx] + p_[idx - sv + sw] - p_[idx - sv]) + density_negv * gravity * (z_[idx + sw] - z_[idx] + z_[idx - sv + sw] - z_[idx - sv]));
        dh_face_posv[V] = (0.5 / dw) * ((p_[idx + sv + sw] - p_[idx + sv] + p_[idx + sw] - p_[idx]) + density_posv * gravity * (z_[idx + sv + sw] - z_[idx + sv] + z_[idx + sw] - z_[idx]));
        // W Face Permutation: (u', v', w') -> (w, u, v)
        // dh_face_negw[U] = 0.0; // does not matter
        // dh_face_posw[U] = 0.0; // does not matter

        // corrections for the edges
        if(is_face_negu)
        {
          // dh_face_negu is given by bc
          // dh_face_posu nothing changes
          // dh_face_negv[W] change
          dh_face_negv[W] = (0.5 / du) * ((p_[idx + su] - p_[idx] + p_[idx - sv + su] - p_[idx - sv]) + density_negv * gravity * (z_[idx + su] - z_[idx] + z_[idx - sv + su] - z_[idx - sv]));
          // dh_face_posv[W] change
          dh_face_posv[W] = (0.5 / du) * ((p_[idx + sv + su] - p_[idx + sv] + p_[idx + su] - p_[idx]) + density_posv * gravity * (z_[idx + sv + su] - z_[idx + sv] + z_[idx + su] - z_[idx]));
        }
        if(is_face_posu)
        {
          // dh_face_negu nothing changes
          // dh_face_posu is given by bc
          // dh_face_negv[W] change
          dh_face_negv[W] = (0.5 / du) * ((p_[idx] - p_[idx - su] + p_[idx - sv] - p_[idx - sv - su]) + density_negv * gravity * (z_[idx] - z_[idx - su] + z_[idx - sv] - z_[idx - sv - su]));
          // dh_face_posv[W] change
          dh_face_posv[W] = (0.5 / du) * ((p_[idx + sv] - p_[idx + sv - su] + p_[idx] - p_[idx - su]) + density_posv * gravity * (z_[idx + sv] - z_[idx + sv - su] + z_[idx] - z_[idx - su]));
        }
        if(is_face_negv)
        {
          // dh_face_negu[V] change
          dh_face_negu[V] = (0.5 / dv) * ((p_[idx + sv] - p_[idx] + p_[idx - su + sv] - p_[idx - su]) + density_negu * gravity * (z_[idx + sv] - z_[idx] + z_[idx - su + sv] - z_[idx - su]));
          // dh_face_posu[V] change
          dh_face_posu[V] = (0.5 / dv) * ((p_[idx + su + sv] - p_[idx + su] + p_[idx + sv] - p_[idx]) + density_posu * gravity * (z_[idx + su + sv] - z_[idx + su] + z_[idx + sv] - z_[idx]));
          // dh_face_negv is given by bc
          // dh_face_posv nothing changes
        }
        if(is_face_posv)
        {
          // dh_face_negu[V] change
          dh_face_negu[V] = (0.5 / dv) * ((p_[idx] - p_[idx - sv] + p_[idx - su] - p_[idx - su - sv]) + density_negu * gravity * (z_[idx] - z_[idx - sv] + z_[idx - su] - z_[idx - su - sv]));
          // dh_face_posu[V] change
          dh_face_posu[V] = (0.5 / dv) * ((p_[idx + su] - p_[idx + su - sv] + p_[idx] - p_[idx - sv]) + density_posu * gravity * (z_[idx + su] - z_[idx + su - sv] + z_[idx] - z_[idx - sv]));
          // dh_face_negv nothing changes
          // dh_face_posv is given by bc
        }
      }),
                           FACE(PosWFace, {
        double norm_posw = 0.5 * (lengthA_w_dat[idx] + lengthB_w_dat[idx]);
        u_fluxbc = dArea_posw * flow_magnitude / norm_posw;
        
        is_patch_posw = TRUE;
        make_correction_except_negw = !(is_face_u || is_face_v);
        make_full_correction = (is_face_u || is_face_v);

        /* Override derivatives in w */
        // U Face Permutation: (u', v', w') -> (u, v, w)
        dh_face_negu[W] = (0.5 / dw) * ((p_[idx] - p_[idx - sw] + p_[idx - su] - p_[idx - su - sw]) + density_negu * gravity * (z_[idx] - z_[idx - sw] + z_[idx - su] - z_[idx - su - sw]));
        dh_face_posu[W] = (0.5 / dw) * ((p_[idx + su] - p_[idx + su - sw] + p_[idx] - p_[idx - sw]) + density_posu * gravity * (z_[idx + su] - z_[idx + su - sw] + z_[idx] - z_[idx - sw]));
        // V Face Permutation: (u', v', w') -> (v, w, u)
        dh_face_negv[V] = (0.5 / dw) * ((p_[idx] - p_[idx - sw] + p_[idx - sv] - p_[idx - sv - sw]) + density_negv * gravity * (z_[idx] - z_[idx - sw] + z_[idx - sv] - z_[idx - sv - sw]));
        dh_face_posv[V] = (0.5 / dw) * ((p_[idx + sv] - p_[idx + sv - sw] + p_[idx] - p_[idx - sw]) + density_posv * gravity * (z_[idx + sv] - z_[idx + sv - sw] + z_[idx] - z_[idx - sw]));
        // W Face Permutation: (u', v', w') -> (w, u, v)
        // dh_face_negw[U] = 0.0; // does not matter
        // dh_face_posw[U] = 0.0; // does not matter

        // corrections for the edges
        if(is_face_negu)
        {
          // dh_face_negu is given by bc
          // dh_face_posu nothing changes
          // dh_face_negv[W] change
          dh_face_negv[W] = (0.5 / du) * ((p_[idx + su] - p_[idx] + p_[idx - sv + su] - p_[idx - sv]) + density_negv * gravity * (z_[idx + su] - z_[idx] + z_[idx - sv + su] - z_[idx - sv]));
          // dh_face_posv[W] change
          dh_face_posv[W] = (0.5 / du) * ((p_[idx + sv + su] - p_[idx + sv] + p_[idx + su] - p_[idx]) + density_posv * gravity * (z_[idx + sv + su] - z_[idx + sv] + z_[idx + su] - z_[idx]));
        }
        if(is_face_posu)
        {
          // dh_face_negu nothing changes
          // dh_face_posu is given by bc
          // dh_face_negv[W] change
          dh_face_negv[W] = (0.5 / du) * ((p_[idx] - p_[idx - su] + p_[idx - sv] - p_[idx - sv - su]) + density_negv * gravity * (z_[idx] - z_[idx - su] + z_[idx - sv] - z_[idx - sv - su]));
          // dh_face_posv[W] change
          dh_face_posv[W] = (0.5 / du) * ((p_[idx + sv] - p_[idx + sv - su] + p_[idx] - p_[idx - su]) + density_posv * gravity * (z_[idx + sv] - z_[idx + sv - su] + z_[idx] - z_[idx - su]));
        }
        if(is_face_negv)
        {
          // dh_face_negu[V] change
          dh_face_negu[V] = (0.5 / dv) * ((p_[idx + sv] - p_[idx] + p_[idx - su + sv] - p_[idx - su]) + density_negu * gravity * (z_[idx + sv] - z_[idx] + z_[idx - su + sv] - z_[idx - su]));
          // dh_face_posu[V] change
          dh_face_posu[V] = (0.5 / dv) * ((p_[idx + su + sv] - p_[idx + su] + p_[idx + sv] - p_[idx]) + density_posu * gravity * (z_[idx + su + sv] - z_[idx + su] + z_[idx + sv] - z_[idx]));
          // dh_face_negv is given by bc
          // dh_face_posv nothing changes
        }
        if(is_face_posv)
        {
          // dh_face_negu[V] change
          dh_face_negu[V] = (0.5 / dv) * ((p_[idx] - p_[idx - sv] + p_[idx - su] - p_[idx - su - sv]) + density_negu * gravity * (z_[idx] - z_[idx - sv] + z_[idx - su] - z_[idx - su - sv]));
          // dh_face_posu[V] change
          dh_face_posu[V] = (0.5 / dv) * ((p_[idx + su] - p_[idx + su - sv] + p_[idx] - p_[idx - sv]) + density_posu * gravity * (z_[idx + su] - z_[idx + su - sv] + z_[idx] - z_[idx - sv]));
          // dh_face_negv nothing changes
          // dh_face_posv is given by bc
        }
      }),
                           CellFinalize(
      {
        // U Face Permutation: (u', v', w') -> (u, v, w)
        double u_negu = FluxAtFace((idx-stride_u), stride_u, dArea_negu,
                                  dh_face_negu, density_dat, rel_perm_dat,
                                  viscosity, K_uu_dat, K_uv_dat, K_uw_dat);

        double u_posu = FluxAtFace(idx, stride_u, dArea_posu, dh_face_posu,
                                  density_dat, rel_perm_dat, viscosity,
                                  K_uu_dat, K_uv_dat, K_uw_dat);

        // V Face Permutation: (u', v', w') -> (v, w, u)
        double u_negv = FluxAtFace((idx-stride_v), stride_v, dArea_negv,
                                  dh_face_negv, density_dat, rel_perm_dat,
                                  viscosity, K_vv_dat, K_vw_dat, K_vu_dat);

        double u_posv = FluxAtFace(idx, stride_v, dArea_posv, dh_face_posv,
                                  density_dat, rel_perm_dat, viscosity,
                                  K_vv_dat, K_vw_dat, K_vu_dat);

        // W Face Permutation: (u', v', w') -> (w, u, v)
        double u_negw = FluxAtFace((idx-stride_w), stride_w, dArea_negw,
                                  dh_face_negw, density_dat, rel_perm_dat,
                                  viscosity, K_ww_dat, K_wu_dat, K_wv_dat);

        double u_posw = FluxAtFace(idx, stride_w, dArea_posw, dh_face_posw,
                                  density_dat, rel_perm_dat, viscosity,
                                  K_ww_dat, K_wu_dat, K_wv_dat);

        // vel_u_dat[vu_negu_idx] = u_negu / dArea_negu;
        // vel_u_dat[vu_posu_idx] = u_posu / dArea_posu;
        // vel_v_dat[vv_negv_idx] = u_negv / dArea_negv;
        // vel_v_dat[vv_posv_idx] = u_posv / dArea_posv;
        // vel_w_dat[vw_negw_idx] = u_negw / dArea_negw;
        // vel_w_dat[vw_posw_idx] = u_posw / dArea_posw;

        if(is_patch_negu) vel_u_dat[vu_negu_idx] = u_fluxbc / dArea_negu;
        if(is_patch_posu) vel_u_dat[vu_posu_idx] = u_fluxbc / dArea_posu;
        if(is_patch_negv) vel_v_dat[vv_negv_idx] = u_fluxbc / dArea_negv;
        if(is_patch_posv) vel_v_dat[vv_posv_idx] = u_fluxbc / dArea_posv;
        if(is_patch_negw) vel_w_dat[vw_negw_idx] = u_fluxbc / dArea_negw;
        if(is_patch_posw) vel_w_dat[vw_posw_idx] = u_fluxbc / dArea_posw;

        u_correction += (is_patch_negu || is_patch_posu) * ((!is_face_w) * 
                  (- (!is_face_negv) * (u_negv) + (!is_face_posv) * (u_posv))
                  + (!is_face_v) *
                  (- (!is_face_negw) * (u_negw) + (!is_face_posw) * (u_posw)));
        u_correction += (is_patch_negv || is_patch_posv) * ((!is_face_w) * 
                    (- (!is_face_negu) * (u_negu) + (!is_face_posu) * (u_posu))
                     - (!is_face_negw) * (u_negw) + (!is_face_posw) * (u_posw));
        u_correction += (is_patch_negw || is_patch_posw) * 
                    (- (!is_face_negu) * (u_negu) + (!is_face_posu) * (u_posu)
                     - (!is_face_negv) * (u_negv) + (!is_face_posv) * (u_posv));
        u_correction += u_fluxbc;

        PlusEquals(fval_dat[idx], dt * u_correction);
      }),
                           AfterAllCells(DoNothing)
                           ); /* End FluxBC */
    }          /* End ipatch loop */
  }            /* End subgrid loop */

  /*
   * Reset values inserted for the DirichletBC back to the decoupled
   * problem used in the inactive cells.
   *
   * See comments above on why this is needed.
   */
  is = 0;
  ForSubgridI(is, GridSubgrids(grid))
  {
    Subgrid *subgrid = GridSubgrid(grid, is);

    Subvector *fval_sub = VectorSubvector(fval, is);
    double *fval_dat = SubvectorData(fval_sub);

    int stride_u = 1;
    int stride_v = SubvectorNX(fval_sub);
    int stride_w = SubvectorNY(fval_sub) * stride_v;

    double *pressure_dat = SubvectorData(VectorSubvector(pressure, is));

    int ipatch = 0;
    ForBCStructNumPatches(ipatch, bc_struct)
    {
      int ival = 0, i = 0, j = 0, k = 0;
      ForPatchCellsPerFace(DirichletBC,
                           BeforeAllCells(DoNothing),
                           LoopVars(i, j, k, ival, bc_struct, ipatch, is),
                           Locals(int press_idx, idx; ),
                           CellSetup({
        press_idx = 0;
        idx = SubvectorEltIndex(fval_sub, i, j, k);
      }),
                           FACE(NegUFace, { press_idx = idx - stride_u; }),
                           FACE(PosUFace, { press_idx = idx + stride_u; }),
                           FACE(NegVFace, { press_idx = idx - stride_v; }),
                           FACE(PosVFace, { press_idx = idx + stride_v; }),
                           FACE(NegWFace, { press_idx = idx - stride_w; }),
                           FACE(PosWFace, { press_idx = idx + stride_w; }),
                           CellFinalize({
        pressure_dat[press_idx] = 0.0; //-FLT_MAX;
      }),
                           AfterAllCells(DoNothing)
                           );

      ival = 0, i = 0, j = 0, k = 0;
      ForPatchCellsPerFace(BC_ALL,
                           BeforeAllCells(DoNothing),
                           LoopVars(i, j, k, ival, bc_struct, ipatch, is),
                           Locals(int fval_idx, idx; ),
                           CellSetup({
        fval_idx = 0;
        idx = SubvectorEltIndex(fval_sub, i, j, k);
      }),
                           FACE(NegUFace, { fval_idx = idx - stride_u; }),
                           FACE(PosUFace, { fval_idx = idx + stride_u; }),
                           FACE(NegVFace, { fval_idx = idx - stride_v; }),
                           FACE(PosVFace, { fval_idx = idx + stride_v; }),
                           FACE(NegWFace, { fval_idx = idx - stride_w; }),
                           FACE(PosWFace, { fval_idx = idx + stride_w; }),
                           CellFinalize({
        fval_dat[fval_idx] = 0.0;
      }),
                           AfterAllCells(DoNothing)
                           );
    }          /* End ipatch loop */
  }            /* End subgrid loop */

  FreeBCStruct(bc_struct);

  PFModuleInvokeType(RichardsBCInternalInvoke, bc_internal,
                     (problem, problem_data, fval, NULL, time,
                      pressure, CALCFCN));

  EndTiming(public_xtra->time_index);

  POP_NVTX

  return;
}


/*--------------------------------------------------------------------------
 * NlFunctionEvalInitInstanceXtra
 *--------------------------------------------------------------------------*/

PFModule    *NlFunctionEvalInitInstanceXtra(Problem *problem,
                                            Grid *   grid,
                                            double * temp_data)

{
  PFModule      *this_module = ThisPFModule;
  InstanceXtra  *instance_xtra;

  (void)grid;
  (void)temp_data;

  if (PFModuleInstanceXtra(this_module) == NULL)
    instance_xtra = ctalloc(InstanceXtra, 1);
  else
    instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);

  if (problem != NULL)
  {
    (instance_xtra->problem) = problem;
  }

  if (PFModuleInstanceXtra(this_module) == NULL)
  {
    (instance_xtra->density_module) =
      PFModuleNewInstance(ProblemPhaseDensity(problem), ());
    (instance_xtra->saturation_module) =
      PFModuleNewInstanceType(SaturationInitInstanceXtraInvoke,
                              ProblemSaturation(problem), (NULL, NULL));
    (instance_xtra->rel_perm_module) =
      PFModuleNewInstanceType(PhaseRelPermInitInstanceXtraInvoke,
                              ProblemPhaseRelPerm(problem), (NULL, NULL));
    (instance_xtra->phase_source) =
      PFModuleNewInstance(ProblemPhaseSource(problem), ());
    (instance_xtra->bc_pressure) =
      PFModuleNewInstanceType(BCPressurePackageInitInstanceXtraInvoke,
                              ProblemBCPressure(problem), (problem));
    (instance_xtra->bc_internal) =
      PFModuleNewInstance(ProblemBCInternal(problem), ());
    (instance_xtra->overlandflow_module) =
      PFModuleNewInstance(ProblemOverlandFlowEval(problem), ());     //DOK
    (instance_xtra->overlandflow_module_diff) =
      PFModuleNewInstance(ProblemOverlandFlowEvalDiff(problem), ());   //@RMM
    (instance_xtra->overlandflow_module_kin) =
      PFModuleNewInstance(ProblemOverlandFlowEvalKin(problem), ());
  }
  else
  {
    PFModuleReNewInstance((instance_xtra->density_module), ());
    PFModuleReNewInstanceType(SaturationInitInstanceXtraInvoke,
                              (instance_xtra->saturation_module),
                              (NULL, NULL));
    PFModuleReNewInstanceType(PhaseRelPermInitInstanceXtraInvoke,
                              (instance_xtra->rel_perm_module),
                              (NULL, NULL));
    PFModuleReNewInstanceType(BCPressurePackageInitInstanceXtraInvoke,
                              (instance_xtra->phase_source), (NULL));
    PFModuleReNewInstanceType(BCPressurePackageInitInstanceXtraInvoke,
                              (instance_xtra->bc_pressure), (problem));
    PFModuleReNewInstance((instance_xtra->bc_internal), ());
    PFModuleReNewInstance((instance_xtra->overlandflow_module), ());     //DOK
    PFModuleReNewInstance((instance_xtra->overlandflow_module_diff), ());      //@RMM
    PFModuleReNewInstance((instance_xtra->overlandflow_module_kin), ());
  }

  PFModuleInstanceXtra(this_module) = instance_xtra;
  return this_module;
}


/*--------------------------------------------------------------------------
 * NlFunctionEvalFreeInstanceXtra
 *--------------------------------------------------------------------------*/

void  NlFunctionEvalFreeInstanceXtra()
{
  PFModule      *this_module = ThisPFModule;
  InstanceXtra  *instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);

  if (instance_xtra)
  {
    PFModuleFreeInstance(instance_xtra->density_module);
    PFModuleFreeInstance(instance_xtra->saturation_module);
    PFModuleFreeInstance(instance_xtra->rel_perm_module);
    PFModuleFreeInstance(instance_xtra->phase_source);
    PFModuleFreeInstance(instance_xtra->bc_pressure);
    PFModuleFreeInstance(instance_xtra->bc_internal);
    PFModuleFreeInstance(instance_xtra->overlandflow_module);     //DOK
    PFModuleFreeInstance(instance_xtra->overlandflow_module_diff);      //@RMM
    PFModuleFreeInstance(instance_xtra->overlandflow_module_kin);

    tfree(instance_xtra);
  }
}


/*--------------------------------------------------------------------------
 * NlFunctionEvalNewPublicXtra
 *--------------------------------------------------------------------------*/

PFModule   *NlFunctionEvalNewPublicXtra(char *name)
{
  PFModule      *this_module = ThisPFModule;
  PublicXtra    *public_xtra;
  char key[IDB_MAX_KEY_LEN];
  char *switch_name;
  int switch_value;
  NameArray upwind_switch_na;


  public_xtra = ctalloc(PublicXtra, 1);

/* These parameters dampen the transition/switching into overland flow to speedup
 * the spinup process. */
  sprintf(key, "OverlandSpinupDampP1");
  public_xtra->SpinupDampP1 = GetDoubleDefault(key, 0.0);
  sprintf(key, "OverlandSpinupDampP2");
  public_xtra->SpinupDampP2 = GetDoubleDefault(key, 0.0);    //NBE

  ///* parameters for upwinding formulation for TFG */
  upwind_switch_na = NA_NewNameArray("Original UpwindSine Upwind");
  sprintf(key, "Solver.TerrainFollowingGrid.SlopeUpwindFormulation");
  switch_name = GetStringDefault(key, "Original");
  switch_value = NA_NameToIndexExitOnError(upwind_switch_na, switch_name, key);
  switch (switch_value)
  {
    case 0:
    {
      public_xtra->tfgupwind = 0;
      break;
    }

    case 1:
    {
      public_xtra->tfgupwind = 1;
      break;
    }

    case 2:
    {
      public_xtra->tfgupwind = 2;
      break;
    }

    default:
    {
      InputError("Invalid switch value <%s> for key <%s>", switch_name, key);
    }
  }
  NA_FreeNameArray(upwind_switch_na);

  (public_xtra->time_index) = RegisterTiming("NL_F_Eval");

  PFModulePublicXtra(this_module) = public_xtra;

  return this_module;
}


/*--------------------------------------------------------------------------
 * NlFunctionEvalFreePublicXtra
 *--------------------------------------------------------------------------*/

void  NlFunctionEvalFreePublicXtra()
{
  PFModule    *this_module = ThisPFModule;
  PublicXtra  *public_xtra = (PublicXtra*)PFModulePublicXtra(this_module);


  if (public_xtra)
  {
    tfree(public_xtra);
  }
}


/*--------------------------------------------------------------------------
 * NlFunctionEvalSizeOfTempData
 *--------------------------------------------------------------------------*/

int  NlFunctionEvalSizeOfTempData()
{
  return 0;
}
