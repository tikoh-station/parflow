/*BHEADER**********************************************************************
*
*  Copyright (c) 1995-2026, Lawrence Livermore National Security,
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
#include "psblas_session.h"

PSBLASSession* NewPSBLASSession()
{
  PSBLASSession *session = (PSBLASSession*)ctalloc(PSBLASSession, 1);

  /* Create new PSBLAS Context */
  PSBLASSessionContext(session) = psb_c_new_ctxt();
  /* Create new PSBLAS Descriptor */
  PSBLASSessionDescriptor(session) = psb_c_new_descriptor();

  return session;
}

void FreePSBLASSession(PSBLASSession *session)
{
  if (session != NULL)
  {
    SUNMatDestroy(PSBLASSessionSUNMatrix(session));

    psb_c_cdfree(PSBLASSessionDescriptor(session));
    psb_c_delete_descriptor(PSBLASSessionDescriptor(session));

    // psb_c_exit_ctxt(*cctxt);
    psb_c_delete_ctxt(PSBLASSessionContext(session));

    tfree(session);
  }

  return;
}

void InitPSBLASSession(PSBLASSession *session, Grid *grid)
{
  psb_i_t info = 0;

  /* Init PSBLAS Context */
  psb_c_init_from_fint(PSBLASSessionContext(session), MPI_Comm_c2f(amps_CommWorld));
  psb_c_set_index_base(0);

  /* Count number of elements in current process. */
  psb_i_t nl = 0;
  int is = 0;
  ForSubgridI(is, GridSubgrids(grid))
  {
    Subgrid *sg = GridSubgrid(grid, is);
    nl += (psb_i_t)(SubgridNX(sg) * SubgridNY(sg) * SubgridNZ(sg));
  }

  /* Set local to glabal index mapping */
  psb_l_t *vl = ctalloc(psb_l_t, nl);
  Subgrid *user_subgrid = GridSubgrid(GlobalsUserGrid, 0);
  ForSubgridI(is, GridSubgrids(grid))
  {
    Subgrid *subgrid = GridSubgrid(grid, is);

    int ix = SubgridIX(subgrid);
    int iy = SubgridIY(subgrid);
    int iz = SubgridIZ(subgrid);
    int nx = SubgridNX(subgrid);
    int ny = SubgridNY(subgrid);
    int nz = SubgridNZ(subgrid);

    int i = 0, j = 0, k = 0;
    BoxLoopI0(i, j, k, ix, iy, iz, nx, ny, nz,
    {
      int local_idx = SubgridEltIndex(subgrid, i, j, k);
      int global_idx = SubgridEltIndex(user_subgrid, i, j, k);
      vl[local_idx] = (psb_l_t)global_idx;
    });
  }

  /* allocate a context descriptor */
  info = psb_c_cdall_vl(nl, vl, *PSBLASSessionContext(session), PSBLASSessionDescriptor(session));
  free(vl);
  if (info != 0) {
    amps_Printf("Error in psb_c_cdall_vl: %d\n", info);
  }

  /* context descriptor is finalized in psb_c_cdasb */
  info = psb_c_cdasb(PSBLASSessionDescriptor(session));
  if (info != 0) {
    amps_Printf("Error in psb_c_cdasb: %d\n", info);
  }

  /* Create PSBLAS SUNMatrix */
  PSBLASSessionSUNMatrix(session) = SUNPSBLASMatrix(
      PSBLASSessionContext(session), 
      PSBLASSessionDescriptor(session)
    );

  if (PSBLASSessionSUNMatrix(session) == NULL) {
    amps_Printf("Error: Failure to create a SUNMatrix\n");
  }

  return;
}

void Set_N_Vector_From_Vector(N_Vector nvec, Vector *vec)
{
  double *psb_data = N_VGetArrayPointer(nvec);

  Grid *grid = VectorGrid(vec);

  int is = 0;
  ForSubgridI(is, GridSubgrids(grid))
  {
    Subgrid *subgrid = GridSubgrid(grid, is);
    Subvector *v_sub = VectorSubvector(vec, is);
    double *v_data = SubvectorData(v_sub);

    int ix = SubgridIX(subgrid);
    int iy = SubgridIY(subgrid);
    int iz = SubgridIZ(subgrid);
    int nx = SubgridNX(subgrid);
    int ny = SubgridNY(subgrid);
    int nz = SubgridNZ(subgrid);
    int nx_v = SubvectorNX(v_sub);
    int ny_v = SubvectorNY(v_sub);
    int nz_v = SubvectorNZ(v_sub);

    int i = 0, j = 0, k = 0, pf_idx = 0;
    BoxLoopI1(i, j, k, ix, iy, iz, nx, ny, nz,
              pf_idx, nx_v, ny_v, nz_v, 1, 1, 1,
    {
      int psb_idx = SubgridEltIndex(subgrid, i, j, k);
      psb_data[psb_idx] = v_data[pf_idx];
    });
  }
  return;
}

void Set_Vector_From_N_Vector(Vector *vec, N_Vector nvec)
{
  double *psb_data = N_VGetArrayPointer(nvec);

  Grid *grid = VectorGrid(vec);

  int is = 0;
  ForSubgridI(is, GridSubgrids(grid))
  {
    Subgrid *subgrid = GridSubgrid(grid, is);
    Subvector *v_sub = VectorSubvector(vec, is);
    double *v_data = SubvectorData(v_sub);

    int ix = SubgridIX(subgrid);
    int iy = SubgridIY(subgrid);
    int iz = SubgridIZ(subgrid);
    int nx = SubgridNX(subgrid);
    int ny = SubgridNY(subgrid);
    int nz = SubgridNZ(subgrid);
    int nx_v = SubvectorNX(v_sub);
    int ny_v = SubvectorNY(v_sub);
    int nz_v = SubvectorNZ(v_sub);

    int i = 0, j = 0, k = 0, pf_idx = 0;
    BoxLoopI1(i, j, k, ix, iy, iz, nx, ny, nz,
              pf_idx, nx_v, ny_v, nz_v, 1, 1, 1,
    {
      int psb_idx = SubgridEltIndex(subgrid, i, j, k);
      v_data[pf_idx] = psb_data[psb_idx];
    });
  }
  return;
}

void Set_SUNMatrix_From_Matrix(SUNMatrix sunmat,
                               Matrix *JB,
                               Matrix *JC,
                               void *current_state)
{
  Subgrid *user_subgrid = GridSubgrid(GlobalsUserGrid, 0);

  Grid *JB_grid = MatrixGrid(JB);
  Stencil *JB_stencil = MatrixStencil(JB);
  int JB_stencil_size = StencilSize(JB_stencil);
  StencilElt *JB_shape = StencilShape(JB_stencil);

  psb_l_t *idx_row = ctalloc(psb_l_t, JB_stencil_size);
  psb_l_t *idx_col = ctalloc(psb_l_t, JB_stencil_size);
  psb_d_t *psb_val = ctalloc(psb_d_t, JB_stencil_size);

  int isubgrid = 0;
  ForSubgridI(isubgrid, GridSubgrids(JB_grid))
  {
    Subgrid *subgrid = SubgridArraySubgrid(GridSubgrids(JB_grid), isubgrid);

    Submatrix *JB_sub = MatrixSubmatrix(JB, isubgrid);

    int ix = SubgridIX(subgrid);
    int iy = SubgridIY(subgrid);
    int iz = SubgridIZ(subgrid);
    int nx = SubgridNX(subgrid);
    int ny = SubgridNY(subgrid);
    int nz = SubgridNZ(subgrid);
    int nx_m = SubmatrixNX(JB_sub);
    int ny_m = SubmatrixNY(JB_sub);
    int nz_m = SubmatrixNZ(JB_sub);

    /* Insert contributions from JB Matrix */
    int i = 0, j = 0, k = 0, pf_idx = 0;
    BoxLoopI1(i, j, k, ix, iy, iz, nx, ny, nz,
              pf_idx, nx_m, ny_m, nz_m, 1, 1, 1,
    {
      int psb_row_idx = SubgridEltIndex(user_subgrid, i, j, k); // row index

      int num_invalid_elements = 0;
      for(int istencil = 0; istencil < JB_stencil_size; ++istencil)
      {
        double *JB_dat = SubmatrixStencilData(JB_sub, istencil);

        int st_i = JB_shape[istencil][0];
        int st_j = JB_shape[istencil][1];
        int st_k = JB_shape[istencil][2];

        if (i + st_i < SubgridIX(user_subgrid) ||
            i + st_i >= SubgridIX(user_subgrid) + SubgridNX(user_subgrid) ||
            j + st_j < SubgridIY(user_subgrid) ||
            j + st_j >= SubgridIY(user_subgrid) + SubgridNY(user_subgrid) ||
            k + st_k < SubgridIZ(user_subgrid) ||
            k + st_k >= SubgridIZ(user_subgrid) + SubgridNZ(user_subgrid))
        {
          ++num_invalid_elements;
        }
        else
        {
          int psb_col_idx = SubgridEltIndex(user_subgrid,
            (i + st_i), (j + st_j), (k + st_k));

          idx_row[istencil - num_invalid_elements] = psb_row_idx; // row index
          idx_col[istencil - num_invalid_elements] = psb_col_idx; // col index
          psb_val[istencil - num_invalid_elements] = JB_dat[pf_idx]; // value
        }
      }

      int num_valid_elements = JB_stencil_size - num_invalid_elements;
      if (num_valid_elements > 0)
      {
        SUNMatIns_PSBLAS(num_valid_elements, idx_row, idx_col, psb_val, sunmat);
      }
    });
  }

  if(JC == NULL)
  {
    free(idx_row);
    free(idx_col);
    free(psb_val);
    return;
  }

  Stencil *JC_stencil = MatrixStencil(JC);
  int JC_stencil_size = StencilSize(JC_stencil);
  StencilElt *JC_shape = StencilShape(JC_stencil);

  ProblemData *problem_data = StateProblemData(((State*)current_state));
  Vector *top = ProblemDataIndexOfDomainTop(problem_data);

  isubgrid = 0;
  ForSubgridI(isubgrid, GridSubgrids(JB_grid))
  {
    Subgrid *subgrid = SubgridArraySubgrid(GridSubgrids(JB_grid), isubgrid);

    Submatrix *JC_sub = MatrixSubmatrix(JC, isubgrid);
    Subvector *top_sub = VectorSubvector(top, isubgrid);
    double *top_dat = SubvectorData(top_sub);

    int ix = SubgridIX(subgrid);
    int iy = SubgridIY(subgrid);
    int iz = SubgridIZ(subgrid);
    int nx = SubgridNX(subgrid);
    int ny = SubgridNY(subgrid);
    int nz = SubgridNZ(subgrid);
    int nx_m = SubmatrixNX(JC_sub);
    int ny_m = SubmatrixNY(JC_sub);
    int nz_m = SubmatrixNZ(JC_sub);

    /* Insert contributions from JC Matrix */
    int i = 0, j = 0, k = 0, pf_idx = 0;
    BoxLoopI1(i, j, k, ix, iy, iz, nx, ny, 1,
              pf_idx, nx_m, ny_m, nz_m, 1, 1, 1,
    {
      int itop = SubvectorEltIndex(top_sub, i, j, 0);
      int k_ = (int)top_dat[itop];

      if (k_ >= 0)
      {
        int psb_row_idx = SubgridEltIndex(user_subgrid, i, j, k_); // row index

        int num_invalid_elements = 0;
        for(int istencil = 0; istencil < JC_stencil_size; ++istencil)
        {
          double *JC_dat = SubmatrixStencilData(JC_sub, istencil);

          int st_i = JC_shape[istencil][0];
          int st_j = JC_shape[istencil][1];

          itop = SubvectorEltIndex(top_sub, (i + st_i), (j + st_j), 0);
          int kk = (int)top_dat[itop];

          if (kk < 0)
          {
            ++num_invalid_elements;
          }
          else
          {
            int psb_col_idx = SubgridEltIndex(user_subgrid,
              (i + st_i), (j + st_j), (kk));

            idx_row[istencil - num_invalid_elements] = psb_row_idx; // row index
            idx_col[istencil - num_invalid_elements] = psb_col_idx; // col index
            psb_val[istencil - num_invalid_elements] = JC_dat[pf_idx]; // value
          }
        }

        int num_valid_elements = JC_stencil_size - num_invalid_elements;
        if (num_valid_elements > 0)
        {
          SUNMatIns_PSBLAS(num_valid_elements, idx_row, idx_col, psb_val, sunmat);
        }
      }
    });
  }

  free(idx_row);
  free(idx_col);
  free(psb_val);

  return;
}