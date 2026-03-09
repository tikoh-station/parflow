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
  /* Init PSBLAS Context */
  psb_c_init_from_fint(PSBLASSessionContext(session), MPI_Comm_c2f(amps_CommWorld));

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
  psb_c_cdall_vl(nl, vl, *PSBLASSessionContext(session), PSBLASSessionDescriptor(session));
  free(vl);

  /* context descriptor is finalized in psb_c_cdasb */
  psb_i_t info = psb_c_cdasb(PSBLASSessionDescriptor(session));
  if (info != 0) {
    amps_Printf("Error in psb_c_cdasb: %d\n", info);
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
