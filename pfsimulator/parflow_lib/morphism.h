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

#ifndef _MORPHISM_HEADER
#define _MORPHISM_HEADER

#include "parflow.h"
#include "v3algebra.h"
#include "v3basis.h"

// structure to support curvilinear coordinates
typedef struct {

  PFModule *morphism_core;

  // transform from curvilinear `zeta` to cartesian `x` coordinates
  v3 (*eval)(PFModule *morphism_core, v3 zeta);

  // transform from cartesian `x` to curvilinear `zeta` coordinates
  v3 (*inverse)(PFModule *morphism_core, v3 x);

  // compute covariant basis vectors
  v3basis (*del)(PFModule *morphism_core, v3 zeta);

  // compute contravariant basis vectors
  v3basis (*del_inverse)(PFModule *morphism_core, v3 zeta);

  // compute jacobian
  double (*jacobian)(PFModule *morphism_core, v3 zeta);


  // compute contravariant components of `vec` from its cartesian representation
  v3 (*to_contravariant)(v3 vec, v3basis basis_contravariant);

  // compute covariant components of `vec` from its cartesian representation
  v3 (*to_covariant)(v3 vec, v3basis basis_covariant);

  // compute cartesian representation of `vec` from its contravariant components
  v3 (*from_contravariant)(v3 vec, v3basis basis_covariant);

  // compute cartesian representation of `vec` from its covariant components
  v3 (*from_covariant)(v3 vec, v3basis basis_contravariant);

} Morphism;



#define MorphismEval(my_morphism, zeta)                                \
        (*(my_morphism->eval))(my_morphism->morphism_core, zeta)

#define MorphismInverse(my_morphism, zeta)                             \
        (*(my_morphism->inverse))(my_morphism->morphism_core, zeta)

#define MorphismDel(my_morphism, zeta)                                 \
        (*(my_morphism->del))(my_morphism->morphism_core, zeta)

#define MorphismDelInverse(my_morphism, zeta)                          \
        (*(my_morphism->del_inverse))(my_morphism->morphism_core, zeta)

#define MorphismJacobian(my_morphism, zeta)                            \
        (*(my_morphism->jacobian))(my_morphism->morphism_core, zeta)



v3 MorphismToContravariant(v3 vec, v3basis basis_contravariant);

v3 MorphismToCovariant(v3 vec, v3basis basis_covariant);

v3 MorphismFromContravariant(v3 vec, v3basis basis_covariant);

v3 MorphismFromCovariant(v3 vec, v3basis basis_contravariant);


#endif // _MORPHISM_HEADER