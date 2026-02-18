#-----------------------------------------------------------------------------
# - Find PSCTOOLKIT includes and libraries.
#
# This module finds if PSCTOOLKIT is installed and determines where the
# include files and libraries are.  This code sets the following variables:
#  PSCTOOLKIT_FOUND         = PSCTOOLKIT was found
#  PSCTOOLKIT_INCLUDE_DIR   = path to where header files can be found
#  PSCTOOLKIT_LIBRARIES     = link libraries for PSCTOOLKIT
#-----------------------------------------------------------------------------

include(FindPackageHandleStandardArgs)

find_package(PSBLAS REQUIRED COMPONENTS psb_cbind psb_util psb_linsolve psb_prec psb_ext psb_base)

find_package(AMG REQUIRED COMPONENTS amg_cbind amg_prec)

find_package(SUNDIALS REQUIRED COMPONENTS sundials_cvode sundials_kinsol sundials_core sundials_nvecpsblas sundials_sunmatrixpsblas sundials_sunlinsolpsblas sundials_sunlinsolspgmr)

set(PSCTOOLKIT_FOUND NO)

if (PSBLAS_FOUND AND AMG_FOUND AND SUNDIALS_FOUND)
  
  set(PSCTOOLKIT_FOUND YES)

  set(PSCTOOLKIT_INCLUDE_DIR ${SUNDIALS_INCLUDE_DIR} ${AMG_INCLUDE_DIR} ${PSBLAS_INCLUDE_DIR} ${PSBLAS_MODULES_DIR} ${AMG_MODULES_DIR})

  set(PSCTOOLKIT_LIBRARIES ${SUNDIALS_LIBRARIES} ${AMG_LIBRARIES} ${PSBLAS_LIBRARIES})

  set(LINKED_LIBRARIES -lgfortran -L/usr/lib/x86_64-linux-gnu/openmpi/lib
    -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lstdc++ -lm)

  set(PSCTOOLKIT_LINKED_LIBRARIES ${LINKED_LIBRARIES} ${LINK_METIS_LIB} ${LINK_MUMPS_LIB} ${LINK_SLU_LIB} ${LINK_SLUDIST_LIB} ${LINK_UMF_LIB} ${LINK_EXTRA_LIB} ${LINK_BLAS})

endif (PSBLAS_FOUND AND AMG_FOUND AND SUNDIALS_FOUND)

find_package_handle_standard_args(PSCTOOLKIT DEFAULT_MSG PSCTOOLKIT_LIBRARIES PSCTOOLKIT_INCLUDE_DIR)
