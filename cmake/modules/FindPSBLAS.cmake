#-----------------------------------------------------------------------------
# - Find PSBLAS includes and libraries.
#
# This module finds if PSBLAS is installed and determines where the
# include files and libraries are.  This code sets the following variables:
#  PSBLAS_FOUND         = PSBLAS was found
#  PSBLAS_INCLUDE_DIR   = path to where header files can be found
#  PSBLAS_LIBRARIES     = link libraries for PSBLAS
#-----------------------------------------------------------------------------

include(FindPackageHandleStandardArgs)


find_path (PSBLAS_DIR include/Make.inc.psblas HINTS ${PSBLAS_ROOT} 
  DOC "PSBLAS Directory")
# find_path (PSBLAS_DIR include/psb_config.h HINTS ${PSBLAS_ROOT} 
  # DOC "PSBLAS Directory")

if (PSBLAS_DIR)

  set(PSBLAS_FOUND YES)

  set(PSBLAS_INCLUDE_DIR ${PSBLAS_DIR}/include)
  set(PSBLAS_LIBRARY_DIR ${PSBLAS_DIR}/lib)
  set(PSBLAS_MODULES_DIR ${PSBLAS_DIR}/modules)

  if(PSBLAS_FIND_COMPONENTS)
    
    foreach(comp ${PSBLAS_FIND_COMPONENTS})

      # Need to make sure variable to search for isn't set
      unset(PSBLAS_LIB CACHE)

      find_library(PSBLAS_LIB
        NAMES ${comp}
        HINTS ${PSBLAS_LIBRARY_DIR}
        NO_DEFAULT_PATH)

      if(PSBLAS_LIB)
        list(APPEND PSBLAS_LIBRARIES ${PSBLAS_LIB})
      else(PSBLAS_LIB)	    
        message(FATAL_ERROR "Could not find required PSBLAS library : ${comp}")
      endif(PSBLAS_LIB)
    
    endforeach(comp)

  endif(PSBLAS_FIND_COMPONENTS)


  # Now we parse the Make.inc file to set the compilation variables
  # We start with the PSBLAS Make.inc.psblas file, this has to be found so we
  # check for nothing
  set(regex "BLAS=.*")
  file(STRINGS ${PSBLAS_INCLUDE_DIR}/Make.inc.psblas LINK_BLAS REGEX "${regex}")
  set(regex "BLAS=")
  string(REGEX REPLACE "${regex}" "" LINK_BLAS "${LINK_BLAS}")
  separate_arguments(LINK_BLAS)

  set(regex "METIS_LIB=.*")
  file(STRINGS ${PSBLAS_INCLUDE_DIR}/Make.inc.psblas LINK_METIS_LIB REGEX "${regex}")
  set(regex "METIS_LIB=")
  string(REGEX REPLACE "${regex}" "" LINK_METIS_LIB "${LINK_METIS_LIB}")
  separate_arguments(LINK_METIS_LIB)

  set(regex "AMD_LIB=.*")
  file(STRINGS ${PSBLAS_INCLUDE_DIR}/Make.inc.psblas LINK_AMD_LIB REGEX "${regex}")
  set(regex "AMD_LIB=")
  string(REGEX REPLACE "${regex}" "" LINK_AMD_LIB "${LINK_AMD_LIB}")
  separate_arguments(LINK_AMD_LIB)

  set(regex "PSBFDEFINES=.*")
  file(STRINGS ${PSBLAS_INCLUDE_DIR}/Make.inc.psblas PSBFDEFINES REGEX "${regex}")
  set(regex "PSBFDEFINES=")
  string(REGEX REPLACE "${regex}" "" PSBFDEFINES "${PSBFDEFINES}")
  set(regex "-D")
  string(REGEX REPLACE "${regex}" "" PSBFDEFINES "${PSBFDEFINES}")
  separate_arguments(PSBFDEFINES)

  set(regex "PSBCDEFINES=.*")
  file(STRINGS ${PSBLAS_INCLUDE_DIR}/Make.inc.psblas PSBCDEFINES REGEX "${regex}")
  set(regex "PSBCDEFINES=")
  string(REGEX REPLACE "${regex}" "" PSBCDEFINES "${PSBCDEFINES}")
  set(regex "-D")
  string(REGEX REPLACE "${regex}" "" PSBCDEFINES "${PSBCDEFINES}")
  separate_arguments(PSBCDEFINES)


else(PSBLAS_DIR)

  set(PSBLAS_FOUND NO)

endif(PSBLAS_DIR)

find_package_handle_standard_args(PSBLAS DEFAULT_MSG PSBLAS_LIBRARIES PSBLAS_INCLUDE_DIR)
