#-----------------------------------------------------------------------------
# - Find AMG includes and libraries.
#
# This module finds if AMG is installed and determines where the
# include files and libraries are.  This code sets the following variables:
#  AMG_FOUND         = AMG was found
#  AMG_INCLUDE_DIR   = path to where header files can be found
#  AMG_LIBRARIES     = link libraries for AMG
#-----------------------------------------------------------------------------

include(FindPackageHandleStandardArgs)

# find_path (AMG_DIR include/amg_config.h HINTS ${AMG_ROOT} 
#   DOC "AMG Directory")
find_path (AMG_DIR include/Make.inc.amg4psblas HINTS ${AMG_ROOT} 
  DOC "AMG Directory")

if (AMG_DIR)

  set(AMG_FOUND YES)

  set(AMG_INCLUDE_DIR ${AMG_DIR}/include)
  set(AMG_LIBRARY_DIR ${AMG_DIR}/lib)
  set(AMG_MODULES_DIR ${AMG_DIR}/modules)

  if(AMG_FIND_COMPONENTS)
    
    foreach(comp ${AMG_FIND_COMPONENTS})

      # Need to make sure variable to search for isn't set
      unset(AMG_LIB CACHE)

      find_library(AMG_LIB
        NAMES ${comp}
        HINTS ${AMG_LIBRARY_DIR}
        NO_DEFAULT_PATH)

      if(AMG_LIB)
        list(APPEND AMG_LIBRARIES ${AMG_LIB})
      else(AMG_LIB)	    
        message(FATAL_ERROR "Could not find required AMG library : ${comp}")
      endif(AMG_LIB)
    
    endforeach(comp)

  endif(AMG_FIND_COMPONENTS)

  # We now parse the AMG4PSBLAS make.inc file, AMG4PSBLAS can be compiled with
  # linking to several other libraries, so we need to collect these information
  set(regex "MUMPSLIBS=.*")
  file(STRINGS ${AMG_INCLUDE_DIR}/Make.inc.amg4psblas LINK_MUMPS_LIB REGEX "${regex}")
  set(regex ".*#;MUMPSLIBS=")
  string(REGEX REPLACE "${regex}" "" LINK_MUMPS_LIB "${LINK_MUMPS_LIB}")
  separate_arguments(LINK_MUMPS_LIB)

  set(regex "MUMPSFLAGS=.*")
  file(STRINGS ${AMG_INCLUDE_DIR}/Make.inc.amg4psblas FLAGS_MUMPS_LIB REGEX "${regex}")
  set(regex ".*#;MUMPSFLAGS=")
  string(REGEX REPLACE "${regex}" "" FLAGS_MUMPS_LIB "${FLAGS_MUMPS_LIB}")
  separate_arguments(FLAGS_MUMPS_LIB)

  set(regex "SLULIBS=.*")
  file(STRINGS ${AMG_INCLUDE_DIR}/Make.inc.amg4psblas LINK_SLU_LIB REGEX "${regex}")
  set(regex ".*#;SLULIBS=")
  string(REGEX REPLACE "${regex}" "" LINK_SLU_LIB "${LINK_SLU_LIB}")
  separate_arguments(LINK_SLU_LIB)

  set(regex "SLUFLAGS=.*")
  file(STRINGS ${AMG_INCLUDE_DIR}/Make.inc.amg4psblas FLAGS_SLU_LIB REGEX "${regex}")
  set(regex ".*#;SLUFLAGS=")
  string(REGEX REPLACE "${regex}" "" FLAGS_SLU_LIB "${FLAGS_SLU_LIB}")
  separate_arguments(FLAGS_SLU_LIB)

  set(regex "SLUDISTLIBS=.*")
  file(STRINGS ${AMG_INCLUDE_DIR}/Make.inc.amg4psblas LINK_SLUDIST_LIB REGEX "${regex}")
  set(regex ".*#;SLUDISTLIBS=")
  string(REGEX REPLACE "${regex}" "" LINK_SLUDIST_LIB "${LINK_SLUDIST_LIB}")
  separate_arguments(LINK_SLUDIST_LIB)

  set(regex "SLUDISTFLAGS=.*")
  file(STRINGS ${AMG_INCLUDE_DIR}/Make.inc.amg4psblas FLAGS_SLUDIST_LIB REGEX "${regex}")
  set(regex ".*#;SLUDISTFLAGS=")
  string(REGEX REPLACE "${regex}" "" FLAGS_SLUDIST_LIB "${FLAGS_SLUDIST_LIB}")
  separate_arguments(FLAGS_SLUDIST_LIB)

  set(regex "UMFLIBS=.*")
  file(STRINGS ${AMG_INCLUDE_DIR}/Make.inc.amg4psblas LINK_UMF_LIB REGEX "${regex}")
  set(regex ".*#;UMFLIBS=")
  string(REGEX REPLACE "${regex}" "" LINK_UMF_LIB "${LINK_UMF_LIB}")
  separate_arguments(LINK_UMF_LIB)

  set(regex "UMFFLAGS=.*")
  file(STRINGS ${AMG_INCLUDE_DIR}/Make.inc.amg4psblas FLAGS_UMF_LIB REGEX "${regex}")
  set(regex ".*#;UMFFLAGS=")
  string(REGEX REPLACE "${regex}" "" FLAGS_UMF_LIB "${FLAGS_UMF_LIB}")
  separate_arguments(FLAGS_UMF_LIB)

  set(regex "EXTRALIBS=.*")
  file(STRINGS ${AMG_INCLUDE_DIR}/Make.inc.amg4psblas LINK_EXTRA_LIB REGEX "${regex}")
  set(regex "EXTRALIBS=")
  string(REGEX REPLACE "${regex}" "" LINK_EXTRA_LIB "${LINK_EXTRA_LIB}")
  separate_arguments(LINK_EXTRA_LIB)

  set(AMGCDEFINES "${FLAGS_MUMPS_LIB} ${FLAGS_SLU_LIB} ${FLAGS_SLUDIST_LIB} ${FLAGS_UMF_LIB}")
  separate_arguments(AMGCDEFINES)

else(AMG_DIR)

  set(AMG_FOUND NO)

endif(AMG_DIR)

find_package_handle_standard_args(AMG DEFAULT_MSG AMG_LIBRARIES AMG_INCLUDE_DIR)
