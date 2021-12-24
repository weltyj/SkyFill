# minimum cmake version
cmake_minimum_required(VERSION 3.10)

if(POLICY CMP0074)
  # use new <PackageName>_ROOT in underlying FindPackages call
  cmake_policy(SET CMP0074 NEW)
endif()

# now start with skyfill
project(skyfill)

# tweak install path for Windows
IF(WIN32)
    # install into place in build-dir
    SET( CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/INSTALL/ CACHE FILEPATH "install prefix" FORCE)
ENDIF()

# getting install paths other systems
include(GNUInstallDirs)

# find all required packages
# currently only libtiff is required
FIND_PACKAGE(TIFF REQUIRED)
INCLUDE_DIRECTORIES(${TIFF_INCLUDE_DIR})
SET(NEEDED_LIBS ${TIFF_LIBRARIES})
 
# check for libm
include(CheckLibraryExists)
# Does explicit linking against libm work?
check_library_exists(m exp "" LIBMNEEDED)
if(LIBMNEEDED)
  list(APPEND NEEDED_LIBS "m")
endif()
 

# file list
SET(SKYFILL_HEADER 
  mstat.h
  amoeba_06.h
)

SET(SKYFILL_SOURCE
  amoeba_06.c
  mstat.c
  skyfill_tif.c
)

# now create the executable
add_executable(skyfill ${SKYFILL_HEADER} ${SKYFILL_SOURCE})
target_link_libraries(skyfill ${NEEDED_LIBS})

install(TARGETS skyfill RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
