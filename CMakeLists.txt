# minimum cmake version
cmake_minimum_required(VERSION 3.10)

#cmake_policy(SET CMP0074 NEW)

# now start with skyfill
project(skyfill)

# tweak install path for Windows
IF(WIN32)
    # install into place in build-dir
    SET( CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/INSTALL/ CACHE FILEPATH "install prefix" FORCE)
ENDIF()

# find all required packages
# currently only libtiff is required
FIND_PACKAGE(TIFF REQUIRED)
INCLUDE_DIRECTORIES(${TIFF_INCLUDE_DIR})


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
target_link_libraries(skyfill ${TIFF_LIBRARIES} m)

#install(TARGETS skyfill)
