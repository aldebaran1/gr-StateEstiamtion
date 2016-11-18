# Copyright 2011,2012 Free Software Foundation, Inc.
#
# This file is part of GNU Radio
#
# GNU Radio is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
#
# GNU Radio is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GNU Radio; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street,
# Boston, MA 02110-1301, USA.


########################################################################
# Project setup
########################################################################
cmake_minimum_required(VERSION 2.6)
project(gr-StateEstimation CXX C)
enable_testing()

#select the release build type by default to get optimization flags
if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE "Release")
   message(STATUS "Build type not specified: defaulting to release.")
endif(NOT CMAKE_BUILD_TYPE)
set(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "")

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules)

########################################################################
# Compiler specific setup
########################################################################
if(CMAKE_COMPILER_IS_GNUCXX AND NOT WIN32)
    #http://gcc.gnu.org/wiki/Visibility
    add_definitions(-fvisibility=hidden)
endif()

########################################################################
# Find boost
########################################################################
if(UNIX AND EXISTS "/usr/lib64")
    list(APPEND BOOST_LIBRARYDIR "/usr/lib64") #fedora 64-bit fix
endif(UNIX AND EXISTS "/usr/lib64")
set(Boost_ADDITIONAL_VERSIONS
    "1.35.0" "1.35" "1.36.0" "1.36" "1.37.0" "1.37" "1.38.0" "1.38" "1.39.0" "1.39"
    "1.40.0" "1.40" "1.41.0" "1.41" "1.42.0" "1.42" "1.43.0" "1.43" "1.44.0" "1.44"
    "1.45.0" "1.45" "1.46.0" "1.46" "1.47.0" "1.47" "1.48.0" "1.48" "1.49.0" "1.49"
    "1.50.0" "1.50" "1.51.0" "1.51" "1.52.0" "1.52" "1.53.0" "1.53" "1.54.0" "1.54"
    "1.55.0" "1.55" "1.56.0" "1.56" "1.57.0" "1.57" "1.58.0" "1.58" "1.59.0" "1.59"
    "1.60.0" "1.60" "1.61.0" "1.61" "1.62.0" "1.62" "1.63.0" "1.63" "1.64.0" "1.64"
    "1.65.0" "1.65" "1.66.0" "1.66" "1.67.0" "1.67" "1.68.0" "1.68" "1.69.0" "1.69"
)
find_package(Boost "1.35" COMPONENTS filesystem system)

if(NOT Boost_FOUND)
    message(FATAL_ERROR "Boost required to compile StateEstimation")
endif()

########################################################################
# Install directories
########################################################################
include(GrPlatform) #define LIB_SUFFIX
set(GR_RUNTIME_DIR      bin)
set(GR_LIBRARY_DIR      lib${LIB_SUFFIX})
set(GR_INCLUDE_DIR      include/StateEstimation)
set(GR_DATA_DIR         share)
set(GR_PKG_DATA_DIR     ${GR_DATA_DIR}/${CMAKE_PROJECT_NAME})
set(GR_DOC_DIR          ${GR_DATA_DIR}/doc)
set(GR_PKG_DOC_DIR      ${GR_DOC_DIR}/${CMAKE_PROJECT_NAME})
set(GR_CONF_DIR         etc)
set(GR_PKG_CONF_DIR     ${GR_CONF_DIR}/${CMAKE_PROJECT_NAME}/conf.d)
set(GR_LIBEXEC_DIR      libexec)
set(GR_PKG_LIBEXEC_DIR  ${GR_LIBEXEC_DIR}/${CMAKE_PROJECT_NAME})
set(GRC_BLOCKS_DIR      ${GR_PKG_DATA_DIR}/grc/blocks)

########################################################################
# Find gnuradio build dependencies
########################################################################
find_package(Gruel)
find_package(GnuradioCore)
find_package(CppUnit)

if(NOT GRUEL_FOUND)
    message(FATAL_ERROR "Gruel required to compile StateEstimation")
endif()

if(NOT GNURADIO_CORE_FOUND)
    message(FATAL_ERROR "GnuRadio Core required to compile StateEstimation")
endif()

if(NOT CPPUNIT_FOUND)
    message(FATAL_ERROR "CppUnit required to compile StateEstimation")
endif()

########################################################################
# Setup the include and linker paths
########################################################################
include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${Boost_INCLUDE_DIRS}
    ${GRUEL_INCLUDE_DIRS}
    ${CPPUNIT_INCLUDE_DIRS}
    ${GNURADIO_CORE_INCLUDE_DIRS}
)

link_directories(
    ${Boost_LIBRARY_DIRS}
    ${GRUEL_LIBRARY_DIRS}
    ${CPPUNIT_LIBRARY_DIRS}
    ${GNURADIO_CORE_LIBRARY_DIRS}
)

# Set component parameters
set(GR_STATEESTIMATION_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include CACHE INTERNAL "" FORCE)
set(GR_STATEESTIMATION_SWIG_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/swig CACHE INTERNAL "" FORCE)

########################################################################
# Create uninstall target
########################################################################
configure_file(
    ${CMAKE_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake
@ONLY)

add_custom_target(uninstall
    ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake
)

################################################################################
# Armadillo - http://arma.sourceforge.net/
################################################################################
if(OS_IS_LINUX)
#############################################################################
# Check that LAPACK is found in the system
# LAPACK is required for matrix decompositions (eg. SVD) and matrix inverse.
#############################################################################
find_library(LAPACK lapack)
if(NOT LAPACK)
message(" The LAPACK library has not been found.")
message(" You can try to install it by typing:")
if(${LINUX_DISTRIBUTION} MATCHES "Fedora" OR ${LINUX_DISTRIBUTION} MATCHES "Red Hat")
message(" sudo yum install lapack-devel")
elseif(${LINUX_DISTRIBUTION} MATCHES "openSUSE")
message(" sudo zypper install lapack-devel")
else(${LINUX_DISTRIBUTION} MATCHES "Fedora" OR ${LINUX_DISTRIBUTION} MATCHES "Red Hat")
message(" sudo apt-get install liblapack-dev")
endif(${LINUX_DISTRIBUTION} MATCHES "Fedora" OR ${LINUX_DISTRIBUTION} MATCHES "Red Hat")
message(FATAL_ERROR "LAPACK is required to build gnss-sdr")
endif(NOT LAPACK)
#############################################################################
# Check that BLAS is found in the system
# BLAS is used for matrix multiplication.
# Without BLAS, matrix multiplication will still work, but might be slower.
#############################################################################
find_library(BLAS blas)
if(NOT BLAS)
message(" The BLAS library has not been found.")
message(" You can try to install it by typing:")
if(${LINUX_DISTRIBUTION} MATCHES "Fedora" OR ${LINUX_DISTRIBUTION} MATCHES "Red Hat")
message(" sudo yum install blas-devel")
else(${LINUX_DISTRIBUTION} MATCHES "Fedora" OR ${LINUX_DISTRIBUTION} MATCHES "Red Hat")
message(" sudo apt-get install libopenblas-dev")
endif(${LINUX_DISTRIBUTION} MATCHES "Fedora" OR ${LINUX_DISTRIBUTION} MATCHES "Red Hat")
message(FATAL_ERROR "BLAS is required to build gnss-sdr")
endif(NOT BLAS)
#############################################
# Check if GFORTRAN is found in the system
#############################################
find_library(GFORTRAN NAMES gfortran
PATHS /usr/lib
/usr/local/lib
/usr/local/lib/i386
/usr/lib/gcc/x86_64-linux-gnu
/usr/lib/gcc/i686-linux-gnu
/usr/lib/gcc/i386-linux-gnu
/usr/lib/gcc/x86_64-linux-gnu/4.6 # Ubuntu 12.04
/usr/lib/gcc/i686-linux-gnu/4.6
/usr/lib/gcc/x86_64-linux-gnu/4.7
/usr/lib/gcc/i686-linux-gnu/4.7
/usr/lib/gcc/x86_64-linux-gnu/4.8
/usr/lib/gcc/i686-linux-gnu/4.8
/usr/lib/gcc/x86_64-linux-gnu/4.9
/usr/lib/gcc/i686-linux-gnu/4.9
/usr/lib/gcc/x86_64-redhat-linux/4.7.2 # Fedora 18
/usr/lib/gcc/i686-redhat-linux/4.7.2
/usr/lib/gcc/x86_64-redhat-linux/4.8.1 # Fedora 19
/usr/lib/gcc/x86_64-redhat-linux/4.8.3 # Fedora 20
/usr/lib/gcc/x86_64-redhat-linux/4.9.1 # Fedora 21
/usr/lib/gcc/i686-redhat-linux/4.8.1
/usr/lib/gcc/i686-redhat-linux/4.8.3
/usr/lib/gcc/i686-redhat-linux/4.9.1
/usr/lib/gcc/x86_64-redhat-linux/4.4.4 # CentOS 6
/usr/lib/gcc/i686-redhat-linux/4.4.4
/usr/lib/gcc/x86_64-redhat-linux/4.8.2
/usr/lib/gcc/i686-redhat-linux/4.8.2
/usr/lib/gcc/i586-suse-linux/4.8 # OpenSUSE 13.1
/usr/lib/gcc/i586-suse-linux/4.9
/usr/lib/gcc/x86_64-suse-linux/4.8
/usr/lib/gcc/x86_64-suse-linux/4.9
/usr/lib/gcc/i486-linux-gnu # Debian 7
/usr/lib/gcc/i486-linux-gnu/4.4
/usr/lib/gcc/i486-linux-gnu/4.6
/usr/lib/gcc/i486-linux-gnu/4.7
/usr/lib/gcc/i486-linux-gnu/4.8
/usr/lib/gcc/i486-linux-gnu/4.9
/usr/lib/gcc/i586-linux-gnu/4.9
/usr/lib/gcc/arm-linux-gnueabihf/4.4 # Debian armhf
/usr/lib/gcc/arm-linux-gnueabihf/4.5
/usr/lib/gcc/arm-linux-gnueabihf/4.6
/usr/lib/gcc/arm-linux-gnueabihf/4.7
/usr/lib/gcc/arm-linux-gnueabihf/4.8
/usr/lib/gcc/arm-linux-gnueabihf/4.9
/usr/lib/gcc/aarch64-linux-gnu/4.9 # Debian arm64
/usr/lib/gcc/arm-linux-gnueabi/4.7 # Debian armel
/usr/lib/gcc/arm-linux-gnueabi/4.9
)
if(NOT GFORTRAN)
message(STATUS "The gfortran library has not been found.")
message(STATUS " You can try to install it by typing:")
if(${LINUX_DISTRIBUTION} MATCHES "Fedora" OR ${LINUX_DISTRIBUTION} MATCHES "Red Hat")
message(STATUS " sudo yum install gcc-fortran")
elseif(${LINUX_DISTRIBUTION} MATCHES "openSUSE")
message(STATUS " sudo zypper install gcc-fortran")
else(${LINUX_DISTRIBUTION} MATCHES "Fedora" OR ${LINUX_DISTRIBUTION} MATCHES "Red Hat")
message(STATUS " sudo apt-get install gfortran")
endif(${LINUX_DISTRIBUTION} MATCHES "Fedora" OR ${LINUX_DISTRIBUTION} MATCHES "Red Hat")
endif(NOT GFORTRAN)
endif(OS_IS_LINUX)
find_package(Armadillo)
if(NOT ARMADILLO_FOUND)
message(STATUS " Armadillo has not been found.")
message(STATUS " Armadillo will be downloaded and built automatically ")
message(STATUS " when doing 'make'. ")
set(armadillo_RELEASE 6.200.4)
set(armadillo_MD5 "20c9de6ad06b4339a7fd3f67d78922e4")
ExternalProject_Add(
armadillo-${armadillo_RELEASE}
PREFIX ${CMAKE_CURRENT_BINARY_DIR}/armadillo-${armadillo_RELEASE}
URL http://sourceforge.net/projects/arma/files/armadillo-${armadillo_RELEASE}.tar.gz
DOWNLOAD_DIR ${CMAKE_CURRENT_BINARY_DIR}/download/armadillo-${armadillo_RELEASE}
URL_MD5 ${armadillo_MD5}
CMAKE_ARGS -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DBUILD_SHARED_LIBS=OFF
BUILD_IN_SOURCE 1
BUILD_COMMAND make
UPDATE_COMMAND ""
INSTALL_COMMAND ""
)
# Set up variables
ExternalProject_Get_Property(armadillo-${armadillo_RELEASE} binary_dir)
set(ARMADILLO_INCLUDE_DIRS ${binary_dir}/include )
find_library(LAPACK NAMES lapack HINTS /usr/lapack /usr/lib /usr/local/lib /usr/lib64)
if(OS_IS_MACOSX)
find_library(BLAS blas)
endif(OS_IS_MACOSX)
find_package(OpenBLAS)
if(OPENBLAS_FOUND)
set(BLAS ${OPENBLAS})
endif(OPENBLAS_FOUND)
if(NOT GFORTRAN)
set(GFORTRAN "")
endif(NOT GFORTRAN)
set(ARMADILLO_LIBRARIES ${BLAS} ${LAPACK} ${GFORTRAN} ${binary_dir}/${CMAKE_FIND_LIBRARY_PREFIXES}armadillo.a)
set(LOCAL_ARMADILLO true CACHE STRING "Armadillo downloaded and built automatically" FORCE)
# Save a copy at the thirdparty folder
file(COPY ${CMAKE_CURRENT_BINARY_DIR}/armadillo-${armadillo_RELEASE}
DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/armadillo
)
else(NOT ARMADILLO_FOUND)
set(armadillo_RELEASE ${ARMADILLO_VERSION_STRING})
add_library(armadillo-${armadillo_RELEASE} UNKNOWN IMPORTED)
set_property(TARGET armadillo-${armadillo_RELEASE} PROPERTY IMPORTED_LOCATION "${ARMADILLO_LIBRARIES}")
endif(NOT ARMADILLO_FOUND)

########################################################################
# Add subdirectories
########################################################################
add_subdirectory(include/StateEstimation)
add_subdirectory(lib)
add_subdirectory(swig)
add_subdirectory(python)
add_subdirectory(grc)
add_subdirectory(apps)
add_subdirectory(docs)
