# Ceres Solver - A fast non-linear least squares minimizer
# Copyright 2013 Google Inc. All rights reserved.
# http://code.google.com/p/ceres-solver/
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice,
#   this list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
# * Neither the name of Google Inc. nor the names of its contributors may be
#   used to endorse or promote products derived from this software without
#   specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
# Author: alexs.mac@gmail.com (Alex Stewart)
#
# FindSuiteSparse.cmake - Find SuiteSparse libraries & dependencies.
#
# This module defines the following variables:
#
# SUITESPARSE_FOUND: TRUE iff SuiteSparse and all dependencies have been found.
# SUITESPARSE_INCLUDE_DIRS: Include directories for all SuiteSparse components.
# SUITESPARSE_LIBRARIES: Libraries for all SuiteSparse component libraries and
#                        dependencies.
# SUITESPARSE_VERSION: Extracted from UFconfig.h (<= v3) or
#                      SuiteSparse_config.h (>= v4).
# SUITESPARSE_MAIN_VERSION: Equal to 4 if SUITESPARSE_VERSION = 4.2.1
# SUITESPARSE_SUB_VERSION: Equal to 2 if SUITESPARSE_VERSION = 4.2.1
# SUITESPARSE_SUBSUB_VERSION: Equal to 1 if SUITESPARSE_VERSION = 4.2.1
#
# SUITESPARSE_IS_BROKEN_SHARED_LINKING_UBUNTU_SYSTEM_VERSION: TRUE iff running
#     on Ubuntu, SUITESPARSE_VERSION is 3.4.0 and found SuiteSparse is a system
#     install, in which case found version of SuiteSparse cannot be used to link
#     a shared library due to a bug (static linking is unaffected).
#
# The following variables control the behaviour of this module:
#
# SUITESPARSE_INCLUDE_DIR_HINTS: List of additional directories in which to
#                                search for SuiteSparse includes,
#                                e.g: /timbuktu/include.
# SUITESPARSE_LIBRARY_DIR_HINTS: List of additional directories in which to
#                                search for SuiteSparse libraries,
#                                e.g: /timbuktu/lib.
#
# The following variables define the presence / includes & libraries for the
# SuiteSparse components searched for, the SUITESPARSE_XX variables are the
# union of the variables for all components.
#
# == Symmetric Approximate Minimum Degree (AMD)
# AMD_FOUND
# AMD_INCLUDE_DIR
# AMD_LIBRARY
#
# == Constrained Approximate Minimum Degree (CAMD)
# CAMD_FOUND
# CAMD_INCLUDE_DIR
# CAMD_LIBRARY
#
# == Column Approximate Minimum Degree (COLAMD)
# COLAMD_FOUND
# COLAMD_INCLUDE_DIR
# COLAMD_LIBRARY
#
# Constrained Column Approximate Minimum Degree (CCOLAMD)
# CCOLAMD_FOUND
# CCOLAMD_INCLUDE_DIR
# CCOLAMD_LIBRARY
#
# == Sparse Supernodal Cholesky Factorization and Update/Downdate (CHOLMOD)
# CHOLMOD_FOUND
# CHOLMOD_INCLUDE_DIR
# CHOLMOD_LIBRARY
#
# == Multifrontal LU Factorization (UMFPACK)
# UMFPACK_FOUND
# UMFPACK_INCLUDE_DIR
# UMFPACK_LIBRARY
#
# == Multifrontal Sparse QR (SuiteSparseQR)
# SUITESPARSEQR_FOUND
# SUITESPARSEQR_INCLUDE_DIR
# SUITESPARSEQR_LIBRARY
#
# == Common configuration for all but CSparse (SuiteSparse version >= 4).
# SUITESPARSE_CONFIG_FOUND
# SUITESPARSE_CONFIG_INCLUDE_DIR
# SUITESPARSE_CONFIG_LIBRARY
#
# == Common configuration for all but CSparse (SuiteSparse version < 4).
# UFCONFIG_FOUND
# UFCONFIG_INCLUDE_DIR
#
# Optional SuiteSparse Dependencies:
#
# == Serial Graph Partitioning and Fill-reducing Matrix Ordering (METIS)
# METIS_FOUND
# METIS_LIBRARY
#
# == Intel Thread Building Blocks (TBB)
# TBB_FOUND
# TBB_LIBRARIES
# Called if we failed to find SuiteSparse or any of it's required dependencies,
# unsets all public (designed to be used externally) variables and reports
# error message at priority depending upon [REQUIRED/QUIET/<NONE>] argument.
macro(SUITESPARSE_REPORT_NOT_FOUND REASON_MSG)
    unset(SUITESPARSE_FOUND)
    unset(SUITESPARSE_INCLUDE_DIRS)
    unset(SUITESPARSE_LIBRARIES)
    unset(SUITESPARSE_VERSION)
    unset(SUITESPARSE_MAIN_VERSION)
    unset(SUITESPARSE_SUB_VERSION)
    unset(SUITESPARSE_SUBSUB_VERSION)
    # Note <package>_FIND_[REQUIRED/QUIETLY] variables defined by FindPackage()
    # use the camelcase library name, not uppercase.
    if(SuiteSparse_FIND_QUIETLY)
        message(STATUS "Failed to find SuiteSparse - " ${REASON_MSG} ${ARGN})
    elseif(SuiteSparse_FIND_REQUIRED)
        message(FATAL_ERROR "Failed to find SuiteSparse - " ${REASON_MSG} ${ARGN})
    else()
        # Neither QUIETLY nor REQUIRED, use WARNING which emits a message
        # but continues configuration and allows generation.
        message(WARNING "Failed to find SuiteSparse - " ${REASON_MSG} ${ARGN})
    endif(SuiteSparse_FIND_QUIETLY)
endmacro(SUITESPARSE_REPORT_NOT_FOUND)
# Specify search directories for include files and libraries (this is the union
# of the search directories for all OSs).  Search user-specified hint
# directories first if supplied, and search user-installed locations first
# so that we prefer user installs to system installs where both exist.
list(
    APPEND SUITESPARSE_CHECK_INCLUDE_DIRS
    ${SUITESPARSE_INCLUDE_DIR_HINTS}
    /opt/local/include
    /opt/local/include/ufsparse # Mac OS X
    /opt/homebrew/include # Mac OS X apple silicon
    /opt/homebrew/include/suitesparse # Mac OS X apple silicon
    /usr/local/homebrew/include # Mac OS X
    /usr/local/include
    /usr/local/include/suitesparse
    /usr/include/suitesparse # Ubuntu
    /usr/include
)
list(
    APPEND SUITESPARSE_CHECK_LIBRARY_DIRS
    ${SUITESPARSE_LIBRARY_DIR_HINTS}
    /opt/local/lib
    /opt/local/lib/ufsparse # Mac OS X
    /opt/homebrew/lib # Mac OS X apple silicon
    /usr/local/homebrew/lib # Mac OS X
    /usr/local/lib
    /usr/local/lib/suitesparse
    /usr/lib/suitesparse # Ubuntu
    /usr/lib
)
# BLAS.
find_package(BLAS QUIET)
if(NOT BLAS_FOUND)
    SUITESPARSE_REPORT_NOT_FOUND("Did not find BLAS library (required for SuiteSparse).")
endif(NOT BLAS_FOUND)
# LAPACK.
find_package(LAPACK QUIET)
if(NOT LAPACK_FOUND)
    SUITESPARSE_REPORT_NOT_FOUND("Did not find LAPACK library (required for SuiteSparse).")
endif(NOT LAPACK_FOUND)
# AMD.
set(AMD_FOUND TRUE)
find_library(AMD_LIBRARY NAMES amd PATHS ${SUITESPARSE_CHECK_LIBRARY_DIRS})
if(EXISTS ${AMD_LIBRARY})
    message(STATUS "Found AMD library: ${AMD_LIBRARY}")
else(EXISTS ${AMD_LIBRARY})
    SUITESPARSE_REPORT_NOT_FOUND("Did not find AMD library.")
    set(AMD_FOUND FALSE)
endif(EXISTS ${AMD_LIBRARY})
mark_as_advanced(AMD_LIBRARY)
find_path(AMD_INCLUDE_DIR NAMES amd.h PATHS ${SUITESPARSE_CHECK_INCLUDE_DIRS})
if(EXISTS ${AMD_INCLUDE_DIR})
    message(STATUS "Found AMD header in: ${AMD_INCLUDE_DIR}")
else(EXISTS ${AMD_INCLUDE_DIR})
    SUITESPARSE_REPORT_NOT_FOUND("Did not find AMD header.")
    set(AMD_FOUND FALSE)
endif(EXISTS ${AMD_INCLUDE_DIR})
mark_as_advanced(AMD_INCLUDE_DIR)
# CAMD.
set(CAMD_FOUND TRUE)
find_library(CAMD_LIBRARY NAMES camd PATHS ${SUITESPARSE_CHECK_LIBRARY_DIRS})
if(EXISTS ${CAMD_LIBRARY})
    message(STATUS "Found CAMD library: ${CAMD_LIBRARY}")
else(EXISTS ${CAMD_LIBRARY})
    SUITESPARSE_REPORT_NOT_FOUND("Did not find CAMD library.")
    set(CAMD_FOUND FALSE)
endif(EXISTS ${CAMD_LIBRARY})
mark_as_advanced(CAMD_LIBRARY)
find_path(CAMD_INCLUDE_DIR NAMES camd.h PATHS ${SUITESPARSE_CHECK_INCLUDE_DIRS})
if(EXISTS ${CAMD_INCLUDE_DIR})
    message(STATUS "Found CAMD header in: ${CAMD_INCLUDE_DIR}")
else(EXISTS ${CAMD_INCLUDE_DIR})
    SUITESPARSE_REPORT_NOT_FOUND("Did not find CAMD header.")
    set(CAMD_FOUND FALSE)
endif(EXISTS ${CAMD_INCLUDE_DIR})
mark_as_advanced(CAMD_INCLUDE_DIR)
# COLAMD.
set(COLAMD_FOUND TRUE)
find_library(COLAMD_LIBRARY NAMES colamd PATHS ${SUITESPARSE_CHECK_LIBRARY_DIRS})
if(EXISTS ${COLAMD_LIBRARY})
    message(STATUS "Found COLAMD library: ${COLAMD_LIBRARY}")
else(EXISTS ${COLAMD_LIBRARY})
    SUITESPARSE_REPORT_NOT_FOUND("Did not find COLAMD library.")
    set(COLAMD_FOUND FALSE)
endif(EXISTS ${COLAMD_LIBRARY})
mark_as_advanced(COLAMD_LIBRARY)
find_path(COLAMD_INCLUDE_DIR NAMES colamd.h PATHS ${SUITESPARSE_CHECK_INCLUDE_DIRS})
if(EXISTS ${COLAMD_INCLUDE_DIR})
    message(STATUS "Found COLAMD header in: ${COLAMD_INCLUDE_DIR}")
else(EXISTS ${COLAMD_INCLUDE_DIR})
    SUITESPARSE_REPORT_NOT_FOUND("Did not find COLAMD header.")
    set(COLAMD_FOUND FALSE)
endif(EXISTS ${COLAMD_INCLUDE_DIR})
mark_as_advanced(COLAMD_INCLUDE_DIR)
# CCOLAMD.
set(CCOLAMD_FOUND TRUE)
find_library(CCOLAMD_LIBRARY NAMES ccolamd PATHS ${SUITESPARSE_CHECK_LIBRARY_DIRS})
if(EXISTS ${CCOLAMD_LIBRARY})
    message(STATUS "Found CCOLAMD library: ${CCOLAMD_LIBRARY}")
else(EXISTS ${CCOLAMD_LIBRARY})
    SUITESPARSE_REPORT_NOT_FOUND("Did not find CCOLAMD library.")
    set(CCOLAMD_FOUND FALSE)
endif(EXISTS ${CCOLAMD_LIBRARY})
mark_as_advanced(CCOLAMD_LIBRARY)
find_path(CCOLAMD_INCLUDE_DIR NAMES ccolamd.h PATHS ${SUITESPARSE_CHECK_INCLUDE_DIRS})
if(EXISTS ${CCOLAMD_INCLUDE_DIR})
    message(STATUS "Found CCOLAMD header in: ${CCOLAMD_INCLUDE_DIR}")
else(EXISTS ${CCOLAMD_INCLUDE_DIR})
    SUITESPARSE_REPORT_NOT_FOUND("Did not find CCOLAMD header.")
    set(CCOLAMD_FOUND FALSE)
endif(EXISTS ${CCOLAMD_INCLUDE_DIR})
mark_as_advanced(CCOLAMD_INCLUDE_DIR)
# CHOLMOD.
set(CHOLMOD_FOUND TRUE)
find_library(CHOLMOD_LIBRARY NAMES cholmod PATHS ${SUITESPARSE_CHECK_LIBRARY_DIRS})
if(EXISTS ${CHOLMOD_LIBRARY})
    message(STATUS "Found CHOLMOD library: ${CHOLMOD_LIBRARY}")
else(EXISTS ${CHOLMOD_LIBRARY})
    SUITESPARSE_REPORT_NOT_FOUND("Did not find CHOLMOD library.")
    set(CHOLMOD_FOUND FALSE)
endif(EXISTS ${CHOLMOD_LIBRARY})
mark_as_advanced(CHOLMOD_LIBRARY)
find_path(CHOLMOD_INCLUDE_DIR NAMES cholmod.h PATHS ${SUITESPARSE_CHECK_INCLUDE_DIRS})
if(EXISTS ${CHOLMOD_INCLUDE_DIR})
    message(STATUS "Found CHOLMOD header in: ${CHOLMOD_INCLUDE_DIR}")
else(EXISTS ${CHOLMOD_INCLUDE_DIR})
    SUITESPARSE_REPORT_NOT_FOUND("Did not find CHOLMOD header.")
    set(CHOLMOD_FOUND FALSE)
endif(EXISTS ${CHOLMOD_INCLUDE_DIR})
mark_as_advanced(CHOLMOD_INCLUDE_DIR)
# UMFPACK.
set(UMFPACK_FOUND TRUE)
find_library(UMFPACK_LIBRARY NAMES umfpack PATHS ${SUITESPARSE_CHECK_LIBRARY_DIRS})
if(EXISTS ${UMFPACK_LIBRARY})
    message(STATUS "Found UMFPACK library: ${UMFPACK_LIBRARY}")
else(EXISTS ${UMFPACK_LIBRARY})
    SUITESPARSE_REPORT_NOT_FOUND("Did not find UMFPACK library.")
    set(UMFPACK_FOUND FALSE)
endif(EXISTS ${UMFPACK_LIBRARY})
mark_as_advanced(UMFPACK_LIBRARY)
find_path(UMFPACK_INCLUDE_DIR NAMES umfpack.h PATHS ${SUITESPARSE_CHECK_INCLUDE_DIRS})
if(EXISTS ${UMFPACK_INCLUDE_DIR})
    message(STATUS "Found UMFPACK header in: ${UMFPACK_INCLUDE_DIR}")
else(EXISTS ${UMFPACK_INCLUDE_DIR})
    SUITESPARSE_REPORT_NOT_FOUND("Did not find UMFPACK header.")
    set(UMFPACK_FOUND FALSE)
endif(EXISTS ${UMFPACK_INCLUDE_DIR})
mark_as_advanced(UMFPACK_INCLUDE_DIR)
# SuiteSparseQR.
set(SUITESPARSEQR_FOUND TRUE)
find_library(SUITESPARSEQR_LIBRARY NAMES spqr PATHS ${SUITESPARSE_CHECK_LIBRARY_DIRS})
if(EXISTS ${SUITESPARSEQR_LIBRARY})
    message(STATUS "Found SuiteSparseQR library: ${SUITESPARSEQR_LIBRARY}")
else(EXISTS ${SUITESPARSEQR_LIBRARY})
    SUITESPARSE_REPORT_NOT_FOUND("Did not find SUITESPARSEQR library.")
    set(SUITESPARSEQR_FOUND FALSE)
endif(EXISTS ${SUITESPARSEQR_LIBRARY})
mark_as_advanced(SUITESPARSEQR_LIBRARY)
find_path(SUITESPARSEQR_INCLUDE_DIR NAMES SuiteSparseQR.hpp PATHS ${SUITESPARSE_CHECK_INCLUDE_DIRS})
if(EXISTS ${SUITESPARSEQR_INCLUDE_DIR})
    message(STATUS "Found SuiteSparseQR header in: ${SUITESPARSEQR_INCLUDE_DIR}")
else(EXISTS ${SUITESPARSEQR_INCLUDE_DIR})
    SUITESPARSE_REPORT_NOT_FOUND("Did not find SUITESPARSEQR header.")
    set(SUITESPARSEQR_FOUND FALSE)
endif(EXISTS ${SUITESPARSEQR_INCLUDE_DIR})
mark_as_advanced(SUITESPARSEQR_INCLUDE_DIR)
if(SUITESPARSEQR_FOUND)
    # SuiteSparseQR may be compiled with Intel Threading Building Blocks,
    # we assume that if TBB is installed, SuiteSparseQR was compiled with
    # support for it, this will do no harm if it wasn't.
    set(TBB_FOUND TRUE)
    find_library(TBB_LIBRARIES NAMES tbb PATHS ${SUITESPARSE_CHECK_LIBRARY_DIRS})
    if(EXISTS ${TBB_LIBRARIES})
        message(
            STATUS
            "Found Intel Thread Building Blocks (TBB) library: "
            "${TBB_LIBRARIES}, assuming SuiteSparseQR was compiled with TBB."
        )
    else(EXISTS ${TBB_LIBRARIES})
        message(
            STATUS
            "Did not find Intel TBB library, assuming SuiteSparseQR was "
            "not compiled with TBB."
        )
        set(TBB_FOUND FALSE)
    endif(EXISTS ${TBB_LIBRARIES})
    mark_as_advanced(TBB_LIBRARIES)
    if(TBB_FOUND)
        find_library(TBB_MALLOC_LIB NAMES tbbmalloc PATHS ${SUITESPARSE_CHECK_LIBRARY_DIRS})
        if(EXISTS ${TBB_MALLOC_LIB})
            message(
                STATUS
                "Found Intel Thread Building Blocks (TBB) Malloc library: "
                "${TBB_MALLOC_LIB}"
            )
            # Append TBB malloc library to TBB libraries list whilst retaining
            # any CMake generated help string (cache variable).
            list(APPEND TBB_LIBRARIES ${TBB_MALLOC_LIB})
            get_property(HELP_STRING CACHE TBB_LIBRARIES PROPERTY HELPSTRING)
            set(TBB_LIBRARIES "${TBB_LIBRARIES}" CACHE STRING ${HELP_STRING})
            # Add the TBB libraries to the SuiteSparseQR libraries (the only
            # libraries to optionally depend on TBB).
            list(APPEND SUITESPARSEQR_LIBRARY ${TBB_LIBRARIES})
        else(EXISTS ${TBB_MALLOC_LIB})
            # If we cannot find all required TBB components do not include it as
            # a dependency.
            message(
                STATUS
                "Did not find Intel Thread Building Blocks (TBB) Malloc "
                "Library, discarding TBB as a dependency."
            )
            set(TBB_FOUND FALSE)
        endif(EXISTS ${TBB_MALLOC_LIB})
        mark_as_advanced(TBB_MALLOC_LIB)
    endif(TBB_FOUND)
endif(SUITESPARSEQR_FOUND)
# UFconfig / SuiteSparse_config.
#
# If SuiteSparse version is >= 4 then SuiteSparse_config is required.
# For SuiteSparse 3, UFconfig.h is required.
set(SUITESPARSE_CONFIG_FOUND TRUE)
set(UFCONFIG_FOUND TRUE)
find_library(
    SUITESPARSE_CONFIG_LIBRARY
    NAMES suitesparseconfig
    PATHS ${SUITESPARSE_CHECK_LIBRARY_DIRS}
)
if(EXISTS ${SUITESPARSE_CONFIG_LIBRARY})
    message(STATUS "Found SuiteSparse_config library: " "${SUITESPARSE_CONFIG_LIBRARY}")
else(EXISTS ${SUITESPARSE_CONFIG_LIBRARY})
    message(
        STATUS
        "Did not find SuiteSparse_config library (should be present "
        "in SuiteSparse >= v4 installs)."
    )
endif(EXISTS ${SUITESPARSE_CONFIG_LIBRARY})
mark_as_advanced(SUITESPARSE_CONFIG_LIBRARY)
find_path(
    SUITESPARSE_CONFIG_INCLUDE_DIR
    NAMES SuiteSparse_config.h
    PATHS ${SUITESPARSE_CHECK_INCLUDE_DIRS}
)
if(EXISTS ${SUITESPARSE_CONFIG_INCLUDE_DIR})
    message(STATUS "Found SuiteSparse_config header in: " "${SUITESPARSE_CONFIG_INCLUDE_DIR}")
    set(UFCONFIG_FOUND FALSE)
else(EXISTS ${SUITESPARSE_CONFIG_INCLUDE_DIR})
    message(
        STATUS
        "Did not find SuiteSparse_config header (should be present "
        "in SuiteSparse >= v4 installs)."
    )
endif(EXISTS ${SUITESPARSE_CONFIG_INCLUDE_DIR})
mark_as_advanced(SUITESPARSE_CONFIG_INCLUDE_DIR)
if(EXISTS ${SUITESPARSE_CONFIG_LIBRARY} AND EXISTS ${SUITESPARSE_CONFIG_INCLUDE_DIR})
    # SuiteSparse_config (SuiteSparse version >= 4) requires librt library for
    # timing by default when compiled on Linux or Unix, but not on OSX (which
    # does not have librt).
    if(CMAKE_SYSTEM_NAME MATCHES "Linux" OR UNIX AND NOT APPLE)
        find_library(LIBRT_LIBRARY NAMES rt PATHS ${SUITESPARSE_CHECK_LIBRARY_DIRS})
        if(LIBRT_LIBRARY)
            message(
                STATUS
                "Adding librt: ${LIBRT_LIBRARY} to "
                "SuiteSparse_config libraries (required on Linux & Unix [not OSX] if "
                "SuiteSparse is compiled with timing)."
            )
        else(LIBRT_LIBRARY)
            message(
                STATUS
                "Could not find librt, but found SuiteSparse_config, "
                "assuming that SuiteSparse was compiled without timing."
            )
        endif(LIBRT_LIBRARY)
        mark_as_advanced(LIBRT_LIBRARY)
        list(APPEND SUITESPARSE_CONFIG_LIBRARY ${LIBRT_LIBRARY})
    endif(CMAKE_SYSTEM_NAME MATCHES "Linux" OR UNIX AND NOT APPLE)
else(EXISTS ${SUITESPARSE_CONFIG_LIBRARY} AND EXISTS ${SUITESPARSE_CONFIG_INCLUDE_DIR})
    set(SUITESPARSE_CONFIG_FOUND FALSE)
    find_path(UFCONFIG_INCLUDE_DIR NAMES UFconfig.h PATHS ${SUITESPARSE_CHECK_INCLUDE_DIRS})
    if(EXISTS ${UFCONFIG_INCLUDE_DIR})
        message(STATUS "Found UFconfig header in: ${UFCONFIG_INCLUDE_DIR}")
    else(EXISTS ${UFCONFIG_INCLUDE_DIR})
        message(
            STATUS
            "Did not find UFconfig header (should be present "
            "in SuiteSparse < v4 installs)"
        )
        set(UFCONFIG_FOUND FALSE)
    endif(EXISTS ${UFCONFIG_INCLUDE_DIR})
    mark_as_advanced(UFCONFIG_INCLUDE_DIR)
endif(EXISTS ${SUITESPARSE_CONFIG_LIBRARY} AND EXISTS ${SUITESPARSE_CONFIG_INCLUDE_DIR})
# Extract the SuiteSparse version from the appropriate header (UFconfig.h for
# <= v3, SuiteSparse_config.h for >= v4).
if(UFCONFIG_FOUND)
    # SuiteSparse version <= 3.
    set(SUITESPARSE_VERSION_FILE ${UFCONFIG_INCLUDE_DIR}/UFconfig.h)
    if(NOT EXISTS ${SUITESPARSE_VERSION_FILE})
        SUITESPARSE_REPORT_NOT_FOUND(
            "Could not find file: ${SUITESPARSE_VERSION_FILE} containing version "
            "information for <= v3 SuiteSparse installs, but UFconfig was found "
            "(only present in <= v3 installs)."
        )
    else(NOT EXISTS ${SUITESPARSE_VERSION_FILE})
        file(READ ${SUITESPARSE_VERSION_FILE} UFCONFIG_CONTENTS)
        string(
            REGEX MATCH "#define SUITESPARSE_MAIN_VERSION [0-9]+"
            SUITESPARSE_MAIN_VERSION
            "${UFCONFIG_CONTENTS}"
        )
        string(
            REGEX REPLACE "#define SUITESPARSE_MAIN_VERSION ([0-9]+)"
            "\\1"
            SUITESPARSE_MAIN_VERSION
            "${SUITESPARSE_MAIN_VERSION}"
        )
        string(
            REGEX MATCH "#define SUITESPARSE_SUB_VERSION [0-9]+"
            SUITESPARSE_SUB_VERSION
            "${UFCONFIG_CONTENTS}"
        )
        string(
            REGEX REPLACE "#define SUITESPARSE_SUB_VERSION ([0-9]+)"
            "\\1"
            SUITESPARSE_SUB_VERSION
            "${SUITESPARSE_SUB_VERSION}"
        )
        string(
            REGEX MATCH "#define SUITESPARSE_SUBSUB_VERSION [0-9]+"
            SUITESPARSE_SUBSUB_VERSION
            "${UFCONFIG_CONTENTS}"
        )
        string(
            REGEX REPLACE "#define SUITESPARSE_SUBSUB_VERSION ([0-9]+)"
            "\\1"
            SUITESPARSE_SUBSUB_VERSION
            "${SUITESPARSE_SUBSUB_VERSION}"
        )
        # This is on a single line s/t CMake does not interpret it as a list of
        # elements and insert ';' separators which would result in 4.;2.;1 nonsense.
        set(SUITESPARSE_VERSION
            "${SUITESPARSE_MAIN_VERSION}.${SUITESPARSE_SUB_VERSION}.${SUITESPARSE_SUBSUB_VERSION}"
        )
    endif(NOT EXISTS ${SUITESPARSE_VERSION_FILE})
endif(UFCONFIG_FOUND)
if(SUITESPARSE_CONFIG_FOUND)
    # SuiteSparse version >= 4.
    set(SUITESPARSE_VERSION_FILE ${SUITESPARSE_CONFIG_INCLUDE_DIR}/SuiteSparse_config.h)
    if(NOT EXISTS ${SUITESPARSE_VERSION_FILE})
        SUITESPARSE_REPORT_NOT_FOUND(
            "Could not find file: ${SUITESPARSE_VERSION_FILE} containing version "
            "information for >= v4 SuiteSparse installs, but SuiteSparse_config was "
            "found (only present in >= v4 installs)."
        )
    else(NOT EXISTS ${SUITESPARSE_VERSION_FILE})
        file(READ ${SUITESPARSE_VERSION_FILE} SUITESPARSE_CONFIG_CONTENTS)
        string(
            REGEX MATCH "#define SUITESPARSE_MAIN_VERSION [0-9]+"
            SUITESPARSE_MAIN_VERSION
            "${SUITESPARSE_CONFIG_CONTENTS}"
        )
        string(
            REGEX REPLACE "#define SUITESPARSE_MAIN_VERSION ([0-9]+)"
            "\\1"
            SUITESPARSE_MAIN_VERSION
            "${SUITESPARSE_MAIN_VERSION}"
        )
        string(
            REGEX MATCH "#define SUITESPARSE_SUB_VERSION [0-9]+"
            SUITESPARSE_SUB_VERSION
            "${SUITESPARSE_CONFIG_CONTENTS}"
        )
        string(
            REGEX REPLACE "#define SUITESPARSE_SUB_VERSION ([0-9]+)"
            "\\1"
            SUITESPARSE_SUB_VERSION
            "${SUITESPARSE_SUB_VERSION}"
        )
        string(
            REGEX MATCH "#define SUITESPARSE_SUBSUB_VERSION [0-9]+"
            SUITESPARSE_SUBSUB_VERSION
            "${SUITESPARSE_CONFIG_CONTENTS}"
        )
        string(
            REGEX REPLACE "#define SUITESPARSE_SUBSUB_VERSION ([0-9]+)"
            "\\1"
            SUITESPARSE_SUBSUB_VERSION
            "${SUITESPARSE_SUBSUB_VERSION}"
        )
        # This is on a single line s/t CMake does not interpret it as a list of
        # elements and insert ';' separators which would result in 4.;2.;1 nonsense.
        set(SUITESPARSE_VERSION
            "${SUITESPARSE_MAIN_VERSION}.${SUITESPARSE_SUB_VERSION}.${SUITESPARSE_SUBSUB_VERSION}"
        )
    endif(NOT EXISTS ${SUITESPARSE_VERSION_FILE})
endif(SUITESPARSE_CONFIG_FOUND)
# METIS (Optional dependency).
find_library(METIS_LIBRARY NAMES libmetis.so.5 PATHS ${SUITESPARSE_CHECK_LIBRARY_DIRS})
if(EXISTS ${METIS_LIBRARY})
    message(STATUS "Found METIS library: ${METIS_LIBRARY}.")
else(EXISTS ${METIS_LIBRARY})
    message(STATUS "Did not find METIS library (optional SuiteSparse dependency)")
endif(EXISTS ${METIS_LIBRARY})
mark_as_advanced(METIS_LIBRARY)
# Only mark SuiteSparse as found if all required dependencies have been found.
set(SUITESPARSE_FOUND FALSE)
if(
    AMD_FOUND
    AND CAMD_FOUND
    AND COLAMD_FOUND
    AND CCOLAMD_FOUND
    AND CHOLMOD_FOUND
    AND UMFPACK_FOUND
    AND SUITESPARSEQR_FOUND
    AND (SUITESPARSE_CONFIG_FOUND OR UFCONFIG_FOUND)
    AND BLAS_FOUND
    AND LAPACK_FOUND
)
    set(SUITESPARSE_FOUND TRUE)
    list(
        APPEND SUITESPARSE_INCLUDE_DIRS
        ${AMD_INCLUDE_DIR}
        ${CAMD_INCLUDE_DIR}
        ${COLAMD_INCLUDE_DIR}
        ${CCOLAMD_INCLUDE_DIR}
        ${CHOLMOD_INCLUDE_DIR}
        ${UMFPACK_INCLUDE_DIR}
        ${SUITESPARSEQR_INCLUDE_DIR}
    )
    # Handle config separately, as otherwise at least one of them will be set
    # to NOTFOUND which would cause any check on SUITESPARSE_INCLUDE_DIRS to fail.
    if(SUITESPARSE_CONFIG_FOUND)
        list(APPEND SUITESPARSE_INCLUDE_DIRS ${SUITESPARSE_CONFIG_INCLUDE_DIR})
    endif(SUITESPARSE_CONFIG_FOUND)
    if(UFCONFIG_FOUND)
        list(APPEND SUITESPARSE_INCLUDE_DIRS ${UFCONFIG_INCLUDE_DIR})
    endif(UFCONFIG_FOUND)
    # As SuiteSparse includes are often all in the same directory, remove any
    # repetitions.
    list(REMOVE_DUPLICATES SUITESPARSE_INCLUDE_DIRS)
    # Important: The ordering of these libraries is *NOT* arbitrary, as these
    # could potentially be static libraries their link ordering is important.
    list(
        APPEND SUITESPARSE_LIBRARIES
        ${SUITESPARSEQR_LIBRARY}
        ${UMFPACK_LIBRARY}
        ${CHOLMOD_LIBRARY}
        ${CCOLAMD_LIBRARY}
        ${CAMD_LIBRARY}
        ${COLAMD_LIBRARY}
        ${AMD_LIBRARY}
    )
    if(SUITESPARSE_CONFIG_FOUND)
        list(APPEND SUITESPARSE_LIBRARIES ${SUITESPARSE_CONFIG_LIBRARY})
    endif(SUITESPARSE_CONFIG_FOUND)
    if(METIS_FOUND)
        list(APPEND SUITESPARSE_LIBRARIES ${METIS_LIBRARY})
    endif(METIS_FOUND)
else()
    SUITESPARSE_REPORT_NOT_FOUND("Failed to find some/all required components of SuiteSparse.")
endif()
# Determine if we are running on Ubuntu with the package install of SuiteSparse
# which is broken and does not support linking a shared library.
set(SUITESPARSE_IS_BROKEN_SHARED_LINKING_UBUNTU_SYSTEM_VERSION FALSE)
if(CMAKE_SYSTEM_NAME MATCHES "Linux" AND SUITESPARSE_VERSION VERSION_EQUAL 3.4.0)
    find_program(LSB_RELEASE_EXECUTABLE lsb_release)
    if(LSB_RELEASE_EXECUTABLE)
        # Any even moderately recent Ubuntu release (likely to be affected by
        # this bug) should have lsb_release, if it isn't present we are likely
        # on a different Linux distribution (should be fine).
        execute_process(
            COMMAND ${LSB_RELEASE_EXECUTABLE} -si
            OUTPUT_VARIABLE LSB_DISTRIBUTOR_ID
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        if(LSB_DISTRIBUTOR_ID MATCHES "Ubuntu" AND SUITESPARSE_LIBRARIES MATCHES "/usr/lib/libamd")
            # We are on Ubuntu, and the SuiteSparse version matches the broken
            # system install version and is a system install.
            set(SUITESPARSE_IS_BROKEN_SHARED_LINKING_UBUNTU_SYSTEM_VERSION TRUE)
            message(
                STATUS
                "Found system install of SuiteSparse "
                "${SUITESPARSE_VERSION} running on Ubuntu, which has a known bug "
                "preventing linking of shared libraries (static linking unaffected)."
            )
        endif(
            LSB_DISTRIBUTOR_ID MATCHES "Ubuntu"
            AND SUITESPARSE_LIBRARIES MATCHES "/usr/lib/libamd"
        )
    endif(LSB_RELEASE_EXECUTABLE)
endif(CMAKE_SYSTEM_NAME MATCHES "Linux" AND SUITESPARSE_VERSION VERSION_EQUAL 3.4.0)
# Handle REQUIRED and QUIET arguments to FIND_PACKAGE
include(FindPackageHandleStandardArgs)
# A change to CMake after release 2.8.10.2 means that
# FindPackageHandleStandardArgs() unsets <LibraryName>_FOUND without checking
# if it is one of the variables passed whose existence & validity is verified
# by FindPackageHandleStandardArgs() in conjunction with handling the REQUIRED
# and QUIET optional arguments, as such we use an intermediary variable.
set(SUITESPARSE_FOUND_COPY ${SUITESPARSE_FOUND})
find_package_handle_standard_args(
    SuiteSparse
    REQUIRED_VARS SUITESPARSE_INCLUDE_DIRS SUITESPARSE_LIBRARIES SUITESPARSE_FOUND_COPY
    VERSION_VAR SUITESPARSE_VERSION
)
