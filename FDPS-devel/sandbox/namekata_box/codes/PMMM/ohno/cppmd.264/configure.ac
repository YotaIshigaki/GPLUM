# Copyright (c) 2008, 2009, 2010, 2011, 2012 RIKEN. All Rights Reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
# 
#    1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#    2. Redistributions in binary form must reproduce the above
#    copyright notice, this list of conditions and the following
#    disclaimer in the documentation and/or other materials provided
#    with the distribution.
#    3. The name of the author may not be used to endorse or promote
#    products derived from this software without specific prior written
#    permission.
# 
# 
# THIS SOFTWARE IS PROVIDED BY RIKEN ``AS IS'' AND ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
# IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

AC_INIT([cppmd], [0.3.0], [cppmd-maintainer@riken.jp], [cppmd])
AC_PREREQ([2.68])
AC_COPYRIGHT([Copyright (c) 2008, 2009, 2010, 2011, 2012 RIKEN. All Rights Reserved.])
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_AUX_DIR([build-aux])
AC_CANONICAL_BUILD
AC_CANONICAL_HOST
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([1.11 -Wall subdir-objects foreign])
AM_MAINTAINER_MODE


# Checks for programs. #################################################
AC_PROG_CXX([mpic++ mpicxx mpiCC mpiicpc mpiclang++ mpiFCCpx mpiFCC])
AC_LANG(C++)
AC_GNU_SOURCE
AC_PROG_CXXCPP
AC_PROG_RANLIB
AM_MISSING_PROG(AUTOM4TE, [autom4te])
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
#AC_ARG_VAR(AM_PROG_AR, [archiver command])
#AC_SUBST(AR, [${AR:-ar}])
#AC_ARG_VAR(ARFLAGS, [archiver flags])
#AC_SUBST(ARFLAGS, [${ARFLAGS:-cru}])

# Conditionals. ########################################################
CPPMD_FEATURE_PMMM
CPPMD_FEATURE_MR3EXAFMM
CPPMD_FEATURE_PME
CPPMD_FEATURE_SIMPLE_FFT
CPPMD_FEATURE_OLDPME
CPPMD_FEATURE_EWALD
CPPMD_FEATURE_LARGEMODEL
CPPMD_FEATURE_FORTRAN_KERNEL
CPPMD_FEATURE_TABLE_FOR_EWALD_REAL
CPPMD_FEATURE_FFT_INTERFACE
CPPMD_FEATURE_HDF

# Checks for libraries. ################################################
AX_BOOST_BASE([1.36.0])

# Checks for header files. #############################################
AC_CHECK_HEADERS([array tr1/array boost/array.hpp])
AC_CHECK_HEADERS([unordered_map tr1/unordered_map boost/unordered_map.hpp])
AC_CHECK_HEADERS([boost/multi_array.hpp])

# Checks for typedefs, structures, and compiler characteristics. #######
AC_C_BIGENDIAN
AX_CHECK_DATA_MODEL
AC_CHECK_SIZEOF(long long)
AC_CHECK_SIZEOF(size_t)
AC_CHECK_SIZEOF(off_t)
AC_CHECK_SIZEOF(ptrdiff_t)
AC_CHECK_SIZEOF(float)
AC_CHECK_SIZEOF(double)
AC_CHECK_SIZEOF(long double)
AX_CREATE_STDINT_H([cppmd/cppmd-int.h])

# Checks for library functions. ########################################
AC_CHECK_FUNCS([getline])

# Checks for system services. ##########################################
AC_SYS_LARGEFILE

# Output ###############################################################
AC_CONFIG_FILES(Makefile)
AC_CONFIG_FILES(cppmd/Makefile)
AC_OUTPUT

dnl vim:se ft=config:
