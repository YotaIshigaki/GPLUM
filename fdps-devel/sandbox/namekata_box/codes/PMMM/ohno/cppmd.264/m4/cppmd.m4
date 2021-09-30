# -*- mode: c; -*-
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

# CPPMD_FEATURE_PMMM
# --------------------------------------
AC_DEFUN([CPPMD_FEATURE_PMMM],
         [AC_ARG_ENABLE([pmmm],
                        [AS_HELP_STRING([--enable-pmmm], [enable PMMM long range interaction calculation])],
                        [AS_IF([test "x$enableval" != xyes && test "x$enableval" != xno], AC_MSG_ERROR([bad value ${enableval} for --enable-pmmm]))],
			[enable_pmmm=no])
          AS_IF([test "x$enable_pmmm" != xno],
                [AC_DEFINE([CPPMD_ENABLE_PMMM], [1], [Define if you enable PMMM for long range interaction calculation.])
                 AC_DEFINE([CPPMD_ENABLE_LONGRANGE], [1], [Define if you enable long range interaction calculation.])])
          AM_CONDITIONAL([ENABLE_PMMM], [test "x$enable_PMmm" != xno])
          ])dnl CPPMD_FEATURE_PMMM

# CPPMD_FEATURE_MR3EXAFMM
# --------------------------------------
AC_DEFUN([CPPMD_FEATURE_MR3EXAFMM],
         [AC_ARG_ENABLE([mr3exafmm],
                        [AS_HELP_STRING([--enable-mr3exafmm], [enable MR3EXAFMM long range interaction calculation])],
                        [AS_IF([test "x$enableval" != xyes && test "x$enableval" != xno], AC_MSG_ERROR([bad value ${enableval} for --enable-mr3exafmm]))],
                        [enable_mr3exafmm=no])
          AS_IF([test "x$enable_mr3exafmm" != xno],
                [AC_DEFINE([CPPMD_ENABLE_MR3EXAFMM], [1], [Define if you enable MR3EXAFMM for long range interaction calculation.])
                 AC_DEFINE([CPPMD_ENABLE_FMM], [1], [Define if you enable FMM for long range interaction calculation.])
                 AC_DEFINE([CPPMD_ENABLE_LONGRANGE], [1], [Define if you enable long range interaction calculation.])])
          AM_CONDITIONAL([ENABLE_MR3EXAFMM], [test "x$enable_mr3exafmm" != xno])
          ])dnl CPPMD_FEATURE_MR3EXAFMM

# CPPMD_FEATURE_PME
# --------------------------------------
AC_DEFUN([CPPMD_FEATURE_PME],
         [AC_ARG_ENABLE([pme],
                        [AS_HELP_STRING([--enable-pme], [enable PME long range interaction calculation])],
                        [AS_IF([test "x$enableval" != xyes && test "x$enableval" != xno], AC_MSG_ERROR([bad value ${enableval} for --enable-pme]))],
                        [enable_pme=no])
          AS_IF([test "x$enable_pme" != xno],
                [AC_DEFINE([CPPMD_ENABLE_PME], [1], [Define if you enable Original PME for long range interaction calculation.])
                 AC_DEFINE([CPPMD_ENABLE_LONGRANGE], [1], [Define if you enable long range interaction calculation.])
                 ])
          AM_CONDITIONAL([ENABLE_PME], [test "x$enable_pme" != xno])
          ])dnl CPPMD_FEATURE_PME

# CPPMD_FEATURE_OLDPME
# --------------------------------------
AC_DEFUN([CPPMD_FEATURE_OLDPME],
         [AC_ARG_ENABLE([oldpme],
                        [AS_HELP_STRING([--enable-oldpme], [enable OLDPME long range interaction calculation])],
                        [AS_IF([test "x$enableval" != xyes && test "x$enableval" != xno], AC_MSG_ERROR([bad value ${enableval} for --enable-oldpme]))],
                        [enable_oldpme=no])
          AS_IF([test "x$enable_oldpme" != xno],
                [AC_DEFINE([CPPMD_ENABLE_OLDPME], [1], [Define if you enable Original OLDPME for long range interaction calculation.])
                 AC_DEFINE([CPPMD_ENABLE_LONGRANGE], [1], [Define if you enable long range interaction calculation.])
                 ])
          AM_CONDITIONAL([ENABLE_OLDPME], [test "x$enable_oldpme" != xno])
          ])dnl CPPMD_FEATURE_OLDPME

# CPPMD_FEATURE_SIMPLE_FFT
# --------------------------------------
AC_DEFUN([CPPMD_FEATURE_SIMPLE_FFT],
         [AC_ARG_ENABLE([simple-fft],
                        [AS_HELP_STRING([--enable-simple-fft], [enable SIMPLE_FFT long range interaction calculation])],
                        [AS_IF([test "x$enableval" != xyes && test "x$enableval" != xno], AC_MSG_ERROR([bad value ${enableval} for --enable-simple-fft]))],
                        [enable_simple_fft=no])
          AS_IF([test "x$enable_simple_fft" != xno],
                [AC_DEFINE([CPPMD_ENABLE_SIMPLE_FFT], [1], [Define if you enable SIMPLE_FFT for long range interaction calculation.])
                 AC_DEFINE([CPPMD_ENABLE_LONGRANGE], [1], [Define if you enable long range interaction calculation.])
                 ])
          AM_CONDITIONAL([ENABLE_SIMPLE_FFT], [test "x$enable_simple_fft" != xno])
          ])dnl CPPMD_FEATURE_SIMPLE_FFT

# CPPMD_FEATURE_EWALD
# --------------------------------------
AC_DEFUN([CPPMD_FEATURE_EWALD],
         [AC_ARG_ENABLE([ewald],
                        [AS_HELP_STRING([--enable-ewald], [enable Ewald long range interaction calculation])],
                        [AS_IF([test "x$enableval" != xyes && test "x$enableval" != xno], AC_MSG_ERROR([bad value ${enableval} for --enable-ewald]))],
                        [enable_ewald=no])
          AS_IF([test "x$enable_ewald" != xno],
                [AC_DEFINE([CPPMD_ENABLE_EWALD], [1], [Define if you enable Original Ewald for long range interaction calculation.])
                 AC_DEFINE([CPPMD_ENABLE_LONGRANGE], [1], [Define if you enable long range interaction calculation.])])
          AM_CONDITIONAL([ENABLE_EWALD], [test "x$enable_ewald" != xno])
          ])dnl CPPMD_FEATURE_EWALD

## CPPMD_FEATURE_STMEWALD
## --------------------------------------
#AC_DEFUN([CPPMD_FEATURE_STMEWALD],
#         [AC_ARG_ENABLE([stmewald],
#                        [AS_HELP_STRING([--enable-stmewald], [enable STMEWALD long range interaction calculation])],
#                        [AS_IF([test "x$enableval" != xyes && test "x$enableval" != xno], AC_MSG_ERROR([bad value ${enableval} for --enable-stmewald]))],
#                        [enable_stmewald=no])
#          AS_IF([test "x$enable_stmewald" != xno],
#                [AC_DEFINE([CPPMD_ENABLE_STMEWALD], [1], [Define if you enable STMEWALD for long range interaction calculation.])
#                 AC_DEFINE([CPPMD_ENABLE_LONGRANGE], [1], [Define if you enable long range interaction calculation.])
#                 AC_DEFINE([STMEWALD], [1], [Enable STMEWALD])])
#          AM_CONDITIONAL([ENABLE_STMEWALD], [test "x$enable_stmewald" != xno])
#          ])dnl CPPMD_FEATURE_STMEWALD

## CPPMD_FEATURE_GFMM
## --------------------------------------
#AC_DEFUN([CPPMD_FEATURE_GFMM],
#         [AC_ARG_ENABLE([gfmm],
#         [AS_HELP_STRING([--enable-gfmm], [enable GFMM long range interaction calculation])],
#         [AS_IF([test "x$enableval" != xyes && test "x$enableval" != xno], AC_MSG_ERROR([bad value ${enableval} for --enable-gfmm]))],
#         [enable_gfmm=no])
#         AS_IF([test "x$enable_gfmm" != xno],
#         [AC_DEFINE([CPPMD_ENABLE_GFMM], [1], [Define if you enable GFMM for long range interaction calculation.])
#                 AC_DEFINE([CPPMD_ENABLE_FMM], [1], [Define if you enable FMM for long range interaction calculation.])
#                 AC_DEFINE([CPPMD_ENABLE_LONGRANGE], [1], [Define if you enable long range interaction calculation.])
#                 AC_DEFINE([FMM], [1], [Enable FMM])])
#         AM_CONDITIONAL([ENABLE_GFMM], [test "x$enable_gfmm" != xno])
#         ])dnl CPPMD_FEATURE_GFMM

# CPPMD_FEATURE_HDF
# --------------------------------------
AC_DEFUN([CPPMD_FEATURE_HDF],
         [
	AC_MSG_CHECKING(for HDF)
	AC_ARG_WITH([hdf],
                        [AS_HELP_STRING([--with-hdf@<:@=DIR@:>@], [use HDF (default is no)])],
			[AS_CASE([$withval],
				[no], [with_hdf=no],
				[yes],[with_hdf=yes hdf_path=""],
				[with_hdf=path hdf_path=$withval]
				)],
                        [with_hdf=no,
			AC_MSG_RESULT(shared no)
			])
	AC_ARG_WITH([static-hdf],
		    [AS_HELP_STRING([--with-static-hdf@<:@=DIR@:>@], [use static HDF (default is no)])],
		    [AS_CASE([$withval],
                                [no], [with_static_hdf=no],
                                [yes],[with_static_hdf=yes hdf_path=""],
                                [with_static_hdf=path hdf_path=$withval]
                                )],
                        [with_static_hdf=no,
                        AC_MSG_RESULT(static no)
                        ])
	AS_CASE([$with_hdf],
		[yes], [AC_DEFINE([USE_HDF], [1], [Define if you use HDF])
			HDF_LDFLAGS="-lhdf5_cpp -lhdf5 -lz"
			HDFINC=""
	AC_MSG_RESULT(shared yes $HDF_LDFLAGS)
			AC_SUBST(HDF_LDFLAGS)
			],
		[path], [AC_DEFINE([USE_HDF], [1], [Define if you use HDF at PATH])
			HDFINC="-I$hdf_path/include"
			HDF_LDFLAGS="-L$hdf_path/lib -lhdf5_cpp -lhdf5 -lz"
	AC_MSG_RESULT($HDFINC $HDF_LDFLAGS)
			AC_SUBST(HDF_LDFLAGS)
			])
	AS_CASE([$with_static_hdf],
		[yes], [AC_DEFINE([USE_HDF], [1], [Define if you use HDF])
			AS_IF([test -d /usr/lib64],
                              [lib_arch_name=64])
			HDF_LDADD="/usr/lib$lib_arch_name/libhdf5_cpp.a /usr/lib$lib_arch_name/libhdf5.a -lz"
			HDFINC=""
	AC_MSG_RESULT(static yes $HDF_LDADD)
			AC_SUBST(HDF_LDADD)
			],
		[path], [AC_DEFINE([USE_HDF], [1], [Define if you use HDF at PATH])
			HDFINC="-I$hdf_path/include"
			HDF_LDADD="$hdf_path/lib/libhdf5_cpp.a $hdf_path/lib/libhdf5.a -lz"
	AC_MSG_RESULT($HDFINC $HDF_LDADD)
			AC_SUBST(HDF_LDADD)
			])
	AM_CONDITIONAL([USE_HDF], [test "x$with_hdf" != xno])
	])dnl CPPMD_FEATURE_HDF

# CPPMD_FEATURE_LARGEMODEL
# --------------------------------------
AC_DEFUN([CPPMD_FEATURE_LARGEMODEL],
         [AC_ARG_ENABLE([largemodel],
                        [AS_HELP_STRING([--enable-largemodel], [use 64-bit integer in Amber file I/O])],
                        [AS_IF([test "x$enableval" != xyes && test "x$enableval" != xno], AC_MSG_ERROR([bad value ${enableval} for --enable-largemodel]))],
                        [enable_largemodel=no])
          AS_IF([test "x$enable_largemodel" != xno], [AC_DEFINE([CPPMD_ENABLE_LARGEMODEL], [1], [Define if you use 64-bit integer in Amber file I/O.])])
          AM_CONDITIONAL([ENABLE_LARGEMODEL], [test "x$enable_largemodel" != xno])
          ])dnl CPPMD_FEATURE_LARGEMODEL

# CPPMD_FEATURE_TABLE_FOR_EWALD_REAL
# --------------------------------------
AC_DEFUN([CPPMD_FEATURE_TABLE_FOR_EWALD_REAL],
         [AC_ARG_ENABLE([table-for-ewald-real],
                        [AS_HELP_STRING([--disable-table-for-ewald-real], [disable fast evaluation of Ewald real space])],
                        [AS_IF([test "x$enableval" != xyes && test "x$enableval" != xno], AC_MSG_ERROR([bad value ${enableval} for --disable-table-for-ewald-real]))],
                        [enable_table_for_ewald_real=yes])
          AS_IF([test "x$enable_table_for_ewald_real" != xno], [AC_DEFINE([CPPMD_ENABLE_TABLE_FOR_EWALD_REAL], [1], [Define if you use fast evaluation of Ewald real space])])
          AM_CONDITIONAL([ENABLE_TABLE_FOR_EWALD_REAL], [test "x$enable_table_for_ewald_real" != xno])
          ])dnl CPPMD_FEATURE_TABLE_FOR_EWALD_REAL

# CPPMD_FEATURE_FORTRAN_KERNEL
# --------------------------------------
AC_DEFUN([CPPMD_FEATURE_FORTRAN_KERNEL],
         [AC_ARG_ENABLE([fortran-kernel],
                        [AS_HELP_STRING([--enable-fortran-kernel], [use kernel code written in Fortran])],
                        [AS_IF([test "x$enableval" != xyes &&  test "x$enableval" != xno], AC_MSG_ERROR([bad value ${enableval} for --enable-fortran-kernel]))],
                        [enable_fortran_kernel=no])
          AS_IF([test "x$enable_fortran_kernel" != xno], [AC_DEFINE([CPPMD_ENABLE_FORTRAN_KERNEL], [1], [Define if you use Fortran kernel code])
                                                          AC_PROG_FC])
          AM_CONDITIONAL([ENABLE_FORTRAN_KERNEL], [test "x$enable_fortran_kernel" != xno])
          ])dnl CPPMD_FEATURE_FORTRAN_KERNEL

# CPPMD_FEATURE_FFT_INTERFACE
# --------------------------------------
AC_DEFUN([CPPMD_FEATURE_FFT_INTERFACE],
         [AC_ARG_WITH([fft],
                      [AS_HELP_STRING([--with-fft=@<:@ffte|fftw2|fftw3l|mkl@:>@], [select fft library interface])],
                      [AS_CASE([$withval],
                               [ffte], [with_fft=ffte],
                               [fftw2], [with_fft=fftw2],
                               [fftw3], [with_fft=fftw3],
                               [mkl], [with_fft=mkl],
                               [with_fft=no])],
                      [with_fft=no])
          AS_CASE([$with_fft],
                  [ffte], [AC_DEFINE([USE_FFTE], [1], [Define if you use FFTE library])],
                  [fftw2], [AC_DEFINE([USE_FFTW2], [1], [Define if you use FFTW2 libray])
dnl                            LIBS="-lrfftw_mpi -lfftw_mpi -lrfftw -lfftw $LIBS"
                            ],
                  [fftw3], [AC_DEFINE([USE_FFTW3], [1], [Define if you use FFTW3 libray])
dnl                            LIBS="-lfftw3_mpi -lfftw3 $LIBS"
                            ],
                  [mkl], [AC_DEFINE([USE_MKL], [1], [Define if you use MKL library])])
          ])dnl CPPMD_FEATURE_FFT_INTERFACE
