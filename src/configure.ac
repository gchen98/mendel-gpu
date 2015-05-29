AC_INIT([mendel-gpu], [1.1])

AC_CONFIG_SRCDIR([configure.ac])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_FILES([Makefile])

AM_INIT_AUTOMAKE([foreign])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

AC_LANG([C++])
AC_PROG_CXX

AX_LIB_SAMTOOLS()
AX_LIB_OPENCL()
AC_CHECK_HEADERS([boost/property_tree/ptree.hpp], [],
                 [AC_MSG_ERROR(You need the Boost libraries)])

AC_OUTPUT
