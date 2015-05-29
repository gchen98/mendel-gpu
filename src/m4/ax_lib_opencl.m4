# SYNOPSIS
#
#   AX_LIB_OPENCL()
#
# DESCRIPTION
#
#   This macro checks if the OpenCL headers and library are available and
#   reachable. A non-standard installation path can be provided by using
#   the --with-opencl=<dir> option, in which case headers will be searched
#   in <dir>/includes and the library in <dir>/lib64.
#
#   If the checks are successful, it sets OPENCL_INCLUDE and OPENCL_LIB
#   and defines USE_GPU.
#
#   These checks can be disabled altogether by passing --with-opencl=no or
#   --without-opencl to configure.
#
# LICENSE
#
#   Copyright (c) 2015 Marco Colombo <m.colombo@ed.ac.uk>
#
#   This program is free software; you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation; either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.

AC_DEFUN([AX_LIB_OPENCL], [
    AC_REQUIRE([AC_PROG_CC])

    AC_ARG_WITH([opencl],
                [AC_HELP_STRING([--with-opencl=<dir>],
                                [search for OpenCL below the given directory])])

    # inspect the argument of --with-opencl
    case $with_opencl in
        no|"") ;;
        yes) AC_MSG_WARN([OpenCL directory missing in --with-opencl]) ;;
        *) OPENCL_DIR=$with_opencl ;;
    esac

    if test x$with_opencl != xno ; then

        # save the preexisting variables
        old_CPPFLAGS=$CPPFLAGS
        old_LDFLAGS=$LDFLAGS
        old_LIBS=$LIBS

        # check the presence of the header file
        CPPFLAGS="$CPPFLAGS -I$OPENCL_DIR/include"
        AC_CHECK_HEADERS([CL/cl.hpp], [opencl_hpp=yes], [opencl_hpp=no])

        if test x$opencl_hpp = xyes ; then
           OPENCL_INCLUDE="-I$OPENCL_DIR/include"
        fi

        # check the presence of the library
        LDFLAGS="$LDFLAGS -L$OPENCL_DIR/lib64"
        LIBS="$LIBS -lOpenCL"
        AC_MSG_CHECKING(for presence of libOpenCL)
        AC_LINK_IFELSE(
                [AC_LANG_PROGRAM([#include<CL/cl.hpp>],
                                 [cl::Platform p;])],
                                 [AC_MSG_RESULT(ok)
                                  opencl_lib=yes],
                                 [AC_MSG_RESULT(failed)])

        # restore the original variables
        CPPFLAGS=$old_CPPFLAGS
        LDFLAGS=$old_LDFLAGS
        LIBS=$old_LIBS
    fi

    # substitute makefile variables
    if test x$opencl_hpp = xyes -a x$opencl_lib = xyes ; then
       AC_SUBST(OPENCL_INCLUDE, $OPENCL_INCLUDE)
       AC_SUBST(OPENCL_LIB, -lOpenCL)
       AC_SUBST(USE_GPU, -DUSEGPU)
    else
        AC_MSG_WARN(Building without GPU support)
    fi
])
