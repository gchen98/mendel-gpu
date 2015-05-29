# SYNOPSIS
#
#   AX_LIB_OPENCL()
#
# DESCRIPTION
#
#   This macro checks if the OpenCL header and library are available and
#   reachable.
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

    #
    # check the presence of the header file
    #
    AC_CHECK_HEADERS([CL/cl.hpp], [],
                     [AC_MSG_ERROR(Could not find the OpenCL header. You could specify its location using CPPFLAGS.)])

    #
    # check the presence of the library
    #
    AC_MSG_CHECKING(for presence of libOpenCL)
    old_LIBS=$LIBS
    LIBS="$LIBS -lOpenCL"
    AC_LINK_IFELSE(
        [AC_LANG_PROGRAM([#include<CL/cl.hpp>],
                         [cl::Platform p;])],
        [AC_MSG_RESULT(ok)],
        [AC_MSG_RESULT(failed)
         AC_MSG_ERROR(Could not find the OpenCL library. You could specify its location using LDFLAGS.)])
    LIBS=$oldLIBS
])
