#***********************************************************************
#                   GNU Lesser General Public License
#
# This file is part of the GFDL FRE NetCDF tools package (FRE-NCTools).
#
# FRE-NCTools is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.
#
# FRE-NCTools is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
# for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with FRE-NCTools (LICENSE.md).  If not, see
# <http://www.gnu.org/licenses/>.
#***********************************************************************

# Prelude
AC_PREREQ([2.63])
AC_INIT([Ocean Model Grid Generator],
        [0.0.1],
	[oar.gfdl.help@noaa.gov],
	[omgg],
	[https://github.com/NOAA-GFDL/ocean_model_grid_generator])
# Place the git description informatin in the configure script
m4_define([git_revision],
  [m4_esyscmd_s([git describe --always --dirty 2> /dev/null || echo "unknown version"])])
AC_REVISION([$git_revision])

AC_CONFIG_SRCDIR([ocean_grid_generator.py])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([foreign -Wall -Werror -Woverride tar-ustar])

AC_ARG_ENABLE([venv],
  [AS_HELP_STRING([--enable-venv],
    [Build a venv to contain the script and all python dependencies])],
    [],
    [enable_venv=no])
AM_CONDITIONAL([ENABLE_VENV], [test x"$enable_venv" = x"yes"])

# Find a python version >=3.6
AC_PATH_PROGS([PYTHON], [python3 python])
AX_PROG_PYTHON_VERSION([3.6], ,AC_MSG_ERROR([unable to find a python version >=3.6 in PATH]))

# Check for required python modules
# pytest is needed for tests to run
AX_PYTHON_MODULE([pytest])
AC_PATH_PROG([PYTEST], [pytest])

AS_IF([test ${enable_venv:-no} = "yes"], [
# venv is needed if installing in a virtural environment
AX_PYTHON_MODULE([venv], [fatal])],[
# The following are only needed if not using a venv to install
AX_PYTHON_MODULE([numpypi], [fatal])
AX_PYTHON_MODULE([netCDF4], [fatal])])
AM_CONDITIONAL([RUN_PYTESTS], [test ${HAVE_PYMOD_PYTEST:-no} = "yes" || test ${enable_venv:-no} = "yes"])


AC_CONFIG_FILES([
  Makefile
  t/Makefile
])
AC_OUTPUT
