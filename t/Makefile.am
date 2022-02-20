#***********************************************************************
#                   GNU Lesser General Public License
#
# This file is part of the ocean grid generator package.
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
#
# Copyright 2021 - @underwoo

TEST_EXTENSIONS = .py
AM_PY_LOG_FALGS = -v
if ENABLE_VENV
  AM_TESTS_ENVIRONMENT = . $(abs_builddir)/venv-test/bin/activate;
  PY_LOG_COMPILER = $(builddir)/venv-test/bin/pytest
else
  PY_LOG_COMPILER = $(PYTEST)
endif

if RUN_PYTESTS
TESTS = \
  test_funcs.py \
  test_ocean_grid_gen.py
endif

if ENABLE_VENV
venv-test:
	$(PYTHON) -m venv $(builddir)/venv-test
	. $(builddir)/venv-test/bin/activate && \
	$(PYTHON) -m pip install pytest && \
	$(PYTHON) -m pip install git+https://github.com/adcroft/numpypi.git && \
	$(PYTHON) -m pip install $(top_srcdir)

test_funcs.py: venv-test

test_ocean_grid_gen.py: venv-test
endif

clean-local:
	-rm -rf venv-test
