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

SUBDIRS = t

bin_SCRIPTS = ocean_grid_generator.py

# Need to use python3 for the two `pip install`s to ensure the packages are
# installed in the virtual env.
if ENABLE_VENV
install-exec-local:
	$(PYTHON) -m venv $(DESTDIR)$(pkglibexecdir)
	. $(DESTDIR)$(pkglibexecdir)/bin/activate && \
	python3 -m pip install --upgrade pip && \
	python3 -m pip install git+https://github.com/adcroft/numpypi.git && \
	python3 -m pip install $(srcdir)

install-exec-hook:
	$(SED) -i -e "s|^#!/usr/bin/env python$$|#!$(pkglibexecdir)/bin/python|" $(DESTDIR)$(bindir)/ocean_grid_generator.py
else
install-exec-hook:
	$(SED) -i -e "s|^#!/usr/bin/env python$$|#!$(PYTHON)|" $(DESTDIR)$(bindir)/ocean_grid_generator.py
endif
