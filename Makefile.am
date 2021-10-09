bin_SCRIPTS = ocean_grid_generator.py

if ENABLE_VENV
install-exec-local:
	python3 -m venv $(DESTDIR)$(pkglibexecdir)
	. $(DESTDIR)$(pkglibexecdir)/bin/activate && \
	python3 -m pip install --use-feature=in-tree-build git+https://github.com/underwoo/numpypi.git@pip.installable && \
	python3 -m pip install --use-feature=in-tree-build $(srcdir)

install-exec-hook:
	sed -i -e "s|^#!/usr/bin/env python$$|#!$(pkglibexecdir)/bin/python|" $(DESTDIR)$(bindir)/ocean_grid_generator.py
endif
