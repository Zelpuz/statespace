# Build the documentation: `make -C docs html` from the repository root.
# The shared library is built first because autodoc imports ssfortran.
ROOT := $(abspath ..)
VENV ?= $(ROOT)/.venv/bin
# CMake finds ninja (and Sphinx its tools) in the same environment.
export PATH := $(VENV):$(PATH)
OUT := $(ROOT)/build/docs

.PHONY: html doctest lib clean

lib:
	$(VENV)/cmake -S $(ROOT) -B $(ROOT)/build/cmake -G Ninja >/dev/null
	$(VENV)/cmake --build $(ROOT)/build/cmake >/dev/null

html: lib
	$(VENV)/sphinx-build -W --keep-going -b html . $(OUT)/html

doctest: lib
	$(VENV)/sphinx-build -W -b doctest . $(OUT)/doctest

clean:
	rm -rf $(OUT) reference/python/generated
