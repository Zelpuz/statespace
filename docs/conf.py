"""Sphinx configuration for the statespace / ssfortran documentation.

Build from the repository root with ``make -C docs html`` (builds the
shared library first; autodoc imports ``ssfortran``).
"""

import doctest
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python" / "src"))
os.environ.setdefault(
    "SSFORTRAN_LIB", str(ROOT / "build" / "cmake" / "libstatespace.so")
)

project = "statespace"
author = "Taylor D. Edwards"
copyright = "2026, Taylor D. Edwards"
release = "0.1.4"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.doctest",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "numpydoc",
    "myst_parser",
]

source_suffix = {".rst": "restructuredtext", ".md": "markdown"}
exclude_patterns = ["_build", "README.md"]
templates_path = ["_templates"]

html_theme = "pydata_sphinx_theme"
html_title = f"statespace {release}"
html_baseurl = "https://zelpuz.github.io/statespace/"
html_theme_options = {
    "navigation_depth": 3,
    "show_toc_level": 2,
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/Zelpuz/statespace",
            "icon": "fa-brands fa-github",
        },
        {
            "name": "PyPI",
            "url": "https://pypi.org/project/ssfortran/",
            "icon": "fa-brands fa-python",
        },
    ],
}

autosummary_generate = True
autodoc_default_options = {"members": True}
autodoc_typehints = "none"

numpydoc_show_class_members = False
numpydoc_class_members_toctree = False
numpydoc_xref_param_type = False
# Structure checks, applied to every public docstring.
numpydoc_validation_checks = {
    "all",
    "GL01",
    "GL02",
    "ES01",
    "SA01",
    "EX01",
    "RT02",
    "PR09",
    "RT05",
    "SS06",
    "GL08",
}
# Private names; and the results dataclasses, whose fields are documented as
# Attributes rather than constructor Parameters.
numpydoc_validation_exclude = {
    r"\._",
    r"^ssfortran\._lib",
    r"^ssfortran\.(Filter|Smoother|Augmented|Fit|Forecast)Results$",
}

intersphinx_mapping = {
    "numpy": ("https://numpy.org/doc/stable/", None),
    "python": ("https://docs.python.org/3", None),
}

myst_enable_extensions = ["dollarmath"]
myst_heading_anchors = 3

doctest_default_flags = (
    doctest.NORMALIZE_WHITESPACE | doctest.ELLIPSIS | doctest.DONT_ACCEPT_TRUE_FOR_1
)

doctest_global_setup = f"""
import os
import numpy as np
import ssfortran as ss
os.chdir(r"{ROOT}")   # examples read data/ relative to the repository root
"""
