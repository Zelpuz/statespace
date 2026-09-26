"""pytest setup for the docstring examples (pytest --doctest-modules
python/ssfortran): numpy and ssfortran are available as np and ss."""

import numpy as np
import pytest

import ssfortran as ss


@pytest.fixture(autouse=True)
def _doctest_namespace(doctest_namespace):
    doctest_namespace["np"] = np
    doctest_namespace["ss"] = ss
