"""The version agrees everywhere it is written."""

import re
from pathlib import Path

import ssfortran as ss

ROOT = Path(__file__).resolve().parents[2]


def source_version():
    text = (ROOT / "src" / "statespace_capi.f90").read_text()
    return re.search(r'parameter :: version = "([^"]+)"', text).group(1)


def test_versions_agree():
    v = source_version()
    assert ss.__version__ == v
    assert (
        re.search(r'^version = "([^"]+)"', (ROOT / "fpm.toml").read_text(), re.M).group(
            1
        )
        == v
    )
    cmake = (ROOT / "CMakeLists.txt").read_text()
    assert re.search(r"project\(statespace VERSION ([0-9.]+)", cmake).group(1) == v
    conf = (ROOT / "docs" / "conf.py").read_text()
    assert re.search(r'^release = "([^"]+)"', conf, re.M).group(1) == v
