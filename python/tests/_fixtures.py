"""Reading the statsmodels fixtures of the Fortran tests (test/fixtures)."""

from pathlib import Path

import numpy as np

import ssfortran as ss

FIXTURES = Path(__file__).resolve().parents[2] / "test" / "fixtures"
RTOL = 1e-9


def read_fixture(name):
    """Named arrays: a header "name ndim shape..." then values in F order."""
    out = {}
    with open(FIXTURES / f"{name}.txt") as f:
        lines = f.read().split("\n")
    i = 0
    while i < len(lines) and lines[i]:
        head = lines[i].split()
        key, ndim = head[0], int(head[1])
        shape = tuple(int(s) for s in head[2:2 + ndim])
        size = int(np.prod(shape))
        vals = np.array([float(v) for v in lines[i + 1:i + 1 + size]])
        out[key] = vals.reshape(shape, order="F")
        i += 1 + size
    return out


def rep_from_fixture(fx):
    """As test/fixture_io.f90: matrices, then init_blocks, a1/P1 or stationary."""
    y = fx["y"]
    rep = ss.Representation(y.T, k_states=fx["T"].shape[0], k_posdef=fx["Q"].shape[0])
    for name, key in [("design", "Z"), ("obs_cov", "H"), ("transition", "T"),
                      ("selection", "R"), ("state_cov", "Q"), ("state_intercept", "c"),
                      ("obs_intercept", "d")]:
        rep[name] = fx[key]
    if "init_blocks" in fx:
        for first, last, kind in np.rint(fx["init_blocks"]).astype(int):
            rep.initialize_block(first - 1, last, kind)
    elif "a1" in fx:
        rep.initialize_known(fx["a1"], fx["P1"])
    else:
        rep.initialize_stationary()
    return rep


def close(actual, expected, mask=None):
    actual, expected = np.asarray(actual), np.asarray(expected)
    assert actual.shape == expected.shape
    use = ~np.isnan(expected) if mask is None else mask & ~np.isnan(expected)
    err = np.max(np.abs(actual[use] - expected[use]) / np.maximum(1.0, np.abs(expected[use])))
    assert err <= RTOL, f"relative error {err:.3e}"
