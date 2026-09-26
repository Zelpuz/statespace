"""Results conveniences: fitted values, forecasts, components, summaries."""

import numpy as np
import pytest

import ssfortran as ss
from _fixtures import close, read_fixture

NILE = read_fixture("nile_llevel_exact")["y"].ravel()


def test_fitted_values_and_forecasts():
    res = ss.StructuralModel(NILE, [ss.Irregular(), ss.Level()]).fit()
    fv = res.fittedvalues
    assert fv.shape == (NILE.size,)
    close(res.resid, NILE - fv)
    fc = res.get_forecast(5)
    lo, hi = fc.conf_int(0.1)
    assert np.all(lo < fc.predicted_mean) and np.all(fc.predicted_mean < hi)
    mean, cov = res.model.representation(res.params).forecast(5)
    close(fc.predicted_mean, mean[0])
    close(fc.se_mean, np.sqrt(cov[0, 0]))
    close(res.forecast(5), mean[0])


def test_multivariate_shapes():
    y = read_fixture("mv_invariant")["y"].T
    res = ss.StructuralModel(y, [ss.Irregular(cov="full"), ss.Level(cov="full")]).fit()
    assert res.fittedvalues.shape == y.shape
    fc = res.get_forecast(2)
    assert fc.predicted_mean.shape == (2, 2) and fc.var_pred_mean.shape == (2, 2, 2)


def test_components_add_up():
    fx = read_fixture("uc_bsm_dummy")
    mod = ss.StructuralModel(fx["y"].T, [ss.Irregular(), ss.Trend(), ss.Seasonal(4)])
    comps, var = mod.components(fx["params"], variance=True)
    assert list(comps) == ["irregular", "trend", "seasonal"]
    close(comps["trend"], fx["level"])
    close(comps["seasonal"], fx["seasonal"])
    y = fx["y"].ravel()
    close(comps["irregular"] + comps["trend"] + comps["seasonal"], y)
    assert np.all(var["trend"] >= 0)


def test_summary_and_diagnostics():
    res = ss.StructuralModel(NILE, [ss.Irregular(), ss.Level()]).fit()
    text = res.summary()
    for key in ("P>|z|", "Ljung-Box", "Jarque-Bera", "sigma2.level", "Diffuse periods"):
        assert key in text
    d = res.diagnostics()
    e = res.filter().standardized_forecasts_error[0, 1:]
    jb = ss.diagnostics.jarque_bera(e)
    assert d["jarque_bera"] == pytest.approx(jb[:2])
