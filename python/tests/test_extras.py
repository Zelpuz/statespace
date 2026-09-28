"""The DK ch. 4-7 extras through the Python layer: against statsmodels
fixtures where they exist, otherwise against identities between methods."""

import numpy as np
import pytest
import ssfortran as ss
from _fixtures import close, read_fixture, rep_from_fixture


@pytest.fixture
def mv():
    return rep_from_fixture(read_fixture("mv_invariant"))


def test_smoother_variants_agree(mv):
    ref = mv.smooth()
    a, V = ref.smoothed_state, ref.smoothed_state_cov
    close(mv.fast_smoother(), a)
    for alt in (mv.classical_smoother(), mv.two_filter_smoother()):
        close(alt[0], a)
        close(alt[1], V)
    # The Whittle recursion runs backwards through T^-1 and loses accuracy
    # with the length of the series (DK 4.6.3); as in test_recursions.f90,
    # 25 periods of the local linear trend.
    fx = read_fixture("nile_lltrend_exact")
    short = rep_from_fixture(dict(fx, y=fx["y"][:, :25]))
    w = short.whittle_smoother()
    exact = short.smooth().smoothed_state
    assert np.max(np.abs(w - exact) / np.maximum(1.0, np.abs(exact))) < 1e-6
    sq = mv.smooth(method="sqrt")
    close(sq.smoothed_state, a)
    close(sq.smoothed_state_cov, V)
    close(mv.filter(method="sqrt").predicted_state_cov, ref.filter.predicted_state_cov)


def test_fixed_point_and_lag(mv):
    n = mv.nobs
    f = mv.filter()
    ref = mv.smooth()
    pa, pV = mv.fixed_point_smoother(10)
    close(pa[:, 0], f.filtered_state[:, 10])  # given y up to period 10
    close(pa[:, -1], ref.smoothed_state[:, 10])  # given all data
    a, V = mv.fixed_lag_smoother(n - 1)
    close(a[:, 0], ref.smoothed_state[:, 0])
    assert np.all(np.isnan(a[:, 1:]))


def test_update_smoothed():
    fx = read_fixture("mv_invariant")
    rep = rep_from_fixture(fx)
    head = rep_from_fixture(dict(fx, y=fx["y"][:, :25])).smooth()
    a, V = rep.update_smoothed(head.smoothed_state, head.smoothed_state_cov, 25)
    ref = rep.smooth()
    close(a, ref.smoothed_state)
    close(V, ref.smoothed_state_cov)


def test_autocov_cov_between_and_weights(mv):
    ex = read_fixture("extras_mv_invariant")
    close(mv.smoothed_state_autocov(), ex["autocov"])
    n = mv.nobs
    cov3 = np.stack(
        [mv.smoothed_state_cov_between(t, t + 3) for t in range(n - 3)], axis=2
    )
    close(cov3, ex["cov_shift3"])
    W, C, A = mv.smoothed_state_weights()
    close(W, ex["weights"])
    close(C, ex["c_weights"])
    close(A, ex["prior_weights"])
    Wa, Watt = mv.filtered_state_weights()
    assert Wa.shape == (mv.k_states, mv.k_endog, n, n)
    assert np.allclose(Wa[:, :, 5, 5:], 0.0)  # a_t uses y_1..y_t-1 only
    L = mv.innovation_transition(3)
    assert L.shape == (mv.k_states, mv.k_states)


def test_augmented_filter():
    rep = rep_from_fixture(read_fixture("nile_llevel_exact"))
    aug = rep.augmented()
    exact = rep.smooth()
    assert aug.llf == pytest.approx(exact.filter.llf, rel=1e-10)
    assert aug.delta.shape == (1,) and aug.delta_cov.shape == (1, 1)
    close(aug.smoothed_state, exact.smoothed_state)
    rep.marginal_likelihood = True
    assert rep.loglike() == pytest.approx(aug.llf_marginal, rel=1e-10)


def test_em_reaches_mle():
    rep = rep_from_fixture(read_fixture("nile_llevel_exact"))
    v = np.nanvar(rep.endog)
    rep["obs_cov"] = [[v / 2]]
    rep["state_cov"] = [[v / 2]]
    llf, niter, path = rep.em(maxiter=2000, tol=1e-10)
    assert np.all(np.diff(path) >= -1e-8)  # EM never decreases llf
    assert rep["obs_cov"][0, 0] == pytest.approx(15099, rel=1e-3)
    assert rep["state_cov"][0, 0] == pytest.approx(1469.1, rel=1e-2)


def test_collapse_and_restrictions(mv):
    # collapsing needs Z of full column rank: p = 4 series on m = 2 states
    rng = np.random.default_rng(2)
    wide = ss.Representation(rng.standard_normal((80, 4)), k_states=2)
    wide["design"] = rng.standard_normal((4, 2))
    wide["obs_cov"] = np.diag([1.0, 2.0, 0.5, 1.5])
    wide["transition"] = np.diag([0.9, 0.5])
    wide["selection"] = np.eye(2)
    wide["state_cov"] = np.eye(2)
    wide.initialize_stationary()
    crep, adj = wide.collapse()
    assert crep.k_endog == wide.k_states
    assert crep.loglike() + adj.sum() == pytest.approx(wide.loglike(), rel=1e-10)
    close(crep.smooth().smoothed_state, wide.smooth().smoothed_state)
    # alpha_1 + alpha_2 = 0 at every period
    R = np.zeros((1, mv.k_states))
    R[0, :2] = 1.0
    rr = mv.add_state_restrictions(R, np.zeros((1, mv.nobs)))
    rr.filter_method = ss.FILTER_UNIVARIATE
    a = rr.smooth().smoothed_state
    assert np.max(np.abs(a[0] + a[1])) < 1e-8


def test_residual_diagnostics():
    rep = rep_from_fixture(read_fixture("nile_llevel_exact"))
    e, u = rep.auxiliary_residuals()
    ev, uv = rep.auxiliary_residuals(vector=True)
    close(ev[:, 1:], e[:, 1:])  # p = 1: the same
    r_stat, e_stat = rep.de_jong_penzer()
    assert r_stat.shape == (1, rep.nobs) and e_stat.shape == (1, rep.nobs)
    r2 = rep.r2_diffuse()
    assert np.isfinite(r2) and r2 < 1.0


def test_least_squares_residuals():
    rng = np.random.default_rng(0)
    n, k = 60, 2
    X = np.column_stack([np.ones(n), rng.standard_normal(n)])
    y = X @ [1.0, 2.0] + rng.standard_normal(n)
    rep = ss.Representation(y, k_states=k, k_posdef=k)
    rep["design"] = X.T[None, :, :]
    rep["obs_cov"] = [[1.0]]
    rep["transition"] = np.eye(k)
    rep["selection"] = np.eye(k)
    rep["state_cov"] = np.zeros((k, k))
    rep.initialize_diffuse()
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    close(rep.least_squares_residuals(0, 2)[0], y - X @ beta)


def test_djs_simulation_smoother(mv):
    a1, e1, u1 = mv.djs_simulation_smoother(rng=3)
    a2, _, _ = mv.djs_simulation_smoother(rng=3)
    assert np.array_equal(a1, a2)
    # mean over draws approaches the smoothed state
    ref = mv.smooth().smoothed_state
    draws = np.mean([mv.djs_simulation_smoother(rng=s)[0] for s in range(400)], axis=0)
    assert np.max(np.abs(draws - ref)) < 0.5 * np.sqrt(
        mv.smooth().smoothed_state_cov.max()
    )


def test_estimation_bias():
    y = read_fixture("nile_llevel_exact")["y"].ravel()
    res = ss.StructuralModel(y, [ss.Irregular(), ss.Level()]).fit()
    b1 = res.estimation_bias(ndraw=40, seed=5)
    b2 = res.estimation_bias(ndraw=40, seed=5)
    assert b1.shape == (1, y.size) and np.array_equal(b1, b2)
    bias, bV = res.estimation_bias(ndraw=40, seed=5, variance=True)
    assert bV.shape == (1, 1, y.size)
    sd = np.sqrt(res.smooth().smoothed_state_cov[0, 0])
    assert np.max(np.abs(bias[0]) / sd) < 1.0
