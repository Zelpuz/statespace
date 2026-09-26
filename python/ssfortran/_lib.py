"""Loading the Fortran library and declaring its C interface (ctypes)."""

import ctypes
import os
import sys
from pathlib import Path

import numpy as np

# Status codes (statespace_kinds)
SS_OK = 0
_ERRORS = {
    1: "inconsistent array dimensions",
    2: "matrix not positive definite",
    3: "missing or invalid initialization",
    4: "not supported for this model",
    5: "singular linear system",
    6: "transition matrix not stationary",
    7: "iteration did not converge",
}


class StateSpaceError(RuntimeError):
    """An error reported by the Fortran library.

    Parameters
    ----------
    code : int
        Status code: 1 inconsistent dimensions, 2 matrix not positive
        definite, 3 missing or invalid initialization, 4 not supported for
        this model, 5 singular linear system, 6 transition matrix not
        stationary, 7 iteration did not converge.
    where : str
        The library routine that reported it.

    Attributes
    ----------
    code : int
        The status code.
    """

    def __init__(self, code, where):
        self.code = code
        super().__init__(f"{where}: {_ERRORS.get(code, 'error')} (code {code})")


def _library_names():
    if sys.platform == "darwin":
        return ["libstatespace.dylib"]
    if sys.platform == "win32":
        return ["statespace.dll", "libstatespace.dll"]
    return ["libstatespace.so"]


def _find_library():
    """The library: $SSFORTRAN_LIB, next to this package (installed wheel),
    or the repository's CMake build directory (development)."""
    env = os.environ.get("SSFORTRAN_LIB")
    if env:
        return env
    here = Path(__file__).resolve().parent
    candidates = [here] + [here.parents[1] / "build" / "cmake"]
    for d in candidates:
        for name in _library_names():
            if (d / name).exists():
                return str(d / name)
    raise OSError("libstatespace not found; build it with CMake (see README) "
                  "or set SSFORTRAN_LIB")


lib = ctypes.CDLL(_find_library())

_i = ctypes.c_int
_d = ctypes.c_double
_p = ctypes.c_void_p
_pp = ctypes.POINTER(ctypes.c_void_p)
_pd = ctypes.POINTER(ctypes.c_double)
_pi = ctypes.POINTER(ctypes.c_int)

_SIGNATURES = {
    "ss_version": [ctypes.c_char_p, _i],
    "ss_rep_new": [_i, _i, _i, _i, _pp],
    "ss_rep_free": [_p],
    "ss_rep_set": [_p, _i, _i, _p],
    "ss_rep_set_int": [_p, _i, _i],
    "ss_rep_set_real": [_p, _i, _d],
    "ss_rep_init_known": [_p, _p, _p],
    "ss_rep_init_diffuse": [_p],
    "ss_rep_init_approx_diffuse": [_p, _d],
    "ss_rep_init_stationary": [_p],
    "ss_rep_init_general": [_p, _p, _p, _p],
    "ss_rep_init_block": [_p, _i, _i, _i, _p, _p, _d],
    "ss_loglike": [_p, _pd],
    "ss_loglike_concentrated": [_p, _pd, _pd],
    "ss_filter": [_p, _pp],
    "ss_filter_free": [_p],
    "ss_filter_scalars": [_p, _pd, _pi, _pi, _pi],
    "ss_filter_get": [_p, _i, _p],
    "ss_smooth": [_p, _p, _pp],
    "ss_smoother_free": [_p],
    "ss_smoother_get": [_p, _i, _p],
    "ss_rep_info": [_p, _pi, _pi, _pi, _pi, _p],
    "ss_rep_get": [_p, _i, _p],
    # models
    "ss_struct_new": [_i, _i, _p, _pp],
    "ss_struct_free": [_p],
    "ss_struct_add_irregular": [_p, _i, _p],
    "ss_struct_add_level": [_p, _i, _i],
    "ss_struct_add_trend": [_p, _i, _i, _i],
    "ss_struct_add_seasonal": [_p, _i, _i, _i, _i],
    "ss_struct_add_cycle": [_p, _i, _i, _d, _d, _i],
    "ss_struct_add_regression": [_p, _i, _p, _p, _i, _i],
    "ss_struct_add_arima": [_p, _i, _i, _i, _i, _i, _i, _i, _i, _i, _i, _i],
    "ss_struct_add_continuous": [_p, _i, _p, _i],
    "ss_struct_build": [_p, _i, _p, _p, _pp],
    "ss_model_free": [_p],
    "ss_model_info": [_p, _pi, _pi, _pi, _pi, _pi],
    "ss_model_param_names": [_p, _i, ctypes.c_char_p],
    "ss_model_start_params": [_p, _p],
    "ss_model_transform": [_p, _p, _p],
    "ss_model_untransform": [_p, _p, _p],
    "ss_model_loglike": [_p, _p, _pd],
    "ss_model_rep": [_p, _pp],
    "ss_model_rep_at": [_p, _p, _pp],
    "ss_model_set_concentrate": [_p, _i],
    "ss_model_scale": [_p, _pd],
    "ss_model_fit": [_p, _p, _i, _i, _d, _d, _i, _i, _pp],
    "ss_fit_many": [_i, _p, _i, _i, _d, _d, _i, _i, _p, _p],
    "ss_fit_free": [_p],
    "ss_fit_scalars": [_p, _pd, _pd, _pd, _pd, _pi, _pi, _pi, _pi, _pi],
    "ss_fit_get": [_p, _i, _p],
    "ss_fit_message": [_p, _i, ctypes.c_char_p],
    "ss_mapped_new": [_p, _i, _pp],
    "ss_mapped_entry": [_p, _i, _i, _i, _i, _i, _d],
    "ss_mapped_block": [_p, _i, _i, _i, _i],
    "ss_mapped_group": [_p, _i, _i, _i, _d, _d],
    "ss_mapped_start": [_p, _p],
    "ss_model_set_names": [_p, _i, ctypes.c_char_p],
    "ss_callback_new": [_p, _i, _p, _p, _p, _pp],
    "ss_model_failures": [_p, _pi],
    "ss_forecast": [_p, _p, _i, _p, _p],
    "ss_simulate": [_p, _p, _p, _p, _p, _p, _p, _p],
    "ss_simulation_smoother": [_p, _p, _p, _p, _p, _p, _p],
    "ss_steady_state": [_p, _p, _p],
    "ss_standardized_residuals": [_p, _p],
    "ss_diagnostic_start": [_p, _p, _pi],
    "ss_ljung_box": [_i, _p, _i, _i, _p, _p],
    "ss_jarque_bera": [_i, _p, _pd, _pd, _pd, _pd],
    "ss_breakvar": [_i, _p, _i, _pd, _pd],
    # extras
    "ss_fast_smoother": [_p, _p, _p],
    "ss_classical_smoother": [_p, _p, _p, _p],
    "ss_two_filter_smoother": [_p, _p, _p, _p],
    "ss_whittle_smoother": [_p, _p, _p],
    "ss_fixed_point_smoother": [_p, _p, _i, _p, _p],
    "ss_fixed_lag_smoother": [_p, _p, _i, _p, _p],
    "ss_update_smoothed": [_p, _p, _i, _p, _p],
    "ss_smoothed_state_autocov": [_p, _p, _p, _p],
    "ss_smoothed_state_cov_between": [_p, _p, _p, _i, _i, _p],
    "ss_filtered_state_weights": [_p, _p, _p, _p],
    "ss_smoothed_state_weights": [_p, _p, _p, _p],
    "ss_innovation_transition": [_p, _p, _i, _p],
    "ss_sqrt_filter": [_p, _pp],
    "ss_sqrt_smoother": [_p, _p, _pp],
    "ss_augmented_filter": [_p, _pp],
    "ss_augmented_free": [_p],
    "ss_augmented_info": [_p, _pi, _pd, _pd, _pd],
    "ss_augmented_get": [_p, _i, _p],
    "ss_augmented_filter_result": [_p, _pp],
    "ss_augmented_smoother": [_p, _p, _p, _p],
    "ss_em": [_p, _i, _d, _i, _i, _pd, _pi, _p],
    "ss_collapse": [_p, _pp, _p],
    "ss_add_restrictions": [_p, _i, _i, _p, _p, _pp],
    "ss_auxiliary_residuals": [_p, _p, _i, _p, _p],
    "ss_de_jong_penzer": [_p, _p, _p, _p, _p],
    "ss_least_squares_residuals": [_p, _i, _i, _p],
    "ss_r2_diffuse": [_p, _p, _i, _pd],
    "ss_djs_simulation_smoother": [_p, _p, _p, _p, _p, _p, _p, _p],
    "ss_estimation_bias": [_p, _p, _p, _i, _i, _i, _p, _p, _pi],
    "ss_model_components": [_p, _i, _pi, _p, _p, _p],
}

#: C types of the callbacks of ss_callback_new.
UPDATE_CB = ctypes.CFUNCTYPE(_i, _i, _pd, _p)
TRANSFORM_CB = ctypes.CFUNCTYPE(_i, _i, _pd, _pd)
for _name, _args in _SIGNATURES.items():
    _f = getattr(lib, _name)
    _f.argtypes = _args
    _f.restype = _i


def call(name, *args):
    """Call a library routine and raise StateSpaceError on a nonzero status."""
    code = getattr(lib, name)(*args)
    if code != SS_OK:
        raise StateSpaceError(code, name)


def farray(x, shape=None):
    """x as a Fortran-ordered float64 array (a copy if needed)."""
    a = np.asarray(x, dtype=np.float64)
    if shape is not None:
        a = a.reshape(shape, order="F")
    return np.asfortranarray(a)


def ptr(a):
    """The data pointer of a numpy array, for ctypes."""
    return a.ctypes.data_as(ctypes.c_void_p)


def version():
    """Return the version of the Fortran library.

    Returns
    -------
    str
        The version, e.g. ``"0.1.0"``.
    """
    buf = ctypes.create_string_buffer(32)
    call("ss_version", buf, len(buf))
    return buf.value.decode()
