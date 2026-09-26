C interface
===========

The modules ``statespace_capi``, ``statespace_capi_models`` and
``statespace_capi_extras`` export the library through ``bind(C)`` routines
in the shared library ``libstatespace``. The Python package uses them
through ctypes; any language with a C foreign function interface can.

Conventions
-----------

* **Handles.** Objects (representations, filter and smoother results,
  models, fits) are opaque ``void *`` handles. A routine that creates an
  object returns its handle through an output argument; release it with the
  matching ``*_free``. A handle obtained from ``ss_model_rep`` or
  ``ss_augmented_filter_result`` is borrowed: do not free it.
* **Arrays.** Arrays are ``double`` in Fortran (column-major) order with time
  as the last dimension. Inputs are copied in. Outputs are copied into
  buffers the caller allocates; their sizes follow from the dimensions
  (``ss_rep_info``, ``ss_model_info``) and are given below.
* **Status.** Every routine returns ``int``: 0 on success or a status code
  of :doc:`fortran/statespace_kinds`.
* **Indices.** Periods, states and parameters passed as arguments are
  1-based.
* **Missing values** are NaN in ``y``.
* **Threads.** Different handles can be used from different threads.
  ``ss_fit_many`` runs its own OpenMP threads and rejects callback models.

Example
-------

The local level model in C:

.. code-block:: c

   void *rep, *fres;
   double one = 1.0, h = 15099.0, q = 1469.1, llf;
   ss_rep_new(1, 1, 1, n, &rep);               /* p, m, r, n */
   ss_rep_set(rep, 1, n, y);                   /* SS_ARR_Y */
   ss_rep_set(rep, 2, 1, &one);                /* Z */
   ss_rep_set(rep, 4, 1, &one);                /* T */
   ss_rep_set(rep, 3, 1, &h);                  /* H */
   ss_rep_set(rep, 6, 1, &q);                  /* Q */
   ss_rep_init_diffuse(rep);
   ss_loglike(rep, &llf);
   ss_filter(rep, &fres);
   /* ... ss_filter_get(fres, SS_F_A, buffer) ... */
   ss_filter_free(fres);
   ss_rep_free(rep);

Codes
-----

Arrays (``ss_rep_set``, ``ss_rep_get``, ``ss_mapped_entry``): 1 ``y``
(p, n), 2 Z (p, m, nt), 3 H (p, p, nt), 4 T (m, m, nt), 5 R (m, r, nt),
6 Q (r, r, nt), 7 c (m, nt), 8 d (p, nt); nt is 1 or n.

Options: integer (``ss_rep_set_int``) 1 filter method, 2 diffuse method,
3 log likelihood burn-in, 4 marginal likelihood (0/1); real
(``ss_rep_set_real``) 1 ``tol_diffuse``, 2 ``tol_steady``.

Filter outputs (``ss_filter_get``): 1 a (m, n+1), 2 P (m, m, n+1), 3 P∞
(m, m, n+1), 4 a_t|t (m, n), 5 P_t|t (m, m, n), 6 ŷ (p, n), 7 v (p, n),
8 F (p, p, n), 9 F∞ (p, p, n), 10 F⁻¹ (p, p, n), 11 K (m, p, n),
12 log likelihood contributions (n).

Smoother outputs (``ss_smoother_get``): 1 α̂ (m, n), 2 V (m, m, n), 3 r
(m, n+1), 4 N (m, m, n+1), 5 ε̂ (p, n), 6 Var(ε|Y) (p, p, n), 7 η̂ (r, n),
8 Var(η|Y) (r, r, n).

Fit outputs (``ss_fit_get``): 1 parameters (k), 2 standard errors (k),
3 covariance (k, k).

Representations
---------------

.. code-block:: c

   int ss_version(char *buf, int len);
   int ss_rep_new(int p, int m, int r, int n, void **rep);
   int ss_rep_free(void *rep);
   int ss_rep_set(void *rep, int code, int nt, const double *data);
   int ss_rep_get(void *rep, int code, double *out);
   int ss_rep_info(void *rep, int *p, int *m, int *r, int *n, int nts[8]);
   int ss_rep_set_int(void *rep, int code, int value);
   int ss_rep_set_real(void *rep, int code, double value);
   int ss_rep_init_known(void *rep, const double *a1, const double *P1);
   int ss_rep_init_diffuse(void *rep);
   int ss_rep_init_approx_diffuse(void *rep, double kappa);
   int ss_rep_init_stationary(void *rep);
   int ss_rep_init_general(void *rep, const double *a1, const double *Pstar,
                           const double *Pinf);
   int ss_rep_init_block(void *rep, int first, int last, int kind,
                         const double *a1, const double *P1, double kappa);

``ss_rep_new`` creates a time-invariant representation with zero matrices
(R the leading m × r identity) and y = 0. ``ss_rep_info`` returns the
dimensions and the time dimension of each array. In ``ss_rep_init_block``
pass NULL for ``a1`` and ``P1`` unless ``kind`` is ``INIT_KNOWN``.

Filtering and smoothing
-----------------------

.. code-block:: c

   int ss_loglike(void *rep, double *llf);
   int ss_loglike_concentrated(void *rep, double *llf, double *scale);
   int ss_filter(void *rep, void **fres);
   int ss_sqrt_filter(void *rep, void **fres);
   int ss_filter_free(void *fres);
   int ss_filter_scalars(void *fres, double *llf, int *nobs_diffuse,
                         int *k_diffuse, int *t_steady);
   int ss_filter_get(void *fres, int code, double *out);
   int ss_smooth(void *rep, void *fres, void **sres);
   int ss_sqrt_smoother(void *rep, void *fres, void **sres);
   int ss_smoother_free(void *sres);
   int ss_smoother_get(void *sres, int code, double *out);
   int ss_steady_state(void *rep, double *P, double *F);

The smoothing extras of DK ch. 4 (:doc:`fortran/statespace_smoothing`)
take the representation and a filter result:

.. code-block:: c

   int ss_fast_smoother(void *rep, void *fres, double *alphahat);
   int ss_classical_smoother(void *rep, void *fres, double *alphahat, double *V);
   int ss_two_filter_smoother(void *rep, void *fres, double *alphahat, double *V);
   int ss_whittle_smoother(void *rep, void *fres, double *alphahat);
   int ss_fixed_point_smoother(void *rep, void *fres, int t,
                               double *path_a, double *path_V);   /* (m, n-t+1) */
   int ss_fixed_lag_smoother(void *rep, void *fres, int j, double *alphahat, double *V);
   int ss_update_smoothed(void *rep, void *fres, int n0, double *alphahat, double *V);
   int ss_smoothed_state_autocov(void *rep, void *fres, void *sres,
                                 double *acov);                    /* (m, m, n-1) */
   int ss_smoothed_state_cov_between(void *rep, void *fres, void *sres,
                                     int t, int j, double *C);
   int ss_filtered_state_weights(void *rep, void *fres, double *Wa,
                                 double *Watt);                    /* (m, p, n, n) */
   int ss_smoothed_state_weights(void *rep, double *W, double *C, double *A);
   int ss_innovation_transition(void *rep, void *fres, int t, double *L);

Augmented filter, EM, collapsing, restrictions
----------------------------------------------

.. code-block:: c

   int ss_augmented_filter(void *rep, void **aug);
   int ss_augmented_free(void *aug);
   int ss_augmented_info(void *aug, int *k, double *llf, double *llf_fixed,
                         double *llf_marginal);
   int ss_augmented_get(void *aug, int code, double *out);   /* 1 delta (k), 2 cov (k, k) */
   int ss_augmented_filter_result(void *aug, void **fres);   /* borrowed */
   int ss_augmented_smoother(void *rep, void *aug, double *alphahat, double *V);
   int ss_em(void *rep, int maxiter, double tol, int diagonal_H, int diagonal_Q,
             double *llf, int *niter, double *path);         /* path (maxiter) or NULL */
   int ss_collapse(void *rep, void **collapsed, double *llf_adjust);
   int ss_add_restrictions(void *rep, int q, int nt, const double *Rmat,
                           const double *rval, void **restricted);

Forecasting, simulation, diagnostics
------------------------------------

.. code-block:: c

   int ss_forecast(void *rep, void *fres, int h, double *mean, double *cov);
   int ss_simulate(void *rep, const double *u_init, const double *u_eps,
                   const double *u_eta, double *y, double *alpha, double *eps,
                   double *eta);
   int ss_simulation_smoother(void *rep, const double *u_init, const double *u_eps,
                              const double *u_eta, double *state, double *eps,
                              double *eta);
   int ss_djs_simulation_smoother(void *rep, void *fres, const double *u_eps,
                                  const double *u_eta, const double *u_init,
                                  double *eps, double *eta, double *alpha);
   int ss_standardized_residuals(void *fres, double *out);
   int ss_diagnostic_start(void *rep, void *fres, int *t0);
   int ss_auxiliary_residuals(void *rep, void *sres, int vector, double *eps_std,
                              double *eta_std);
   int ss_de_jong_penzer(void *rep, void *fres, void *sres, double *r_stat,
                         double *e_stat);
   int ss_least_squares_residuals(void *rep, int first, int last, double *vplus);
   int ss_r2_diffuse(void *rep, void *fres, int i, double *r2);
   int ss_ljung_box(int n, const double *x, int lags, int model_df, double *stat,
                    double *pvalue);
   int ss_jarque_bera(int n, const double *x, double *jb, double *pvalue,
                      double *skew, double *kurtosis);
   int ss_breakvar(int n, const double *x, int h, double *stat, double *pvalue);

The variates of the simulation routines are standard normal: ``u_init``
(m), ``u_eps`` (p, n), ``u_eta`` (r, n).

Structural and ARIMA models
---------------------------

A builder collects components; ``ss_struct_build`` makes the model.

.. code-block:: c

   int ss_struct_new(int p, int n, const double *y, void **builder);
   int ss_struct_free(void *builder);
   int ss_struct_add_irregular(void *b, int cov, const double *weights);  /* or NULL */
   int ss_struct_add_level(void *b, int cov, int at_obs);
   int ss_struct_add_trend(void *b, int cov_level, int cov_slope, int at_obs);
   int ss_struct_add_seasonal(void *b, int period, int form, int cov, int at_obs);
   int ss_struct_add_cycle(void *b, int cov, int damped, double period_min,
                           double period_max, int at_obs);
   int ss_struct_add_regression(void *b, int kx, const double *x,
                                const int *random_walk, int series, int at_obs);
   int ss_struct_add_arima(void *b, int p, int d, int q, int P, int D, int Q, int s,
                           int series, int enforce_stationarity,
                           int enforce_invertibility, int at_obs);
   int ss_struct_add_continuous(void *b, int kind, const double *times, int at_obs);
   int ss_struct_build(void *b, int psig, const double *loading,
                       const int *loading_free, void **model);

``cov`` is 0 none, 1 diagonal, 2 full; ``form`` 1 dummy, 2 trigonometric,
3 Harrison-Stevens; ``kind`` 1 continuous level, 2 continuous trend.
``psig`` = 0 for no loadings.

Mapped and callback models
--------------------------

.. code-block:: c

   int ss_mapped_new(void *rep, int k, void **model);
   int ss_mapped_entry(void *model, int param, int array, int i, int j, int t,
                       double coef);
   int ss_mapped_block(void *model, int array, int offset, int dim, int first);
   int ss_mapped_group(void *model, int kind, int first, int last, double lo,
                       double hi);
   int ss_mapped_start(void *model, const double *start);
   int ss_model_set_names(void *model, int width, const char *names);

   typedef int (*ss_update_cb)(int k, const double *params, void *rep);
   typedef int (*ss_transform_cb)(int k, const double *x, double *out);
   int ss_callback_new(void *rep, int k, ss_update_cb update,
                       ss_transform_cb transform, ss_transform_cb untransform,
                       void **model);
   int ss_model_failures(void *model, int *n);

See :doc:`fortran/statespace_mapped` and :doc:`fortran/statespace_callback`.
``ss_mapped_start`` also sets the start values of a callback model.

Models and estimation
---------------------

.. code-block:: c

   int ss_model_free(void *model);
   int ss_model_info(void *model, int *k_params, int *p, int *m, int *r, int *n);
   int ss_model_param_names(void *model, int width, char *buf);   /* k * width chars */
   int ss_model_start_params(void *model, double *out);
   int ss_model_transform(void *model, const double *x, double *out);
   int ss_model_untransform(void *model, const double *params, double *out);
   int ss_model_loglike(void *model, const double *params, double *llf);
   int ss_model_rep(void *model, void **rep);                     /* borrowed */
   int ss_model_rep_at(void *model, const double *params, void **rep);
   int ss_model_set_concentrate(void *model, int flag);
   int ss_model_scale(void *model, double *scale);
   int ss_model_components(void *model, int maxcomp, int *ncomp, int *first,
                           int *size, int *at_obs);
   int ss_model_fit(void *model, const double *start, int maxiter, int m,
                    double factr, double pgtol, int compute_cov, int gradient,
                    void **fit);
   int ss_fit_many(int nm, void *const *models, int maxiter, int m, double factr,
                   double pgtol, int compute_cov, int gradient, void **fits,
                   int *infos);
   int ss_fit_free(void *fit);
   int ss_fit_scalars(void *fit, double *llf, double *scale, double *aic,
                      double *bic, int *niter, int *nfev, int *converged,
                      int *analytic_gradient, int *has_cov);
   int ss_fit_get(void *fit, int code, double *out);
   int ss_fit_message(void *fit, int len, char *buf);
   int ss_estimation_bias(void *model, const double *params, const double *cov,
                          int ndraw, int antithetic, int seed, double *bias_alpha,
                          double *bias_V, int *failed);

``ss_model_fit`` returns a fit handle also when the status is 2 (estimates
without a covariance matrix). ``ss_fit_many`` returns 0 when every handle is
valid and sets each fit's status in ``infos``; ``fits[i]`` is NULL where a
fit failed.
