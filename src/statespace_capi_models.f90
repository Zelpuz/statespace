!> C interface to models and estimation: the structural and ARIMA builder,
!> generic model operations (parameters, log likelihood, the representation
!> at given parameters), `fit` and `fit_many`. Conventions as in
!> statespace_capi. Model handles come from ss_struct_build (and the other
!> model constructors) and are released with ss_model_free.
module statespace_capi_models
  use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_char, c_null_ptr, c_loc, &
                                         c_f_pointer, c_associated, c_null_char
  use statespace_kinds, only: dp, SS_OK, SS_ERR_DIM, SS_ERR_NOT_PD, SS_ERR_UNSUPPORTED
  use statespace_rep, only: ssm_rep_t
  use statespace_model, only: ssm_model_t
  use statespace_components, only: component_t, component_holder_t, structural_model_t, &
                                   structural_model
  use statespace_structural, only: irregular_t, level_t, trend_t, seasonal_t, cycle_t, &
                                   regression_t, continuous_level_t, continuous_trend_t
  use statespace_arima, only: arima_t
  use statespace_mle, only: fit_options_t, fit_result_t, fit
  use statespace_mapped, only: mapped_model_t, mapped_model
  use statespace_callback, only: callback_model_t
  use, intrinsic :: iso_c_binding, only: c_funptr, c_null_funptr
  implicit none
  private

  !> A model of any type (structural, mapped, callback).
  type, public :: model_box
    class(ssm_model_t), allocatable :: model
  end type model_box

  type :: builder_box
    real(dp), allocatable :: y(:, :)
    type(component_holder_t), allocatable :: comps(:)
  end type builder_box

  type :: fit_box
    type(fit_result_t) :: res
    integer :: info = SS_OK
  end type fit_box

  public :: get_model
  public :: ss_struct_new, ss_struct_free, ss_struct_add_irregular, ss_struct_add_level
  public :: ss_struct_add_trend, ss_struct_add_seasonal, ss_struct_add_cycle
  public :: ss_struct_add_regression, ss_struct_add_arima, ss_struct_add_continuous
  public :: ss_struct_build
  public :: ss_model_free, ss_model_info, ss_model_param_names, ss_model_start_params
  public :: ss_model_transform, ss_model_untransform, ss_model_loglike, ss_model_rep
  public :: ss_model_rep_at, ss_model_set_concentrate, ss_model_scale
  public :: ss_model_fit, ss_fit_many, ss_fit_free, ss_fit_scalars, ss_fit_get, ss_fit_message
  public :: ss_mapped_new, ss_mapped_entry, ss_mapped_block, ss_mapped_group, ss_mapped_start
  public :: ss_model_set_names, ss_callback_new, ss_model_failures

contains

  ! ------------------------------------------------------------ builder

  function get_builder(handle) result(b)
    type(c_ptr), intent(in) :: handle
    type(builder_box), pointer :: b

    b => null()
    if (c_associated(handle)) call c_f_pointer(handle, b)
  end function get_builder

  !> Start a structural model for data y (p, n); add components, then build.
  integer(c_int) function ss_struct_new(p, n, y, bhandle) bind(C, name="ss_struct_new") &
    result(info)
    integer(c_int), value :: p, n
    type(c_ptr), value :: y
    type(c_ptr), intent(out) :: bhandle
    type(builder_box), pointer :: b
    real(dp), pointer :: yp(:, :)

    bhandle = c_null_ptr
    info = SS_ERR_DIM
    if (p < 1 .or. n < 1 .or. .not. c_associated(y)) return
    allocate (b)
    call c_f_pointer(y, yp, [p, n])
    b%y = yp
    allocate (b%comps(0))
    bhandle = c_loc(b)
    info = SS_OK
  end function ss_struct_new

  integer(c_int) function ss_struct_free(bhandle) bind(C, name="ss_struct_free") result(info)
    type(c_ptr), value :: bhandle
    type(builder_box), pointer :: b

    b => get_builder(bhandle)
    if (associated(b)) deallocate (b)
    info = SS_OK
  end function ss_struct_free

  !> Append a component (moved in).
  subroutine append(b, c, at_obs)
    type(builder_box), intent(inout) :: b
    class(component_t), allocatable, intent(inout) :: c
    integer(c_int), intent(in) :: at_obs
    type(component_holder_t), allocatable :: tmp(:)
    integer :: i, k

    c%at_observations = at_obs /= 0
    k = size(b%comps)
    allocate (tmp(k + 1))
    do i = 1, k
      call move_alloc(b%comps(i)%c, tmp(i)%c)
    end do
    call move_alloc(c, tmp(k + 1)%c)
    call move_alloc(tmp, b%comps)
  end subroutine append

  !> Irregular with covariance type cov (COV_*); weights (n) or NULL.
  integer(c_int) function ss_struct_add_irregular(bhandle, cov, weights) &
    bind(C, name="ss_struct_add_irregular") result(info)
    type(c_ptr), value :: bhandle, weights
    integer(c_int), value :: cov
    type(builder_box), pointer :: b
    class(component_t), allocatable :: c
    real(dp), pointer :: w(:)

    info = SS_ERR_DIM
    b => get_builder(bhandle)
    if (.not. associated(b)) return
    if (c_associated(weights)) then
      call c_f_pointer(weights, w, [size(b%y, 2)])
      c = irregular_t(cov=cov, weights=w)
    else
      c = irregular_t(cov=cov)
    end if
    call append(b, c, 1_c_int)
    info = SS_OK
  end function ss_struct_add_irregular

  integer(c_int) function ss_struct_add_level(bhandle, cov, at_obs) &
    bind(C, name="ss_struct_add_level") result(info)
    type(c_ptr), value :: bhandle
    integer(c_int), value :: cov, at_obs
    type(builder_box), pointer :: b
    class(component_t), allocatable :: c

    info = SS_ERR_DIM
    b => get_builder(bhandle)
    if (.not. associated(b)) return
    c = level_t(cov=cov)
    call append(b, c, at_obs)
    info = SS_OK
  end function ss_struct_add_level

  integer(c_int) function ss_struct_add_trend(bhandle, cov_level, cov_slope, at_obs) &
    bind(C, name="ss_struct_add_trend") result(info)
    type(c_ptr), value :: bhandle
    integer(c_int), value :: cov_level, cov_slope, at_obs
    type(builder_box), pointer :: b
    class(component_t), allocatable :: c

    info = SS_ERR_DIM
    b => get_builder(bhandle)
    if (.not. associated(b)) return
    c = trend_t(cov_level=cov_level, cov_slope=cov_slope)
    call append(b, c, at_obs)
    info = SS_OK
  end function ss_struct_add_trend

  !> Seasonal of period `period`, form SEASONAL_DUMMY/TRIG/HS.
  integer(c_int) function ss_struct_add_seasonal(bhandle, period, form, cov, at_obs) &
    bind(C, name="ss_struct_add_seasonal") result(info)
    type(c_ptr), value :: bhandle
    integer(c_int), value :: period, form, cov, at_obs
    type(builder_box), pointer :: b
    class(component_t), allocatable :: c

    info = SS_ERR_DIM
    b => get_builder(bhandle)
    if (.not. associated(b) .or. period < 2) return
    c = seasonal_t(period=period, form=form, cov=cov)
    call append(b, c, at_obs)
    info = SS_OK
  end function ss_struct_add_seasonal

  !> Cycle; period_max <= 0 means the number of observations.
  integer(c_int) function ss_struct_add_cycle(bhandle, cov, damped, period_min, period_max, &
                                              at_obs) bind(C, name="ss_struct_add_cycle") result(info)
    type(c_ptr), value :: bhandle
    integer(c_int), value :: cov, damped, at_obs
    real(c_double), value :: period_min, period_max
    type(builder_box), pointer :: b
    class(component_t), allocatable :: c

    info = SS_ERR_DIM
    b => get_builder(bhandle)
    if (.not. associated(b)) return
    c = cycle_t(cov=cov, damped=damped /= 0, period_min=period_min, period_max=period_max)
    call append(b, c, at_obs)
    info = SS_OK
  end function ss_struct_add_cycle

  !> Regression effects on series `series` (1-based): x (n, kx);
  !> random_walk (kx, 0/1) or NULL for fixed coefficients.
  integer(c_int) function ss_struct_add_regression(bhandle, kx, x, random_walk, series, at_obs) &
    bind(C, name="ss_struct_add_regression") result(info)
    type(c_ptr), value :: bhandle, x, random_walk
    integer(c_int), value :: kx, series, at_obs
    type(builder_box), pointer :: b
    class(component_t), allocatable :: c
    real(dp), pointer :: xp(:, :)
    integer(c_int), pointer :: rw(:)

    info = SS_ERR_DIM
    b => get_builder(bhandle)
    if (.not. associated(b) .or. kx < 1 .or. .not. c_associated(x)) return
    call c_f_pointer(x, xp, [size(b%y, 2), kx])
    if (c_associated(random_walk)) then
      call c_f_pointer(random_walk, rw, [kx])
      c = regression_t(x=xp, random_walk=rw /= 0, series=series)
    else
      c = regression_t(x=xp, series=series)
    end if
    call append(b, c, at_obs)
    info = SS_OK
  end function ss_struct_add_regression

  !> ARIMA(p, d, q)(P, D, Q)_s on series `series`.
  integer(c_int) function ss_struct_add_arima(bhandle, ar, d, ma, sar, sd, sma, s, series, &
                                              enforce_stationarity, enforce_invertibility, at_obs) &
    bind(C, name="ss_struct_add_arima") result(info)
    type(c_ptr), value :: bhandle
    integer(c_int), value :: ar, d, ma, sar, sd, sma, s, series, enforce_stationarity, &
                             enforce_invertibility, at_obs
    type(builder_box), pointer :: b
    class(component_t), allocatable :: c

    info = SS_ERR_DIM
    b => get_builder(bhandle)
    if (.not. associated(b) .or. sd > 1 .or. min(ar, d, ma, sar, sd, sma, s) < 0) return
    c = arima_t(ar=ar, d=d, ma=ma, sar=sar, sd=sd, sma=sma, s=s, series=series, &
                enforce_stationarity=enforce_stationarity /= 0, &
                enforce_invertibility=enforce_invertibility /= 0)
    call append(b, c, at_obs)
    info = SS_OK
  end function ss_struct_add_arima

  !> Continuous-time level (kind 1) or smooth trend (kind 2) at times (n).
  integer(c_int) function ss_struct_add_continuous(bhandle, kind, times, at_obs) &
    bind(C, name="ss_struct_add_continuous") result(info)
    type(c_ptr), value :: bhandle, times
    integer(c_int), value :: kind, at_obs
    type(builder_box), pointer :: b
    class(component_t), allocatable :: c
    real(dp), pointer :: tp(:)

    info = SS_ERR_DIM
    b => get_builder(bhandle)
    if (.not. associated(b) .or. .not. c_associated(times)) return
    call c_f_pointer(times, tp, [size(b%y, 2)])
    select case (kind)
    case (1)
      c = continuous_level_t(times=tp)
    case (2)
      c = continuous_trend_t(times=tp)
    case default
      return
    end select
    call append(b, c, at_obs)
    info = SS_OK
  end function ss_struct_add_continuous

  !> Build the model. psig = 0 for no signal loadings; otherwise loading
  !> (p, psig) and loading_free (p, psig, 0/1, or NULL).
  integer(c_int) function ss_struct_build(bhandle, psig, loading, loading_free, mhandle) &
    bind(C, name="ss_struct_build") result(info)
    type(c_ptr), value :: bhandle, loading, loading_free
    integer(c_int), value :: psig
    type(c_ptr), intent(out) :: mhandle
    type(builder_box), pointer :: b
    type(model_box), pointer :: mb
    real(dp), pointer :: L(:, :)
    integer(c_int), pointer :: Lf(:, :)
    type(structural_model_t) :: sm
    integer :: stat, p

    mhandle = c_null_ptr
    info = SS_ERR_DIM
    b => get_builder(bhandle)
    if (.not. associated(b)) return
    if (size(b%comps) == 0) return
    p = size(b%y, 1)
    if (psig > 0) then
      if (.not. c_associated(loading)) return
      call c_f_pointer(loading, L, [p, int(psig)])
      if (c_associated(loading_free)) then
        call c_f_pointer(loading_free, Lf, [p, int(psig)])
        sm = structural_model(b%y, b%comps, stat, loading=L, loading_free=Lf /= 0)
      else
        sm = structural_model(b%y, b%comps, stat, loading=L)
      end if
    else
      sm = structural_model(b%y, b%comps, stat)
    end if
    info = stat
    if (stat /= SS_OK) return
    allocate (mb)
    allocate (mb%model, source=sm)
    mhandle = c_loc(mb)
  end function ss_struct_build

  ! ----------------------------------------------- mapped and callback

  !> A mapped model (statespace_mapped) with k parameters over a copy of the
  !> representation rhandle (fixed entries and initialization).
  integer(c_int) function ss_mapped_new(rhandle, k, mhandle) bind(C, name="ss_mapped_new") &
    result(info)
    use statespace_capi, only: get_rep
    type(c_ptr), value :: rhandle
    integer(c_int), value :: k
    type(c_ptr), intent(out) :: mhandle
    type(ssm_rep_t), pointer :: rp
    type(model_box), pointer :: mb

    mhandle = c_null_ptr
    info = SS_ERR_DIM
    rp => get_rep(rhandle)
    if (.not. associated(rp) .or. k < 1) return
    allocate (mb)
    allocate (mb%model, source=mapped_model(rp, int(k)))
    mhandle = c_loc(mb)
    info = SS_OK
  end function ss_mapped_new

  function get_mapped(handle) result(mm)
    type(c_ptr), intent(in) :: handle
    type(mapped_model_t), pointer :: mm
    type(model_box), pointer :: mb

    mm => null()
    mb => get_model(handle)
    if (.not. associated(mb)) return
    select type (m => mb%model)
    type is (mapped_model_t)
      mm => m
    end select
  end function get_mapped

  !> Parameter `param` (1-based) sets X(i, j, t) = coef * param for X the
  !> array `array` (SS_ARR_Z..SS_ARR_D); t = 0 for all time slices; for c
  !> and d, j is ignored.
  integer(c_int) function ss_mapped_entry(mhandle, param, array, i, j, t, coef) &
    bind(C, name="ss_mapped_entry") result(info)
    type(c_ptr), value :: mhandle
    integer(c_int), value :: param, array, i, j, t
    real(c_double), value :: coef
    type(mapped_model_t), pointer :: mm

    info = SS_ERR_DIM
    mm => get_mapped(mhandle)
    if (.not. associated(mm)) return
    if (param < 1 .or. param > mm%k_params .or. array < 2 .or. array > 8) return
    if (.not. in_bounds(mm%rep, array, i, j, t)) return
    call mm%add_entry(int(param), int(array), int(i), int(j), int(t), real(coef, dp))
    info = SS_OK
  end function ss_mapped_entry

  logical function in_bounds(r, array, i, j, t) result(ok)
    type(ssm_rep_t), intent(in) :: r
    integer(c_int), intent(in) :: array, i, j, t
    integer :: s1, s2, st

    select case (array)
    case (2); s1 = size(r%Z, 1); s2 = size(r%Z, 2); st = size(r%Z, 3)
    case (3); s1 = size(r%H, 1); s2 = size(r%H, 2); st = size(r%H, 3)
    case (4); s1 = size(r%T, 1); s2 = size(r%T, 2); st = size(r%T, 3)
    case (5); s1 = size(r%R, 1); s2 = size(r%R, 2); st = size(r%R, 3)
    case (6); s1 = size(r%Q, 1); s2 = size(r%Q, 2); st = size(r%Q, 3)
    case (7); s1 = size(r%c, 1); s2 = 1; st = size(r%c, 2)
    case default; s1 = size(r%d, 1); s2 = 1; st = size(r%d, 2)
    end select
    ok = i >= 1 .and. i <= s1 .and. (array >= 7 .or. (j >= 1 .and. j <= s2)) .and. &
         t >= 0 .and. t <= st
  end function in_bounds

  !> Parameters first..first+q-1 are the lower triangle of a dim x dim
  !> covariance block of H (array 3) or Q (6) at rows/columns offset..
  integer(c_int) function ss_mapped_block(mhandle, array, offset, dim, first) &
    bind(C, name="ss_mapped_block") result(info)
    type(c_ptr), value :: mhandle
    integer(c_int), value :: array, offset, dim, first
    type(mapped_model_t), pointer :: mm
    integer :: n

    info = SS_ERR_DIM
    mm => get_mapped(mhandle)
    if (.not. associated(mm)) return
    if (array /= 3 .and. array /= 6) return
    n = merge(size(mm%rep%H, 1), size(mm%rep%Q, 1), array == 3)
    if (dim < 1 .or. offset < 1 .or. offset + dim - 1 > n) return
    if (first < 1 .or. first + dim * (dim + 1) / 2 - 1 > mm%k_params) return
    call mm%add_block(int(array), int(offset), int(dim), int(first))
    info = SS_OK
  end function ss_mapped_block

  !> A transform group (MAP_POSITIVE, MAP_INTERVAL with lo/hi,
  !> MAP_STATIONARY, MAP_INVERTIBLE) over parameters first..last.
  integer(c_int) function ss_mapped_group(mhandle, kind, first, last, lo, hi) &
    bind(C, name="ss_mapped_group") result(info)
    type(c_ptr), value :: mhandle
    integer(c_int), value :: kind, first, last
    real(c_double), value :: lo, hi
    type(mapped_model_t), pointer :: mm

    info = SS_ERR_DIM
    mm => get_mapped(mhandle)
    if (.not. associated(mm)) return
    if (first < 1 .or. last < first .or. last > mm%k_params .or. kind < 1 .or. kind > 4) return
    call mm%add_group(int(kind), int(first), int(last), real(lo, dp), real(hi, dp))
    info = SS_OK
  end function ss_mapped_group

  !> Start values (constrained) for a mapped or callback model.
  integer(c_int) function ss_mapped_start(mhandle, start) bind(C, name="ss_mapped_start") &
    result(info)
    type(c_ptr), value :: mhandle, start
    type(model_box), pointer :: mb
    real(dp), pointer :: s(:)

    info = SS_ERR_DIM
    mb => get_model(mhandle)
    if (.not. associated(mb) .or. .not. c_associated(start)) return
    call c_f_pointer(start, s, [mb%model%k_params])
    select type (m => mb%model)
    type is (mapped_model_t)
      m%start = s
    type is (callback_model_t)
      m%start = s
    class default
      return
    end select
    info = SS_OK
  end function ss_mapped_start

  !> Parameter names for a mapped model, `width` characters each.
  integer(c_int) function ss_model_set_names(mhandle, width, buf) &
    bind(C, name="ss_model_set_names") result(info)
    type(c_ptr), value :: mhandle
    integer(c_int), value :: width
    character(kind=c_char), intent(in) :: buf(*)
    type(mapped_model_t), pointer :: mm
    integer :: i, j

    info = SS_ERR_DIM
    mm => get_mapped(mhandle)
    if (.not. associated(mm) .or. width < 1) return
    do i = 1, mm%k_params
      mm%names(i) = ''
      do j = 1, min(int(width), 32)
        mm%names(i)(j:j) = buf((i - 1) * width + j)
      end do
    end do
    info = SS_OK
  end function ss_model_set_names

  !> A callback model (statespace_callback) with k parameters over a copy of
  !> the representation rhandle; transform_cb and untransform_cb may be NULL.
  integer(c_int) function ss_callback_new(rhandle, k, update_cb, transform_cb, untransform_cb, &
                                          mhandle) bind(C, name="ss_callback_new") result(info)
    use statespace_capi, only: get_rep
    type(c_ptr), value :: rhandle
    integer(c_int), value :: k
    type(c_funptr), value :: update_cb, transform_cb, untransform_cb
    type(c_ptr), intent(out) :: mhandle
    type(ssm_rep_t), pointer :: rp
    type(model_box), pointer :: mb

    mhandle = c_null_ptr
    info = SS_ERR_DIM
    rp => get_rep(rhandle)
    if (.not. associated(rp) .or. k < 1) return
    allocate (mb)
    allocate (callback_model_t :: mb%model)
    select type (m => mb%model)
    type is (callback_model_t)
      m%rep = rp
      m%k_params = k
      m%update_cb = update_cb
      m%transform_cb = transform_cb
      m%untransform_cb = untransform_cb
      allocate (m%start(k), source=0.1_dp)
      m%rep_handle = c_loc(m%rep)
    end select
    mhandle = c_loc(mb)
    info = SS_OK
  end function ss_callback_new

  !> Number of failed callbacks of a callback model (0 for other models).
  integer(c_int) function ss_model_failures(mhandle, n) bind(C, name="ss_model_failures") &
    result(info)
    type(c_ptr), value :: mhandle
    integer(c_int), intent(out) :: n
    type(model_box), pointer :: mb

    n = 0
    info = SS_ERR_DIM
    mb => get_model(mhandle)
    if (.not. associated(mb)) return
    select type (m => mb%model)
    type is (callback_model_t)
      n = m%failures
    end select
    info = SS_OK
  end function ss_model_failures

  ! ------------------------------------------------------------- models

  function get_model(handle) result(mb)
    type(c_ptr), intent(in) :: handle
    type(model_box), pointer :: mb

    mb => null()
    if (c_associated(handle)) call c_f_pointer(handle, mb)
    if (associated(mb)) then
      if (.not. allocated(mb%model)) mb => null()
    end if
  end function get_model

  integer(c_int) function ss_model_free(mhandle) bind(C, name="ss_model_free") result(info)
    type(c_ptr), value :: mhandle
    type(model_box), pointer :: mb

    if (c_associated(mhandle)) then
      call c_f_pointer(mhandle, mb)
      deallocate (mb)
    end if
    info = SS_OK
  end function ss_model_free

  !> Number of parameters and the dimensions of the model's representation.
  integer(c_int) function ss_model_info(mhandle, k_params, p, m, r, n) &
    bind(C, name="ss_model_info") result(info)
    type(c_ptr), value :: mhandle
    integer(c_int), intent(out) :: k_params, p, m, r, n
    type(model_box), pointer :: mb

    info = SS_ERR_DIM
    mb => get_model(mhandle)
    if (.not. associated(mb)) return
    k_params = mb%model%k_params
    p = mb%model%rep%k_endog; m = mb%model%rep%k_states
    r = mb%model%rep%k_posdef; n = mb%model%rep%nobs
    info = SS_OK
  end function ss_model_info

  !> Parameter names, each in a field of `width` characters (blank padded).
  integer(c_int) function ss_model_param_names(mhandle, width, buf) &
    bind(C, name="ss_model_param_names") result(info)
    type(c_ptr), value :: mhandle
    integer(c_int), value :: width
    character(kind=c_char), intent(out) :: buf(*)
    type(model_box), pointer :: mb
    character(len=32), allocatable :: names(:)
    integer :: i, j

    info = SS_ERR_DIM
    mb => get_model(mhandle)
    if (.not. associated(mb) .or. width < 1) return
    call model_names(mb%model, names)
    do i = 1, size(names)
      do j = 1, width
        if (j <= len(names(i))) then
          buf((i - 1) * width + j) = names(i)(j:j)
        else
          buf((i - 1) * width + j) = ' '
        end if
      end do
    end do
    info = SS_OK
  end function ss_model_param_names

  !> (A separate routine: gfortran 15 -O3 has an internal compiler error on
  !> `names = mb%model%param_names()` in place.)
  subroutine model_names(model, names)
    class(ssm_model_t), intent(in) :: model
    character(len=32), allocatable, intent(out) :: names(:)

    names = model%param_names()
  end subroutine model_names

  integer(c_int) function ss_model_start_params(mhandle, out) &
    bind(C, name="ss_model_start_params") result(info)
    type(c_ptr), value :: mhandle, out
    type(model_box), pointer :: mb
    real(dp), pointer :: o(:)

    info = SS_ERR_DIM
    mb => get_model(mhandle)
    if (.not. associated(mb) .or. .not. c_associated(out)) return
    call c_f_pointer(out, o, [mb%model%k_params])
    o = mb%model%start_params()
    info = SS_OK
  end function ss_model_start_params

  !> Unconstrained x -> constrained params.
  integer(c_int) function ss_model_transform(mhandle, x, out) bind(C, name="ss_model_transform") &
    result(info)
    type(c_ptr), value :: mhandle, x, out
    type(model_box), pointer :: mb
    real(dp), pointer :: xi(:), o(:)

    info = SS_ERR_DIM
    mb => get_model(mhandle)
    if (.not. associated(mb) .or. .not. c_associated(x) .or. .not. c_associated(out)) return
    call c_f_pointer(x, xi, [mb%model%k_params])
    call c_f_pointer(out, o, [mb%model%k_params])
    o = mb%model%transform_params(xi)
    info = SS_OK
  end function ss_model_transform

  !> Constrained params -> unconstrained x.
  integer(c_int) function ss_model_untransform(mhandle, params, out) &
    bind(C, name="ss_model_untransform") result(info)
    type(c_ptr), value :: mhandle, params, out
    type(model_box), pointer :: mb
    real(dp), pointer :: pi(:), o(:)

    info = SS_ERR_DIM
    mb => get_model(mhandle)
    if (.not. associated(mb) .or. .not. c_associated(params) .or. .not. c_associated(out)) return
    call c_f_pointer(params, pi, [mb%model%k_params])
    call c_f_pointer(out, o, [mb%model%k_params])
    o = mb%model%untransform_params(pi)
    info = SS_OK
  end function ss_model_untransform

  !> Log likelihood at constrained params (the model is left at them).
  integer(c_int) function ss_model_loglike(mhandle, params, llf) bind(C, name="ss_model_loglike") &
    result(info)
    type(c_ptr), value :: mhandle, params
    real(c_double), intent(out) :: llf
    type(model_box), pointer :: mb
    real(dp), pointer :: pi(:)
    integer :: stat

    llf = 0.0_dp
    info = SS_ERR_DIM
    mb => get_model(mhandle)
    if (.not. associated(mb) .or. .not. c_associated(params)) return
    call c_f_pointer(params, pi, [mb%model%k_params])
    llf = mb%model%loglike(pi, stat)
    info = stat
  end function ss_model_loglike

  !> A handle to the model's own representation (borrowed: do not free it;
  !> valid while the model exists). Use it to set filter options.
  integer(c_int) function ss_model_rep(mhandle, rhandle) bind(C, name="ss_model_rep") result(info)
    type(c_ptr), value :: mhandle
    type(c_ptr), intent(out) :: rhandle
    type(model_box), pointer :: mb

    rhandle = c_null_ptr
    info = SS_ERR_DIM
    mb => get_model(mhandle)
    if (.not. associated(mb)) return
    rhandle = c_loc(mb%model%rep)
    info = SS_OK
  end function ss_model_rep

  !> A new representation (owned: free with ss_rep_free) at constrained
  !> params, on the data's scale when the scale is concentrated out.
  integer(c_int) function ss_model_rep_at(mhandle, params, rhandle) &
    bind(C, name="ss_model_rep_at") result(info)
    type(c_ptr), value :: mhandle, params
    type(c_ptr), intent(out) :: rhandle
    type(model_box), pointer :: mb
    type(ssm_rep_t), pointer :: rp
    real(dp), pointer :: pi(:)
    integer :: stat

    rhandle = c_null_ptr
    info = SS_ERR_DIM
    mb => get_model(mhandle)
    if (.not. associated(mb) .or. .not. c_associated(params)) return
    call c_f_pointer(params, pi, [mb%model%k_params])
    allocate (rp)
    call mb%model%rep_at(pi, rp, stat)
    info = stat
    if (stat /= SS_OK) then
      deallocate (rp)
      return
    end if
    rhandle = c_loc(rp)
  end function ss_model_rep_at

  integer(c_int) function ss_model_set_concentrate(mhandle, flag) &
    bind(C, name="ss_model_set_concentrate") result(info)
    type(c_ptr), value :: mhandle
    integer(c_int), value :: flag
    type(model_box), pointer :: mb

    info = SS_ERR_DIM
    mb => get_model(mhandle)
    if (.not. associated(mb)) return
    mb%model%concentrate_scale = flag /= 0
    info = SS_OK
  end function ss_model_set_concentrate

  !> The concentrated scale at the last log likelihood evaluation.
  integer(c_int) function ss_model_scale(mhandle, scale) bind(C, name="ss_model_scale") &
    result(info)
    type(c_ptr), value :: mhandle
    real(c_double), intent(out) :: scale
    type(model_box), pointer :: mb

    info = SS_ERR_DIM
    mb => get_model(mhandle)
    if (.not. associated(mb)) return
    scale = mb%model%scale
    info = SS_OK
  end function ss_model_scale

  ! ---------------------------------------------------------- estimation

  function options(maxiter, m, factr, pgtol, compute_cov, gradient) result(opts)
    integer(c_int), intent(in) :: maxiter, m, compute_cov, gradient
    real(c_double), intent(in) :: factr, pgtol
    type(fit_options_t) :: opts

    opts%maxiter = maxiter
    opts%m = m
    opts%factr = factr
    opts%pgtol = pgtol
    opts%compute_cov = compute_cov /= 0
    opts%gradient = gradient
  end function options

  !> Maximum likelihood (statespace_mle's fit). start is the constrained
  !> starting point or NULL. fhandle receives the result whenever the
  !> estimates exist (also with SS_ERR_NOT_PD: no covariance); free it with
  !> ss_fit_free.
  integer(c_int) function ss_model_fit(mhandle, start, maxiter, m, factr, pgtol, compute_cov, &
                                       gradient, fhandle) bind(C, name="ss_model_fit") result(info)
    type(c_ptr), value :: mhandle, start
    integer(c_int), value :: maxiter, m, compute_cov, gradient
    real(c_double), value :: factr, pgtol
    type(c_ptr), intent(out) :: fhandle
    type(model_box), pointer :: mb
    type(fit_box), pointer :: fb
    real(dp), pointer :: s(:)
    integer :: stat

    fhandle = c_null_ptr
    info = SS_ERR_DIM
    mb => get_model(mhandle)
    if (.not. associated(mb)) return
    allocate (fb)
    if (c_associated(start)) then
      call c_f_pointer(start, s, [mb%model%k_params])
      call fit(mb%model, fb%res, start_params=s, &
               options=options(maxiter, m, factr, pgtol, compute_cov, gradient), info=stat)
    else
      call fit(mb%model, fb%res, options=options(maxiter, m, factr, pgtol, compute_cov, gradient), &
               info=stat)
    end if
    info = stat
    fb%info = stat
    if (stat /= SS_OK .and. stat /= SS_ERR_NOT_PD) then
      deallocate (fb)
      return
    end if
    fhandle = c_loc(fb)
  end function ss_model_fit

  !> Fit nm models (handles in mhandles) in parallel with OpenMP. Results go
  !> to fhandles and status codes to infos, as for ss_model_fit. Returns
  !> SS_OK if every handle was valid.
  integer(c_int) function ss_fit_many(nm, mhandles, maxiter, m, factr, pgtol, compute_cov, &
                                      gradient, fhandles, infos) bind(C, name="ss_fit_many") &
    result(info)
    integer(c_int), value :: nm, maxiter, m, compute_cov, gradient
    real(c_double), value :: factr, pgtol
    type(c_ptr), intent(in) :: mhandles(nm)
    type(c_ptr), intent(out) :: fhandles(nm)
    integer(c_int), intent(out) :: infos(nm)
    type(fit_options_t) :: opts
    type(model_box), pointer :: mb
    integer :: i

    fhandles = c_null_ptr
    infos = SS_ERR_DIM
    info = SS_ERR_DIM
    do i = 1, nm
      mb => get_model(mhandles(i))
      if (.not. associated(mb)) return
      ! Callbacks into another runtime cannot run on OpenMP threads.
      select type (mdl => mb%model)
      type is (callback_model_t)
        info = SS_ERR_UNSUPPORTED
        return
      end select
    end do
    opts = options(maxiter, m, factr, pgtol, compute_cov, gradient)
    !$omp parallel do schedule(dynamic)
    do i = 1, nm
      call fit_one(i)
    end do
    !$omp end parallel do
    info = SS_OK

  contains

    subroutine fit_one(i)
      integer, intent(in) :: i
      type(model_box), pointer :: mb
      type(fit_box), pointer :: fb
      integer :: stat

      mb => get_model(mhandles(i))
      allocate (fb)
      call fit(mb%model, fb%res, options=opts, info=stat)
      infos(i) = stat
      fb%info = stat
      if (stat /= SS_OK .and. stat /= SS_ERR_NOT_PD) then
        deallocate (fb)
      else
        fhandles(i) = c_loc(fb)
      end if
    end subroutine fit_one
  end function ss_fit_many

  function get_fit(handle) result(fb)
    type(c_ptr), intent(in) :: handle
    type(fit_box), pointer :: fb

    fb => null()
    if (c_associated(handle)) call c_f_pointer(handle, fb)
  end function get_fit

  integer(c_int) function ss_fit_free(fhandle) bind(C, name="ss_fit_free") result(info)
    type(c_ptr), value :: fhandle
    type(fit_box), pointer :: fb

    fb => get_fit(fhandle)
    if (associated(fb)) deallocate (fb)
    info = SS_OK
  end function ss_fit_free

  !> Scalars of a fit; has_cov = 1 when cov_params and bse are available.
  integer(c_int) function ss_fit_scalars(fhandle, llf, scale, aic, bic, niter, nfev, converged, &
                                         analytic_gradient, has_cov) &
    bind(C, name="ss_fit_scalars") result(info)
    type(c_ptr), value :: fhandle
    real(c_double), intent(out) :: llf, scale, aic, bic
    integer(c_int), intent(out) :: niter, nfev, converged, analytic_gradient, has_cov
    type(fit_box), pointer :: fb

    info = SS_ERR_DIM
    fb => get_fit(fhandle)
    if (.not. associated(fb)) return
    llf = fb%res%llf; scale = fb%res%scale; aic = fb%res%aic; bic = fb%res%bic
    niter = fb%res%niter; nfev = fb%res%nfev
    converged = merge(1, 0, fb%res%converged)
    analytic_gradient = merge(1, 0, fb%res%analytic_gradient)
    has_cov = merge(1, 0, allocated(fb%res%bse))
    info = SS_OK
  end function ss_fit_scalars

  !> Copy params (1), bse (2) or cov_params (3, k x k) into out.
  integer(c_int) function ss_fit_get(fhandle, code, out) bind(C, name="ss_fit_get") result(info)
    type(c_ptr), value :: fhandle, out
    integer(c_int), value :: code
    type(fit_box), pointer :: fb
    real(dp), pointer :: o(:)
    integer :: k

    info = SS_ERR_DIM
    fb => get_fit(fhandle)
    if (.not. associated(fb) .or. .not. c_associated(out)) return
    k = size(fb%res%params)
    select case (code)
    case (1)
      call c_f_pointer(out, o, [k]); o = fb%res%params
    case (2)
      if (.not. allocated(fb%res%bse)) return
      call c_f_pointer(out, o, [k]); o = fb%res%bse
    case (3)
      if (.not. allocated(fb%res%cov_params)) return
      call c_f_pointer(out, o, [k * k]); o = reshape(fb%res%cov_params, [k * k])
    case default
      return
    end select
    info = SS_OK
  end function ss_fit_get

  !> The optimizer's final message, NUL-terminated, into buf(len).
  integer(c_int) function ss_fit_message(fhandle, len, buf) bind(C, name="ss_fit_message") &
    result(info)
    type(c_ptr), value :: fhandle
    integer(c_int), value :: len
    character(kind=c_char), intent(out) :: buf(len)
    type(fit_box), pointer :: fb
    integer :: i, k

    info = SS_ERR_DIM
    fb => get_fit(fhandle)
    if (.not. associated(fb) .or. len < 1) return
    k = min(len - 1, len_trim(fb%res%message))
    do i = 1, k
      buf(i) = fb%res%message(i:i)
    end do
    buf(k + 1) = c_null_char
    info = SS_OK
  end function ss_fit_message
end module statespace_capi_models
