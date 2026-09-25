!> Model components and their assembly into a state space model (DK ch. 3).
!>
!> A component owns a block of the state vector: its columns of Z, its
!> diagonal blocks of T, R and Q, and optionally a contribution to H (the
!> irregular). `structural_model` stacks components as DK do in (3.10):
!>
!>     Z_t = (Z[1], Z[2], ...),   T_t = diag(T[1], T[2], ...),
!>     R_t = diag(R[1], R[2], ...),   Q_t = diag(Q[1], Q[2], ...),
!>
!> and initializes each block as the component asks (diffuse for
!> nonstationary components, stationary for ARMA and damped cycles, DK 5.6).
!> The result is an ordinary `ssm_model_t`, so fitting, smoothing and the
!> diagnostics apply unchanged.
!>
!> Signal loadings (DK 3.3.2, 3.3.3, 3.7): with `loading` (p x p_sig), the
!> components describe p_sig signals theta_t and the observations load on
!> them, Z_t = Lambda Z_sig,t; `loading_free` marks entries of Lambda that
!> are parameters (appended after the components' parameters). Components
!> whose `observation_level()` is true (the irregular, per-series
!> intercepts) act on the observations directly instead.
!>
!> To write a new component, extend `component_t` and implement `setup`
!> (set m, r, k and which matrices vary over time) and `fill` (write the
!> blocks from constrained parameters); override `transform_params`,
!> `untransform_params`, `start_params`, `param_names` and `init_blocks` as
!> needed. Parameters are ordered component by component.
module statespace_components
  use statespace_kinds, only: dp, SS_OK, SS_ERR_DIM
  use statespace_rep, only: ssm_rep_t, ssm_rep, INIT_DIFFUSE, INIT_STATIONARY, INIT_KNOWN
  use statespace_model, only: ssm_model_t, constrain_positive, unconstrain_positive
  implicit none
  private

  public :: component_t, component_holder_t, structural_model_t, structural_model
  public :: constrain_interval, unconstrain_interval, diff_variance

  type, abstract :: component_t
    integer :: p = 0        !< observed series
    integer :: n = 0        !< observations
    integer :: m = 0        !< states
    integer :: r = 0        !< state disturbances
    integer :: k = 0        !< parameters
    logical :: tv_Z = .false., tv_H = .false., tv_T = .false., tv_R = .false., tv_Q = .false.
    !> Place this component at the observations rather than the signals.
    logical :: at_observations = .false.
  contains
    procedure :: observation_level => comp_observation_level
    procedure(setup_iface), deferred :: setup
    procedure(fill_iface), deferred :: fill
    procedure :: transform_params => comp_transform
    procedure :: untransform_params => comp_untransform
    procedure :: start_params => comp_start
    procedure :: param_names => comp_names
    procedure :: init_blocks => comp_init_blocks
  end type component_t

  abstract interface
    !> Set m, r, k and the time-variation flags for p series of n observations.
    subroutine setup_iface(self, y)
      import :: component_t, dp
      class(component_t), intent(inout) :: self
      real(dp), intent(in) :: y(:, :)
    end subroutine setup_iface

    !> Write this component's blocks from constrained parameters. Z, T, R, Q
    !> are the component's own blocks (third dimension 1 or n as flagged);
    !> H is the full observation covariance, to be added to.
    subroutine fill_iface(self, params, Z, H, T, R, Q)
      import :: component_t, dp
      class(component_t), intent(in) :: self
      real(dp), intent(in) :: params(:)
      real(dp), intent(inout) :: Z(:, :, :), H(:, :, :), T(:, :, :), R(:, :, :), Q(:, :, :)
    end subroutine fill_iface
  end interface

  type :: component_holder_t
    class(component_t), allocatable :: c
  end type component_holder_t

  !> A model assembled from components.
  type, extends(ssm_model_t) :: structural_model_t
    type(component_holder_t), allocatable :: comps(:)
    integer, allocatable :: s0(:), e0(:), k0(:)   !< offsets of states, disturbances, params
    real(dp), allocatable :: loading(:, :)        !< (p, p_sig) Lambda
    logical, allocatable :: loading_free(:, :)    !< entries of Lambda that are parameters
    integer :: k_comp = 0                         !< parameters of the components
    real(dp), allocatable :: Zsig(:, :, :)        !< (p_sig, m, 1|n) signal design
  contains
    procedure :: update => sm_update
    procedure :: start_params => sm_start
    procedure :: transform_params => sm_transform
    procedure :: untransform_params => sm_untransform
    procedure :: param_names => sm_names
  end type structural_model_t

contains

  !> Assemble the components into a model for data y (p, n).
  function structural_model(y, comps, info, loading, loading_free) result(mod)
    real(dp), intent(in) :: y(:, :)
    type(component_holder_t), intent(in) :: comps(:)
    integer, intent(out) :: info
    real(dp), intent(in), optional :: loading(:, :)       !< (p, p_sig)
    logical, intent(in), optional :: loading_free(:, :)   !< (p, p_sig)
    type(structural_model_t) :: mod
    real(dp), allocatable :: ysig(:, :)
    integer :: i, j, m, r, k, n, nz, nh, nt, nr, nq, p, psig
    integer, allocatable :: blocks(:, :)

    n = size(y, 2)
    p = size(y, 1)
    info = SS_ERR_DIM
    if (present(loading)) then
      if (size(loading, 1) /= p) return
      mod%loading = loading
    else
      mod%loading = identity(p)
    end if
    psig = size(mod%loading, 2)
    allocate (mod%loading_free(p, psig), source=.false.)
    if (present(loading_free)) mod%loading_free = loading_free
    ! Signals get the leading observed series for start values.
    allocate (ysig(psig, n), source=0.0_dp)
    ysig(1:min(p, psig), :) = y(1:min(p, psig), :)
    mod%comps = comps
    allocate (mod%s0(size(comps)), mod%e0(size(comps)), mod%k0(size(comps)))
    m = 0; r = 0; k = 0
    nz = 1; nh = 1; nt = 1; nr = 1; nq = 1
    do i = 1, size(comps)
      associate (c => mod%comps(i)%c)
        if (c%observation_level()) then
          call c%setup(y)
        else
          call c%setup(ysig)
        end if
        mod%s0(i) = m; mod%e0(i) = r; mod%k0(i) = k
        m = m + c%m; r = r + c%r; k = k + c%k
        if (c%tv_Z) nz = n
        if (c%tv_H) nh = n
        if (c%tv_T) nt = n
        if (c%tv_R) nr = n
        if (c%tv_Q) nq = n
      end associate
    end do
    if (m == 0) return
    mod%k_comp = k
    mod%k_params = k + count(mod%loading_free)
    allocate (mod%Zsig(psig, m, nz), source=0.0_dp)
    mod%rep = ssm_rep(y, m, max(r, 1))
    mod%rep%R = 0.0_dp
    if (nz > 1) mod%rep%Z = spread(mod%rep%Z(:, :, 1), 3, nz)
    if (nh > 1) mod%rep%H = spread(mod%rep%H(:, :, 1), 3, nh)
    if (nt > 1) mod%rep%T = spread(mod%rep%T(:, :, 1), 3, nt)
    if (nr > 1) mod%rep%R = spread(mod%rep%R(:, :, 1), 3, nr)
    if (nq > 1) mod%rep%Q = spread(mod%rep%Q(:, :, 1), 3, nq)
    do i = 1, size(comps)
      associate (c => mod%comps(i)%c)
        if (c%m == 0) cycle
        blocks = c%init_blocks()
        do j = 1, size(blocks, 1)
          call mod%rep%initialize_block(mod%s0(i) + blocks(j, 1), mod%s0(i) + blocks(j, 2), &
                                        blocks(j, 3))
        end do
      end associate
    end do
    call mod%update(mod%start_params())
    info = SS_OK
  end function structural_model

  subroutine sm_update(self, params)
    class(structural_model_t), intent(inout) :: self
    real(dp), intent(in) :: params(:)
    real(dp), allocatable :: Zb(:, :, :), Tb(:, :, :), Rb(:, :, :), Qb(:, :, :)
    integer :: i, s, e, k, t

    self%rep%H = 0.0_dp
    do i = 1, size(self%comps)
      associate (c => self%comps(i)%c)
        s = self%s0(i); e = self%e0(i); k = self%k0(i)
        ! Work on copies of the blocks: sections with a time dimension of 1
        ! or n, as the component expects.
        if (c%observation_level()) then
          Zb = self%rep%Z(:, s + 1:s + c%m, 1:merge(size(self%rep%Z, 3), 1, c%tv_Z))
        else
          Zb = self%Zsig(:, s + 1:s + c%m, 1:merge(size(self%Zsig, 3), 1, c%tv_Z))
        end if
        Tb = self%rep%T(s + 1:s + c%m, s + 1:s + c%m, 1:merge(size(self%rep%T, 3), 1, c%tv_T))
        Rb = self%rep%R(s + 1:s + c%m, e + 1:e + c%r, 1:merge(size(self%rep%R, 3), 1, c%tv_R))
        Qb = self%rep%Q(e + 1:e + c%r, e + 1:e + c%r, 1:merge(size(self%rep%Q, 3), 1, c%tv_Q))
        Zb = 0.0_dp; Tb = 0.0_dp; Rb = 0.0_dp; Qb = 0.0_dp
        call c%fill(params(k + 1:k + c%k), Zb, self%rep%H, Tb, Rb, Qb)
        if (c%observation_level()) then
          self%rep%Z(:, s + 1:s + c%m, :) = broadcast(Zb, size(self%rep%Z, 3))
          self%Zsig(:, s + 1:s + c%m, :) = 0.0_dp
        else
          self%Zsig(:, s + 1:s + c%m, :) = broadcast(Zb, size(self%Zsig, 3))
        end if
        self%rep%T(s + 1:s + c%m, s + 1:s + c%m, :) = broadcast(Tb, size(self%rep%T, 3))
        self%rep%R(s + 1:s + c%m, e + 1:e + c%r, :) = broadcast(Rb, size(self%rep%R, 3))
        self%rep%Q(e + 1:e + c%r, e + 1:e + c%r, :) = broadcast(Qb, size(self%rep%Q, 3))
      end associate
    end do
    ! Z = Lambda Z_sig for the signal columns (observation-level columns of
    ! Z_sig are zero, so they are unaffected).
    if (any(self%loading_free)) self%loading = unpack(params(self%k_comp + 1:), &
                                                      self%loading_free, self%loading)
    do t = 1, size(self%rep%Z, 3)
      do i = 1, size(self%comps)
        if (self%comps(i)%c%observation_level()) cycle
        s = self%s0(i)
        self%rep%Z(:, s + 1:s + self%comps(i)%c%m, t) = &
          matmul(self%loading, self%Zsig(:, s + 1:s + self%comps(i)%c%m, t))
      end do
    end do
  contains
    !> A block with time dimension 1 copied to all nt slices, or as is.
    function broadcast(B, nt) result(Bn)
      real(dp), intent(in) :: B(:, :, :)
      integer, intent(in) :: nt
      real(dp), allocatable :: Bn(:, :, :)

      if (size(B, 3) == nt) then
        Bn = B
      else
        Bn = spread(B(:, :, 1), 3, nt)
      end if
    end function broadcast
  end subroutine sm_update

  function sm_start(self) result(params)
    class(structural_model_t), intent(in) :: self
    real(dp), allocatable :: params(:)
    integer :: i

    allocate (params(0))
    do i = 1, size(self%comps)
      params = [params, self%comps(i)%c%start_params(signal_data(self, i))]
    end do
    params = [params, pack(self%loading, self%loading_free)]
  end function sm_start

  !> The data a component was set up with: the observations, or the leading
  !> series standing in for the signals.
  function signal_data(self, i) result(y)
    class(structural_model_t), intent(in) :: self
    integer, intent(in) :: i
    real(dp), allocatable :: y(:, :)
    integer :: psig, p

    if (self%comps(i)%c%observation_level()) then
      y = self%rep%y
    else
      p = self%rep%k_endog
      psig = size(self%loading, 2)
      allocate (y(psig, self%rep%nobs), source=0.0_dp)
      y(1:min(p, psig), :) = self%rep%y(1:min(p, psig), :)
    end if
  end function signal_data

  function sm_transform(self, unconstrained) result(constrained)
    class(structural_model_t), intent(in) :: self
    real(dp), intent(in) :: unconstrained(:)
    real(dp), allocatable :: constrained(:)
    integer :: i, k

    constrained = unconstrained
    do i = 1, size(self%comps)
      k = self%k0(i)
      associate (c => self%comps(i)%c)
        if (c%k > 0) constrained(k + 1:k + c%k) = c%transform_params(unconstrained(k + 1:k + c%k))
      end associate
    end do
  end function sm_transform

  function sm_untransform(self, constrained) result(unconstrained)
    class(structural_model_t), intent(in) :: self
    real(dp), intent(in) :: constrained(:)
    real(dp), allocatable :: unconstrained(:)
    integer :: i, k

    unconstrained = constrained
    do i = 1, size(self%comps)
      k = self%k0(i)
      associate (c => self%comps(i)%c)
        if (c%k > 0) unconstrained(k + 1:k + c%k) = c%untransform_params(constrained(k + 1:k + c%k))
      end associate
    end do
  end function sm_untransform

  function sm_names(self) result(names)
    class(structural_model_t), intent(in) :: self
    character(len=32), allocatable :: names(:)
    integer :: i

    character(len=32) :: buf
    integer :: j, l

    allocate (names(0))
    do i = 1, size(self%comps)
      names = [names, self%comps(i)%c%param_names()]
    end do
    do l = 1, size(self%loading, 2)
      do j = 1, size(self%loading, 1)
        if (.not. self%loading_free(j, l)) cycle
        write (buf, '("loading.", i0, ".", i0)') j, l
        names = [names, buf]
      end do
    end do
  end function sm_names

  ! Defaults: every parameter is a variance (sigma^2 = exp(2 psi), DK 7.3.2),
  ! started at a share of the variance of the differenced data.

  function comp_transform(self, unconstrained) result(constrained)
    class(component_t), intent(in) :: self
    real(dp), intent(in) :: unconstrained(:)
    real(dp), allocatable :: constrained(:)

    constrained = constrain_positive(unconstrained)
  end function comp_transform

  function comp_untransform(self, constrained) result(unconstrained)
    class(component_t), intent(in) :: self
    real(dp), intent(in) :: constrained(:)
    real(dp), allocatable :: unconstrained(:)

    unconstrained = unconstrain_positive(constrained)
  end function comp_untransform

  function comp_start(self, y) result(params)
    class(component_t), intent(in) :: self
    real(dp), intent(in) :: y(:, :)
    real(dp), allocatable :: params(:)

    allocate (params(self%k), source=diff_variance(y) / 4.0_dp)
  end function comp_start

  function comp_names(self) result(names)
    class(component_t), intent(in) :: self
    character(len=32), allocatable :: names(:)
    integer :: i

    allocate (names(self%k))
    do i = 1, self%k
      write (names(i), '("param", i0)') i
    end do
  end function comp_names

  !> Whether the component acts on the observations (not the signals).
  logical function comp_observation_level(self)
    class(component_t), intent(in) :: self

    comp_observation_level = self%at_observations
  end function comp_observation_level

  pure function identity(n) result(A)
    integer, intent(in) :: n
    real(dp) :: A(n, n)
    integer :: i

    A = 0.0_dp
    do i = 1, n
      A(i, i) = 1.0_dp
    end do
  end function identity

  !> Default initialization: the whole block diffuse.
  function comp_init_blocks(self) result(blocks)
    class(component_t), intent(in) :: self
    integer, allocatable :: blocks(:, :)

    blocks = reshape([1, self%m, INIT_DIFFUSE], [1, 3])
  end function comp_init_blocks

  !> Sample variance of the first differences of the first series (NaNs skipped).
  real(dp) function diff_variance(y) result(v)
    real(dp), intent(in) :: y(:, :)
    real(dp), allocatable :: d(:)
    integer :: n

    n = size(y, 2)
    if (n < 3) then
      v = 1.0_dp
      return
    end if
    d = y(1, 2:) - y(1, :n - 1)
    d = pack(d, d == d)
    if (size(d) < 2) then
      v = 1.0_dp
      return
    end if
    v = sum((d - sum(d) / size(d))**2) / (size(d) - 1)
    if (.not. v > 0.0_dp) v = 1.0_dp
  end function diff_variance

  !> DK's bounded transform (7.3.2), shifted to an interval:
  !> chi = mid + half psi / sqrt(1 + psi^2) maps the real line onto (lo, hi).
  elemental real(dp) function constrain_interval(psi, lo, hi) result(chi)
    real(dp), intent(in) :: psi, lo, hi

    chi = 0.5_dp * (lo + hi) + 0.5_dp * (hi - lo) * psi / sqrt(1.0_dp + psi**2)
  end function constrain_interval

  !> Inverse of `constrain_interval`.
  elemental real(dp) function unconstrain_interval(chi, lo, hi) result(psi)
    real(dp), intent(in) :: chi, lo, hi
    real(dp) :: u

    u = (2.0_dp * chi - lo - hi) / (hi - lo)
    psi = u / sqrt(1.0_dp - u**2)
  end function unconstrain_interval
end module statespace_components
