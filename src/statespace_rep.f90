!> Linear Gaussian state space representation (DK 3.1):
!>
!>     y_t       = d_t + Z_t alpha_t + eps_t,        eps_t ~ N(0, H_t)
!>     alpha_t+1 = c_t + T_t alpha_t + R_t eta_t,    eta_t ~ N(0, Q_t)
!>
!> with y_t (p), alpha_t (m), eta_t (r), t = 1..n.
!>
!> Arrays are column-major with time as the last dimension. That dimension is
!> either 1 (time-invariant) or n (time-varying), and `tidx` maps t to the
!> right slice. The components are public: build a representation with
!> `ssm_rep(y, m, r)`, then assign the system matrices, e.g.
!> `rep%T(:,:,1) = ...` or `rep%Z = Z_time_varying` (reallocates).
!>
!> Missing observations are NaN in y and may be any subset of y_t (DK 4.10).
!>
!> Initialization (DK 5.1): alpha_1 = a + A delta + R_0 eta_0 with delta
!> diffuse, i.e. alpha_1 ~ N(a, P_star + kappa P_inf) with kappa -> infinity,
!> P_star = R_0 Q_0 R_0' and P_inf = A A'. It is specified either for the
!> whole state vector (`initialize_known`, `_approximate_diffuse`,
!> `_stationary`, `_diffuse`, `_general`) or block by block with
!> `initialize_block`, e.g. a diffuse trend alongside a stationary AR block.
module statespace_rep
  use statespace_kinds, only: dp, SS_OK, SS_ERR_DIM, SS_ERR_INIT
  use statespace_linalg, only: eye, gemm, symmetrize, solve, solve_lyapunov, ldl_psd
  implicit none
  private

  public :: ssm_rep_t, ssm_rep, tidx

  !> Initialization kinds (DK ch. 5).
  integer, parameter, public :: INIT_NONE = 0
  integer, parameter, public :: INIT_KNOWN = 1           !< N(a1, P1)
  integer, parameter, public :: INIT_APPROX_DIFFUSE = 2  !< N(a1, kappa I)
  integer, parameter, public :: INIT_STATIONARY = 3      !< unconditional distribution
  integer, parameter, public :: INIT_DIFFUSE = 4         !< exact diffuse, P_inf = I
  integer, parameter, public :: INIT_GENERAL = 5         !< given a1, P_star, P_inf

  !> Filter methods. The conventional filter processes y_t as a vector (DK
  !> 4.3); periods with a diffuse state always use the univariate treatment
  !> (DK 5.2 with 6.4). FILTER_UNIVARIATE uses the univariate treatment
  !> for every period.
  integer, parameter, public :: FILTER_CONVENTIONAL = 0
  integer, parameter, public :: FILTER_UNIVARIATE = 1

  !> Treatment of the diffuse periods when filter_method is conventional:
  !> the univariate exact initial filter (DK 5.2 with 6.4), or the
  !> multivariate one (DK 5.2) where F_inf is zero or nonsingular, falling
  !> back to the univariate step for periods where it is singular but nonzero.
  integer, parameter, public :: DIFFUSE_UNIVARIATE = 0
  integer, parameter, public :: DIFFUSE_MULTIVARIATE = 1

  type :: ssm_rep_t
    integer :: k_endog = 0, k_states = 0, k_posdef = 0, nobs = 0
    real(dp), allocatable :: y(:, :)        !< (p, n)
    real(dp), allocatable :: Z(:, :, :)     !< (p, m, 1|n)
    real(dp), allocatable :: H(:, :, :)     !< (p, p, 1|n)
    real(dp), allocatable :: T(:, :, :)     !< (m, m, 1|n)
    real(dp), allocatable :: R(:, :, :)     !< (m, r, 1|n)
    real(dp), allocatable :: Q(:, :, :)     !< (r, r, 1|n)
    real(dp), allocatable :: c(:, :)        !< (m, 1|n)
    real(dp), allocatable :: d(:, :)        !< (p, 1|n)

    !> Initialization blocks: states blk_first(k)..blk_last(k) have kind
    !> blk_kind(k). User-supplied means and variances live in a1, Pstar1,
    !> Pinf1 at the block's positions.
    integer, allocatable :: blk_first(:), blk_last(:), blk_kind(:)
    real(dp), allocatable :: a1(:)          !< (m)
    real(dp), allocatable :: Pstar1(:, :)   !< (m, m)
    real(dp), allocatable :: Pinf1(:, :)    !< (m, m)

    integer :: filter_method = FILTER_CONVENTIONAL
    integer :: diffuse_method = DIFFUSE_UNIVARIATE
    !> Threshold on F_inf (and on ||P_inf||_F^2) for diffuse updates, as in
    !> statsmodels.
    real(dp) :: tol_diffuse = 1.0e-10_dp

    !> Number of leading observations excluded from the total log likelihood.
    integer :: loglikelihood_burn = 0

    !> Use the marginal likelihood of Francke, Koopman and de Vos (2010)
    !> instead of the diffuse likelihood (DK 7.2.6). They differ by a term
    !> that depends on Z, T and the diffuse directions A only, so estimates
    !> differ only when parameters enter Z or T (e.g. factor loadings).
    logical :: marginal_likelihood = .false.
  contains
    procedure :: initialize_known
    procedure :: initialize_approximate_diffuse
    procedure :: initialize_stationary
    procedure :: initialize_diffuse
    procedure :: initialize_general
    procedure :: initialize_block
    procedure :: initial_state
    procedure :: k_diffuse
    procedure :: scale_by
    procedure :: validate
  end type ssm_rep_t

contains

  !> Create a time-invariant representation for data y(p, n) with m states
  !> and r state disturbances. All system matrices start at zero except R,
  !> which starts as the leading m x r identity.
  function ssm_rep(y, m, r) result(rep)
    real(dp), intent(in) :: y(:, :)
    integer, intent(in) :: m, r
    type(ssm_rep_t) :: rep

    rep%k_endog = size(y, 1)
    rep%nobs = size(y, 2)
    rep%k_states = m
    rep%k_posdef = r
    rep%y = y
    allocate (rep%Z(rep%k_endog, m, 1), rep%H(rep%k_endog, rep%k_endog, 1), rep%T(m, m, 1), &
              rep%R(m, r, 1), rep%Q(r, r, 1), rep%c(m, 1), rep%d(rep%k_endog, 1), source=0.0_dp)
    rep%R(:, :, 1) = eye(m, r)
  end function ssm_rep

  !> Slice index for time t into an array whose time dimension has length nt.
  pure integer function tidx(nt, t)
    integer, intent(in) :: nt, t

    tidx = min(t, nt)
  end function tidx

  !> alpha_1 ~ N(a1, P1) with known a1 and P1.
  subroutine initialize_known(self, a1, P1)
    class(ssm_rep_t), intent(inout) :: self
    real(dp), intent(in) :: a1(:), P1(:, :)

    call reset_init(self)
    call self%initialize_block(1, self%k_states, INIT_KNOWN, a1, P1)
  end subroutine initialize_known

  !> alpha_1 ~ N(a1, kappa I), a1 = 0 unless given (DK 5.1).
  subroutine initialize_approximate_diffuse(self, kappa, a1)
    class(ssm_rep_t), intent(inout) :: self
    real(dp), intent(in), optional :: kappa
    real(dp), intent(in), optional :: a1(:)

    call reset_init(self)
    call self%initialize_block(1, self%k_states, INIT_APPROX_DIFFUSE, a1=a1, kappa=kappa)
  end subroutine initialize_approximate_diffuse

  !> alpha_1 from the unconditional distribution of a stationary state:
  !> a1 = (I - T)^-1 c and P1 solves P1 = T P1 T' + R Q R', using the system
  !> matrices at t = 1. Computed at the start of each filter run, so it
  !> follows parameter updates.
  subroutine initialize_stationary(self)
    class(ssm_rep_t), intent(inout) :: self

    call reset_init(self)
    call self%initialize_block(1, self%k_states, INIT_STATIONARY)
  end subroutine initialize_stationary

  !> Exact diffuse alpha_1: P_inf = I, P_star = 0, mean a1 (default 0) (DK 5.1).
  !> The diffuse log likelihood (DK 7.2.2) already accounts for the diffuse
  !> states, so `loglikelihood_burn` should normally be 0.
  subroutine initialize_diffuse(self, a1)
    class(ssm_rep_t), intent(inout) :: self
    real(dp), intent(in), optional :: a1(:)

    call reset_init(self)
    call self%initialize_block(1, self%k_states, INIT_DIFFUSE, a1=a1)
  end subroutine initialize_diffuse

  !> DK 5.1 in full: alpha_1 ~ N(a1, Pstar + kappa Pinf), kappa -> infinity.
  subroutine initialize_general(self, a1, Pstar, Pinf)
    class(ssm_rep_t), intent(inout) :: self
    real(dp), intent(in) :: a1(:), Pstar(:, :), Pinf(:, :)

    call reset_init(self)
    call add_block(self, 1, self%k_states, INIT_GENERAL)
    self%a1 = a1
    self%Pstar1 = Pstar
    self%Pinf1 = Pinf
  end subroutine initialize_general

  !> Initialize states first..last with `kind`:
  !>   INIT_KNOWN           needs a1 and P1
  !>   INIT_APPROX_DIFFUSE  optional a1 (default 0) and kappa (default 1e6)
  !>   INIT_STATIONARY      from T, R, Q, c restricted to the block
  !>   INIT_DIFFUSE         optional a1 (default 0)
  !> Blocks must not overlap and together must cover every state. A stationary
  !> block must not depend on states outside it through T.
  subroutine initialize_block(self, first, last, kind, a1, P1, kappa)
    class(ssm_rep_t), intent(inout) :: self
    integer, intent(in) :: first, last, kind
    real(dp), intent(in), optional :: a1(:), P1(:, :), kappa

    if (.not. allocated(self%blk_kind)) call reset_init(self)
    if (any(self%blk_kind == INIT_GENERAL)) call reset_init(self)
    call add_block(self, first, last, kind)
    if (present(a1)) self%a1(first:last) = a1
    select case (kind)
    case (INIT_KNOWN)
      if (present(P1)) self%Pstar1(first:last, first:last) = P1
    case (INIT_APPROX_DIFFUSE)
      self%Pstar1(first:last, first:last) = eye(last - first + 1)
      if (present(kappa)) then
        self%Pstar1(first:last, first:last) = kappa * self%Pstar1(first:last, first:last)
      else
        self%Pstar1(first:last, first:last) = 1.0e6_dp * self%Pstar1(first:last, first:last)
      end if
    end select
  end subroutine initialize_block

  subroutine reset_init(self)
    class(ssm_rep_t), intent(inout) :: self
    integer :: m

    m = self%k_states
    self%blk_first = [integer ::]
    self%blk_last = [integer ::]
    self%blk_kind = [integer ::]
    self%a1 = spread(0.0_dp, 1, m)
    self%Pstar1 = 0.0_dp * eye(m)
    self%Pinf1 = self%Pstar1
  end subroutine reset_init

  subroutine add_block(self, first, last, kind)
    class(ssm_rep_t), intent(inout) :: self
    integer, intent(in) :: first, last, kind

    self%blk_first = [self%blk_first, first]
    self%blk_last = [self%blk_last, last]
    self%blk_kind = [self%blk_kind, kind]
  end subroutine add_block

  !> Mean a1 and the variance parts Pstar, Pinf of alpha_1.
  subroutine initial_state(self, a1, Pstar, Pinf, info)
    class(ssm_rep_t), intent(in) :: self
    real(dp), intent(out) :: a1(:), Pstar(:, :), Pinf(:, :)
    integer, intent(out) :: info
    integer :: k, i, j

    info = SS_OK
    a1 = self%a1
    Pstar = self%Pstar1
    Pinf = self%Pinf1
    do k = 1, size(self%blk_kind)
      i = self%blk_first(k)
      j = self%blk_last(k)
      select case (self%blk_kind(k))
      case (INIT_DIFFUSE)
        Pinf(i:j, i:j) = eye(j - i + 1)
      case (INIT_STATIONARY)
        call stationary_block(self, i, j, a1(i:j), Pstar(i:j, i:j), info)
        if (info /= SS_OK) return
      end select
    end do
  end subroutine initial_state

  !> Unconditional mean and variance of states i..j.
  subroutine stationary_block(self, i, j, a, P, info)
    class(ssm_rep_t), intent(in) :: self
    integer, intent(in) :: i, j
    real(dp), intent(out) :: a(:), P(:, :)
    integer, intent(out) :: info
    real(dp), allocatable :: RQ(:, :), V(:, :), IminusT(:, :), rhs(:, :)
    integer :: nb

    nb = j - i + 1
    allocate (RQ(nb, self%k_posdef), V(nb, nb))
    call gemm('N', 'N', 1.0_dp, self%R(i:j, :, 1), self%Q(:, :, 1), 0.0_dp, RQ)
    call gemm('N', 'T', 1.0_dp, RQ, self%R(i:j, :, 1), 0.0_dp, V)
    call symmetrize(V)
    call solve_lyapunov(self%T(i:j, i:j, 1), V, P, info)
    if (info /= SS_OK) return
    IminusT = eye(nb) - self%T(i:j, i:j, 1)
    rhs = reshape(self%c(i:j, 1), [nb, 1])
    call solve(IminusT, rhs, info)
    if (info /= SS_OK) return
    a = rhs(:, 1)
  end subroutine stationary_block

  !> Number of diffuse elements in alpha_1, the rank of P_inf.
  integer function k_diffuse(self)
    class(ssm_rep_t), intent(in) :: self
    real(dp), allocatable :: L(:, :), D(:)
    integer :: k

    k_diffuse = 0
    if (.not. allocated(self%blk_kind)) return
    do k = 1, size(self%blk_kind)
      select case (self%blk_kind(k))
      case (INIT_DIFFUSE)
        k_diffuse = k_diffuse + self%blk_last(k) - self%blk_first(k) + 1
      case (INIT_GENERAL)
        allocate (L(self%k_states, self%k_states), D(self%k_states))
        call ldl_psd(self%Pinf1, L, D)
        k_diffuse = count(D > 0.0_dp)
      end select
    end do
  end function k_diffuse

  !> Multiply H, Q and the P_star part of the initialization by s, the scale
  !> of a model whose scale was concentrated out (`loglike_concentrated`).
  !> Stationary blocks scale through Q; P_inf is not affected.
  subroutine scale_by(self, s)
    class(ssm_rep_t), intent(inout) :: self
    real(dp), intent(in) :: s

    self%H = s * self%H
    self%Q = s * self%Q
    if (allocated(self%Pstar1)) self%Pstar1 = s * self%Pstar1
  end subroutine scale_by

  !> Check that dimensions are consistent and that the model is initialized.
  subroutine validate(self, info)
    class(ssm_rep_t), intent(in) :: self
    integer, intent(out) :: info
    integer :: p, m, r, n, k
    integer, allocatable :: covered(:)

    p = self%k_endog; m = self%k_states; r = self%k_posdef; n = self%nobs
    info = SS_ERR_DIM
    if (.not. allocated(self%y)) return
    if (any(shape(self%y) /= [p, n])) return
    if (.not. ok3(self%Z, p, m)) return
    if (.not. ok3(self%H, p, p)) return
    if (.not. ok3(self%T, m, m)) return
    if (.not. ok3(self%R, m, r)) return
    if (.not. ok3(self%Q, r, r)) return
    if (.not. ok2(self%c, m)) return
    if (.not. ok2(self%d, p)) return

    ! Every state is initialized exactly once.
    info = SS_ERR_INIT
    if (.not. allocated(self%blk_kind)) return
    if (size(self%blk_kind) == 0) return
    if (size(self%a1) /= m .or. any(shape(self%Pstar1) /= [m, m]) .or. &
        any(shape(self%Pinf1) /= [m, m])) return
    allocate (covered(m), source=0)
    do k = 1, size(self%blk_kind)
      if (self%blk_first(k) < 1 .or. self%blk_last(k) > m .or. &
          self%blk_first(k) > self%blk_last(k)) return
      if (self%blk_kind(k) < INIT_KNOWN .or. self%blk_kind(k) > INIT_GENERAL) return
      covered(self%blk_first(k):self%blk_last(k)) = covered(self%blk_first(k):self%blk_last(k)) + 1
    end do
    if (any(covered /= 1)) return

    info = SS_OK
  contains
    logical function ok3(A, d1, d2)
      real(dp), allocatable, intent(in) :: A(:, :, :)
      integer, intent(in) :: d1, d2

      ok3 = .false.
      if (.not. allocated(A)) return
      ok3 = size(A, 1) == d1 .and. size(A, 2) == d2 .and. &
            (size(A, 3) == 1 .or. size(A, 3) == n)
    end function ok3

    logical function ok2(A, d1)
      real(dp), allocatable, intent(in) :: A(:, :)
      integer, intent(in) :: d1

      ok2 = .false.
      if (.not. allocated(A)) return
      ok2 = size(A, 1) == d1 .and. (size(A, 2) == 1 .or. size(A, 2) == n)
    end function ok2
  end subroutine validate
end module statespace_rep
