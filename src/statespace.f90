!> Umbrella module: `use statespace` gives access to the whole public API.
module statespace
  use statespace_kinds
  use statespace_linalg, only: eye
  use statespace_rep
  use statespace_filter
  use statespace_smoother
  use statespace_forecast
  use statespace_simsmooth
  use statespace_smoothing
  use statespace_augmented
  use statespace_sqrt
  use statespace_collapse
  use statespace_restrict
  use statespace_model
  use statespace_mle
  use statespace_score
  use statespace_em
  use statespace_special
  use statespace_diagnostics
  use statespace_dense
  use statespace_components
  use statespace_structural
  use statespace_arima
  implicit none
  public
end module statespace
