# References

This page lists the works cited in the source code and docs, the datasets, and the
software this project uses. Where a project asks users to cite it in a particular way,
that is the citation given here.

## Methods

The library follows the first of these, and section numbers in the source ("DK 4.3.4")
refer to it.

- Durbin, J. and Koopman, S. J. (2012). *Time Series Analysis by State Space Methods*,
  2nd ed. Oxford University Press. ISBN 978-0-19-964117-8.
  [doi:10.1093/acprof:oso/9780199641178.001.0001](https://doi.org/10.1093/acprof:oso/9780199641178.001.0001)
- Anderson, B. D. O. and Moore, J. B. (1979). *Optimal Filtering*. Prentice-Hall.
  (Steady state of the Kalman filter.)
- Chu, E. K.-W., Fan, H.-Y., Lin, W.-W. and Wang, C.-S. (2004). Structure-preserving
  algorithms for periodic discrete-time algebraic Riccati equations. *International
  Journal of Control* 77(8), 767–788.
  [doi:10.1080/00207170410001714988](https://doi.org/10.1080/00207170410001714988)
  (The doubling algorithm in `steady_state`.)
- de Jong, P. and Penzer, J. (1998). Diagnosing shocks in time series. *Journal of the
  American Statistical Association* 93(442), 796–806.
  [doi:10.1080/01621459.1998.10473731](https://doi.org/10.1080/01621459.1998.10473731)
- de Jong, P. and Shephard, N. (1995). The simulation smoother for time series models.
  *Biometrika* 82(2), 339–350.
  [doi:10.1093/biomet/82.2.339](https://doi.org/10.1093/biomet/82.2.339)
- Durbin, J. and Koopman, S. J. (2002). A simple and efficient simulation smoother for
  state space time series analysis. *Biometrika* 89(3), 603–615.
  [doi:10.1093/biomet/89.3.603](https://doi.org/10.1093/biomet/89.3.603)
- Francke, M. K., Koopman, S. J. and de Vos, A. F. (2010). Likelihood functions for state
  space models with diffuse initial conditions. *Journal of Time Series Analysis* 31(6),
  407–414. [doi:10.1111/j.1467-9892.2010.00673.x](https://doi.org/10.1111/j.1467-9892.2010.00673.x)
- Harvey, A. C. (1989). *Forecasting, Structural Time Series Models and the Kalman
  Filter*. Cambridge University Press.
  [doi:10.1017/CBO9781107049994](https://doi.org/10.1017/CBO9781107049994)
- Jungbacker, B. and Koopman, S. J. (2015). Likelihood-based dynamic factor analysis for
  measurement and forecasting. *The Econometrics Journal* 18(2), C1–C21.
  [doi:10.1111/ectj.12029](https://doi.org/10.1111/ectj.12029)
  (Collapsing the observation vector, DK 6.5. DK cite the 2008 working paper version.)
- Koopman, S. J. and Shephard, N. (1992). Exact score for time series models in state
  space form. *Biometrika* 79(4), 823–826.
  [doi:10.1093/biomet/79.4.823](https://doi.org/10.1093/biomet/79.4.823)
- Monahan, J. F. (1984). A note on enforcing stationarity in autoregressive-moving
  average models. *Biometrika* 71(2), 403–404.
  [doi:10.1093/biomet/71.2.403](https://doi.org/10.1093/biomet/71.2.403)
- Nelson, C. R. and Siegel, A. F. (1987). Parsimonious modeling of yield curves.
  *Journal of Business* 60(4), 473–489. [doi:10.1086/296409](https://doi.org/10.1086/296409)
- Thompson, I. J. and Barnett, A. R. (1986). Coulomb and Bessel functions of complex
  arguments and order. *Journal of Computational Physics* 64(2), 490–509.
  [doi:10.1016/0021-9991(86)90046-X](https://doi.org/10.1016/0021-9991(86)90046-X)
  (The modified Lentz method for the continued fractions in `statespace_special`.)
- Wahba, G. (1978). Improper priors, spline smoothing and the problem of guarding
  against model errors in regression. *Journal of the Royal Statistical Society,
  Series B* 40(3), 364–372.
  [doi:10.1111/j.2517-6161.1978.tb01050.x](https://doi.org/10.1111/j.2517-6161.1978.tb01050.x)

## Data

See `data/README.md` for the files themselves.

- **Nile.** Cobb, G. W. (1978). The problem of the Nile: conditional solution to a
  changepoint problem. *Biometrika* 65(2), 243–251.
  [doi:10.1093/biomet/65.2.243](https://doi.org/10.1093/biomet/65.2.243)
- **Seat belts (DK 8.2–8.3).** Harvey, A. C. and Durbin, J. (1986). The effects of seat
  belt legislation on British road casualties: a case study in structural time series
  modelling. *Journal of the Royal Statistical Society, Series A* 149(3), 187–227.
  [doi:10.2307/2981553](https://doi.org/10.2307/2981553)
  Distributed as `Seatbelts` in R's `datasets` package.
- **Internet users (DK 8.4).** Makridakis, S., Wheelwright, S. C. and Hyndman, R. J.
  (1998). *Forecasting: Methods and Applications*, 3rd ed. Wiley. Distributed as
  `WWWusage` in R's `datasets` package.
- **Motorcycle acceleration (DK 8.5).** Silverman, B. W. (1985). Some aspects of the
  spline smoothing approach to non-parametric regression curve fitting. *Journal of the
  Royal Statistical Society, Series B* 47(1), 1–52.
  [doi:10.1111/j.2517-6161.1985.tb01327.x](https://doi.org/10.1111/j.2517-6161.1985.tb01327.x)
  Distributed as `mcycle` in R's `MASS` package:
  Venables, W. N. and Ripley, B. D. (2002). *Modern Applied Statistics with S*, 4th ed.
  Springer. [doi:10.1007/978-0-387-21706-2](https://doi.org/10.1007/978-0-387-21706-2)
- **US Treasury yields (DK 8.6 substitute).** FRED asks for this citation for each
  series: Board of Governors of the Federal Reserve System (US), Market Yield on U.S.
  Treasury Securities at 3-Month, 6-Month, 1-Year, 2-Year, 3-Year, 5-Year, 7-Year and
  10-Year Constant Maturity, Quoted on an Investment Basis [GS3M, GS6M, GS1, GS2, GS3,
  GS5, GS7, GS10], retrieved from FRED, Federal Reserve Bank of St. Louis,
  <https://fred.stlouisfed.org/>. The data DK use are from:
  Diebold, F. X. and Li, C. (2006). Forecasting the term structure of government bond
  yields. *Journal of Econometrics* 130(2), 337–364.
  [doi:10.1016/j.jeconom.2005.03.005](https://doi.org/10.1016/j.jeconom.2005.03.005)
- **R and Rdatasets.** R Core Team. *R: A Language and Environment for Statistical
  Computing*. R Foundation for Statistical Computing, Vienna. <https://www.R-project.org/>.
  Arel-Bundock, V. *Rdatasets: A collection of datasets originally distributed in R
  packages*. <https://github.com/vincentarelbundock/Rdatasets>.

## Software

### Build dependencies

- **L-BFGS-B**, via [jacobwilliams/lbfgsb](https://github.com/jacobwilliams/lbfgsb):
  - Byrd, R. H., Lu, P., Nocedal, J. and Zhu, C. (1995). A limited memory algorithm for
    bound constrained optimization. *SIAM Journal on Scientific Computing* 16(5),
    1190–1208. [doi:10.1137/0916069](https://doi.org/10.1137/0916069)
  - Zhu, C., Byrd, R. H., Lu, P. and Nocedal, J. (1997). Algorithm 778: L-BFGS-B:
    Fortran subroutines for large-scale bound-constrained optimization. *ACM
    Transactions on Mathematical Software* 23(4), 550–560.
    [doi:10.1145/279232.279236](https://doi.org/10.1145/279232.279236)
  - Morales, J. L. and Nocedal, J. (2011). Remark on "Algorithm 778: L-BFGS-B: Fortran
    subroutines for large-scale bound constrained optimization". *ACM Transactions on
    Mathematical Software* 38(1), Article 7.
    [doi:10.1145/2049662.2049669](https://doi.org/10.1145/2049662.2049669)
- **LAPACK**, which asks that proper credit be given to its authors:
  Anderson, E. et al. (1999). *LAPACK Users' Guide*, 3rd ed. SIAM.
  [doi:10.1137/1.9780898719604](https://doi.org/10.1137/1.9780898719604)
- **BLAS**: any implementation, e.g. [OpenBLAS](https://www.openmathlib.org/OpenBLAS/)
  or the reference BLAS.
- **test-drive** (tests only): <https://github.com/fortran-lang/test-drive>.

### Reference implementation and fixture generation

The test fixtures are generated by `test/fixtures/make_fixtures.py`, using:

- **statsmodels**, used as the reference implementation. Its requested citation:
  Seabold, S. and Perktold, J. (2010). statsmodels: Econometric and statistical
  modeling with Python. *Proceedings of the 9th Python in Science Conference*.
  [doi:10.25080/Majora-92bf1922-011](https://doi.org/10.25080/Majora-92bf1922-011)
  Software archive: [doi:10.5281/zenodo.593847](https://doi.org/10.5281/zenodo.593847)
- **SciPy** (the smoothing-spline reference):
  Virtanen, P. et al. (2020). SciPy 1.0: fundamental algorithms for scientific computing
  in Python. *Nature Methods* 17, 261–272.
  [doi:10.1038/s41592-019-0686-2](https://doi.org/10.1038/s41592-019-0686-2)
- **NumPy**: Harris, C. R. et al. (2020). Array programming with NumPy. *Nature* 585,
  357–362. [doi:10.1038/s41586-020-2649-2](https://doi.org/10.1038/s41586-020-2649-2)
- **pandas**: The pandas development team. *pandas-dev/pandas: Pandas*. Zenodo.
  [doi:10.5281/zenodo.3509134](https://doi.org/10.5281/zenodo.3509134);
  McKinney, W. (2010). Data structures for statistical computing in Python.
  *Proceedings of the 9th Python in Science Conference*, 56–61.
  [doi:10.25080/Majora-92bf1922-00a](https://doi.org/10.25080/Majora-92bf1922-00a)
