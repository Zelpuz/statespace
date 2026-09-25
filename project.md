# State space time series models: project notes

## Requirements

- Full capability of modelling linear Gaussian state space time series models _a la_ Part I of Durbin and Koopman's _Time Series Analysis by State Space Methods, Second Edition_.
- Support for MLE via L-BFGS-B for all variants of likelihood function outlined in Durbin and Koopman (ch. 7).
- User-defined state and transition matrices, alongside other aspects of model structure, with the program handling initialization, filtering, smoothing, MLE, etc.
- The user defines the model structure by providing the matrices. Model structure not supported implicitly (i.e. atypical matrices, custom cycle components, etc.) will require the user to provide short FORTRAN module(s) describing the necessary component(s). This will require recompilation, naturally.
- Common or popular structural models (local level, trend, seasonal, cycle, regression components, AR terms, etc.) should be provided as examples.
