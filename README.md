# spec-state-trans: Spectral state transitions project

## Dependencies - clone next to this repository:

* [proektlab/Ethan_misc](https://github.com/proektlab/Ethan_misc)
  * Add `Ethan_misc` to the MATLAB path.
  * In `Ethan_misc/prepSR.m`, change `box_dir` to be a local path to the contents of the "Blackwood" data folder on Box.

* [proektlab/SpectralAnalysis](https://github.com/proektlab/SpectralAnalysis)

* [proektlab/NonNegativeMatrixFactorization](https://github.com/proektlab/NonNegativeMatrixFactorization)

* [probml/pmtk4](https://github.com/probml/pmtk3)

## Code organization

* The root folder contains common functions for doing each step of the analysis, as well as scripts that use data from all recordings.

* Each individual dataset folder contains `multitaper_layers.m`, which performs a multitaper time-frequency analysis.

* `nnmf_all_byday.m` preprocesses multitaper data and reduces dimensionality using non-negative matrix factorization.

