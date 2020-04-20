# spec-state-trans: Spectral state transitions project

**Goal:** Identify times in each LFP channel when major state transitions occur (based on time-frequency spectral analysis) and test whether temporally matched transitions between channels tend to happen in a particular order.*

**Hypothesis:** State transitions that occur in both analyzed regions (V1 and motor cortex) tend to consistently happen in one region before the other.

## Dependencies - clone next to this repository:

* [proektlab/Ethan_misc](https://github.com/proektlab/Ethan_misc)
  * Add `Ethan_misc` to the MATLAB path.
  * In `Ethan_misc/prepSR.m`, change `box_dir` to be a local path to the contents of the "Blackwood" data folder on Box.

* [proektlab/SpectralAnalysis](https://github.com/proektlab/SpectralAnalysis)

* [proektlab/NonNegativeMatrixFactorization](https://github.com/proektlab/NonNegativeMatrixFactorization)

## Code organization

* The root folder contains common functions for doing each step of the analysis, as well as the final analysis script that uses data from all recordings.

* Each individual dataset folder contains two scripts for analyzing that dataset (and saving results in Box):
  
  * `multitaper.m` performs a multitaper time-frequency analysis.

  * `nnmf_clustering.m`, in each day's folder, does the rest of the analysis - preprocessing, non-negative matrix factorization, and clustering using a GMM.

* Finally, once each result file has a matrix of GMM state transition times from `nnmf_clustering.m`, `gmm_transition_comp.m` combines them and makes a histogram.
