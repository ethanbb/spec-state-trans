# spec-state-trans: Spectral state transitions project

**Goal:** Identify times in each LFP channel when major state transitions occur (based on time-frequency spectral analysis) and test whether temporally matched transitions between channels tend to happen in a particular order.*

**Hypothesis:** State transitions that occur in both analyzed regions (V1 and motor cortex) tend to consistently happen in one region before the other.

## Dependencies - clone next to this repository:

* [proektlab/Ethan_misc](https://github.com/proektlab/Ethan_misc)
  * Add `Ethan_misc` to the MATLAB path.
  * In `Ethan_misc/prepSR.m`, change `box_dir` to be a local path to the contents of the "Blackwood" data folder on Box.

* [proektlab/SpectralAnalysis](https://github.com/proektlab/SpectralAnalysis)

## Code organization

* The root folder contains common functions for doing each step of the analysis, as well as the final analysis script that uses data from all recordings.

* Each individual dataset folder contains two scripts for analyzing that dataset (and saving results in Box):
  
  * `multitaper.m` performs a multitaper time-frequency analysis.

  * `change_compare_chans.m` does the rest of the analysis - PCA, extracting change velocity, finding peaks and troughs and matching them.

* Finally, once each result file has matched peaks and troughs from `change_compare_chans`, `change_extrema_comp.m` combines them and makes a histogram.
