# spec-state-trans: Spectral state transitions project

**Goal:** Identify times in each LFP channel when major state transitions occur (based on time-frequency spectral analysis) and test whether temporally matched transitions between channels tend to happen in a particular order.*

**Hypothesis:** State transitions that occur in both analyzed regions (V1 and motor cortex) tend to consistently happen in one region before the other.

## Code organization

* Depends on utilities in ethanbb/proektlab-matlab-tools. The tools repository should be cloned next to this one (i.e. they should be in sibling directories).
  * In `proektlab-matlab-tools/prepSR.m`, change `box_dir` to be a local path to the contents of the "Blackwood" data folder on Box.

* The `multitaper.m` script in each <date>/<time> dataset directory performs the initial time-frequency analysis on the corresponding pre-processed data on Box.

* The script `2020-02-06/16-01-00/change_compare_chans.m` contains the rest of the analysis for the 2/6 4:01PM dataset.