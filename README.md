# spec-state-trans: Spectral state transitions project

## Dependencies - clone next to this repository:

* [proektlab/Ethan_misc](https://github.com/proektlab/Ethan_misc)
  * Add `Ethan_misc` to the MATLAB path.
  * In `Ethan_misc/prepSR.m`, change `box_dir` to be a local path to the contents of the "Blackwood" data folder on Box.

* [proektlab/SpectralAnalysis](https://github.com/proektlab/SpectralAnalysis)

* [proektlab/NonNegativeMatrixFactorization](https://github.com/proektlab/NonNegativeMatrixFactorization)

## Current analysis pipeline

TODO: clean up m-file organization. For now, many folders contain an "old" folder which has obsolete scripts that were at that location in the file tree.

1. Compute CSDs for each region on each recording day using `plot_csd.m`

    * **Output:** `<date>/csd_{V1,MC}.mat`

2. Use CSDs to assign channel numbers to layers (2/3, 4, 5) in the CSD matfiles using `pick_channels_from_csd.m`
  
    * **Output:** `<date>/csd_{V1,MC}.mat`

3. Do multitaper analyses of each recording's selected channels, with artifact filtering etc., using individual `multitaper_layers.m` files in each date/time folder

    * **Output:** `<date>/<time>/mt_res_layers.mat`

4. Do NMF and extract classes and transitions for each channel (based on dominant component per time) using `concat_and_nmf.m`

    * **Output:** `<date>/nmf_res.mat`

5. Generate "null model" score data for each channel (shuffled using Markov chain - negative control for connectivity) using `gen_null_model_data.m`

    * **Output:** `<date>/nmf_res.mat`

6. Compute KL divergence between aligned NMF scores of each directed pair of channels (and plot results) with `kl_divergence_analysis.m`

    * **Output:** `<date>/nmf_res.mat`
