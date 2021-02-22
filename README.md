# spec-state-trans: Spectral state transitions project

## Dependencies:

* [proektlab/Ethan_misc](https://github.com/proektlab/Ethan_misc) - clone next to this repository
  * Add `Ethan_misc` to the MATLAB path.
  * Edit `Ethan_misc/prepSR.m` to point to local folders, as described in the comments.

* [proektlab/SpectralAnalysis](https://github.com/proektlab/SpectralAnalysis)

* [proektlab/NonNegativeMatrixFactorization](https://github.com/proektlab/NonNegativeMatrixFactorization)

## Current analysis pipeline

See `mt_new_2020_recs.m` followed by `full_analysis_byday.m` for an example.

1. Compute CSDs for each region on each recording day using `plot_csd.m`

    * **Output:** `<date>/csd_<region>.mat`

2. Use CSDs to assign channel numbers to layers in the CSD matfiles using `pick_channels_from_csd.m`
  
    * **Output:** `<date>/csd_<region>.mat`

3. Find burst suppression using `find_likely_bs.mat` and mark any other artifacts manually in the `artifacts` option to `multitaper_analysis`

4. Do multitaper analyses of each recording's selected channels, with artifact filtering etc. using `multitaper_analysis.m`

    * **Output:** `<date>/<time>/mt_res_layers.mat`

5. Do NMF and extract classes and transitions for each channel (based on dominant component per time) using `concat_and_nmf.m`

    * **Output:** `<date>/nmf_res.mat`

6. Generate "null model" score data for each channel (shuffled using Markov chain - negative control for connectivity) using `gen_null_model_data.m`

    * **Output:** `<date>/nmf_res.mat`

7. Compute Euclidean distance between aligned NMF scores of each directed pair of channels (and plot results) with `score_dist_analysis.m`

    * **Output:** `<date>/nmf_res.mat`

8. Get normalized mutual information between discrete classes of each channel using `class_mut_info.m` from `Ethan_misc` (see `full_analysis_byday.m`)

9. Compute a measure of synchrony for each transition and compare to bootstrap with transitions shuffled using `get_state_transitions.m` and `calc_transition_synchrony.m`
