# bnbwordfish

Replication materials for Däubler, Thomas and Kenneth Benoit (2002) "[Scaling hand-coded political texts to learn more about left-right policy content](/pdfs/daubler-benoit-2021-scaling-hand-coded-political-texts-to-learn-more-about-left-right-policy-content.pdf)."  _Party Politics_. 28(5, September): 834–844. [10.1177/13540688211026076](https://doi.org/10.1177/13540688211026076)

## Files for replication

### Code

- `PP_preparedata.R`: this prepares the data for analysis

- `PP_runsimul.R`: this does the HMC sampling of the measurement models

- `PP_analysis.R`: this analyzes the results and prepares the graphs

- `HMC_R_functions_PP.R`: functions that the above scripts use

- `HMC_Stan_code_PP.R`: the STAN code PP_runsimul uses through one of the R_functions.

- `resultsPP`: the results from the simulations for the CMP data.

- `resultsCAPPP`: the results from the simulations for the CAP data.  This file saves results that are too large for GitHub but that is available here: [`resultsCAPPP.RData`](https://www.dropbox.com/scl/fi/lf68t80k5q1bphla881aj/resultsCAPPP.RData?rlkey=gjnd8wae2b9blstxv4um524rn&dl=0)https://www.dropbox.com/scl/fi/lf68t80k5q1bphla881aj/resultsCAPPP.RData?rlkey=gjnd8wae2b9blstxv4um524rn&dl=0)

- `cmpcatlab- : some category labels for the CMP data

### Data

- `CMP2018b.RData`: 2018b version of the [Manifesto Project](https://manifestoproject.wzb.eu) data

- `2014_CHES_dataset_means.dta` and `CHES_means_2017.dta`: [Chapel Hill](https://www.chesdata.eu) expert survey data

- `MPDataset_MPDS2014b_plusexperts.dta`: an earlier CMP dataset with expert survey information added

- `partyfacts-external-parties.dta`: party link file from [PartyFacts](https://partyfacts.herokuapp.com).  Too large for GitHub, but is available [here](https://www.dropbox.com/scl/fi/2og1dnsdnybrk9tofw4jr/partyfacts-external-parties.dta?rlkey=muf5d0d5uu7yrxb28wno4wzyp&dl=0).

- `qog_tomerge_2019.dta`: Quality of Government data on country-level democracy indicators

- `PP_CAPdata.RData` and `PP_data.RData`: the files created by the script `PP_preparedata.R`

### Additional files too large for GitHub

- `resultsPP.RData` (available [here](https://www.dropbox.com/scl/fi/b5zc0wc6c62b6z1mpg1qh/resultsPP.RData?rlkey=zmzap26jg9f4c5bwofe7j8lui&dl=0))
