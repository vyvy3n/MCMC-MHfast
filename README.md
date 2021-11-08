# R Re-implementation: Code for Log-concave sampling 
R implementation of the paper: [Log-concave sampling: Metropolis-Hastings algorithms are fast!](https://arxiv.org/abs/1801.02309)

The original implementation in Python by the author can be founded [here](https://github.com/yuachen/mala_public).
## Description
`main.R`: to run repeated experiments, by taking distribution of starting point as proposed in Table 2 on page 8 of the paper. Output the estimation L1 error over iterations.

`mainSingle.R`: to run experiment for only once, provided with more plot analysis; including ACF (Auto-Correlation Function), CCF  (Cross-Correlation Function), Trace Plot, Estimation Density.

`another-version-plots/logConcave.R`: an early version, with other plots.

`samplers.R`: MCMC sampling methods discussed in the paper, including MRW, ULA, MALA.

`by-stan/`: MCMC for the same posterior estimation problem by stan.

`simple-starts/`: simple examples of MH, MCMC algorihtms.

## Example results
The examples of results of `main.R` can be found in `eg-results/`;

The examples of results of `mainSingle.R` can be found in `eg-results-single/`;

More related plots can be found in `another-version-plots/eg-results/`.
## Notes 
Please note that ULA stands for unadjusted and hence there is no accept-reject step. 

MRW stands for  "Metropolized random walk", therefore the accept-reject step (the same as "MH step") is considered.
