

The files in this dataverse repository will replicated the findings in Blackwell and Glynn, "How to Make Causal Inferences with Time-Series Cross-Sectional Data under Selection on Observables" in the American Political Science Review. `causal-tscs-sims.R` produces the simulations from the main paper and the online appendix, along with the figures. `burgoon.R` produces the reanalysis of the Burgoon (2006) paper and the figures therein. `swank-steinmo.R` contains the additional replication in the online appendix and `panel-utils.R` provides some support functions for the main code. These files are independent and so can be run in any order. The following information summarizes the computing environment used to conduct the analyses for the paper:

- Operating system: MacOS 13.4
- R version: 3.4.3
- R packages: stargazer, xtable, reshape, rms, foreign, lmtest, sandwich, msm, coefplot, ggplot2, aod.
- CPU cores required: 1 (though we use 12 cores to speed up the simulations)
- Running time on iMac Pro (2017, 3Ghz Xeon W): `causal-tscs-sims.R` (~10 mins with 12 cores), `burgoon.R` (~5 mins with single core), `swank-steinmo.R` (~3 secs on a )
