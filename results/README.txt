Notes on Results:

## qqplots.pdf
This document includes quantile-quantile plots of parameter distributions between Stan (500 warmup iterations) and JAGS.
For a sample to be considered as being drawn from the same distribution, the regression should lie as close to the line
y = x as possible. If the distributions are linearly related, then the points should lie along a line but not necessarily
along y = x.

## jags_ggmcmc_output.pdf, stan_ggmcmc_output_500.pdf, stan_ggmcmc_output_20.pdf
These documents contain various plots based on the MCMC chains from the ggmcmc package.

- Histograms
  aggregate histogram counts over all three chains
- Density plots
  estimated probability density curves per chain
- Traceplots
  time-series traceplots of values per chain, to visually gauge chain mixing
- Running Means
  overall mean value at different points of each chain
- Partial Chain Densities
  estimated probability density curves from the last 10% entries of each chain against complete chain density curves
- Autocorrelation
  autocorrelation values at different lags for each chain, bounded between -1 and 1. Low autocorrelation is associated
  with short peaks and the absenceo of a discernible pattern. Moderate autocorrelation involves several comparable 
  peaks at adjacent lags
- Cross-correlation
  shows correlations between all parameters
- Potential Scale Reduction Factor
  a plot of the Gelman-Rubin convergence diagnostic function (R_hat) for each parameter. A value close to 1.0 implies
  good convergence
- Geweke Diagnostic Function
  a plot of Geweke's diagnostic function for each parameter - a standardized Z-score of the means over the first 10%
  and last 50% of the chain.
- Caterpillar Plots
  plots of the Highest Posterior Densities (HPD) of each parameter. The dot represents the mean, the thick line 
  represents the 90% HPD and the thin line depicts the 95% HPD.
