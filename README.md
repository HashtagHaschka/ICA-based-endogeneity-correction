The function implements the ICA-based endogeneity correction by Haschka & Dost (2025). It corrects omitted variable bias without instruments: the response and the endogenous regressors are decomposed by Independent Component Analysis, the one recovered source that is Gaussian estimates the omitted confounder, and that source enters the structural equation as a control function.

The required arguments should be specified as follows:

- formula should be depvar ~ endog_var1 + endog_var2 + ... | exog_var1 + exog_var2 + ...
- depvar is the dependent variable, endog_var1, etc., are the continuous endogenous variables, exog_var1, etc., are the exogenous variables
- the exogenous block may be omitted entirely
- formula accepts -1 to remove intercept, which must be placed in the first part
- transformations are allowed on both sides, e.g. log(sales) ~ log(price) | promotion
- dummy variables can be modelled using as.factor(exog_var1), etc.
- a variable must not appear in both blocks
- data should be a data.frame
- Example: IcaReg(formula = Y ~ endog1 + endog2 + endog3 | exog1 + exog2 + as.factor(exog3), data = data1)
- The file 1EXAMPLES.R contains examples with freely-available datasets and loads the function automatically

Optional arguments:

- method is the ICA algorithm passed to ica::ica() and can be selected via method = c("jade", "fast", "imax"). The default "jade" separates sources using fourth-order cumulants and matches the identification argument of the paper
- CF = TRUE additionally reports the coefficient on the control function, which is a Hausman-type test statistic for endogeneity: a coefficient indistinguishable from zero indicates no detectable endogeneity. All other coefficients are identical to CF = FALSE
- select is the rule identifying the Gaussian source and can be selected via select = c("kurtosis", "ks"). The default "kurtosis" picks the source with the smallest absolute excess kurtosis, which is the quantity identification turns on; "ks" picks the smallest Kolmogorov-Smirnov distance to normality. The two agree whenever identification is clear-cut, so comparing them is a cheap robustness check
- kurt_tol is the threshold above which the diagnostic warning is issued (default 1)
- parallel = TRUE computes the jackknife replicates in parallel, with ncores defaulting to detectCores() - 1. Recommended for n beyond a few thousand

The function returns a list with two elements: a matrix of coefficients and delete-one jackknife standard errors, and the estimated control function. Every jackknife replicate re-runs the complete pipeline, so the sampling variability of the first-stage projection and of the estimated control function is accounted for automatically.

IDENTIFICATION

The excess kurtosis of every recovered source is reported. Identification requires exactly one source close to zero while the others are clearly away from it. Under normality the sampling standard deviation of the sample excess kurtosis is about sqrt(24 / n), which gives a scale for judging the readout. For example,

    Excess kurtosis of the ICA sources: +2.103  +53.669  -0.101

at n = 250 is clear-cut: one source sits within sampling noise of zero and the nearest competitor is several standard deviations away. Two sources close to zero mean the Gaussian source is not pinned down and the estimates should not be trusted. The value reported for the selected source is a minimum over all sources and is therefore biased towards zero, so the warning is conservative.

Further points to note:

- endogenous regressors must be continuous; ties are reported because discreteness destabilises the higher-order moments the method relies on
- identification is weak when the exogenous variation in the endogenous regressor is close to normal, which shows up as large standard errors
- the jackknife assumes independent observations. With clustered or panel data the standard errors will be too small

REFERENCES

- Haschka, R. E. and F. Dost (2025). ICA at the cocktail party: Casting instrument-free omitted variable bias correction as a blind source separation problem. Available at SSRN 5361801.
