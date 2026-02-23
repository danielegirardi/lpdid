# Local Projections Difference-in-Differences (LP-DiD)

LP-DiD is a convenient and flexible regression-based framework for implementing Difference-in-Differences, introduced in <a href="https://doi.org/10.1002/jae.70000" target="_blank" rel="noopener noreferrer">Dube, Girardi, Jorda' and Taylor (2025)</a> - DGJT hereafter.

LP-DiD uses local projections to estimate dynamic effects, while restricting the estimation sample to units entering treatment and 'clean' controls, thus avoiding the 'negative-weighting' bias of TWFE estimators. 

A baseline LP-DiD specification for staggered and absorbing treatment settings can be written as follows.[^1] Assume panel data for $N$ units (indexed by $i$) and $T$ periods (indexed by $t$). Let $y$ denote the outcome variable of interest, while $D$ is a binary treatment indicator, $\delta$ are time-specific intercepts, and $e$ is an error term. Let $H$ denote the length of the post-treatment window for estimating effects, and $Q$ a possible number of pre-treatment periods for assessing pre-treatment trends. Then, the baseline LP-DiD specification is

```math
y_{i,t+h} - y_{i,t-1} 
 = \,\,\beta^{LP-DiD}_h \Delta D_{it} + \delta^h_t + e^h_{it} \quad \text{for } h = -Q, ..., 0, ..., H\,,
```

restricting the estimation sample to observations that are either
```math
\begin{cases}
\text{newly treated:} \qquad \Delta D_{it}=1 \,, \\
\text{or clean control:} \quad \, D_{i, t+h} = 0 \,.

\end{cases}
```

This 'clean control condition' ensures that estimates are obtained from a set of clean DiD comparisons between newly treated units and not-yet treated ones, thus avoiding the unclean (or 'forbidden') comparisons that can introduce negative weighting bias in conventional TWFE specifications.

Under the usual DiD assumptions, OLS estimation of this baseline LP-DiD specification yields a variance-weighted average treatment effect on the treated (VWATT), giving more weight to more precisely estimated cohort-specific effects. A reweighed version of the LP-DiD specification, easily implemented through a weighted regression or regression adjustment, estimates an equally-weighted ATT, giving equal weight to all treated observations.  

By appropriately modifying the outcome variable, one can obtain a 'pooled' DiD estimate giving an overall average effect over the post-treatment window instead of (or in addition to) dynamic event study estimates. Moreover, one can choose whether to compare post-treatment outcomes to the last pre-treatment period (as in the baseline LP-DiD specification above) or to an average of several pre-treatment periods (as in the 'pre-mean differenced' (PMD) LP-DiD specification).

The LP-DiD approach allows for inclusion of covariates. Under conditional parallel trends, a regression-adjustment LP-DiD specification with covariates yields an unbiased estimate of the ATT. Under an additional assumption that treatment effects are independent of covariates, adding covariates directly to the OLS LP-DiD specification yields a variance-weighted effect (VWATT) with same weights as in the baseline OLS version.

LP-DiD can accommodate non-absorbing treatment, meaning that units can enter and exit treatment multiple times, through appropriate modification of the clean control condition. In particular, DGJT derive (i) a LP-DiD estimator for the effect of first-time treatment entry, and (ii) a LP-DiD estimator for the average effect of a treatment event under an additional effect stabilization assumption.

See <a href="https://www.nber.org/papers/w31184" target="_blank" rel="noopener noreferrer">DGJT</a> for a detailed exposition of the LP-DiD method.

[^1]: Here 'staggered' means that different units can enter treatment at different time periods, and 'absorbing' means that treatment is permanent (once a unit gets treated, it stays treated forever).

***

# LP-DiD Stata command

A Stata command implementing LP-DiD (lpdid, by Daniele Girardi and Alexander Busch) is available and can be installed from within Stata by typing 
```stata
ssc install lpdid, replace
```
You can then type 
```stata
help lpdid 
```
to read the help file, which explains the syntax and working of the command, and provides examples using simulated datasets. More details can be found on <a href="https://econpapers.repec.org/software/bocbocode/S459273.htm" target="_blank" rel="noopener noreferrer">the package listing</a>, and you can ask questions (or signal possible problems/bugs) by sending us an email or by replying to the corresponding <a href="https://www.statalist.org/forums/forum/general-stata-discussion/general/1736005-lpdid-new-module-implementing-local-projections-difference-in-differences" target="_blank" rel="noopener noreferrer">Statalist discussion</a>.

***

# Example codes for "manual" implementation of LP-DiD in Stata

While the lpdid Stata command is available, it is also easy to implement the LP-DiD estimator "manually", in the sense of writing your own Stata code for implementing LP-DiD. You might want to do this either because your application requires some bespoke adjustment, or to make sure you understand how the sausage is made. "Manually" implementing LP-DiD is easy, because the method essentially consists in estimating a simple regression (or a regression-adjustment specification) in an estimation sample defined by a 'clean control' condition.

This repository contains six Stata do files that implement the LP-DiD estimator in simulated datasets. 

Four examples illustrate the case of binary, staggered and absorbing treatment (where "absorbing" means that once a unit gets treated, it stays treated forever):

- [**LP_DiD_STATA_example.do**](LP_DiD_STATA_example.do) uses a simulated dataset somewhat similar to the Montecarlo simulations presented in DGJT, and shows how to implement the baseline variance-weighted version and the reweighted version to obtain event-study estimates or pooled estimates. 

- [**LP_DiD_examplefile.do**](LP_DiD_examplefile.do) is similar, but builds the simulated dataset from scratch instead of uploading it (allowing the user to adjust some features of the dataset if they wish so), and includes the pre-mean differenced (PMD) version of LP-DiD. 

- [**lpdid_test.do**](lpdid_test.do) applies LP-DiD in a simulated dataset from Borusyak (2021).

- [**lpdid_regression_adjustment.do**](lpdid_regression_adjustment.do) focuses on the Regression Adjustment LP-DiD estimator with covariates (<a href="https://doi.org/10.1002/jae.70000" target="_blank" rel="noopener noreferrer">Dube, Girardi, Jorda' and Taylor 2025</a>, Sec 4.1.1) and illustrates four alternative ways to implement it, using the Stata official commands 'teffects ra' or 'margins' or the user-written alternatives 'kmatch ra' or 'listreg'. They produce identical results but 'margins', 'kmatch ra' and 'listreg' are faster than 'teffects ra'. This can be especially useful in settings where 'teffects ra' is slow.

The other two examples illustrate the case of binary non-absorbing treatment, meaning that units can enter and exit treatment multiple times:

- [**LPDiD_nonabsorbing_example.do**](LPDiD_nonabsorbing_example.do) illustrates the 'persistent treatment' setting:  after a unit enters treatment, its treatment status persists (ie, the treatment variable remains equal to 1) until a possible exit or reversal. This is the case assumed in Dube, Girardi, Jorda' and Taylor (2023).  An example of this type of 'persistent' treatment is democracy:  after democratization, a polity remains a democracy until a possible reversal. 

- [**lpdid_oneoff_nonabs_example.do**](lpdid_oneoff_nonabs_example.do) deals with binary, non-absorbing and 'one-off' (or 'shock') treatments. Here "one-off" means that the treatment lasts only for 1 period by construction (although its effects can still be dynamic and persistent, and a unit can still receive treatment multiple times). Formally, we have $D_{it}=1$ if unit i experiences an event at time t, and $D_{it}=0$ in all other periods. An example of this type of "one-off" treatment is hurricanes: the treatment indicator equals 1 if the unit is hit by a hurricane at time t, and 0 in all other periods.
***

# LP-DiD example scripts for R

This repository contains scripts demonstrating how to implement the LP-DiD estimator in R:

- [**LP_DiD_R_example_VW.R**](LP_DiD_R_example_VW.R) illustrates the baseline variance-weighted version of LP-DiD in R, for both event-study estimates and pooled estimates. 

*(Additional R scripts covering the reweighted LP-DiD and non-absorbing treatments are coming soon.)*

***

Related paper: 
Dube, Girardi, Jorda' and Taylor (2025) "A Local Projections Approach to Difference-in-Differences", Journal of Applied Econometrics 40 (7), <a href="https://doi.org/10.1002/jae.70000" target="_blank" rel="noopener noreferrer">https://doi.org/10.1002/jae.70000</a>.
