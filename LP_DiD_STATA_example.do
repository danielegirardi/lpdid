/*
	This do file applies the Local Projections Difference-in-Differences (LP-DiD) estimator in a simulated example dataset
	Related paper: Dube, Girardi, Jorda' and Taylor (2025) "A Local Projections Approach to Difference-in-Differences", Journal of Applied Econometrics 40 (7), https://doi.org/10.1002/jae.70000
	This and other example files can be downloaded at https://github.com/danielegirardi/lpdid/ 
    Author: Daniele Girardi (King's College London), daniele.girardi@kcl.ac.uk 
	20 Feb 2026
*/

version 19
clear*
set scheme plotplainblind

********************************************
***                                      ***
***  (1) Upload dataset and preparation  ***
***                                      ***
********************************************

*** Upload simulated dataset with staggered binary treatment. 
* Note that treatment here is absorbing: once a unit gets treated, it stays treated forever.
* (Other example files in this repository demonstrate the use of LP-DiD with non-absorbing treatment, where units can enter and exit treatment multiple times.)
use http://fmwww.bc.edu/repec/bocode/l/lpdidtestdata1.dta

*** Set estimation window
local post_window 	5
local pre_window 	5

*** Generate long differences of the outcome, to be used on the left side of the LP-DiD regressions
forval h = 0/`post_window' {
	qui gen D`h'y = F`h'.Y - L.Y
	}

forval h = 2/`pre_window' {
	qui gen Dm`h'y = L`h'.Y - L.Y
}

** Generate a variable storing the time horizon of the estimate
gen horizon = _n - `pre_window' - 1
replace horizon=. if horizon>`post_window'

*** Compute & store the true (equally-weighted) average effect
gen true_att =.
gen min_true_att = .
gen max_true_att = .

forval h = 0/`post_window'{
	qui sum effect if event_time==`h', det
	qui replace true_att = r(mean) 	  if horizon==`h'
	qui replace min_true_att = r(min) if horizon==`h'
	qui replace max_true_att = r(max) if horizon==`h'

}

qui replace true_att =0 		if horizon<0
qui replace min_true_att = 0  	if horizon<0
qui replace max_true_att =0 	if horizon<0


********************************************
***                                      ***
***  (2) LP-DiD event study estimates    ***
***                                      ***
********************************************

*** Perform LP-DiD estimates - baseline version (estimates a variance-weighted ATT with strictly positive weights)

* Generate a variable where variance-weighted estimates and standard errors will be stored
qui gen b_lpdid_vw = .
qui gen se_lpdid_vw = .

* Run LP-DiD regressions
forval h = 0/`post_window' {
	qui reghdfe D`h'y 	    						///
			D.treat    							  	///   /* treatment indicator */
	if 		D.treat==1 | F`h'.treat==0,				///   /* clean control condition */
			absorb(time) vce(cluster unit)				  /* time indicators */

	qui replace b_lpdid = _b[D.treat] if horizon==`h'
	qui replace se_lpdid = _se[D.treat] if horizon==`h'
	}

forval h = 2/`pre_window' {
		qui reghdfe Dm`h'y  						///
			D.treat 								///   /* treatment indicator */
		if 	D.treat==1 | treat==0,					///   /* clean controls condition */
			absorb(time) vce(cluster unit)			  	  /* time indicators */

	qui replace b_lpdid_vw = _b[D.treat] if horizon==-`h'
	qui replace se_lpdid = _se[D.treat] if horizon==-`h'
}

replace b_lpdid = 0 if horizon==-1
replace se_lpdid = . if horizon==-1


*** Perform LP-DiD estimates - reweighted version (estimates the equally-weighted ATT)

* Generate a variable where reweighted estimates and standard errors will be stored
qui gen b_rw_lpdid=.
qui gen se_rw_lpdid=.

* compute the weight assigned to each observation in the variance-weighted version, and take the inverse  of it
forval h = 0/`post_window' {
	qui gen group_h`h'=.
	qui replace group_h`h'=time if (D.treat==1|F`h'.treat==0)
	qui reghdfe D.treat if (D.treat==1|F`h'.treat==0), absorb(time) residuals(num_weights_`h') nosample
	qui replace num_weights_`h'=. if D.treat!=1
	qui egen den_weights_`h' = total(num_weights_`h')
	qui gen weight_`h' = num_weights_`h'/den_weights_`h'
	qui bysort group_h`h': egen gweight_`h'=max(weight_`h')
	qui replace weight_`h'=gweight_`h' if weight_`h'==.
	qui replace weight_`h'=round(weight_`h',0.000000001)
	qui gen reweight_`h'=1/weight_`h'
	qui sort unit time
}
drop weight_* den_weights* group_h* gweight* num_weights*
order reweight*

* run re-weighted LP-DiD regressions
forval h = 0/`post_window' {
	qui reghdfe D`h'y      							  ///
			D.treat    							  ///   		/* treatment indicator */
			if 		(D.treat==1 | F`h'.treat==0)  ///   		/* clean controls condition */
			[pweight=reweight_`h'],				  ///		    /* re-weighting to get equally-weighted ATT */
			absorb(time) vce(cluster unit)						/* time indicators */

	qui replace b_rw_lpdid = _b[D.treat] if horizon==`h'
	qui replace se_rw_lpdid = _se[D.treat] if horizon==`h'
}

forval h = 2/`pre_window' {
	qui reghdfe Dm`h'y  							///
			D.treat 							///   /* treatment indicator */
			if 		(D.treat==1 | treat==0)		///   /* clean controls condition */
			[pweight=reweight_0],				///	  /* re-weighting to get equally-weighted ATT */
			absorb(time) vce(cluster unit)			  /* time indicators */

	qui replace b_rw_lpdid = _b[D.treat] if horizon==-`h'
	qui replace se_rw_lpdid = _se[D.treat] if horizon==-`h'
}

qui replace b_rw_lpdid = 0 if horizon==-1
qui replace se_rw_lpdid = . if horizon==-1


*** Perform LP-DiD Estimates - alternative (and equivalent) way to implement reweighted LP-DiD, using regression adjustment
qui gen b_lpdid_ra=.
qui gen se_lpdid_ra=.
qui gen dtreat=D.treat
forval h = 0/`post_window' {
	dis as text "Estimating regression adjustment at horizon " `h'
	qui teffects ra (D`h'y  i.time) (dtreat)	if 		D.treat==1 | F`h'.treat==0, atet iterate(0) vce(cluster unit)
	qui replace b_lpdid_ra = r(table)[1,1] if horizon==`h'
	qui replace se_lpdid_ra = r(table)[2,1] if horizon==`h'
}

forval h = 2/`pre_window' {
	dis as text "Estimating regression adjustment at horizon minus " `h'
	qui teffects ra (Dm`h'y  i.time) (dtreat)	if 		D.treat==1 | treat==0, atet iterate(0) vce(cluster unit)
	qui replace b_lpdid_ra = r(table)[1,1] if horizon==-`h'
	qui replace se_lpdid_ra = r(table)[2,1] if horizon==-`h'
}

qui replace b_lpdid_ra = 0 if horizon==-1
qui replace se_lpdid_ra = . if horizon==-1


********************************************
***                                      ***
***  (3) Pooled estimates                ***
***                                      ***
********************************************

*** Perform pooled LP-DiD Estimates (average effect over the post-event and pre-event windows)

*** Generate pooled dependent variables: average Y over the relevant window, then subtract L.Y (value at t-1)
* In filter(), negative lags denote leads; "normalize" divides by the number of terms

* Post pooled: average of Y in periods 0 to +post_window, minus Y at t-1
qui egen aveFY = filter(Y), lags(-`post_window'/0) normalize
qui gen pooled_post_y = aveFY - L.Y
drop aveFY

* Pre pooled: average of Y in periods -pre_window to -2, minus Y at t-1 (h=-1 is the reference period)
qui egen aveFY = filter(Y), lags(`pre_window'/2) normalize
qui gen pooled_pre_y = aveFY - L.Y
drop aveFY

*** Pooled estimates - variance-weighted version
* Post: clean control condition is the most restrictive (units clean through end of the post window)
qui reghdfe pooled_post_y D.treat							///
	if (D.treat==1 | F`post_window'.treat==0),				///
	absorb(time) vce(cluster unit)
scalar b_lpdid_vw_post_pool  = _b[D.treat]
scalar se_lpdid_vw_post_pool = _se[D.treat]
* Pre: clean control condition same as pre event-study
qui reghdfe pooled_pre_y D.treat							///
	if (D.treat==1 | treat==0),								///
	absorb(time) vce(cluster unit)
scalar b_lpdid_vw_pre_pool  = _b[D.treat]
scalar se_lpdid_vw_pre_pool = _se[D.treat]

*** Pooled estimates - reweighted version (equally-weighted ATT), implemented with weighted regression
* Post: use weights for h=post_window 
qui reghdfe pooled_post_y D.treat							///
	if (D.treat==1 | F`post_window'.treat==0)				///
	[pweight=reweight_`post_window'],						///
	absorb(time) vce(cluster unit)
scalar b_rw_lpdid_post_pool  = _b[D.treat]
scalar se_rw_lpdid_post_pool = _se[D.treat]
* Pre: for absorbing treatment, pre-period weights equal reweight_0
qui reghdfe pooled_pre_y D.treat							///
	if (D.treat==1 | treat==0)								///
	[pweight=reweight_0],									///
	absorb(time) vce(cluster unit)
scalar b_rw_lpdid_pre_pool  = _b[D.treat]
scalar se_rw_lpdid_pre_pool = _se[D.treat]

*** Pooled estimates - reweighted version (equally-weighted ATT), implemented with regression adjustment 
* Post
qui teffects ra (pooled_post_y i.time) (dtreat)				///
	if D.treat==1 | F`post_window'.treat==0,				///
	atet iterate(0) vce(cluster unit)
scalar b_lpdid_ra_post_pool  = r(table)[1,1]
scalar se_lpdid_ra_post_pool = r(table)[2,1]
* Pre
qui teffects ra (pooled_pre_y i.time) (dtreat)				///
	if D.treat==1 | treat==0,								///
	atet iterate(0) vce(cluster unit)
scalar b_lpdid_ra_pre_pool  = r(table)[1,1]
scalar se_lpdid_ra_pre_pool = r(table)[2,1]

* Compute true pooled effect: average of the effect variable over the post-treatment window
* (true pooled pre effect is zero since the no anticipation assumption holds in this data)
qui sum effect if event_time >= 0 & event_time <= `post_window'
scalar true_pooled_post = r(mean)

* Collect pooled results into a matrix and display
* Row 1 is the true ATT; rows 2-4 are the three LP-DiD estimators (actually two, since the two versions of reweighted LP-DiD are numerically equivalent)
mat pooled = J(4, 4, .)
mat rownames pooled = "True_ATE" "VW_LP-DiD" "RW_LP-DiD" "RW_LP-DiD_RA"
mat colnames pooled = "Post-Coef" "Post-SE" "Pre-Coef" "Pre-SE"
mat pooled[1,1] = true_pooled_post
mat pooled[1,3] = 0
mat pooled[2,1] = b_lpdid_vw_post_pool
mat pooled[2,2] = se_lpdid_vw_post_pool
mat pooled[2,3] = b_lpdid_vw_pre_pool
mat pooled[2,4] = se_lpdid_vw_pre_pool
mat pooled[3,1] = b_rw_lpdid_post_pool
mat pooled[3,2] = se_rw_lpdid_post_pool
mat pooled[3,3] = b_rw_lpdid_pre_pool
mat pooled[3,4] = se_rw_lpdid_pre_pool
mat pooled[4,1] = b_lpdid_ra_post_pool
mat pooled[4,2] = se_lpdid_ra_post_pool
mat pooled[4,3] = b_lpdid_ra_pre_pool
mat pooled[4,4] = se_lpdid_ra_pre_pool


********************************************
***                                      ***
***  (4) Display results                 ***
***                                      ***
********************************************

*** Event study estimates table
format true* b_* se_* %9.2f
list horizon true* b_lpdid_vw se_lpdid_vw b_rw_lpdid se_rw_lpdid b_lpdid_ra se_lpdid_ra if horizon!=.

*** Event study graph
tw rarea  min_true_att max_true_att     horizon if horizon<=`post_window', fcolor(blue%20)  lcolor(white) 				 ||   ///
   line true_att 		   				horizon if horizon<=`post_window', lcolor(gray%50)  lpat(solid) lwidth(thick) 	 ||   ///
   line b_lpdid_vw 		   				horizon if horizon<=`post_window', lcolor(green%80) lpat(longdash) 				 ||   ///
   line b_rw_lpdid		   				horizon if horizon<=`post_window', lcolor(orange)   lpat(solid)	lwidth(thin)	 ||   ///
   line b_lpdid_ra		   				horizon if horizon<=`post_window', lcolor(brown)    lpat(longdash)	lwidth(thin) 	  ///
   name(effects_1) 		 ///
   legend(order( 2 "True Equally Weighted ATE" 1 "Full Range of Treatment Effects" 3 "Variance-weighted LP-DiD"  4 "Reweighted LP-DiD" 5 "Reweighted LP-DiD (RA implementation)") pos(6) col(2)) ///
   xtitle(Event Time) ytitle(Treatment Effect)  ///
   title("Actual and Estimated Treatment Effects")

*** Pooled estimates table
matlist pooled, format(%9.2f) title("LP-DiD Pooled Estimates")

**********************************************************
***                                      			   ***
***  (5) Do the same using the lpdid STATA command     ***
***                                      			   ***
**********************************************************

* Install the lpdid command (if not installed yet)
ssc install lpdid, replace

* baseline variance-weighted LP-DiD estimates
lpdid Y, unit(unit) time(time) treat(treat) pre_window(`pre_window') post_window(`post_window')

* reweighted LP-DiD estimates using the lpdid 
lpdid Y, unit(unit) time(time) treat(treat) pre_window(`pre_window') post_window(`post_window') rw
