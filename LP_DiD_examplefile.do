/* 
	This do file applies the LP-DiD estimator (Dube, Girardi, Jorda' and Taylor, 2023) in a simulated example dataset
    Author: Daniele Girardi (University of Massachusetts Amherst), dgirardi@umass.edu (2/24/2023)
	
*/

version 17
clear*

set scheme plotplainblind
set seed 20210522

* Set DGP parameters
local units 2000
local time_periods 50
local shocks_sd 25
local persistence 0.25
local post_window 10
local pre_window 5

**************************
***  1 - THE DGP       ***
**************************

* gen random event dates (missing value for never treated units)
set obs `units'
gen unit = _n
gen event_date = round(runiform(11,30))
gen never_treated = round(runiform(0,0.60))
replace event_date=. if never_treated==1
keep unit event_date
tempfile event_dates
save `event_dates'

* gen time-specific fixed effects 
clear
qui set obs `time_periods'
qui gen time=_n
qui gen time_fe = rnormal(0, `shocks_sd')
tempfile time_fe
qui save `time_fe'

* gen unit-specific fixed effects
clear
qui set obs `units'
qui gen unit_fe = rnormal(0, `shocks_sd')
qui gen unit = _n
keep unit unit_fe 
tempfile unit_fe
qui save `unit_fe'

clear

* gen observations
local tot_obs = `units'* `time_periods'
set obs `tot_obs'
egen unit = seq(), b(`time_periods')
egen time = seq(), f(1) t(`time_periods')

xtset unit time

* match observations with their fixed effects and event dates
merge m:1 time using `time_fe', nogen
merge m:1 unit using `unit_fe', nogen
merge m:1 unit using `event_dates', nogen
xtset unit time

* compute untreated potential outcomes (Y_{0it})
gen white_noise = rnormal(0, `shocks_sd')
gen y_0=.
replace y_0 = unit_fe + time_fe + white_noise 							if time==1
replace y_0 = unit_fe + time_fe + white_noise + (`persistence' * l.y_0) if time>1

* create event time indicator (=0 at the time of the event, missing if never treated)
gen event_time = time - event_date

* create binary treatment indicator
gen treat= 0
replace treat= 1 if event_time>=0 & event_time!=.

* assign (dynamic & heterogeneous) treatment effects to treated units
gen effect =0
scalar P = 0.05
egen min_ed = min(event_date)
gen denominator = (event_date/min_ed)^(2)
replace effect = 2*(event_time+1) + P*(((event_time+1)^2)) + (1-P)*((event_time+1)^2)/((denominator)) if event_time>=0 & event_time<=20
replace effect = L.effect if event_time > 20 & event_time!=.

* compute observed outcomes (Y_{it})
qui gen Y = y_0 + effect
	 
* keep only relevant variables
keep unit time Y effect treat event_time event_date 
 
 
******************************
***  2 - ESTIMATION        ***
******************************
xtset unit time
local preobs = `pre_window'+1
gen t =_n - `preobs'
replace t=. if t>`post_window'


*** COMPUTE TRUE (EQUALLY WEIGHTED) AVERAGE EFFECT
egen treateffect = mean(effect) , by(event_time)

gen true_att =.
gen min_true_att = .
gen max_true_att = .

forval j = 0/`post_window'{
	
	sum effect if event_time==`j', det
	
	replace true_att = r(mean) 	  if t==`j'
	replace min_true_att = r(min) if t==`j'
	replace max_true_att = r(max) if t==`j'
	
}

replace true_att =0 		if t<0
replace min_true_att = 0  	if t<0
replace max_true_att =0 	if t<0   


*** EVENT-STUDY TWFE ESTIMATES
gen b_twfe =.
reghdfe Y L(-10/20).D.treat, absorb(unit time) cluster(unit)
	
forval j = 2/`pre_window' {
	lincom F`j'D.treat - FD.treat
	replace b_twfe = r(estimate) if t==-`j'
}
	
replace b_twfe = 0 if t==-1
	
lincom D.treat - FD.treat
replace b_twfe = r(estimate) if t==0
	
lincom LD.treat - FD.treat
replace b_twfe = r(estimate) if t==1
	
forval j = 2/`post_window' {
	lincom L`j'D.treat - FD.treat
	replace b_twfe = r(estimate) if t==`j'
	}


*** LP-DiD ESTIMATES
gen b_lpdid = .

* Gen forward changes in outcome, to be used on the left side of the LP-DiD regressions	
forval j = 0/`post_window' {
	gen D`j'y = F`j'.Y - L.Y
	}

forval j = 2/`pre_window' {
	gen Dm`j'y = L`j'.Y - L.Y
}

* Run LP-DiD regressions
forval j = 0/`post_window' {
	reghdfe D`j'y 	    							  ///
			D.treat    							  	 ///   /* treatment indicator */
	if 		D.treat==1 | F`j'.treat==0,				 ///   /* clean control condition */
			absorb(time) vce(cluster unit)				   /* time indicators */
		
		replace b_lpdid = _b[D.treat] if t==`j'
	
		if `j'>1 & `j'<=`pre_window' {
			reghdfe Dm`j'y  							///
					D.treat 							///   /* treatment indicator */
			if 		D.treat==1 | treat==0,				///   /* clean controls condition */
					absorb(time) vce(cluster unit)			  /* time dummies */
			
			replace b_lpdid = _b[D.treat] if t==-`j'
		}
	}
	
	replace b_lpdid = 0 if t==-1
	


*** ALTERNATIVE (AND EQUIVALENT) IMPLEMENTATION: LP-DiD WITH INTERACTION TERMS
gen b_int_lpdid=.

* gen 'unclean observation' indicator
forval j = 0/`post_window' {
gen UC_`j'=1
replace UC_`j'=0 if D.treat==1 | F`j'.treat==0
	}

* estimate LP-DiD - interaction version
forval j = 0/`post_window' {
	reghdfe D`j'y      							     ///
			D.treat    							  	 ///    /* treatment indicator */
			UC_`j',									/// 	/* unclean indicator */
			absorb(time#UC_`j')	 vce(cluster unit)			/* time dummies x unclean indicator */
		
		replace b_int_lpdid = _b[D.treat] if t==`j'
	
		if `j'>1 & `j'<=`pre_window' {
			reghdfe Dm`j'y  								    ///
					D.treat 								    ///   /* treatment indicator */
					UC_0,										///   /* unclean indicator */
					absorb(time#UC_0) vce(cluster unit)				  /* time dummies */
			
			replace b_int_lpdid = _b[D.treat] if t==-`j'
		}
	}
	
	replace b_int_lpdid = 0 if t==-1
	
	
	
*** Re-weighted LP-DiD to get equally-weighted ATE
gen b_rw_lpdid=.

* compute the weight assigned to each observation in the variance-weighted version, and take the inverse  of it
forval j = 0/`post_window' {
	qui gen group_h`j'=.
	qui replace group_h`j'=time if (D.treat==1|F`j'.treat==0)
	qui reghdfe D.treat if (D.treat==1|F`j'.treat==0), absorb(time) residuals(num_weights_`j')
	qui replace num_weights_`j'=. if D.treat!=1
	qui egen den_weights_`j' = total(num_weights_`j')
	qui gen weight_`j' = num_weights_`j'/den_weights_`j'
	qui bysort group_h`j': egen gweight_`j'=max(weight_`j')
	qui replace weight_`j'=gweight_`j' if weight_`j'==.
	qui replace weight_`j'=round(weight_`j',0.00000001)
	qui gen reweight_`j'=1/weight_`j'
	qui sort unit time
} 

* run re-weighted LP-DiD regressions
forval j = 0/`post_window' {
	reghdfe D`j'y      							  ///
			D.treat    							  ///   		/* treatment indicator */
			if 		(D.treat==1 | F`j'.treat==0)  ///   		/* clean controls condition */
			[pweight=reweight_`j'],				  ///		    /* re-weighting to get equally-weighted ATT */
			absorb(time) vce(cluster unit)						/* time indicators */
		
		replace b_rw_lpdid = _b[D.treat] if t==`j' 
	
		if `j'>1 & `j'<=`pre_window' {
			reghdfe Dm`j'y  							///
					D.treat 							///   /* treatment indicator */
					if 		(D.treat==1 | treat==0)		///   /* clean controls condition */
					[pweight=reweight_0],				///	  /* re-weighting to get equally-weighted ATT */
					absorb(time) vce(cluster unit)			  /* time indicators */
		
		replace b_rw_lpdid = _b[D.treat] if t==-`j'
		}
	}
	
replace b_rw_lpdid = 0 if t==-1
	
*** Alternative (and equivalent) way to get the equally-weighted ATT: LP-DiD implemented through regression adjustment	
gen b_lpdid_ra=.
gen dtreat=D.treat
qui tab time, gen(tdum)
drop tdum1 tdum50
forval j = 0/`post_window' {
	teffects ra (D`j'y  tdum*) (dtreat)	if 		D.treat==1 | F`j'.treat==0, atet iterate(0) vce(cluster unit)
	replace b_lpdid_ra = e(b)[1,1] if t==`j'
	
	if `j'>1 & `j'<=`pre_window' {
		teffects ra (Dm`j'y  tdum*) (dtreat)	if 		D.treat==1 | treat==0, atet iterate(0) vce(cluster unit)
		replace b_lpdid_ra = e(b)[1,1] if t==-`j'
		}
		
}

*** Pre-mean-differenced LP-DiD
gen b_pmd_lpdid=.
* gen pre-mean-differenced long differences
bysort unit (time) : gen cumulative_y = sum(Y)
gen aveLY = L.cumulative_y/(time-1)
forval j = 0/`post_window' {
	gen PMD`j'y = F`j'.Y - aveLY
	}
forval j = 2/`pre_window' {
	gen PMDm`j'y = L`j'.Y - aveLY
}
* run PMD LP-DID regressions
forval j = 0/`post_window' {
	reghdfe PMD`j'y 	    							  ///
			D.treat    							  	 ///   /* treatment indicator */
	if 		D.treat==1 | F`j'.treat==0,				 ///   /* clean control condition */
			absorb(time) vce(cluster unit)				   /* time indicators */
		
		replace b_pmd_lpdid = _b[D.treat] if t==`j'
	
		if `j'>1 & `j'<=`pre_window' {
			reghdfe PMDm`j'y  							///
					D.treat 							///   /* treatment indicator */
			if 		D.treat==1 | treat==0,				///   /* clean controls condition */
					absorb(time) vce(cluster unit)			  /* time dummies */
			
			replace b_pmd_lpdid = _b[D.treat] if t==-`j'
		}
	}
	
	replace b_pmd_lpdid = 0 if t==-1

*** Reweighted and pre-mean-differenced LP-DiD	
gen b_rwpmd_lpdid=.
forval j = 0/`post_window' {
	reghdfe PMD`j'y      							  ///
			D.treat    							  ///   		/* treatment indicator */
			if 		(D.treat==1 | F`j'.treat==0)  ///   		/* clean controls condition */
			[pweight=reweight_`j'],				  ///		    /* re-weighting to get equally-weighted ATT */
			absorb(time) vce(cluster unit)						/* time indicators */
		
		replace b_rwpmd_lpdid = _b[D.treat] if t==`j' 
	
		if `j'>1 & `j'<=`pre_window' {
			reghdfe PMDm`j'y  							///
					D.treat 							///   /* treatment indicator */
					if 		(D.treat==1 | treat==0)		///   /* clean controls condition */
					[pweight=reweight_0],				///	  /* re-weighting to get equally-weighted ATT */
					absorb(time) vce(cluster unit)			  /* time indicators */
		
		replace b_rwpmd_lpdid = _b[D.treat] if t==-`j'
		}
	}
	
replace b_rwpmd_lpdid = 0 if t==-1
	
order t true* b_* 

list t true* b_* if t<=`post_window'
 
******************************
***  ESTIMATES GRAPH       ***
******************************

tw rarea  min_true_att max_true_att     t if t<=`post_window', fcolor(blue%20)  lcolor(white) 				 ||   ///
   line true_att 		   				t if t<=`post_window', lcolor(gray%50)  lpat(solid) lwidth(thick) 	 ||   ///
   line b_twfe 		   				    t if t<=`post_window', lcolor(red) 		lpat(dash) 					 ||   ///
   line b_lpdid 		   				t if t<=`post_window', lcolor(green%80) lpat(longdash) 				 ||   ///
   line b_rw_lpdid		   				t if t<=`post_window', lcolor(orange)   lpat(solid)	lwidth(thin)	 ||   ///
   line b_pmd_lpdid		   				t if t<=`post_window', lcolor(violet)   lpat(longdash)	lwidth(thin) ||   ///
   line b_rwpmd_lpdid		   			t if t<=`post_window', lcolor(brown)    lpat(longdash)	lwidth(thin) 	  ///
   name(effects_1) 		 ///
   legend(order( 2 "True Equally Weighted ATE" 1 "Full Range of Treatment Effects" 3 "TWFE Distributed Lags Estimate" 4 "LP-DiD"  5 "Reweighted LP-DiD" 6 "PMD LP-DiD" 7 "Reweighted PMD LP-DiD") pos(6) col(2)) ///
   xtitle(Event Time) ytitle(Treatment Effect) ///
   note("Units: `units', Periods: `time_periods'" "Shock SD: `shocks_sd', Persistence: `persistence'", ring(0) pos(10) size(*0.8) box) ///
   title("Actual and Estimated Treatment Effects") 

