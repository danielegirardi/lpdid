/* 
	This do file applies the LP-DiD estimator (Dube, Girardi, Jorda' and Taylor, 2023) in a simulated example dataset.
	Related paper: Dube, Girardi, Jorda' and Taylor (2023) "A Local Projections Approach to Difference-in-Differences Event Studies", NBER Working Paper 31184, http://www.nber.org/papers/w31184
	This example illustrates the case of non-absorbing treatment (treatment can turn on and off), using a time-window for defining clean controls (see Section 3.2 in DGJT, 2023).
	This and other example files can be downloaded at https://github.com/danielegirardi/lpdid/ 
    Author: Daniele Girardi (King's College London and University of Massachusetts Amherst), dgirardi@umass.edu
	June 25, 2023 (revised Aug 14, 2023)
*/

clear*
version 17

set scheme plotplainblind
set seed 20210522

// Set DGP parameters for the simulated dataset
local units 		5000
local time_periods  60
local shocks_sd 	25
local persistence 	0.10
local post_window 	10
local pre_window 	5

**********************************************************************
***  1 - SIMULATE DATASET WITH BINARY NON-ABSORBING TREATMENT      ***
**********************************************************************

// gen random event dates (missing value for never treated units)
local max_duration = `time_periods' - 10
local eventsend = `time_periods' - 10
qui set obs `units'
qui gen unit = _n
qui gen never_treated = round(runiform(0,0.55))
qui gen first_edate = round(runiform(11,20)) if never_treated==0
qui gen treatment_duration = round(runiform(5,`max_duration')) if never_treated==0
qui gen exit_date = first_edate + treatment_duration
qui gen exit_duration = round(runiform(5,`max_duration')) if never_treated==0
qui gen second_edate = exit_date + exit_duration
qui replace exit_date = .  if exit_date>=`eventsend' 
qui replace second_edate=. if second_edate>=`eventsend' 
qui keep unit first_edate exit_date second_edate never_treated
tempfile event_dates
qui save `event_dates'

// gen time-specific fixed effects 
clear
qui set obs `time_periods'
qui gen time=_n
qui gen time_fe = rnormal(0, `shocks_sd')
tempfile time_fe
qui save `time_fe'

// gen unit-specific fixed effects
clear
qui set obs `units'
qui gen unit_fe = rnormal(0, `shocks_sd')
qui gen unit = _n
keep unit unit_fe 
tempfile unit_fe
qui save `unit_fe'

clear

// gen observations
local tot_obs = `units' * `time_periods'
qui set obs `tot_obs'
qui egen unit = seq(), b(`time_periods')
qui egen time = seq(), f(1) t(`time_periods')

qui xtset unit time

// match observations with their fixed effects and event dates
qui merge m:1 time using `time_fe', nogen
qui merge m:1 unit using `unit_fe', nogen
qui merge m:1 unit using `event_dates', nogen
qui xtset unit time

// compute untreated potential outcomes (Y(0)_{it})
qui gen white_noise = rnormal(0, `shocks_sd')
qui gen y_0=.
qui replace y_0 = unit_fe + time_fe + white_noise 							if time==1
qui replace y_0 = unit_fe + time_fe + white_noise + (`persistence' * l.y_0) if time>1

// create event time indicators (=0 at the time of the event, missing if never treated)
qui gen 	first_etime = time - first_edate
qui replace first_etime = . if (exit_date!=. & time>=exit_date)

qui gen 	exit_etime = time - exit_date
qui replace exit_etime = . if (second_edate!=. & time>=second_edate)

qui gen 	second_etime = time - second_edate

qui gen event_time = .
qui replace event_time = first_etime  if first_etime!=.
qui replace event_time = second_etime if first_etime==. & second_etime!=.

// create binary treatment indicator
qui gen 	treat= 0
qui replace treat= 1 if (first_etime>=0 & first_etime!=.) | (second_etime>=0 & second_etime!=.)

// Assign (dynamic & heterogeneous) treatment effects to treated units.
* Effects of entering and exiting treatment are both gradual and both stabilize 5 periods after treatment. Early adopters experience stronger effects.
qui gen effect =0
qui scalar P = 0.5
qui egen min_ed = min(first_edate)
qui gen denominator1 = (first_edate/min_ed)^(2)
qui gen denominator2 = (second_edate/min_ed)^(2)

qui replace effect = 2*(first_etime+1) + P*((first_etime+1)^2) + (1-P)*(((first_etime+1)^2)/(denominator1)) if first_etime>=0 & first_etime<=5
qui replace effect = L.effect if first_etime > 5 & first_etime!=.
qui replace effect = L.effect/(2*(exit_etime+2)) if exit_etime>=0 & exit_etime<5
qui replace effect = 0 if exit_etime!=. & exit_etime>=5
qui replace effect = 2*(second_etime+1) + P*((second_etime+1)^2) + (1-P)*(((second_etime+1)^2)/(denominator2)) if second_etime>=0 & second_etime<=5
qui replace effect = L.effect if second_etime > 5 & second_etime!=.

// compute observed outcomes (Y_{it})
qui gen Y = y_0 + effect
	 
// keep only relevant variables
qui keep unit time Y effect treat first_etime exit_etime second_etime first_edate exit_date second_edate never_treated event_time
 
 
************************************
***  2 - COMPUTE TRUE ATT        ***
************************************
xtset unit time
gen t =_n - `pre_window' - 1
replace t=. if t>`post_window'


// Compute and store true (equally-weighted) average effect
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

****************************************
***  3 - EVENT-STUDY TWFE ESTIMATION ***
****************************************

// Event-study TWFE estimates
gen b_twfe =.
reghdfe Y L(-10/20).D.treat, absorb(unit time) cluster(unit)
	
forval h = 2/`pre_window' {
	lincom F`h'D.treat - FD.treat
	replace b_twfe = r(estimate) if t==-`h'
}
	
replace b_twfe = 0 if t==-1
	
forval h = 0/`post_window' {
	lincom L`h'D.treat - FD.treat
	replace b_twfe = r(estimate) if t==`h'
	}


**********************************************************
***  4 - LP-DiD USING ONLY NOT-YET-TREATED AS CONTROLS ***
**********************************************************	
* [Note: using only not-yet-treated units as controls is possible in this case because there is a significant number of not yet treated units at all time periods, but might not be possible in other settings]
	
// Create variables that are useful for running LP-DiD regressions

* Gen forward changes in outcome, to be used on the left side of the LP-DiD regressions	
forval h = 0/`post_window' {
	qui gen D`h'y = F`h'.Y - L.Y
	}

forval h = 2/`pre_window' {
	qui gen Dm`h'y = L`h'.Y - L.Y
}

* set time-window for defining clean controls (note: in doing so, you are assuming that effects stabilize after L periods; see Section 3.2 in DGJT, 2023)
local L = 5

* gen indicators for clean control samples at each time horizon 
local string "abs(L.D.treat)!=1"
forv k=2/`L' {
	local string = "`string' & abs(L`k'.D.treat)!=1"
}
disp "`string'"

gen CCS_0 = 0
replace CCS_0 = 1 if `string'

forval j = 1/`post_window' {
	local i = `j'-1 
	gen CCS_`j' = 0
	replace CCS_`j' = 1 if CCS_`i'==1 & abs(F`j'.D.treat)!=1
	}	

gen CCS_m1 = CCS_0

forval j = 2/`pre_window' {
	local i = `j'-1 
	gen CCS_m`j' = 0
	replace CCS_m`j' = 1 if CCS_m`i'==1 & L.CCS_m`i'==1
	}
	
* gen indicators for not-yet treated units
by unit: gen past_events=sum(abs(D.treat))
gen nyt = 0
replace nyt =1 if past_events==0
forval h = 0/`post_window' {
	gen nyt_`h' = 1 if F`h'.past_events==0
	replace nyt_`h' = 0 if nyt_`h' != 1
	}	
	
// 4a - LP-DiD using only not-yet-treated as controls and ruling out any composition effects (see section 2.10 in DGJT, 2023 on ruling out composition effects)
gen b_lpdid_4a=.

forval h = 0/`post_window' {
	qui reghdfe D`h'y 	    							 		 			 			 	 						///
				D.treat    							  	 		 				 			 						///   /* treatment indicator */
	if 			(D.treat==1 & CCS_`post_window'==1 & CCS_m`pre_window'==1) | (D.treat==0 & nyt_`post_window'==1),	///   /* clean control condition */
				absorb(time) vce(cluster unit)				   		   							   							/* time indicators */
		
		qui replace b_lpdid_4a = _b[D.treat] if t==`h'
	
		if `h'>1 & `h'<=`pre_window' {
			qui reghdfe Dm`h'y  																							///
						D.treat 																							///   	/* treatment indicator */
			if 			(D.treat==1 & CCS_`post_window'==1 & CCS_m`pre_window'==1) | (D.treat==0 & nyt_`post_window'==1),	///   	/* clean controls condition */
						absorb(time) vce(cluster unit)			  	  		  					  									/* time dummies */
			
			qui replace b_lpdid_4a = _b[D.treat] if t==-`h'
		}
	}
	
replace b_lpdid_4a = 0 if t==-1
	
// 4b - LP-DiD using only not-yet-treated as controls but without ruling out composition effects
gen b_lpdid_4b=.

forval h = 0/`post_window' {
	qui reghdfe D`h'y 	    							 		 			///
				D.treat    							  	 		 			///   /* treatment indicator */
	if 			(D.treat==1 & CCS_`h'==1) | (D.treat==0 & nyt_`h'==1),	 	///   /* clean control condition */
				absorb(time) vce(cluster unit)				   		   			  /* time indicators */
		
		qui replace b_lpdid_4b = _b[D.treat] if t==`h'
	
		if `h'>1 & `h'<=`pre_window' {
			qui reghdfe Dm`h'y  											///
					D.treat 												///   /* treatment indicator */
			if 		(D.treat==1 & CCS_m`h'==1) | (D.treat==0 & nyt==1),		///   /* clean controls condition */
					absorb(time) vce(cluster unit)			  	  		  		  /* time dummies */
			
			qui replace b_lpdid_4b = _b[D.treat] if t==-`h'
		}
	}
	
replace b_lpdid_4b = 0 if t==-1


*******************************************************
***  5 - LP-DiD USING ALL ADMISSIBLE CLEAN CONTROLS ***
*******************************************************

// 5a - LP-DiD using a time window for clean controls (previously treated units re-enter the control sample L periods after treatment) and ruling out composition effects
gen b_lpdid_5a=.

forval h = 0/`post_window' {
	qui reghdfe D`h'y 	    							 		 			 				///
			D.treat    							  	 		 				 				///   /* treatment indicator */
	if 		(D.treat==1 | D.treat==0)  & (CCS_`post_window'==1 &  CCS_m`pre_window'==1),	///   /* clean control condition */
			absorb(time) vce(cluster unit)				   		   								  /* time indicators */
		
		qui replace b_lpdid_5a = _b[D.treat] if t==`h'
	
		if `h'>1 & `h'<=`pre_window' {
			 qui reghdfe Dm`h'y  																		///
						 D.treat 																		///     /* treatment indicator */
			if 			(D.treat==1 | D.treat==0)  & (CCS_`post_window'==1 & CCS_m`pre_window'==1),		///     /* clean controls condition */
						absorb(time) vce(cluster unit)			  	  		  			  						/* time indicators */
			
			qui replace b_lpdid_5a = _b[D.treat] if t==-`h'
		}
	}
	
replace b_lpdid_5a = 0 if t==-1

// 5b - LP-DiD using a time window for clean controls (without ruling out composition effects)	
gen b_lpdid_5b=.

forval h = 0/`post_window' {
	qui reghdfe D`h'y 	    							 	 ///
				D.treat    							  	 	 ///   /* treatment indicator */
	if 			(D.treat==1 | D.treat==0)  & (CCS_`h'==1),	 ///    /* clean control condition */
				absorb(time) vce(cluster unit)				   		/* time indicators */
		
		qui replace b_lpdid_5b = _b[D.treat] if t==`h'
	
		if `h'>1 & `h'<=`pre_window' {
			qui reghdfe Dm`h'y  								    	///
						D.treat 										///   /* treatment indicator */
			if 			(D.treat==1 | D.treat==0)  & (CCS_m`h'==1),		///   /* clean controls condition */
						absorb(time) vce(cluster unit)			  	  		  /* time dummies */
			
			qui replace b_lpdid_5b = _b[D.treat] if t==-`h'
		}
	}
	
replace b_lpdid_5b = 0 if t==-1
	
// 5c - LP-DiD using a time window for clean controls (without ruling out composition effects) and re-weighting to get equally-weighted effect (see Section 2.7 in DGJT, 2023 about reweighting)
* (Note: here I'm not computing standard errors for brevity; to compute standard errors, one can either use bootstrap, or the 'teffects ra' STATA command - see other LP-DiD example files)
qui gen b_lpdid_5c = .

forval h = 0/`post_window' {
		
	qui reg D`h'y i.time if (D.treat==0 & CCS_`h'==1)
	qui predict yhat_`h' if (D.treat==1 & CCS_`h'==1), xb

	qui sum yhat_`h' if (D.treat==1 & CCS_`h'==1)
	local po = r(mean)

	qui sum D`h'y if (D.treat==1 & CCS_`h'==1)
	local obs = r(mean)
			
	qui replace b_lpdid_5c = (`obs' - `po') if t==`h'
		
	if `h'>1 & `h'<=`pre_window' {
				
		qui reg Dm`h'y i.time if (D.treat==0 & CCS_m`h'==1)
		qui predict yhat_m`h' if (D.treat==1 & CCS_m`h'==1), xb

		qui sum yhat_m`h' if (D.treat==1 & CCS_m`h'==1)
		local po = r(mean)

		qui sum Dm`h'y if (D.treat==1 & CCS_m`h'==1)
		local obs = r(mean)
				
		qui replace b_lpdid_5c = (`obs' - `po') if t==-`h'
				
				}
		}
	
qui replace b_lpdid_5c = 0 if (t==-1)
	
order t true* b_* 

list t true* b_* if t<=`post_window'
 
******************************
***  ESTIMATES GRAPH       ***
******************************

tw rarea  min_true_att max_true_att     t if t<=`post_window', fcolor(blue%20)  lcolor(white) 				 ||   ///
   line true_att 		   				t if t<=`post_window', lcolor(gray%50)  lpat(solid) lwidth(thick) 	 ||   ///
   line b_twfe 		   				    t if t<=`post_window', lcolor(red) 		lpat(dash) 					 ||   ///
   line b_lpdid_4a 		   				t if t<=`post_window', lcolor(green%80) lpat(longdash) 				 ||   ///
   line b_lpdid_5a		   				t if t<=`post_window', lcolor(orange)   lpat(solid)	lwidth(thin)	 ||   ///
   line b_lpdid_5c			   			t if t<=`post_window', lcolor(brown)    lpat(longdash)	lwidth(thin) 	  ///
   name(effects_1) 		 ///
   legend(order(1 "Full Range of Treatment Effects" 2 "True Equally Weighted ATE"  3 "Event-study TWFE" 4 "LP-DiD (4a)"  5 "LP-DiD (5a)" 6 "Rw LP-DiD (5c)") pos(6) col(2)) ///
   xtitle(Event Time) ytitle(Treatment Effect) ///
   note("Units: `units', Periods: `time_periods'" "Shock SD: `shocks_sd', Persistence: `persistence'", ring(0) pos(10) size(*0.8) box) ///
   title("Actual and Estimated Treatment Effects") 

