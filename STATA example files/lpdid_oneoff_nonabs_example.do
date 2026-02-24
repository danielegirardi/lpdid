/* 
	This do file implements the LP-DiD estimator (Dube, Girardi, Jorda' and Taylor, 2023 - DGJT hereafter) in a simulated example dataset.
	We focus here on a setting where treatment is non-absorbing and one-off.
	Here "one-off" means that the treatment lasts only for 1 period by construction, although its effects can still be dynamic and persistent. 
	Formally, we have D_{it}=1 if unit i experiences an event at time t, and D_{it}=0 in all other periods. 
	Examples of this type of "one-off" settings are for example minimum wage increases or hurricanes.
	Note that this definition of the treatment variable is unlike the one assumed in the main derivations in DGJT.
	The specifications and clean control conditions used here are adjusted accordingly.
	
	Related paper: Dube, Girardi, Jorda' and Taylor "A Local Projections Approach to Difference-in-Differences" (https://www.nber.org/papers/w31184)
	This and other example files can be downloaded at https://github.com/danielegirardi/lpdid/ 
	Author: Daniele Girardi (King's College London), daniele.girardi@kcl.ac.uk
	23 May 2025 (revised 24 May 2025)
	
*/

clear*
version 17

set seed 20210522

// Set DGP parameters for the simulated dataset
local units 		2500
local time_periods  50
local shocks_sd 	10
local persistence 	0.10
local post_window 	5
local pre_window 	5
local L 			3

********************************************************************************
***  1 - SIMULATE DATASET WITH BINARY NON-ABSORBING "ONE-OFF" TREATMENT      ***
********************************************************************************

// gen random event dates (missing value for never treated units)
local max_interval = `time_periods' - 20
local eventsend = `time_periods' - 10
qui set obs `units'
qui gen unit = _n
qui gen never_treated = round(runiform(0,0.55))
qui gen edate1 = round(runiform(11,20)) if never_treated==0
qui gen first_treatment_interval = round(runiform(5,`max_interval')) if never_treated==0
qui gen edate2 = edate1 + first_treatment_interval
qui gen second_treatment_interval = round(runiform(5,`max_interval')) if never_treated==0
qui gen edate3 = edate2 + second_treatment_interval
qui replace edate2 = .  if edate2>=`eventsend' 
qui replace edate3=. 	if edate3>=`eventsend' 
qui keep unit edate1 edate2 edate3 never_treated
forvalues i=1/3 {
	count if edate`i'!=. & never_treated==0
	local n`i' = r(N) 
}
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
forvalues i=1/3 {
	qui gen etime`i' = time - edate`i'
	}

// create binary treatment indicator
qui gen 	treat= 0
qui replace treat= 1 if time==edate1|time==edate2|time==edate3

// Assign (dynamic & heterogeneous) treatment effects to treated units.
qui scalar P = 0.5
qui egen min_ed = min(edate1)

forvalues i=1/3 {
	qui gen effect`i' = 0
	qui gen denominator`i' = (edate`i'/min_ed)^(2)
	qui replace effect`i' = 2*(etime`i'+1) + P*((etime`i'+1)^2) + (1-P)*(((etime`i'+1)^2)/(denominator`i')) if etime`i'>=0 & etime`i'<=`L' 
	qui replace effect`i' = L.effect`i' if etime`i' > `L' & etime`i'!=.
	}

// compute overall effect
gen tot_effect = effect1 + effect2 + effect3

// compute observed outcomes (Y_{it})
qui gen Y = y_0 + tot_effect
	 
// keep only relevant variables
qui keep unit time Y effect1 effect2 effect3 treat etime1 etime2 etime3 edate1 edate2 edate3 never_treated tot_effect


************************************
***  2 - COMPUTE TRUE ATT        ***
************************************
xtset unit time
gen horizon =_n - `pre_window' - 1
replace horizon=. if horizon>`post_window'

forvalues i = 1/3 {
	gen true_att`i' =.
	forval j = 0/`post_window'{
		sum effect`i' if etime`i'==`j'
		replace true_att`i' = r(mean) 	  if horizon==`j'
	}
}

gen true_att=.
forval j = 0/`post_window'{
	replace true_att = (true_att1*`n1' + true_att2*`n2' + true_att3*`n3')/(`n1'+`n2'+`n3')
	}

* gen indicators for clean control samples at each time horizon 
local string "L.treat==0"
forv k=2/`L' {
	local string = "`string' & L`k'.treat==0"
}
disp "`string'"

gen CCS_0 = 0
replace CCS_0 = 1 if `string'

forval j = 1/`post_window' {
	local i = `j'-1 
	gen CCS_`j' = 0
	replace CCS_`j' = 1 if CCS_`i'==1 & F`j'.treat==0
	}	

gen CCS_m1 = CCS_0

forval j = 2/`pre_window' {
	local i = `j'-1 
	gen CCS_m`j' = 0
	replace CCS_m`j' = 1 if CCS_m`i'==1 & L.CCS_m`i'==1
	}

***********************************************************
***  3 - LP-DiD ESTIMATES (Variance-weighted ATT)       ***
***********************************************************

* Compute long differences of the outcomes (left-side variable in the LP-DiD regressions)
forval h = 0/`post_window' {
	qui gen D`h'y = F`h'.Y - L.Y
	}

forval h = 2/`pre_window' {
	qui gen Dm`h'y = L`h'.Y - L.Y
}

gen b_lpdid_vw=.

forval h = 0/`post_window' {
	qui reghdfe D`h'y 	    							 	 ///
				treat    							  	 	 ///   /* treatment indicator */
	if 			(treat==1 | treat==0)  & (CCS_`h'==1),	 	 ///    /* clean control condition */
				absorb(time) vce(cluster unit)				   		/* time indicators */
		
	qui replace b_lpdid_vw = _b[treat] if horizon==`h'
}

forval h = 2/`pre_window' {	
		qui reghdfe Dm`h'y  								    ///
					treat 										///   /* treatment indicator */
		if 			(treat==1 | treat==0)  & (CCS_m`h'==1),		///   /* clean controls condition */
						absorb(time) vce(cluster unit)			  	  /* time dummies */
			
		qui replace b_lpdid_vw = _b[treat] if horizon==-`h'
}
	
replace b_lpdid_vw = 0 if horizon==-1

**********************************************************************
***  4 - Re-weighted LP-DiD ESTIMATES (equally-weighted ATT)       ***
**********************************************************************
gen b_lpdid_rw=.

forval h = 0/`post_window' {
	qui reg D`h'y treat##(i.time) if (treat==1 | treat==0)  & (CCS_`h'==1), vce(cluster unit)
	qui margins r.treat if (treat==1 | treat==0)  & (CCS_`h'==1), vce(unconditional) subpop(treat) noestimcheck
	qui replace b_lpdid_rw = r(table)[1,1] if horizon==`h'
	}
	
forval h = 2/`pre_window' {
	qui reg Dm`h'y treat##(i.time) if (treat==1 | treat==0)  & (CCS_m`h'==1), vce(cluster unit)
	qui margins r.treat if (treat==1 | treat==0)  & (CCS_m`h'==1), vce(unconditional) subpop(treat) noestimcheck
	qui replace b_lpdid_rw = r(table)[1,1] if horizon==-`h'
	}
	
replace b_lpdid_rw = 0 if horizon==-1


* Display results against actual effects
list horizon true_att b_lpdid_vw b_lpdid_rw if horizon!=.	
