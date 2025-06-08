/* 
	This do file applies the LP-DiD estimator (Dube, Girardi, Jorda' and Taylor, 2023, https://www.nber.org/papers/w31184 - DGJT hereafter) in a simulated example dataset 
	
	We illustrate here the Regression Adjustment LP-DiD specification with covariates discussed in DGJT Section 4.1.1. 
	We use the lag of the differenced outcome as the covariate of interest in this example, but the same could be done with other covariates of interest, including additional fixed effects.
	We consider a setting with absorbing treatment, using all not yet treated units as controls. 
	
	Three implementation strategies are compared which result in identical estimates, but with varying computational efficiency
	
    Authors: 
		Daniele Girardi (King's College London), daniele.girardi@kcl.ac.uk
		Alexander Busch (Massachusetts Institute of Technology), abusch@mit.edu
		
		6/7/2025
	
*/

version 17
clear*

set seed 20210522

* Set DGP parameters
local units 2000
local time_periods 50
local shocks_sd 25
local persistence 0.25
local post_window 5

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
***  2 - LPDID VARIABLES   ***
******************************
/*
For simplicity & speed: 
	lagged outcome as control
	only estimate post-events 
	only save point estimate and standard error 
*/

xtset unit time
local preobs = 1
gen t =_n - `preobs'
replace t=. if t>`post_window'

* Gen control variable: lag of outcome 
gen ld1=l.d.Y

* Gen forward changes in outcome, to be used on the left side of the LP-DiD regressions	
forval h = 0/`post_window' {
	qui gen D`h'y = F`h'.Y - L.Y
	}

gen dtreat=d.treat

/* Use the following to apply different weights (without specifying this, regression adjustment applies equal weights to all treated units, yielding an equally-weighted ATT)
gen weight = weight_variable
*/	

	
******************************
***  3 - ESTIMATION        ***
******************************

*** Estimating regression adjustment with teffects 
gen b_teffects = . 
gen se_teffects = . 

timer clear 
timer on 1

forval h = 0/`post_window' {
	qui teffects ra (D`h'y i.time c.ld1) (dtreat) if D.treat==1 | F`h'.treat==0 [pweight=weight], iterate(0) atet vce(cluster unit)
	qui replace b_teffects = r(table)[1,1] if t==`h'
	qui replace se_teffects = r(table)[2,1] if t==`h'
	}
	
timer off 1
timer list // display run-time 

*** Estimating regression adjustment with margins
gen b_margins = . 
gen se_margins = . 

timer clear 
timer on 1

forval h = 0/`post_window' {
	qui reg D`h'y i.dtreat##(i.time c.ld1) /// 
		if D.treat==1 | F`h'.treat==0 [pweight=weight], vce(cluster unit)
	qui margins r.dtreat if D.treat==1 | F`h'.treat==0 [pweight=weight], vce(unconditional) subpop(dtreat)
	qui replace b_margins = r(table)[1,1] if t==`h'
	qui replace se_margins = sqrt(r(V)[1,1]*((e(N)-e(df_m)-1)/(e(N)-1))*((e(N_clust)-1)/e(N_clust))) if t==`h' // standard error requires a degrees of freedom adjustment to be equivalent to teffects 
	}
	
timer off 1
timer list // display run-time 



*** Estimating regression adjustment with kmatch 
gen b_kmatch = . 
gen se_kmatch = . 

timer clear 
timer on 1

forval h = 0/`post_window' {
	qui kmatch ra dtreat (D`h'y = i.time c.ld1)	///
		if D.treat==1 | F`h'.treat==0 [pweight=weight], att vce(cluster unit)
	qui replace b_kmatch = _b[ATT] if t==`h'
	qui replace se_kmatch = _se[ATT] * sqrt((e(N_clust)-1) / e(N_clust)) if t==`h' // standard error requires a degrees of freedom adjustment to be equivalent to teffects 
	}
	
timer off 1
timer list // display run-time 

*** list result 	
order t b_* se_*
list t b_* se_* if t<=`post_window'

 
