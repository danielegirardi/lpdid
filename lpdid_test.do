/*
	This do file applies the LP-DiD estimator (Dube, Girardi, Jorda' and Taylor, 2023) in a simulated example dataset from Borusyak et al. 2021 
	Retrieved from https://github.com/danielegirardi/lpdid/

	We thank (without implicating) Thomas Jørgensen (University of Copenhagen) for the idea of applying LP-DiD to this example.
	
	Author of the original example file: Kirill Borusyak (UCL), k.borusyak@ucl.ac.uk (28/05/2021) (https://github.com/borusyak/did_imputation/blob/main/five_estimators_example.do).
		
	Modified by Daniele Girardi (UMass Amherst), dgirardi@umass.edu (30/09/2022): Applied LP-DiD on this example dataset

*/

version 17
clear*
set scheme s2color

* install packages
global do_install 0 /*1: install packages (first time)*/
if ($do_install == 1){
	ssc install event_plot
	ssc install ftools
	ssc install reghdfe	
}

// Generate a complete panel of 5000 units observed in 15 periods
clear all
timer clear
set seed 10
global T = 15
global I = 5000

set obs `=$I*$T'
gen i = int((_n-1)/$T )+1 					// unit id
gen t = mod((_n-1),$T )+1					// calendar period
tsset i t

// Randomly generate treatment rollout years uniformly across Ei=10..16 (note that periods t>=16 would not be useful since all units are treated by then)
gen Ei = ceil(runiform()*7)+$T -6 if t==1	// year when unit is first treated
bys i (t): replace Ei = Ei[1]
gen K = t-Ei 								// "relative time", i.e. the number periods since treated (could be missing if never-treated)
gen D = K>=0 & Ei!=. 						// treatment indicator

// Generate the outcome with parallel trends and heterogeneous treatment effects
gen tau = cond(D==1, (t-12.5), 0) 			// heterogeneous treatment effects (in this case vary over calendar periods)
gen eps = rnormal()							// error term
gen Y = i/$I + 3*t + tau*D + eps 			// the outcome (FEs play no role since all methods control for them)

// True ATT: Construct the vector of true average treatment effects by the number of periods since treatment
matrix btrue = J(1,6,.)
matrix max_btrue = J(1,6,.)
matrix min_btrue = J(1,6,.)
matrix v_btrue = J(1,6,.)
matrix colnames btrue = tau0 tau1 tau2 tau3 tau4 tau5
matrix colnames max_btrue = tau0 tau1 tau2 tau3 tau4 tau5
matrix colnames min_btrue = tau0 tau1 tau2 tau3 tau4 tau5
matrix colnames v_btrue = tau0 tau1 tau2 tau3 tau4 tau5
qui forvalues h = 0/5 {
	sum tau if K==`h'
	
	matrix btrue[1,`h'+1]		=r(mean)
	matrix max_btrue[1,`h'+1]	=r(max)
	matrix min_btrue[1,`h'+1]	=r(min)
	matrix v_btrue[1,`h'+1]		=r(Var)
	
}

// Estimation with event-study TWFE 
forvalues l = 0/5 {
	gen L`l'event = K==`l'
}
forvalues l = 1/14 {
	gen F`l'event = K==-`l'
}
drop F1event // normalize K=-1 (and also K=-15) to zero
reghdfe Y F*event L*event, a(i t) cluster(i)
estimates store ols // saving the estimates for later

// Etimation with variance-weighted LP-DiD (Dube, Girardi, Jorda' and Taylor 2022)
gen treat=D

* Gen forward changes in outcome, to be used on the left side of the LP-DiD regressions	
forval j = 0/5 {
	gen D`j'y = F`j'.Y - L.Y
	}

forval j = 2/5 {
	gen Dm`j'y = L`j'.Y - L.Y
	}
	
* run variance-weighted LP-DiD regressions
matrix blpdid1 = J(1,11,.)
matrix blpdid1var = J(1,11,.)
matrix colnames blpdid1    = pre5 pre4 pre3 pre2 pre1 tau0 tau1 tau2 tau3 tau4 tau5
matrix colnames blpdid1var = pre5 pre4 pre3 pre2 pre1 tau0 tau1 tau2 tau3 tau4 tau5
forval j = 0/5 {
	reghdfe D`j'y      							  	 ///
			D.treat    							  	 ///   /* treatment indicator */
	if 		D.treat==1 | F`j'.treat==0,				 ///   /* clean control condition */
			absorb(t) vce(cluster i)					   /* time dummies */
		
		matrix blpdid1[1,`j'+6]    = _b[D.treat] 
		matrix blpdid1var[1,`j'+6] = e(V)[1,1] 
	
		if `j'>1 & `j'<=5 {
			reghdfe Dm`j'y  						///
					D.treat 						///   /* treatment indicator */
			if 		D.treat==1 | treat==0,			///   /* clean control condition */
					absorb(t) vce(cluster i)			  /* time dummies */
		
		local p = 6-`j'
		matrix blpdid1[1,`p']    = _b[D.treat] 
		matrix blpdid1var[1,`p'] = e(V)[1,1]
		}
	}
	matrix blpdid1[1,5] = 0 
	matrix blpdid1var[1,5] = 0
	

// Estimation with LP-DiD reweighted to get the equally-weighted ATT (Dube, Girardi, Jorda' and Taylor, 2022)
* compute the weight assigned to each observation in the variance-weighted version, and take the inverse  of it
forval j = 0/5 {
	gen group_h`j'=.
	replace group_h`j'=t if (D.treat==1|F`j'.treat==0)
	reghdfe D.treat if (D.treat==1|F`j'.treat==0), absorb(t) residuals(num_weights_`j')
	qui replace num_weights_`j'=. if D.treat!=1
	qui egen den_weights_`j' = total(num_weights_`j')
	qui gen weight_`j' = num_weights_`j'/den_weights_`j'
	bysort group_h`j': egen gweight_`j'=max(weight_`j')
	replace weight_`j'=gweight_`j' if weight_`j'==.
	replace weight_`j'=round(weight_`j',0.00000001)
	gen reweight_`j'=1/weight_`j'
	sort i t
}

* run re-weighted LP-DiD regressions
matrix blpdid2 = J(1,11,.)
matrix blpdid2var = J(1,11,.)
matrix colnames blpdid2 = pre5 pre4 pre3 pre2 pre1 tau0 tau1 tau2 tau3 tau4 tau5
matrix colnames blpdid2var = pre5 pre4 pre3 pre2 pre1 tau0 tau1 tau2 tau3 tau4 tau5
forval j = 0/5 {
	reghdfe D`j'y      							 	 ///
			D.treat    							  	 ///   		/* treatment indicator */
			if 		(D.treat==1 | F`j'.treat==0)     ///   		/* clean controls condition */
			[pweight=reweight_`j'],					 ///		/* re-weighting to get equally-weighted ATT */
			absorb(t) vce(cluster i)							/* time dummies */
		
		matrix blpdid2[1,`j'+6] = _b[D.treat] 
		matrix blpdid2var[1,`j'+6] = e(V)[1,1] 
	
		if `j'>1 & `j'<=5 {
			reghdfe Dm`j'y  							///
					D.treat 							///   /* treatment indicator */
					if 		(D.treat==1 | treat==0)		///   /* clean controls condition */
					[pweight=reweight_0],				///	  /* re-weighting to get equally-weighted ATT */
					absorb(t) vce(cluster i)				  /* time dummies */
		
		local p = 6-`j'
		matrix blpdid2[1,`p'] = _b[D.treat]
		matrix blpdid2var[1,`p'] = e(V)[1,1] 
		}
	}
	matrix blpdid2[1,5] = 0 
	matrix blpdid2var[1,5] = 0 
	
// Alternative (and equivalent) way to get the equally-weighted ATT: LP-DiD implemented through regression adjustment (Dube, Girardi, Jorda' and Taylor, 2022)
matrix blpdid3 = J(1,11,.)
matrix colnames blpdid3 = pre5 pre4 pre3 pre2 pre1 tau0 tau1 tau2 tau3 tau4 tau5
gen dtreat=D.treat
qui tab t, gen(tdum)
drop tdum1
forval j = 0/5 {
	teffects ra (D`j'y  tdum*) (dtreat)	if 		D.treat==1 | F`j'.treat==0, atet
	
	matrix blpdid3[1,`j'+6] = e(b)[1,1]
	
	if `j'>1 & `j'<=5 {
		teffects ra (Dm`j'y  tdum*) (dtreat)	if 		D.treat==1 | treat==0, atet
		local p = 6-`j'
		matrix blpdid3[1,`p'] = e(b)[1,1]
		}
	
	}
	matrix blpdid3[1,5] = 0 

// Pre-mean-differenced LP-DiD
matrix blpdid4 = J(1,11,.)
matrix blpdid4var = J(1,11,.)
matrix colnames blpdid4    = pre5 pre4 pre3 pre2 pre1 tau0 tau1 tau2 tau3 tau4 tau5
matrix colnames blpdid4var = pre5 pre4 pre3 pre2 pre1 tau0 tau1 tau2 tau3 tau4 tau5
* gen pmd forward changes in outcome
bysort i (t) : gen cumulative_y = sum(Y)
gen aveLY = L.cumulative_y/(t-1)
forval j = 0/5 {
	gen PMD`j'y = F`j'.Y - aveLY
	}

forval j = 2/5 {
	gen PMDm`j'y = L`j'.Y - aveLY
	}
forval j = 0/5 {
	reghdfe PMD`j'y      							 ///
			D.treat    							  	 ///   /* treatment indicator */
	if 		D.treat==1 | F`j'.treat==0,				 ///   /* clean control condition */
			absorb(t) vce(cluster i)					   /* time dummies */
		
		matrix blpdid4[1,`j'+6]    = _b[D.treat] 
		matrix blpdid4var[1,`j'+6] = e(V)[1,1] 
	
		if `j'>1 & `j'<=5 {
			reghdfe PMDm`j'y  						///
					D.treat 						///   /* treatment indicator */
			if 		D.treat==1 | treat==0,			///   /* clean control condition */
					absorb(t) vce(cluster i)			  /* time dummies */
		
		local p = 6-`j'
		matrix blpdid4[1,`p']    = _b[D.treat] 
		matrix blpdid4var[1,`p'] = e(V)[1,1]
		}
	}
	matrix blpdid4[1,5] = 0 
	matrix blpdid4var[1,5] = 0	
	
// Pre-mean-differenced and re-weighted LP-DiD (Dube, Girardi, Jorda' and Taylor, 2022)
matrix blpdid5 = J(1,11,.)
matrix blpdid5var = J(1,11,.)
matrix colnames blpdid5    = pre5 pre4 pre3 pre2 pre1 tau0 tau1 tau2 tau3 tau4 tau5
matrix colnames blpdid5var = pre5 pre4 pre3 pre2 pre1 tau0 tau1 tau2 tau3 tau4 tau5
forval j = 0/5 {
	reghdfe PMD`j'y      							 ///
			D.treat    							  	 ///   		/* treatment indicator */
			if 		(D.treat==1 | F`j'.treat==0)     ///   		/* clean controls condition */
			[pweight=reweight_`j'],					 ///		/* re-weighting to get equally-weighted ATT */
			absorb(t) vce(cluster i)							/* time dummies */
		
		matrix blpdid5[1,`j'+6] = _b[D.treat] 
		matrix blpdid5var[1,`j'+6] = e(V)[1,1] 
	
		if `j'>1 & `j'<=5 {
			reghdfe PMDm`j'y  							///
					D.treat 							///   /* treatment indicator */
					if 		(D.treat==1 | treat==0)		///   /* clean controls condition */
					[pweight=reweight_0],				///	  /* re-weighting to get equally-weighted ATT */
					absorb(t) vce(cluster i)				  /* time dummies */
		
		local p = 6-`j'
		matrix blpdid5[1,`p'] = _b[D.treat]
		matrix blpdid5var[1,`p'] = e(V)[1,1] 
		}
	}
	matrix blpdid5[1,5] = 0 
	matrix blpdid5var[1,5] = 0 
	
* True equally-weighted ATT	
mat list btrue,         format(%12.4f)
* Min group-specific ATT
mat list min_btrue,     format(%12.4f)
* Max group-specific ATT
mat list max_btrue,     format(%12.4f)
* Estimates from variance-weigthed LP-DiD
mat list blpdid1, 		format(%12.4f)
* Estimated from equally-weighted LP-DiD (weighted regression approach)
mat list blpdid2,       format(%12.4f)
* Estimates from equally-weighted LP-DiD (regression adjustment approach)
mat list blpdid3,       format(%12.4f)
* Estimates from pre-mean-differenced LP-DiD
mat list blpdid4,       format(%12.4f)
* Estimates from pre-mean-differenced and reweighted LP-DiD
mat list blpdid5,       format(%12.4f)
	
// Display the full range of true treatment effects in this example
event_plot btrue# min_btrue# max_btrue#, default_look stub_lag(tau# tau# tau#)together graph_opt(xtitle("Days since the event") ytitle("Treatment effects") xlabel(0(1)5) ///
	legend(on order(1 "Average treatment effect" 2 "Minimum group-specific effect" 3 "Maximum group-specific effect") rows(2) region(style(none))) ///
	title("True effects: min, max and average") name(true_effects))
	
// Combine all plots using the stored estimates
event_plot btrue# blpdid1#blpdid1var blpdid2#blpdid2var  blpdid4#blpdid4var blpdid5#blpdid5var ols, ///
	stub_lag(tau# tau# tau#  tau# tau# L#event) stub_lead(pre# pre# pre# pre# pre# F#event) plottype(scatter) ciplottype(rcap) ///
	together perturb(-0.325(0.13)0.325) trimlead(5) noautolegend ///
	lag_opt1(msymbol(+) color(red)) lag_ci_opt1(color(red)) ///
	lag_opt2(msymbol(oh) color(green)) lag_ci_opt2(color(green)) ///
	lag_opt3(msymbol(dh) color(blue)) lag_ci_opt3(color(blue)) ///
	lag_opt4(msymbol(sh) color(purple)) lag_ci_opt4(color(purple)) ///
	lag_opt5(msymbol(sh) color(brown)) lag_ci_opt5(color(brown)) ///
	lag_opt6(msymbol(sh) color(orange)) lag_ci_opt6(color(orange)) ///
	graph_opt(title("Event study estimators in a simulated panel (300 units, 15 periods)", size(medlarge)) ///
		xtitle("Periods since the event") ytitle("ATE and estimates") xlabel(-5(1)5) ylabel(0(1)3) ///
		legend(order(1 "True equally-weighted ATT" 2 "Variance-weighted LP-DiD" 4 "Re-weighted LP-DiD" 6 "Variance-weighted PMD LP-DiD" 8 "Re-weighted PMD LP-DiD" 10 "Event-study TWFE") rows(3) region(style(none))) ///
		xline(-0.5, lcolor(gs8) lpattern(dash)) yline(0, lcolor(gs8)) graphregion(color(white)) bgcolor(white) ylabel(, angle(horizontal)) name(estimates) ///
	)
	

