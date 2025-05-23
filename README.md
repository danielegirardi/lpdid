# LP-DiD

Dube, Girardi, Jorda' and Taylor (https://www.nber.org/papers/w31184) propose a local projections approach to estimating difference-in-differences designs - LP-DiD.

LP-DiD STATA COMMAND

A STATA command implementing LP-DiD (lpdid, by Daniele Girardi and Alexander Busch) is available and can be installed from within Stata by typing 
```
ssc install lpdid 
```
You can then type 
```
help lpdid 
```
to read the help file, which explains the syntax and working of the command, and provides examples using simulated datasets. More details can be found on the package listing (https://econpapers.repec.org/software/bocbocode/S459273.htm), and you can ask questions (or signal possible problems/bugs) by sending us an email or by replying to the corresponding Statalist discussion (https://www.statalist.org/forums/forum/general-stata-discussion/general/1736005-lpdid-new-module-implementing-local-projections-difference-in-differences)


EXAMPLE CODES FOR "MANUAL IMPLEMENTATION" IN STATA

While the lpdid STATA command is available, it is also easy to implement the LP-DiD estimator "manually", in the sense of writing your own STATA code for implementing LP-DiD. You might want to do this either because your application requires some bespoke adjustment, or to make sure to understand how the sausage is made.

This repository contains four STATA do files that implement the LP-DiD estimator in simulated datasets. 

Two examples illustrate the case of binary, staggered and absorbing treatment, when only not yet treated units are used as controls:

- "LP_DiD_examplefile.do" uses a simulated dataset similar to the Montecarlo simulations presented in Dube, Girardi, Jorda' and Taylor (2023). 

- "lpdid_test.do" applies LP-DiD in a simulated dataset from Borusyak (2021).

The other two examples illustrate the case of binary non-absorbing treatment, meaning that units can enter and exit treatment multiple times:

- "LPDiD_nonabsorbing_example.do" illustrates the setting where the treatment variable remains turned on after a treatment event, and only turns off if and when there is a 'treatment reversal'. This is the case assumed in Dube, Girardi, Jorda' and Taylor (2023), and it corresponds to applications like "democracy and growth".

- "lpdid_oneoff_nonabs_example.do" deals with binary, non-absorbing and "one-off" treatments. Here "one-off" means that the treatment lasts only for 1 period by construction, although its effects can still be dynamic and persistent. Formally, we have D_{it}=1 if unit i experiences an event at time t, and D_{it}=0 in all other periods. Examples of this type of "one-off" settings are for example minimum wage increases or hurricanes.


***

Related paper: 
Dube, Girardi, Jorda' and Taylor (2023) "A Local Projections Approach to Difference-in-Differences Event Studies", NBER Working Paper 31184, http://www.nber.org/papers/w31184
