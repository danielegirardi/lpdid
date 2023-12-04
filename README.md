# LP-DiD

Dube, Girardi, Jorda' and Taylor (2023) propose a local projections approach to difference-in-differences event studies - LP-DiD.

This repository contains three STATA do files that implement the LP-DiD estimator in simulated datasets. 

The first two examples illustrate the case of binary, staggered and absorbing treatment, when only not yet treated units are used as controls:

- "LP_DiD_examplefile.do" uses a simulated dataset similar to the Montecarlo simulations presented in Dube, Girardi, Jorda' and Taylor (2023). 

- "lpdid_test.do" applies LP-DiD in a simulated dataset from Borusyak (2021).

The third example ("LPDiD_nonabsorbing_example.do") illustrates the case of a binary non-absorbing treatment, meaning that units can enter and exit treatment multiple times.

****** <bf>DEC 2023 UPDATE</bf>: *****
<bf><ul>A canned STATA command implementing LP-DiD (by Daniele Girardi and Alexander Busch) is now available on SSC.</bf></ul> It can be installed from within Stata by typing "ssc install lpdid". You can then type "help lpdid" to read the help file, which explains the syntax and working of the command, and provides examples using simulated datasets. 
More details on the package listing:
https://econpapers.repec.org/software/bocbocode/S459273.htm
Please do not hesitate to reach out if you have any problems using the lpdid command, or if you have any questions or find any bugs, either via email or by replying to the following Statalist discussion
https://www.statalist.org/forums/forum/general-stata-discussion/general/1736005-lpdid-new-module-implementing-local-projections-difference-in-differences
*****************************

Related paper: 
Dube, Girardi, Jorda' and Taylor (2023) "A Local Projections Approach to Difference-in-Differences Event Studies", NBER Working Paper 31184, http://www.nber.org/papers/w31184
