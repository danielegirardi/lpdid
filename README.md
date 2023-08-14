# LP-DiD

Dube, Girardi, Jorda' and Taylor (2023) propose a local projections approach to difference-in-differences event studies - LP-DiD.

This repository contains two STATA do files that implement the LP-DiD estimator in simulated datasets. 

The first two examples illustrate the case of binary, staggered and absorbing treatment, when only not yet treated units are used as controls:

- "LP_DiD_examplefile.do" uses a simulated dataset similar to the Montecarlo simulations presented in Dube, Girardi, Jorda' and Taylor (2023). 

- "lpdid_test.do" applies LP-DiD in a simulated dataset from Borusyak (2021).

The third example ("LPDiD_nonabsorbing_example.do") illustrates the case of a binary non-absorbing treatment, meaning that units can enter and exit treatment multiple times.

Related paper: 
Dube, Girardi, Jorda' and Taylor (2023) "A Local Projections Approach to Difference-in-Differences Event Studies", NBER Working Paper 31184, http://www.nber.org/papers/w31184
