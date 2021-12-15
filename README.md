# Uncertainty-Finance
Created on Dec 15 2021

Author: C. Alexander Grajales  
Professor @ University of Antioquia, CO  
<alexander.grajales@udea.edu.co>


European options: uncertainty and stochastic views on prices and risks.  
* See working paper:  
Uncertainty and stochastic theories on European options valuation and their delta and vega risks

***********************
Folder: data
***********************

### Contents:
uc.txt: call prices under Liu's uncertain stock model  
uc_y.txt: call delta risk under Liu's uncertain stock model  
uc_sigma.txt: call vega risk under Liu's uncertain stock model

sc.txt: call prices under Black and Scholes  
sc_y.txt: call delta risk under Black and Scholes  
sc_sigma.txt: call vega risk under Black and Scholes  
io.mat: Matlab file with all the above outputs

******************************
Folder: matlab_code
******************************

### Contents:
UT_sharedCode.m: main file  
uc_estimation.m: uncertain call price and risk estimation  
sc_estimation.m: stochastic call price and risk estimation  
sparm_RE.m: Relative equality test: call price and risk under stochastic and
uncertain environments

osparcity.m: output relative equality test  
normsm.m: Metrics: Frobenius, norm1, norminf  
ometric.m: output metrics  
flagm.m: Flags (indicators where prices and risk are most appart from each
other under stochastic and uncertainty environment)
