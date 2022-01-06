# Uncertain Finance
European options risks by uncertainty theory  
Created on Jan 1 2022  
Author: Carlos A. Grajales  
Professor at Universidad de Antioquia, Medellín, Colombia  
<alexander.grajales@udea.edu.co>


* See working paper:  
Uncertainty and stochastic theories on European options valuation and their delta and vega risks  
Grajales - Medina - Mongrut, January 1, 2022  
Available at SSRN: https://ssrn.com/abstract=3998282

*********************************
Folder: Data
*********************************

### Contents:
uc.txt: call prices under Liu's uncertain stock model  
uc_y.txt: call delta risk under Liu's uncertain stock model  
uc_sigma.txt: call vega risk under Liu's uncertain stock model

sc.txt: call prices under Black-Scholes-Merton  
sc_y.txt: call delta risk under Black-Scholes-Merton  
sc_sigma.txt: call vega risk under Black-Scholes-Merton
 
o_utst.xlsx: Excel file with summary of the outputs and visualization  
extra_images.zip: 13 extra Matlab figures for put options consisting of  
..... prices, delta and vega risk, and the norm metrics - Frobenius, norm1, norminf

*********************************
Folder: PythonCodes & MatlabCodes
*********************************

### Contents:
UT_sharedCode: main file  
#### For coming:
uc_estimation: uncertain call price and risk estimation  
sc_estimation: stochastic call price and risk estimation

sparm_RE: relative equality test: call price and risk under stochastic and  
..... uncertain environments  
osparcity: output relative equality test, sparsity coefficients

normsm: norm metrics: Frobenius, norm1, norminf  
ometric: output metrics  
flagm: flags (indicators where prices and risk are most apart from each  
..... other under stochastic and uncertainty environment)
