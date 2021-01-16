# SCBmeanfd
 Simultaneous Confidence Bands for the Mean of Functional Data

Statistical methods for estimating and inferring the mean of functional data. The methods include simultaneous confidence bands (parametric/nonparametric bootstrap + Gaussian kinematic formula), local polynomial fitting,  bandwidth selection by plug-in and cross-validation, goodness-of-fit tests for parametric models, equality tests for two-sample problems, and plotting functions. 

To install the package:
``` 
install.packages("boot") # run this line if you don't have package 'boot' installed yet
install.packages("KernSmooth") # run this line if you don't have package 'KernSmooth' installed yet
install.packages("devtools") # run this line if you don't have package 'devtools' installed yet
library(devtools)  
install_github("ddegras/SCBmeanfd") 
```
