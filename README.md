# ExtremeValueTheory_Rainfall_Lombardy_Italy
This code is part of my Master's Thesis that will be discussed on the 20th of October 2023. The project aims to expand an already-existing extreme risk measure called Marginal Expected Shortfall (MES). 
This statistic is generally used in finance to predict the expected value of a random variable given that the same r.v. or another r.v. has exceeded a high threshold with probability p. 
I add to the well-known version of the paper "Estimation of the marginal expected shortfall: the mean when a related variable is extreme" by Cai et al. (2015) the introduction of covariates in the estimation of the Extreme Value index, 
which is then exploited in the estimation of the MES. The goal is to improve the estimate's accuracy by introducing more information. 
Moreover, I develop a model to predict MES in other locations where no data is available to estimate it by using spatial kriging. This idea is very innovative and very valuable as more local data may better guide the development of 
mitigation and recovery plans against the predicted risks. 
The code is the application of the project's theoretical derivation to recent data on rainfall in Lombardy (Italy). Indeed, due to the climatic crisis, the frequency of extreme natural events, like extreme rainfall, has grown dramatically.
Therefore, their prediction has become essential to mitigate their effect or to repair their damages. However, an issue with environmental data is that their measurement requires expensive tools, like weather stations.
For this reason, my Master's Thesis derives a model to predict extreme events in sites where there are no data available. 

The project requires the use of Extreme Value theory and Spatial Modelling theory. For this reason, a great majority of the code is in the R language, since there is wide availability of packages designed for these purposes. 
However, for data cleaning and manipulation, I preferred Python. Also, part of the plotting is done in Python, since I find the seaborn library practical. 

Note: the code was developed in 6 months and I reassembled it when the analysis was finished. Therefore, I tried to align all the variable names. If I missed something, I apologize. 
