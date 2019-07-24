WGCNA_TS_Documentation is description (LaTex) plus code and some step by step implementation.

WGCNA_TS.r is the code that you will actually be using (only contains functions), 
which will import SimData.r for generating CS/AR/... data and TimeSeriesSimilarity.r which will compute the 
similarity between time series. 

To test parameters and so on, you should run WGCNA_TS_test which imports WGCNA_TS.r. This script is only for testing

You have to install the package "rucrdtw" to use fastDTW. Please see the folder "Useful packages" for more details.
