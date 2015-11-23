---
Title: "ReadMe"
Author: "Roger Zoh"
---
# GSA_counts Project
 In this repos, we consider the GSA (Gene set analysis) based on counts data. We propose a semi-parametric model accounting for potential correlation between gene counts. We then compare the power of the proposed test gains that proposed by [Rahmatallah et al (2014)](http://www.biomedcentral.com/1471-2105/15/397#B30). You can also find the supplemental material [here](http://www.biomedcentral.com/content/supplementary/s12859-014-0397-8-s1.pdf). We directly follow their simulation set up.
 
# Mean testing using BNP
We make use of the BNP to test difference between two high dimensonal means. We make use of the idea of clustering.

 + 11/23/15 - Set up simualtion for the case 1 of the RMBPT paper assuming various values of p0 and looking at the case where mu2/sig0 varies for fied values of sig0. I wrote two function staring with Test_based_on_clustering_Power_analysis_wrt_mu2.R. Again all of this is done assuming know covariance matrix!