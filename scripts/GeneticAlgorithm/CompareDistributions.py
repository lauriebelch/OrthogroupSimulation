# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 13:12:42 2026

@author: Biol0216


 KL DIvergence of distributions....
 
 ## get X Y data  (likely just X) 
 Bin into N bins - compare bins/KL divergence..
 #
 
 Freedman–Diaconis rule for no. of bins
 KL for distribution of bins
 """

import numpy
import sys
import scipy
import math
import pandas
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import entropy

##### inputs a numpy array/list for distribution 1 and distribution 2 as input - cna be different lengths
##### output a kl divergence 
# So divergence([1,2,3,4],[1,2,3,4,5])
def divergence(dst1,dst2):
    if type(dst1) == list:
        dst1 = numpy.array(dst1)
    if type(dst2) == list:
        dst2 = numpy.array(dst2)        
    cmbdst =  numpy.concatenate((dst1,dst2))
    Bins = numpy.histogram_bin_edges(cmbdst, bins='fd')
    dist1counts = []
    dist2counts = []    
    for b1,b2 in zip(Bins,Bins[1:]):
        mask1dist1 =  numpy.abs(dst1) <= b2
        dst1_ = dst1[mask1dist1]
        mask2dist1 = b1 <= numpy.abs(dst1_) 
        dst1_res = dst1_[mask2dist1]
        mask1dist2 =  numpy.abs(dst2) <= b2 
        dst2_ = dst2[mask1dist2]
        mask2dist2 = b1 <= numpy.abs(dst2_) 
        dst2_res = dst2_[mask2dist2]      
        ################## needs a better rule
        if len(dst1_res) == 0:
            dist1counts.append(1)
        else:
            dist1counts.append(len(dst1_res))            
        if len(dst2_res) == 0:
            dist2counts.append(1)
        else:            
            dist2counts.append(len(dst2_res))       

        

    kl_res = entropy(dist1counts, dist2counts, base=2)
    if kl_res == numpy.inf:
        return 1
    else:
       return kl_res

    


#################### test ##########################################################################################
def alt_test():
    mu1, sigma1 = 1, 0.1
    mu2, sigma2 = 2, 0.1
    mu3, sigma3 = 1, 0.1
    mu4, sigma4 = 2, 0.1   
    
    dist_1 = numpy.concatenate((numpy.random.normal(mu1, sigma1, 100),numpy.random.normal(mu2, sigma2, 100) ))
    dist_2 = numpy.concatenate((numpy.random.normal(mu3, sigma3, 100),numpy.random.normal(mu4, sigma4, 100) ))   
    
    K_L_Divergence = divergence(dist_1,dist_2)
    
    data = {
      "dist_1": dist_1,
      "dist_2": dist_2
    }
    df = pandas.DataFrame(data)
    
    x1 = df['dist_1']
    x2 = df['dist_2']
    kwargs = dict(hist_kws={'alpha':.6}, kde_kws={'linewidth':2})        
    plt.figure(figsize=(10,7), dpi= 80)
    sns.distplot(x1, color="dodgerblue", label="dist_1", **kwargs)
    sns.distplot(x2, color="orange", label="dist_2", **kwargs)
    plt.suptitle(str(K_L_Divergence), size=16)
    plt.legend();
    plt.show()





def test():  
    mu, sigma = 1, 0.1 
    dist_1 = numpy.random.normal(mu, sigma, 100)
    dist_2 = numpy.random.normal(mu, sigma  + 0.01, 100)
    
    KL_result = []
    Dist = []
    for i in range(0,200):
        for rep in range(0,10):
            dist_1 = numpy.random.normal(1, sigma, 100)
            dist_2 = numpy.random.normal(i/100, sigma, 100)
            K_L_Divergence = divergence(dist_1,dist_2) 
            KL_result.append(K_L_Divergence)
            Dist.append(i)
    plt.scatter(Dist,KL_result)
    plt.show()
    
    
    #################################
    
    data = {
      "dist_1": dist_1,
      "dist_2": dist_2
    }
    df = pandas.DataFrame(data)
    
    x1 = df['dist_1']
    x2 = df['dist_2']
    
    kwargs = dict(hist_kws={'alpha':.6}, kde_kws={'linewidth':2})
    K_L_Divergence = divergence(dist_1,dist_2)
    plt.figure(figsize=(10,7), dpi= 80)
    sns.distplot(x1, color="dodgerblue", label="dist_1", **kwargs)
    sns.distplot(x2, color="orange", label="dist_2", **kwargs)
    plt.suptitle(str(K_L_Divergence), size=16)
    plt.legend();
































