# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.4
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

import numpy as np
import pandas as pd
import math
from ForwardProb import forward_prob, forward_prob_scaling
from BackwardProb import backward_prob, backward_prob_scaling
from TEProb import emission_prob, transition_prob
from calcDistMat import calcNumFlips, calcLambda


def update_c(Hidden, T, obs, a, b, max_iter, total_dist, sampleID, df, dist, SNPs, start_c):
    """
    input: sample ID and coefficients for hmm
    outputs: c the input to the calcDistMat function for calculating lda for the transition function
    """
    i = 0
    converge = False
    c_old = start_c
    likelihood_old = -math.inf
    lda = dist#this is where the initial c comes in, already computed
    while i < max_iter and converge == False:         
        fwd, c_scal = forward_prob_scaling(Hidden,T, obs, b, a, sampleID, df, lda)
        likelihood = -1 * sum(np.log(c_scal))
        if likelihood < likelihood_old: 
            raise ValueError("Likelihood must not be decreasing")
        print("log likelihood is...", likelihood)
        bwd = backward_prob_scaling(Hidden,T, obs, a, b, sampleID, df, lda, c_scal)
    #print(bwd)
    
        eta = np.zeros(shape=(len(Hidden), len(Hidden), T-1))
#i and j are hidden states
        for t in range(T-1):
            for i,sti in enumerate(Hidden): 
                for j, stj in enumerate(Hidden): 
                    #print(j)
                    eta[i,j,t] = fwd[t][sti]*a(sti,stj, lda, t)*bwd[t+1][stj]*b(stj, obs[t + 1], t+1, sampleID, df)
    #sum over t 
        numer = np.sum(eta, axis = 2)
    #print(numer)
    #denominator
        denom = np.sum(numer, axis = 1, keepdims = True)
        
        a_hat = numer/denom
        
        #print(a_hat, np.sum(a_hat), np.sum(np.diag(a_hat)), T-1)
        
        c = ((np.sum(a_hat) - np.sum(np.diag(a_hat))) * (T - 1)) / total_dist

        converge = abs(c - c_old) < 2e-6
        
        #update lda for next run 
        
        lda = calcNumFlips(calcLambda(SNPs, c)[0], len(Hidden))

        print(c, abs(c - c_old), converge)
        c_old = c 
        i += 1
        likelihood_old = likelihood
    #print(c, sampleID)
    return lda

# +
#test update c again
