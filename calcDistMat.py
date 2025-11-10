# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.2
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %%
#import re
#from sklearn.metrics import pairwise_distances
import numpy as np
import math
import time
import pandas as pd
#r=2e-7 #recombination
#K=2 #number of reference panels to be combined
#H=math.comb(K*2, 2) #number of hidden states
#H=2
RECOM_MIN=1e-04


# %%
def extract_int(s):
    return str.split(s, ":")[1] 


# %%
def calcLambda(SNPs, r):
    #points = [int(re.search(pattern, s).group(1)) for s in SNPs] 
    #SNPs = SNPs[min(chunk):max(chunk)]
    extract_int_vec = np.frompyfunc(extract_int, 1, 1)
    points = np.array(extract_int_vec(SNPs), dtype=int)
    #points = np.append(0, points)
    #diffs = np.append(1e-6, np.diff(points))
    diffs = np.diff(points)
    #print(diffs[0:20])
    #return(diffs)
    return(np.maximum(1 - np.exp(-r*diffs), np.array(RECOM_MIN)), diffs) #r/H before



# %%
def calcNumFlips(lda, H):
    #print(H)
    arr = np.zeros((lda.size, 3))
    arr[:,0] = (1 - lda)**2 + (2*(lda - lda**2))/H + lda**2/(H**2) #numflips is 0
    arr[:, 1] = ((1 - lda)*lda)/H + (lda**2)/H**2
    arr[:, 2] = (lda**2)/(H**2)
    return(arr)


# %%
def calcNumFlips_phased(lda, H):
    #print(H)
    arr = np.zeros((lda.size, 3))
    arr[:,0] = (1 - lda)**2 + (2*(lda - lda**2))/H + lda**2/(H**2) #numflips is 0
    arr[:, 1] = ((1 - lda)*lda)/H + 2*(lda**2)/H**2
    return(arr)

# %% [raw]
# GL_path = "ASWGT_testcase.csv"
# gl = pd.read_pickle(GL_path)
# #dist_mat = np.load(dist_path)
# #remove semi-colons
# SNPs = gl.index
# diffs = calcLambda(SNPs)
# test = 1 - np.exp(-r*diffs)
# idx = np.where(test < RECOM_MIN)
# assert np.all((np.maximum(1 - np.exp(-r*diffs), RECOM_MIN) == RECOM_MIN)[idx]) == True

# %% [raw]
# SNPs = np.load("/net/fantasia/home/kiranhk/HMM/2301107allASWtestSNPs.npy", allow_pickle = True)
# #lda_truth = np.load("testlda.npy")
# print(len(SNPs))
#
# lda_test = calcLambda(SNPs)
#
# #np.save("2308093dlda", calcNumFlips(lda_test))
#
# len(lda_test)
