# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.4
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %%
from ForwardProb import forward_prob, forward_prob_scaling
from BackwardProb import backward_prob, backward_prob_scaling
import numpy as np


# %%
def fwd_bwd(Hidden, T, obs, a, b, sampleID, df, dist): 
    fwd = forward_prob(Hidden,T, obs, b, a, sampleID, df, dist)
    #print("fwd", fwd[2000:2020])
    bwd = backward_prob(Hidden,T, obs, a, b, sampleID, df, dist)
    #print("bwd", bwd)
    posterior = []
    unnorm = []
    for i in range(len(fwd)): #fwd and bwd have the same T 
        #create a dictionary for all the states for each time
        #unnormalized values
        unnorm.append({st: fwd[i][st] * bwd[i][st] for s, st in enumerate(Hidden)}) #need to divide by total probability
        #print(unnorm)
        #noramlize
        c = sum(unnorm[i].values())
        #print(c)
        posterior.append({st: fwd[i][st] * bwd[i][st]/c for s,st in enumerate(Hidden)})
    return(posterior)


# %%
def fwd_bwd_scaling(Hidden, T, obs, a, b, sampleID, df, dist): 
    fwd, c = forward_prob_scaling(Hidden,T, obs, b, a, sampleID, df, dist)
    
    bwd = backward_prob_scaling(Hidden,T, obs, a, b, sampleID, df, dist, c)
    posterior = []
    unnorm = []
    for i in range(len(fwd)): #fwd and bwd have the same T 
        #create a dictionary for all the states for each time
        #unnormalized values
        unnorm.append({st: fwd[i][st] * bwd[i][st] for s, st in enumerate(Hidden)}) #need to divide by total probability
        #noramlize
        c_norm = sum(unnorm[i].values())
        posterior.append({st: fwd[i][st] * bwd[i][st]/c_norm for s,st in enumerate(Hidden)})
    return(posterior)

# %% [raw]
# def fwd_bwd(Hidden, T, obs, a, b, sampleID, df, dist):
#     # Assuming forward_prob and backward_prob are predefined functions
#     fwd = forward_prob(Hidden, T, obs, b, a, sampleID, df, dist)
#     #print("fwd", fwd[2000:2020])
#     bwd = backward_prob(Hidden, T, obs, a, b, sampleID, df, dist)
#
#     posterior = []
#     unnorm = []
#     for i in range(len(fwd)): #fwd and bwd have the same T 
#         # Unnormalized values
#         #print("lenfwd", len(fwd))
#         unnorm.append({st: fwd[i][st] * bwd[i][st] for s, st in enumerate(Hidden)})
#
#         # Normalize
#         c = sum(unnorm[i].values())
#         
#         # If no probability mass or normalization factor is zero, this may indicate an issue with the model or inputs.
#         if c == 0:
#             raise ValueError(f"Normalization factor is zero at time {i}, check your model and inputs.")
#
#         posterior_state = {st: fwd[i][st] * bwd[i][st]/c for s, st in enumerate(Hidden)}
#         
#         # Check if posterior probabilities are within the valid range.
#         for state, prob in posterior_state.items():
#             if not (0 <= prob <= 1):
#                 error_msg = (
#                     f"Invalid probability for state {state} at time {i}: {prob}. "
#                     f"fwd: {fwd[i][state]}, bwd: {bwd[i][state]}, posterior: {prob}"
#                 )
#                 raise ValueError(error_msg)
#         
#         posterior.append(posterior_state)
#
#     return posterior


# %% [raw]
# #test function
# from TEProb import emission_prob, transition_prob
# import pickle
# import numpy as np
# Hidden = (((1,1), (1,2)), ((2,1), (2,2)))
# #run TEprob
#
# obs_genotypes = pickle.load(open("ASWGT_testcase.csv", "rb"))
#
#
# og_transformed = obs_genotypes["NA19625"]#[1993:1998]
#
# ad = np.load("230820_ASWallelicdosages_testcase.npy")
# ad = ad[:,:,:, 1994:3000]
#
# dist = np.load("ldatestcase.npy")
# dist = dist[1994:3000,:]

# %% [raw]
# fwd_bwd(Hidden, 3, og_transformed, transition_prob, emission_prob, 0, ad, dist)[1]
