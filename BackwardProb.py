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
JumpThreshold=1e-10
JumpFix=1e10


# %%
def backward_prob(Hidden, T, observations, a, b, sampleID, df, dist):
    '''Input: Hidden states, T number of markers, observed data, a: transition matrix function
    , b: emission matrix function
       Output: Table of backward probabilities 
    '''
   
    
    bwd = []
    #reverse observations to correspond with t: 
    for t in reversed(range(T)): 
        b_curr = {} #save previous backward probability
        #compute emissions probability here to save time
        if t < T-1: #not calculate t-1 
            #T + 1 FOOL
            b_next_value = [b(names_, observations[t+1], t + 1, sampleID, df) for s_, names_ in enumerate(Hidden)]
            #print("marker", t + 1, "emissions for state 1 and 2", b_next_value)
        for s,names in enumerate(Hidden):
            if t==T-1: 
                sum_next = 1 #base case
            else: 
                #for s_, names_ in enumerate(Hidden):
                    #print("from", names_, "to:", names, "marker", t, a(names_,names, dist,t))
                #change this double for loop to a matrix computation in numpy
                sum_next = sum(b_prev[names_]*a(names, names_, dist, t)*b_next_value[s_] for s_,names_ in enumerate(Hidden)) #recursion
    
            b_curr[names] = sum_next
        
        #check for underflow issues
        if b_curr[list(b_curr.keys())[0]] < JumpThreshold: #check first hidden state only
            for key in b_curr: 
                b_curr[key]*=JumpFix
    
        bwd.insert(0,b_curr) #add dictionary to list
        b_prev = b_curr #update for recursion (inside t-loop)
    return(bwd)


# %%
def backward_prob_scaling(Hidden, T, observations, a, b, sampleID, df, dist, c):
    '''Input: Hidden states, T number of markers, observed data, a: transition matrix function
    , b: emission matrix function, sampleID, allelic dosages, lda, c: scaling factors from forward probability
       Output: Table of backward probabilities 
    '''
   
    
    bwd = []
    #reverse observations to correspond with t: 
    for t in reversed(range(T)): 
        b_curr = {} #save previous backward probability
        #compute emissions probability here to save time
        if t < T-1: #not calculate t-1 
            #T + 1 FOOL
            b_next_value = [b(names_, observations[t+1], t + 1, sampleID, df) for s_, names_ in enumerate(Hidden)]
            #print("marker", t + 1, "emissions for state 1 and 2", b_next_value)
        for s,names in enumerate(Hidden):
            if t==T-1: 
                sum_next = 1 #base case
            else: 
                #for s_, names_ in enumerate(Hidden):
                    #print("from", names_, "to:", names, "marker", t, a(names_,names, dist,t))
                #change this double for loop to a matrix computation in numpy
                sum_next = sum(b_prev[names_]*a(names, names_, dist, t)*b_next_value[s_] for s_,names_ in enumerate(Hidden)) #recursion
    
            b_curr[names] = sum_next
        
        b_curr = {state: prob * c[t] for state, prob in b_curr.items()}
        bwd.insert(0,b_curr) #add dictionary to list
        #scale here
        b_prev = b_curr #update for recursion (inside t-loop)
    return(bwd)

# %% [raw]
# from TEProb import emission_prob, transition_prob
# import pickle
# import numpy as np
# Hidden = (((1,1), (1,2)), ((2,1), (2,2)))
#
#
# obs_genotypes = pickle.load(open("ASWGT_testcase.csv", "rb"))
#
#
# og_transformed = obs_genotypes["NA19625"][1993:1998]
#
# ad = np.load("230820_ASWallelicdosages_testcase.npy")
# ad = ad[:,:,:, 1993:1998]
#
# dist = np.load("ldatestcase.npy")
# dist = dist[1993:1998,:]
#

# %% [raw]
# backward_prob(Hidden, 5, og_transformed, transition_prob, emission_prob, 0, ad, dist)
