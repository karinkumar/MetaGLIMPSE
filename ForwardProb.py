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
JumpThreshold = 1e-10
JumpFix = 1e10
import numpy as np


# %%
def forward_prob(Hidden, T, observations, emission_prob, a, sampleID, df, dist):
    '''Input: Tuple of Hidden States, T number of markers, observed data, emission matrix, a: transition matrix function
    , sampleID, dataframe
       Output: Table of forward probabilities, each entry in the list (a dictionary) is a point in time going forwards. 
   So 0 is Time(Marker) 0 and M is Time (Marker) M 
    '''
    fwd = []     #T is the number of Markers, #H is a tuple of hidden states
    
    for t in range(T): ##why T for loop before S for loop??? #because must progress forward in time and calculate
        #each marker j column and THEN progress in time because each state depends on every other state in t-1
        #print(t)
        f_curr = {} #dictionary
        for s, names in enumerate(Hidden):
            b = emission_prob(names, observations[t], t, sampleID, df) #emission_prob(hidden, obs, m, sampleID):
            #print("emission", "state:", s, "marker:", t, b)
            if t==0:
                sum_prev = 1 #base case flat prior
            else: 
                #for s_, names_ in enumerate(Hidden):
                    #print("transition", "from", names_, "to:", names, "marker", t, a(names_,names, dist,t))
                sum_prev = sum(f_prev[names_]*a(names_,names, dist,t) for s_, names_ in enumerate(Hidden))

            f_curr[names] = sum_prev*b
            #if f_curr[names] < 0: 
             #   raise ValueError("fwd prob is negative at", t, "emission_prob", b, "sum_prev", sum_prev, "state", names)
            #print("fwd at state", names, f_curr[names], "previous alpha/base case*transition", sum_prev, "emissions", b)
        #check for underflow issues
        #print(f_curr[list(f_curr.keys())[0]])
        if f_curr[list(f_curr.keys())[0]] < JumpThreshold: #check first hidden state only

            for key in f_curr: 
                f_curr[key]*=JumpFix   
            #check for negative values
            
        fwd.append(f_curr)
        f_prev = f_curr
    return(fwd)


# %%
def forward_prob_scaling(Hidden, T, observations, emission_prob, a, sampleID, df, dist):
    '''Input: Tuple of Hidden States, T number of markers, observed data, emission matrix, a: transition matrix function
    , sampleID, dataframe
       Output: Table of forward probabilities, each entry in the list (a dictionary) is a point in time going forwards. 
   So 0 is Time(Marker) 0 and M is Time (Marker) M and the scaling factor for use in the backward probability
    '''
    fwd = []     #T is the number of Markers, #H is a tuple of hidden states
    
    c = np.zeros(T) #initiate scaling factor
    for t in range(T): ##why T for loop before S for loop??? #because must progress forward in time and calculate
        #each marker j column and THEN progress in time because each state depends on every other state in t-1
        #print(t)
        f_curr = {} #dictionary
        
        for s, names in enumerate(Hidden):
            b = emission_prob(names, observations[t], t, sampleID, df) 
            if t==0:
                sum_prev = 1 #base case flat prior
            else: 
                sum_prev = sum(f_prev[names_]*a(names_,names, dist,t) for s_, names_ in enumerate(Hidden))
                
            f_curr[names] = sum_prev*b
            c[t] += f_curr[names]
            
            
       
        f_curr = {state: prob / c[t] for state, prob in f_curr.items()} #scale
        fwd.append(f_curr)
        f_prev = f_curr
        

    return(fwd, 1/c) 

# %%

# %% [raw]
# #test function
# from TEProb import emission_prob, transition_prob
# import pickle
# import numpy as np
# Hidden = (((1,1), (1,2)), ((2,1), (2,2)))
# T= 1
# #run TEprob
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

# %% [raw]
# forward_prob(Hidden, 2, og_transformed, emission_prob, transition_prob, 0, ad, dist)
