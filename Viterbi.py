# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.17.2
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

import numpy as np
JumpThreshold = 1e-10
JumpFix = 1e10 
def viterbi(Hidden, T, observations, a, b, sampleID, df, dist):
    '''Input: Tuple of Hidden States, T number of markers, observed data, emission matrix, a: transition matrix function
       Output: Backpointers representing the optimal path
   So 0 is Time(Marker) 0 and M is Time (Marker) M 
    '''
#SAME as the fwd algorithm but use max instead of sum + backpointers    
    
    JumpThreshold = 1e-10
    JumpFix = 1e10 
    v = []     #T is the number of Markers, #H is a tuple of hidden states
    backpointers = [] #point to best state at previous time
    
    for t in range(T): 
        v_curr = {} #dictionary
        back_curr = {}
        for s, names in enumerate(Hidden):
            #print(type(s), type(observations[t]))
            e = b(names, observations[t], t, sampleID, df) 
            #emission_prob(hidden, obs, m, sampleID):
            #print("emission", "state:", s, "marker:", t, b)
            if t==0:
                max_prev = 1 #base case flat prior
        
            else: 
                #for s_, names_ in enumerate(Hidden):
                
                #print("transition", "from", names_, "to:", names, "marker", t, a(names_,names, dist,t))
                values = np.array([v_prev[names_] * a(names_,names, dist,t) for s_, names_ in enumerate(Hidden)])

                # Get the maximum value and its corresponding index
                max_prev = np.max(values)
                argmax_index = np.argmax(values)
                back_curr[names] = Hidden[argmax_index]

                
            v_curr[names] = max_prev*e
           
         
        if v_curr[list(v_curr.keys())[0]] < JumpThreshold: #check first hidden state only

            for key in v_curr: 
                v_curr[key]*=JumpFix   
            
            
        if t > 0:  # Only record backpointers starting from t = 1
            backpointers.append(back_curr)
        
        v.append(v_curr)
        v_prev = v_curr
    
    # **Backward Pass: Reconstruct Optimal Path**
    # Step 1: Find the most likely final state
    last_v = v[-1]  # Probabilities at the last time step
    final_state = max(last_v, key=last_v.get)  # State with the highest probability
    
    # Step 2: Backtrack to find the full sequence
    optimal_path = [final_state]  # Start with the final state
    for t in range(T - 1, 0, -1):  # From T-1 to 1
        final_state = backpointers[t - 1][final_state]  # Use backpointers (adjust index for t > 0)
        optimal_path.append(final_state)
    # Step 3: Reverse the path to get the correct order
    optimal_path.reverse()    
        
    return(optimal_path)


