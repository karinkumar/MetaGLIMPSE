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
import numpy as np 


# %%
def calcMetaDosages(posteriors, sample, allelic_dosages):
    '''INPUT: posteriors list of dictionaries, sample number 0 - N, allelic dosages numpy array
       OUTPUT: list of meta genotype dosages one per marker 0 to M
       
       NOTES:multiply across haplotypes that have the SAME posterior weight 1,1 represents reference panel 1 allele 1  
       AND reference panel 2 allele 2 (d_a*0.01 +  d_b*0.01)
    '''
    meta_dosages=list()
    for m in range(len(posteriors)): 
        panels = posteriors[m]
        #print(min(panels), max(panels))
        meta_dosage=0
        for key, value in panels.items():
            a,b = key[0]
            c,d = key[1]
            meta_dosage += allelic_dosages[a-1][b-1][sample][m]*value + allelic_dosages[c-1][d-1][sample][m]*value
        meta_dosages.append(meta_dosage)
        meta_dosage=0 #reset weighted sum
    #print(max(meta_dosages), min(meta_dosages))
    if min(np.round(np.array(meta_dosages),3))>= 0 and max(np.round(np.array(meta_dosages),3)) <= 2: 
        return(meta_dosages)
    else: 
        print(min(np.round(np.array(meta_dosages),3)), max(np.round(np.array(meta_dosages),3)) )
        raise ValueError("Meta Dosages Not Between 0 and 2")


# %%
def calcMetaDosages_nan(posteriors, sample, allelic_dosages):
    '''INPUT: posteriors list of dictionaries, sample number 0 - N, allelic dosages numpy array
       OUTPUT: list of meta genotype dosages one per marker 0 to M
       
       NOTES:multiply across haplotypes that have the SAME posterior weight 1,1 represents reference panel 1 allele 1  
       AND reference panel 2 allele 2 (d_a*0.01 +  d_b*0.01)
    '''
    meta_dosages=list()
    for m in range(len(posteriors)):
        #print(m)
        panels = posteriors[m]
        has_nan = False
        meta_dosage = 0

        # Check for any NaN values first
        for key in panels:
            a, b = key[0]
            c, d = key[1]
            if np.isnan(allelic_dosages[a-1][b-1][sample][m]) or np.isnan(allelic_dosages[c-1][d-1][sample][m]):
                has_nan = True
                break
        
        for key, value in panels.items():
            a, b = key[0]
            c, d = key[1]
            
            # Perform calculations based on whether any NaN was found
            if has_nan:
                if not (np.isnan(allelic_dosages[a-1][b-1][sample][m]) or np.isnan(allelic_dosages[c-1][d-1][sample][m])):
                    meta_dosage += (allelic_dosages[a-1][b-1][sample][m] + allelic_dosages[c-1][d-1][sample][m])
            else:
                meta_dosage += (allelic_dosages[a-1][b-1][sample][m] * value +
                                allelic_dosages[c-1][d-1][sample][m] * value)
        #print(meta_dosage)
        meta_dosages.append(meta_dosage)
        meta_dosage=0 #reset weighted sum
    
    #print(max(meta_dosages), min(meta_dosages))
    if min(np.round(np.array(meta_dosages),3))>= 0 and max(np.round(np.array(meta_dosages),3)) <= 2: 
        return(meta_dosages)
    else: 
        print(min(np.round(np.array(meta_dosages),3)), max(np.round(np.array(meta_dosages),3)) )
        raise ValueError("Meta Dosages Not Between 0 and 2")


# %%
def calcViterbiDosage(opt_path, sample, allelic_dosages):
    meta_dosages=list()
    for m in range(len(opt_path)): 
        allele_1, allele_2 = opt_path[m]
        a,b = allele_1
        c,d = allele_2
        mg = allelic_dosages[a-1][b-1][sample][m] + allelic_dosages[c-1][d-1][sample][m]
        meta_dosages.append(mg)
    if min(np.round(np.array(meta_dosages),3))>= 0 and max(np.round(np.array(meta_dosages),3)) <= 2: 
        return(meta_dosages)
    else: 
        print(min(np.round(np.array(meta_dosages),3)), max(np.round(np.array(meta_dosages),3)) )
        raise ValueError("Meta Dosages Not Between 0 and 2")
    return meta_dosages


# %%
#test case missing values um gpt inspired

posteriors0 = [
    {((1, 1), (1, 2)): 0.333333333333333, ((2, 1), (2, 2)): 0.3333333333, ((1,1), (2,2)): 0.333333333},
]
sample = 0
allelic_dosages = np.array([
    [
        # Reference Panel 1
        [[0.1], [0.6]]  # Alleles dosage for the first marker (e.g., from panel 1)
    ],
    [
        # Reference Panel 2
        [[np.nan], [np.nan]]  # Alleles dosage for the first marker (e.g., from panel 2)
    ]
])
allelic_dosages.shape = (2,2,1,1)
calcMetaDosages_nan(posteriors0, sample, allelic_dosages)

# %%
#test case no missing values
posteriors = [
    {((1, 1), (1, 2)): 0.5, ((2, 1), (2, 2)): 0.5},
]
sample = 0
allelic_dosages = np.array([
    [
        # Reference Panel 1
        [[0.1], [0.6]]  # Alleles dosage for the first marker (e.g., from panel 1)
    ],
    [
        # Reference Panel 2
        [[0.2], [0.2]]  # Alleles dosage for the first marker (e.g., from panel 2)
    ]
])
allelic_dosages.shape = (2,2,1,1)
calcMetaDosages_nan(posteriors0, sample, allelic_dosages)
