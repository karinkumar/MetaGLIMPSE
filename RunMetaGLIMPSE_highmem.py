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
###command line interface###

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--dosages", required=True, nargs ='+', help = "Specify one or more imputed files")
parser.add_argument("--gl", required = True)
parser.add_argument("--out", required = True)
parser.add_argument("--haploid", required = False, dest = 'haploid', action = 'store_true')
parser.add_argument("--nomixedstate", required = False, dest = 'nomixedstate', action = 'store_true')
parser.add_argument("--pickle", required = False, dest = 'pickle', action = 'store_true')
parser.add_argument("--inner", required = False, dest = 'outer', action= 'store_false')
args = parser.parse_args()
parser.set_defaults(haploid=False)
parser.set_defaults(nomixedstate=False)
parser.set_defaults(pickle=False)
haploid = args.haploid #False
mixed_states = not args.nomixedstate #False
if haploid and mixed_states: #sanity check override user 
    raise ValueError("Cannot have mixed states for haploid data")
print("mixed states are...", mixed_states)

GL = args.gl #"/net/fantasia/home/kiranhk/1kg30xEAS/genogvcfs1x.vcf.gz"

DS_list = args.dosages
K=len(DS_list)
if K < 2 or K > 6: 
    raise ValueError("Must have 2 reference panels to meta impute and cannot meta impute more than 6 panels")

#DS_afr = args.afr "/net/fantasia/home/kiranhk/software/GLIMPSE2_for_kiran_kumar/GLIMPSE_ligate/AFR_EASdiploid_chr20_ligated.bcf"

#DS_eur = args.eur #"/net/fantasia/home/kiranhk/software/GLIMPSE2_for_kiran_kumar/GLIMPSE_ligate/EUR_EASdiploid_chr20_ligated.bcf"

# %%
import pickle
import pandas as pd
import numpy as np
from itertools import chain
import time
#import cProfile
from PosteriorProb import fwd_bwd
if haploid:
    from MetaMinimac import emission_prob, transition_prob, calcNumFlips, calcMetaDosages
else: 
    from TEProb import transition_prob, emission_prob
    from calcMetaDosages import calcMetaDosages
if haploid: 
    from calcDistMat import extract_int, calcLambda
else: 
    from calcDistMat import extract_int, calcLambda, calcNumFlips
from IO import write_vcf, ds_gt_map, read_vcfs_genK
from HiddenStates import generate_hidden
Hidden = generate_hidden(K, mixed_states, haploid)    
def sample_map(sampleID):
    return(dicto[sampleID] - 1) #index starts at 0


# %%
print("Reading vcfs ...")
 #start timing
start = time.time()
SNPs, dicto, gl, ad = read_vcfs_genK(GL, *DS_list, outer = args.outer)

# %%
print("Checking vcfs...")
assert ad.size/(K*2*len(dicto)) == len(gl) == len(SNPs) #check file size is consistent indicating markers are lined up
assert len(np.unique(SNPs))==len(SNPs) #check SNPs are unique
assert len(dicto) == gl.shape[1] #check sample names are equivalent

M=len(SNPs); L=30000
chunks = [np.arange(0+L*k,min(L*(k+1) + 1, M + 1)) for k in range(0, M//L + 1)]
print("Number of Chunks is ..", len(chunks))
#chunks

# %%
samples = {}
weights = {}

lda = calcNumFlips(calcLambda(SNPs), len(Hidden)) #do this once and then subset

for sample in dicto.keys(): 
    mdosages = []
    weightsc = []
    for c in chunks:
        print("Meta Imputing sample ...", sample, "from", gl.iloc[min(c),:].name, "to", gl.iloc[max(c) - 1,:].name)
        print("Chunk size is...", max(c)-min(c))
    #subset data structures
        og_transformed = gl[min(c):max(c)][sample]
        adc = ad[:, :, :, min(c):max(c)]
        ldac = lda[min(c):max(c),:]
 

    #calculate posteriors
        # #%timeit 
        pst = fwd_bwd(Hidden, og_transformed.size, list(og_transformed), transition_prob, emission_prob, sample_map(sample), adc, ldac)
    #cProfile.run('pst = fwd_bwd(Hidden, og_transformed.size, list(og_transformed), transition_prob, emission_prob, sample_map(sample), adc, ldac)')
    #calculate meta dosages
        #np.save("231016posteriors_test", pst)
        mdosages.append(calcMetaDosages(pst, sample_map(sample), adc))
        weightsc.append(pst)
    #add to samples
    samples[sample] = list(chain.from_iterable(mdosages))
    weights[sample] = list(chain.from_iterable(weightsc))
   


# %%
if args.pickle: #must use pickle to perserve dict
    #pickle.dump(samples, open(args.out + '.p', 'wb')) #in case write_vcf fails 
    pickle.dump(weights, open(args.out + "weights" + '.p', 'wb'))
print("writing out vcf...")
write_vcf(samples, SNPs, args.out)
end = time.time ()
print("total time is", end - start)


