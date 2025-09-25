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
import numpy as np
parser = argparse.ArgumentParser()
parser.add_argument("--dosages", required=True, nargs ='+', help = "Specify one or more imputed files")
parser.add_argument("--gl", required = True, help = "Specify low coverage sequence data genotype likelihoods in phred scale")
parser.add_argument("--out", required = True, help = "Specify prefix for output files")
parser.add_argument("--haploid", required = False, help = "Option for haploid data", dest = 'haploid', action = 'store_true')
parser.add_argument("--nomixedstate", required = False, help = "Option for MetaGLIMPSE-plain",  dest = 'nomixedstate', action = 'store_true')
#parser.add_argument("--inner", required = False, dest = 'outer', action= 'store_false')
parser.add_argument("--chunks", required = True)
parser.add_argument("--region", required = True, type = int)
parser.add_argument("--recomb_rate", required = False, dest = 'c_start', type = float)
parser.add_argument("--baum_welch", required = False, dest = 'bw', action = "store_true")
parser.add_argument("--viterbi", required = False, dest = 'vbi', action = "store_true")
parser.add_argument("--samedosage", required = False, dest = 'same', action = "store_true")
parser.add_argument("--zerodosage", required = False, dest = 'zero', action = "store_true")
parser.set_defaults(haploid=False)
parser.set_defaults(nomixedstate=False)
parser.set_defaults(pickle=False)
parser.set_defaults(c_start = 2e-7)
parser.set_defaults(bw=False)
parser.set_defaults(same_dosage=False)
args = parser.parse_args()
haploid = args.haploid #False
mixed_states = not args.nomixedstate #False
start_c = args.c_start
bw = args.bw
vbi = args.vbi
idx = args.region
regions = np.genfromtxt(args.chunks, dtype = "str")
#print(start_c)
if haploid and mixed_states: #sanity check override user 
    raise ValueError("Cannot have mixed states for haploid data")
print("mixed states are...", mixed_states)

if vbi and bw: 
    raise ValueError("Cannot use both Baum-Welch and Viterbi algorithms. Must choose 1")

if vbi: 
    print("Viterbi algorithm in use")
    
GL = args.gl 

DS_list = args.dosages
K=len(DS_list)



# %%
import pickle
import pandas as pd
from itertools import chain
import time
#import cProfile
from PosteriorProb import fwd_bwd
if haploid:
    from MetaMinimac import emission_prob, transition_prob, calcNumFlips, calcMetaDosages
else: 
    from TEProb import transition_prob, emission_prob
    from calcMetaDosages import calcMetaDosages, calcMetaDosages_nan, calcViterbiDosage
if haploid: 
    from calcDistMat import extract_int, calcLambda
else: 
    from calcDistMat import extract_int, calcLambda, calcNumFlips
from IO import write_vcf, ds_gt_map, read_vcfs_genK_region, check_sample_names, get_region_list
from HiddenStates import generate_hidden
from Viterbi import viterbi
Hidden = generate_hidden(K, mixed_states, haploid)    
def sample_map(sampleID):
    return(dicto[sampleID] - 1) #index starts at 0
if args.same: 
    missing = "same"
elif args.zero: 
    missing = "zero"
else:
    missing = "flat" #default correct 

# %%
n_iter = 10
start_c = 0
print(idx)
r = regions[idx] 
chunk_num = idx
#print(r, chunk_num)


start = time.time()
SNPs, dicto, gl, ad = read_vcfs_genK_region(GL, *DS_list, region = r, outer = True, missing_handling = missing) 
assert ad.size/(K*2*len(dicto)) == len(gl) == len(SNPs) #check file size is consistent indicating markers are lined up
assert len(np.unique(SNPs))==len(SNPs) #check SNPs are unique
assert len(dicto) == gl.shape[1] #check sample names are equivalent

samples = {}
#weights = {}

lmbda, diffs = calcLambda(SNPs, start_c)
total_distance = np.sum(diffs)

lda = calcNumFlips(lmbda, len(Hidden)) #initial lda for all samples
for sample in dicto.keys(): 
    mdosages = []
    #weightsc = []
    print("Meta Imputing sample ...", sample, "in region", r)
#subset data structures
    og_transformed = gl[sample]

#calculate posteriors
    # #%timeit 
    if bw:
        lda_c = update_c(Hidden, og_transformed.size, list(og_transformed), transition_prob, emission_prob, n_iter, total_distance, sample_map(sample), ad, lda, SNPs, start_c)

    elif vbi: 
        opt_path = viterbi(Hidden, og_transformed.size, list(og_transformed), transition_prob, emission_prob, sample_map(sample), ad, lda)
        mdosages.append(calcViterbiDosage(opt_path, sample_map(sample), ad))
        #print(opt_path)
        #print(mdosages)
    else: #plain fwd-bwd as in paper 1
    #use jumpfix forward backward function--scaling only required for baum welch.
        pst = fwd_bwd(Hidden, og_transformed.size, list(og_transformed), transition_prob, emission_prob, sample_map(sample), ad, lda)
        mdosages.append(calcMetaDosages(pst, sample_map(sample), ad))
#add to samples
    samples[sample] = list(chain.from_iterable(mdosages))

write_vcf(samples, SNPs, args.out)

#print("total time for", len(dicto.keys()), "samples is", time.time() - start)
end = time.time ()
print("total time is", end - start)
