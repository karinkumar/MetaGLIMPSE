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
###command line interface###

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--dosages", required=True, nargs ='+', help = "Specify one or more imputed files")
parser.add_argument("--gl", required = True, help = "Specify low coverage sequence data genotype likelihoods in phred scale")
parser.add_argument("--out", required = True, help = "Specify prefix for output files")
parser.add_argument("--haploid", required = False, help = "Option for haploid data", dest = 'haploid', action = 'store_true')
parser.add_argument("--nomixedstate", required = False, help = "Option for MetaGLIMPSE-plain",  dest = 'nomixedstate', action = 'store_true')
#parser.add_argument("--pickle", required = False, dest = 'pickle', action = 'store_true')
#parser.add_argument("--inner", required = False, dest = 'outer', action= 'store_false')
parser.add_argument("--chunk_size", required = False, dest = 'L')
parser.add_argument("--recomb_rate", required = False, dest = 'c_start', type = float)
parser.add_argument("--baum_welch", required = False, dest = 'bw', action = "store_true")
parser.add_argument("--viterbi", required = False, dest = 'vbi', action = "store_true")
parser.add_argument("--samedosage", required = False, dest = 'same', action = "store_true")
parser.add_argument("--zerodosage", required = False, dest = 'zero', action = "store_true")
parser.add_argument("--chr", required = True, dest = 'chr')
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
if K < 2 or K > 6: 
    raise ValueError("Must have 2 reference panels to meta impute and cannot meta impute more than 6 panels")


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
    missing = "zero" #default correct 

# %%
start = time.time()
print("Checking vcfs...")
assert check_sample_names(GL, *DS_list)
print("Passed checks .. Chunking vcfs ...")
L=30000
regions = get_region_list(*DS_list, chunk_size = L, CHR = args.chr)

# %%
for num, r in enumerate(regions):
    SNPs, dicto, gl, ad = read_vcfs_genK_region(GL, *DS_list, region = r, outer = True, missing_handling = missing) 
  
    assert ad.size/(K*2*len(dicto)) == len(gl) == len(SNPs) #check file size is consistent indicating markers are lined up
    assert len(np.unique(SNPs))==len(SNPs) #check SNPs are unique
    assert len(dicto) == gl.shape[1] #check sample names are equivalent

    samples = {}
    lmbda, diffs = calcLambda(SNPs, start_c)
    total_distance = np.sum(diffs)
    lda = calcNumFlips(lmbda, len(Hidden))
    #lda = calcNumFlips(calcLambda(SNPs, c_start)[0], len(Hidden)) #do this once and then subset
    for sample in dicto.keys(): 
        mdosages = []
        print("Meta Imputing sample ...", sample, "in region", r)
    #subset data structures
        og_transformed = gl[sample]
        if bw: #baum welch option
            lda_c = update_c(Hidden, og_transformed.size, list(og_transformed), transition_prob, emission_prob, n_iter, total_distance, sample_map(sample), ad, lda, SNPs, start_c)
            pst = fwd_bwd(Hidden, og_transformed.size, list(og_transformed), transition_prob, emission_prob, sample_map(sample), ad, lda_c)
        
        elif vbi:
            print("viterbi used")
            opt_path = viterbi(Hidden, og_transformed.size, list(og_transformed), transition_prob, emission_prob, sample_map(sample), ad, lda)
            mdosages.append(calcViterbiDosage(opt_path, sample_map(sample), ad))
            #print("viterbi used")
        else: #fwdbwd 
    #calculate posteriors
        # #%timeit 
            pst = fwd_bwd(Hidden, og_transformed.size, list(og_transformed), transition_prob, emission_prob, sample_map(sample), ad, lda)
    #cProfile.run('pst = fwd_bwd(Hidden, og_transformed.size, list(og_transformed), transition_prob, emission_prob, sample_map(sample), adc, ldac)')
    #calculate meta dosages
        
            if missing=="flat": 
                #mdosages.append(calcMetaDosages_nan(pst, sample_map(sample), ad))
                raise ValueError("must specifiy either --zerodosage or --samedosage")
            else:
                mdosages.append(calcMetaDosages(pst, sample_map(sample), ad))
                #print("should not be here for viterbi")
    #add to samples
        samples[sample] = list(chain.from_iterable(mdosages))

    write_vcf(samples, SNPs, args.out + str(num), args.chr)

end = time.time ()
print("total time is", end - start)
