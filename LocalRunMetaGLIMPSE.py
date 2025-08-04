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
###run through notebook interface###
outer = True
haploid = False
mixed_states = True
if haploid and mixed_states: #sanity check override user 
    raise ValueError("Cannot have mixed states for haploid data")


#GL = "/net/fantasia/home/kiranhk/1kg30xEUR/gl/"

GL = "/net/fantasia/home/kiranhk/Samoans/gl/bcftoolsgenogvcfs2x.vcf.gz"

DS_list = ["/net/fantasia/home/kiranhk/software/GLIMPSE2_for_kiran_kumar/GLIMPSE_ligate/Samoan_samoanpanel_2xchr20.vcf.gz", 
          "/net/fantasia/home/kiranhk/software/GLIMPSE2_for_kiran_kumar/GLIMPSE_ligate/Samoan_topmednosamoans_2xchr20.vcf.gz"]
K = len(DS_list)
print("mixed states are...", mixed_states, "... with", K, "reference panels", "outer join is..", outer)

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
    from calcMetaDosages import calcMetaDosages, calcViterbiDosage
if haploid: 
    from calcDistMat import extract_int, calcLambda
else: 
    from calcDistMat import extract_int, calcLambda, calcNumFlips
from IO import write_vcf, ds_gt_map, read_vcfs_genK_region, check_sample_names, get_region_list
from HiddenStates import generate_hidden
Hidden = generate_hidden(K, mixed_states, haploid)    
def sample_map(sampleID):
    return(dicto[sampleID] - 1) #index starts at 0
from BaumWelch import update_c
from NelderMead import objective, evaluate
from Viterbi import viterbi
EM=False
Nelder=False
Vbi=True

# %%
print("Checking vcfs...")
assert check_sample_names(GL, *DS_list)
print("Passed checks .. Chunking vcfs ...")
L=30000
regions = get_region_list(*DS_list, chunk_size = L)

# %%
#np.save("/net/fantasia/home/kiranhk/HMM/samoan_test_regions.npy", regions)
regions = np.load("/net/fantasia/home/kiranhk/HMM/samoan_test_regions.npy")

# %%
#start timing
n_iter = 10
start_c = 0
r = regions[0] #r = args.region
chunk_num = 0 
#for num, r in enumerate(regions):
  
     
start = time.time()
SNPs, dicto, gl, ad = read_vcfs_genK_region(GL, *DS_list, region = r, outer = True, missing_handling = "same") 
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
    if EM:
        lda_c = update_c(Hidden, og_transformed.size, list(og_transformed), transition_prob, emission_prob, n_iter, total_distance, sample_map(sample), ad, lda, SNPs, start_c)

    elif Nelder:  
        evaluate(Hidden, og_transformed.size, list(og_transformed), emission_prob, transition_prob, sample_map(sample), ad, SNPs)

    elif Vbi: 
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

write_vcf(samples, SNPs, 'results/chunks/samoan_viterbi_samedosage_test' + str(chunk_num))

#print("total time for", len(dicto.keys()), "samples is", time.time() - start)
end = time.time ()
print("total time is", end - start)

