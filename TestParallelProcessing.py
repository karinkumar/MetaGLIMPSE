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

# %% [raw]
# ###command line interface###
#
# import argparse
# parser = argparse.ArgumentParser()
# parser.add_argument("--afr", required=True)
# parser.add_argument("--eur", required = True)
# parser.add_argument("--gl", required = True)
# parser.add_argument("--out", required = True)
# parser.add_argument("--haploid", required = False, dest = 'haploid', action = 'store_true')
# parser.add_argument("--nomixedstate", required = False, dest = 'nomixedstate', action = 'store_true')
# parser.add_argument("--pickle", required = False, dest = 'pickle', action = 'store_true')
# args = parser.parse_args()
# parser.set_defaults(haploid=False)
# parser.set_defaults(nomixedstate=False)
# parser.set_defaults(pickle=False)
# K=2 #count number of args
# haploid = args.haploid #False
# mixed_states = not args.nomixedstate #False
# if haploid and mixed_states: #sanity check override user 
#     raise ValueError("Cannot have mixed states for haploid data")
# print("mixed states are...", mixed_states)
#
# GL = args.gl #"/net/fantasia/home/kiranhk/1kg30xEAS/genogvcfs1x.vcf.gz"
#
# DS_afr = args.afr #"/net/fantasia/home/kiranhk/software/GLIMPSE2_for_kiran_kumar/GLIMPSE_ligate/AFR_EASdiploid_chr20_ligated.bcf"
#
# DS_eur = args.eur #"/net/fantasia/home/kiranhk/software/GLIMPSE2_for_kiran_kumar/GLIMPSE_ligate/EUR_EASdiploid_chr20_ligated.bcf"

# %%
###run through notebook interface###
outer = True
haploid = False
mixed_states = False
if haploid and mixed_states: #sanity check override user 
    raise ValueError("Cannot have mixed states for haploid data")


GL ="/net/fantasia/home/kiranhk/1kg30xASW/gl/bcftoolsgenogvcfs1x.vcf.gz"

DS_list = ["/net/fantasia/home/kiranhk/software/GLIMPSE2_for_kiran_kumar/GLIMPSE_ligate/ASWbcftoolsAFRdiploid_1xchr20.vcf.gz", 
           "/net/fantasia/home/kiranhk/software/GLIMPSE2_for_kiran_kumar/GLIMPSE_ligate/ASWbcftoolsEURdiploid_1xchr20.vcf.gz"]
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
    from calcMetaDosages import calcMetaDosages
if haploid: 
    from calcDistMat import extract_int, calcLambda
else: 
    from calcDistMat import extract_int, calcLambda, calcNumFlips
from IO import write_vcf, ds_gt_map, read_vcfs_genK_region, check_sample_names, get_region_list
from HiddenStates import generate_hidden
from BaumWelch import update_c
Hidden = generate_hidden(K, mixed_states, haploid)    
def sample_map(sampleID, dicto):
    return(dicto[sampleID] - 1) #index starts at 0



# %%
print("Checking vcfs...")
assert check_sample_names(GL, *DS_list)
print("Passed checks .. Chunking vcfs ...")
L=30000
regions = get_region_list(*DS_list, chunk_size = L)


# %%
from multiprocessing import Pool
from itertools import chain
import time
start_c = 2e-7
bw = False
n_iter = 10
# Assuming all functions like calcNumFlips, calcLambda, fwd_bwd, calcMetaDosages, sample_map, write_vcf, etc. are defined elsewhere
from multiprocessing import Pool
from itertools import chain
import time

# Assuming all functions like calcNumFlips, calcLambda, fwd_bwd, calcMetaDosages, sample_map, write_vcf, etc. are defined elsewhere

def process_sample(args):
    """Process an individual sample"""
    sample, dicto, gl, ad, Hidden, transition_prob, emission_prob, lda, SNPs, region, bw, total_distance, start_c = args
    mdosages = []
    log_output = []
    log_output.append(f"Meta Imputing sample ... {sample} in region {region}")

    # Subset data structures for this sample
    og_transformed = gl[sample]
    
    # Baum-Welch (bw) option
    if bw:
        lda_c = update_c(Hidden, og_transformed.size, list(og_transformed), transition_prob, emission_prob, n_iter, total_distance, sample_map(sample, dicto), ad, lda, SNPs, start_c)
        pst = fwd_bwd(Hidden, og_transformed.size, list(og_transformed), transition_prob, emission_prob, sample_map(sample, dicto), ad, lda_c)
    else:
        # Calculate posteriors without Baum-Welch
        pst = fwd_bwd(Hidden, og_transformed.size, list(og_transformed), transition_prob, emission_prob, sample_map(sample, dicto), ad, lda)

    # Calculate meta dosages
    mdosages.append(calcMetaDosages(pst, sample_map(sample, dicto), ad))

    # Return the sample dosages and the log output
    return sample, list(chain.from_iterable(mdosages)), log_output


def parallel_process_samples(region_data, num_processes=None):
    """Parallelize processing of samples in a given region"""
    sample_list, dicto, gl, ad, Hidden, transition_prob, emission_prob, lda, SNPs, region, bw, total_distance, start_c = region_data
    
    #print("length is", len(dicto.keys()))
    # Prepare the arguments for each sample
    sample_args = [(sample, dicto, gl, ad, Hidden, transition_prob, emission_prob, lda, SNPs, region, bw, total_distance, start_c)
                   for sample in sample_list]

    # Create a pool of worker processes
    with Pool(processes=num_processes) as pool:
        # Process each sample in parallel
        results = pool.map(process_sample, sample_args)
    
    # Collect results into a dictionary
    samples = {sample: dosages for sample, dosages, log_output in results}

    # Collect and print all log messages after processing each region
    for _, _, log_output in results:
        for line in log_output:
            print(line)

    return samples


def main():
    start = time.time()

    # Process regions **sequentially** to maintain order
    for num, r in enumerate(reversed(regions)):
        print(f"Processing region {r} (#{num+1} out of {len(regions)})")
        
        # Read the data for the current region
        SNPs, dicto, gl, ad = read_vcfs_genK_region(GL, *DS_list, region=r, outer=True)
        
        # Ensure data consistency
        assert ad.size / (K * 2 * len(dicto)) == len(gl) == len(SNPs)  # Check file size is consistent
        assert len(np.unique(SNPs)) == len(SNPs)  # Check SNPs are unique
        assert len(dicto) == gl.shape[1]  # Check sample names are equivalent

        # Calculate lambda and num flips
        lmbda, diffs = calcLambda(SNPs, start_c)
        total_distance = np.sum(diffs)
        lda = calcNumFlips(lmbda, len(Hidden))  # Perform this once per region

        # List of samples to process in parallel
        sample_list = list(dicto.keys())

        # Ensure all 61 samples are dispatched
        print(f"Dispatching {len(sample_list)} samples for region {r}")
        
        # Prepare data for parallel processing
        region_data = (sample_list, dicto, gl, ad, Hidden, transition_prob, emission_prob, lda, SNPs, r, bw, total_distance, start_c)

        # Parallelize the processing of samples in this region
        samples = parallel_process_samples(region_data, num_processes=None)  # Adjust `num_processes` as needed

        # Write VCF for the region
        write_vcf(samples, SNPs, "results/chunks/parallelsamples" + str(r))

    end = time.time()
    print("Total time is", end - start)



# %%
if __name__ == "__main__":
    main()

# %% [raw]
# from multiprocessing import Pool
#
# # Assuming calcNumFlips, calcLambda, fwd_bwd, calcMetaDosages, sample_map, and other necessary functions are defined elsewhere
#
# def process_region(r):
#     output = []
#    
#     # Iterate over chunks for the given sample
#     SNPs, dicto, gl, ad = read_vcfs_genK_region(GL, *DS_list, region = r, outer = True) 
#     assert ad.size/(K*2*len(dicto)) == len(gl) == len(SNPs) #check file size is consistent indicating markers are lined up
#     assert len(np.unique(SNPs))==len(SNPs) #check SNPs are unique
#     assert len(dicto) == gl.shape[1] #check sample names are equivalent
#
#     if dicto is None: 
#         dicto = dicto
#
#     samples = {}
#  
#     lda = calcNumFlips(calcLambda(SNPs)[0], len(Hidden)) #do this once and then subset
#
#     for sample in dicto.keys(): 
#         mdosages = []
#         output.append(f"Meta Imputing sample {sample} in region {r}")
#     #subset data structures
#         og_transformed = gl[sample]
#
#     #calculate posteriors 
#         pst = fwd_bwd(Hidden, og_transformed.size, list(og_transformed), transition_prob, emission_prob, sample_map(sample, dicto), ad, lda)
#     #calculate meta dosages
#         mdosages.append(calcMetaDosages(pst, sample_map(sample, dicto), ad))
#     #add to samples
#         samples[sample] = list(chain.from_iterable(mdosages))
#    
#     #write vcf
#     write_vcf(samples, SNPs, "results/paralleltest" + str(num))
#
#     
#     # Return the flattened list of dosages for the sample
#     return '\n'.join(output + '\t')
#
#
# def parallel_process_regions(region_list, num_processes=None):
#     # Create a pool of worker processes
#     with Pool(processes=num_processes) as pool:
#         # Process each sample in parallel
#         results = pool.map(process_region, region_list)
#         for result in results:
#             print(result)
#
# # Here, processed_samples will be a dictionary with the same keys as `samples`
# # but with the processed data as values.

# %% [raw]
# parallel_process_regions(regions, num_processes = 16)
