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

# +
from IO import ds_gt_map, check_sample_names, get_region_list
import numpy as np
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--dosages", required=True, nargs ='+', help = "Specify one or more imputed files")
parser.add_argument("--gl", required = True, help = "Specify low coverage sequence data genotype likelihoods in phred scale")
parser.add_argument("--outname", required = True, help = "Specify prefix for output file")
parser.add_argument("--chunk_size", required = False, dest = 'L')
parser.add_argument("--chr", required = True, dest = 'chr')
args = parser.parse_args()

GL = args.gl 
DS_list = args.dosages
print("Checking vcfs...")
assert check_sample_names(GL, *DS_list)
print("Passed checks .. Chunking vcfs ...")
K = len(DS_list)
if K < 2 or K > 6: 
    raise ValueError("Must have 2 reference panels to meta impute and cannot meta impute more than 6 panels")
L=30000
regions = get_region_list(*DS_list, chunk_size = L, CHR = args.chr)
#write out as text file 


np.savetxt(args.outname, regions, fmt="%s")

# + active=""
# GL = "/net/fantasia/home/kiranhk/Samoans/gl/bcftoolsgenogvcfs2x.vcf.gz"
#
# DS_list = ["/net/fantasia/home/kiranhk/software/GLIMPSE2_for_kiran_kumar/GLIMPSE_ligate/Samoan_samoanpanel_2xchr20.vcf.gz", 
#           "/net/fantasia/home/kiranhk/software/GLIMPSE2_for_kiran_kumar/GLIMPSE_ligate/Samoan_topmednosamoans_2xchr20.vcf.gz"]
# K = len(DS_list)
# #print("mixed states are...", mixed_states, "... with", K, "reference panels", "outer join is..", outer)
