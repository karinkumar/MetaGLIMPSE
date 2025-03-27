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
import pysam
from pysam import VariantFile
import pandas as pd
import subprocess
from io import StringIO
import pickle
import json
import numpy as np
from functools import reduce
min_value = 1e-10  # Define the minimum value for the elements in GLs same as GLIMPSE2
epsilon = 1e-5


# %%
def ds_gt_map(ds): 
    ds = round(ds)
    if ds==0: 
        return ((0,0))
    if ds==1: 
        return ((0,1))
    if ds==2: 
        return((1,1))
    else: 
        raise ValueError("DS cannot be other than 0,1,2")


# %%
def write_vcf(samples, SNPs, outname, haploid=False): 
#write vcf
    vcfh = pysam.VariantHeader()
# Add a sample named "ahstram" to our VCF header
    for s in samples.keys():
        vcfh.add_sample(s)
# Add a contig (chromosome 20) to our VCF header
    vcfh.add_meta('contig', items=[('ID','chr20')])
    
# Add GT andDS to FORMAT in our VCF header
    vcfh.add_meta('FORMAT', items=[('ID',"GT"), ('Number',1), ('Type','String'),
    ('Description','Best Guess Genotype')])
    vcfh.add_meta('FORMAT', items=[('ID',"DS"), ('Number','A'), ('Type','Float'),
    ('Description','Meta Imputed Genotype Dosage')])
    vcfh.add_line("##FPLOIDY=1")

    vcfh.add_meta('FORMAT', items=[('ID',"GP"), ('Number',3), ('Type','Float'),
    ('Description','Meta Imputed Genotype Dosage')])
    
# Open a file, "example.vcf" for writing, with our "vcfh" header
    vcf = pysam.VariantFile(outname + ".vcf.gz", "w", header=vcfh)

# Create a record at chr1:1000 A/T which failed filter due to "RF"
# The 'start' value is 0-based, 'stop' is 1-based
    for s,val in enumerate(SNPs): 
        ID = str.split(val, ":")
        r=vcf.new_record(contig='chr20', start=int(ID[1]) - 1 , stop=int(ID[1]),
            alleles=[ID[2],ID[3]])
        for sample in samples.keys():
        # Set dosage
            r.samples[sample]['DS'] = round(samples[sample][s], 3) #round to account for underflow clipping and to match GLIMPSE DS. 
            if haploid: 
                r.samples[sample]['GT'] = round(samples[sample][s])
            else: 
                r.samples[sample]['GT'] = ds_gt_map(samples[sample][s])
            #r.samples[sample]['GP'] = (1,0,0)

        # Write this record to our VCF file
        vcf.write(r)

    # Close the VCF file
    vcf.close()


# %%
def get_sample_names(vcf_file):
    """Extract sample names from a VCF/BCF file and return a dictionary."""
    # The command to list sample names
    cmd = ["bcftools", "query", "-l", vcf_file]

    # Execute the command and get the result
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    # Check for errors
    if result.returncode != 0:
        print(f"Error running bcftools: {result.stderr}")
        return None

    # Split the result by lines and form the dictionary
    samples = result.stdout.strip().split("\n")
    sample_dict = {sample: i+1 for i, sample in enumerate(samples)}

    return sample_dict

# Test the function
#vcf_path = "/net/fantasia/home/kiranhk/software/GLIMPSE2_for_kiran_kumar/GLIMPSE_ligate/AFR_chr20_ligated.bcf"
#dicto = get_sample_names(vcf_path)


# %%
def query_bcftools(vcf_file, dicto, query_field):
    """Run bcftools and return the output as a pandas DataFrame."""
    cmd = [
        "bcftools", "query",
        "-f", f"%CHROM:%POS:%REF:%ALT\t[%{query_field}\t]\n", vcf_file]
    #print(cmd)
    # Execute the command and get the result
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    
    # Check for errors
    if result.returncode != 0:
        print(f"Error running bcftools: {result.stderr}")
        return None

    # Convert the result to a DataFrame
    df = pd.read_csv(StringIO(result.stdout), sep="\t", header=None)
    df.columns = ["ID"] + list(dicto) + ["drop"]

    return df


# %%
def query_bcftools_region(vcf_file, dicto, query_field, region):
    """Run bcftools and return the output as a pandas DataFrame."""
    cmd = [
        "bcftools", "query",
        "-f", f"%CHROM:%POS:%REF:%ALT\t[%{query_field}\t]\n", "-r", region, vcf_file]
    #print(cmd)
    # Execute the command and get the result
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    
    # Check for errors
    if result.returncode != 0:
        print(f"Error running bcftools: {result.stderr}")
        return None

    # Convert the result to a DataFrame
    df = pd.read_csv(StringIO(result.stdout), sep="\t", header=None)
    df.columns = ["ID"] + list(dicto) + ["drop"]

    return df


# %%
def unphred(GL_tuple):
    # Create your list comprehension, incorporating the minimum value check
    GLs = [max(pow(10., -i/10.), min_value) for i in GL_tuple]
    
    return tuple(GLs)  # Convert the list to a tuple before returning


# %%
def check_sample_names(GL_path, *DS_paths): 
 #file checks sample names are consistent
    sample_names_list = [get_sample_names(path) for path in [GL_path] + list(DS_paths)]
    return all(names == sample_names_list[0] for names in sample_names_list)
    


# %%
def get_region_list(*DS_paths, chunk_size):
    # Join the file paths into a single string for the command
    files = " ".join(DS_paths)
    cmd = f"bcftools isec -n +1 {files}"

    # Run the command and capture the output
    result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, text=True)
    if result.returncode != 0:
        print(f"Error running bcftools: {result.stderr}")
        return None

    # Split the output into lines (positions)
    file = result.stdout.strip().split('\n')

    # Split the positions into chunks
    pos = [line.split("\t")[1] for line in file]
    
    M = len(pos)
    L = chunk_size
    chunks = [np.arange(0+L*k,min(L*(k+1) + 1, M + 1)) for k in range(0, M//L + 1)]
    print("Number of Chunks is ..", len(chunks))
    print("Chunk size is...", L)
    
    regions = ['chr20' + ':' + pos[c[0]] + '-' + pos[c[len(c) - 2]] for c in chunks] #-2 because -1 gets to the last value in the python range but python only subsets upto the last value in the range
    
    return regions
    


# %%
def read_vcfs_genK_region(GL_path, *DS_paths, region, outer = True, missing_handling = "flat"): 
    '''
    @GL_path is the path for the genotype likelihoods
    @DS_paths are an unspecified number of paths for the dosages of each of the K panels 
    @region is the chromosomal region
    @outer means the union set of variants from the panels
    @missing handling handles variants that are not in all reference panels: 
        flat = emission probability is 1 regardless of hidden/observed state/data, 
        zero = setting missing variants to zero
        same = setting missing values to the same dosage as the non-missing value(only works for 2 panels)
    '''
    #get PL --> unphred ---> pandas df #get DS ---> create allelic dosage structure numpy
    ## K reference panel
    
 
 #use bcftools to extract data
    dicto = get_sample_names(DS_paths[0])
    gl = query_bcftools_region(GL_path, dicto, "PL", region)
    
    
    #get ploidy
    #ploidy = len(json.loads("[" + gl.iloc[1,1] + "]")) - 1 #assuming gl all have the same ploidy
    for entry in gl.iloc[:, 1]:  # Adjust the column index if necessary
        if not pd.isna(entry) and entry != ".":
            try:
            # Attempt to parse the first suitable JSON string found
                ploidy = len(json.loads("[" + entry + "]")) - 1
                break  # Exit the loop after finding the first valid entry
            except Exception as e:
            # Optional: Print error if JSON parsing fails
                print(f"Error parsing JSON from entry '{entry}': {e}")
              # Skip to the next entry if parsing fails
                continue  # Skip to the next entry if parsing fails
    print("ploidy is ...", ploidy)
    

    
    if ploidy ==1:
        dosage_list = [query_bcftools_region(DS_paths[k], dicto, "DS", region) for k in range(len(DS_paths))]
    else: 
        dosage_list = [query_bcftools_region(DS_paths[k], dicto, "AP", region) for k in range(len(DS_paths))]
    
    #embedded function that knows what dicto is, otherwise its undefined
    def sample_map(sampleID):
        return(dicto[sampleID] - 1)     
        
    for df in dosage_list + [gl]: 
        df.pop("drop")
        
    if outer:
        all_dosage = reduce(lambda  left,right: pd.merge(left,right,on="ID",
                                            how='outer'), dosage_list)  
    else:
        all_dosage = reduce(lambda  left,right: pd.merge(left,right,on="ID",
                                            how='inner'), dosage_list) 
    
    
    pos2 = np.array([int(st.split(":")[1])for st in all_dosage.ID])
    all_dosage["pos"] = pos2
    all_dosage["ref"] = np.array([st.split(":")[2] for st in all_dosage.ID])
    all_dosage["alt"] = np.array([st.split(":")[3] for st in all_dosage.ID])
    all_dosage = all_dosage.sort_values(["pos", "ref", "alt"]) #sort values for write to vcf (#NEED to sort of ref:ALT as well)
    all_dosage.drop(["pos", "ref", "alt"], axis = 1, inplace = True)
    #print(all_dosage)
#OR 
 #inner join is imputed and outer - inner is saved to a vcf file which is then (merged?) to the one generated by the output function
#write a function for this and test it    
    
    
#left merge GL 
    all_GL = pd.merge(all_dosage["ID"], gl, how = 'left', on = "ID")
    all_GL.index = all_GL.pop("ID") #should really be merge GL
    all_dosage.index = all_dosage.pop("ID")
    #print(all_GL.describe())

#all dosage str --> tuple 
    if ploidy==2: #for 0 dosage case only applies in outer merge case but still need to convert to tuple
        if missing_handling == "flat":
            all_dosage = all_dosage.applymap(lambda x: (np.nan, np.nan) if x == '.' or pd.isna(x) else tuple(json.loads("[" + x + "]")))
            print("flat emission probability for missing variants")
        elif missing_handling == "zero": 
            all_dosage = all_dosage.applymap(lambda x: (0,0) if x == '.' or pd.isna(x) else tuple(json.loads("[" + x + "]")))
            print("missing dosages are set equal to 0")
        elif missing_handling == "same" and len(DS_paths) == 2: #only for 2 reference panels
            print("missing dosages are set equal to the dosage in the other reference panel")
        
            for col in all_dosage.columns:
                if col.endswith('_x'):
                    paired_col = col.replace('_x', '_y')
        
                    if paired_col in all_dosage.columns: #sanity check
                        mask = all_dosage[col].apply(lambda x: pd.isna(x)) #find missing values
                        all_dosage.loc[mask, col] = all_dosage.loc[mask, paired_col] #replace missing values
                elif col.endswith('_y'):
                    paired_col = col.replace('_y', '_x')
        
                    if paired_col in all_dosage.columns:
                        mask = all_dosage[col].apply(lambda x: pd.isna(x))
                        all_dosage.loc[mask, col] = all_dosage.loc[mask, paired_col]
    #if its not in both reference panels, it shouldn't be here!
            if all_dosage.isnull().values.any(): 
                raise ValueError("Some marker is not in both reference panels")
            all_dosage = all_dosage.applymap(lambda x: tuple(json.loads("[" + x + "]")))
    
#GL tuple, unphred
    if ploidy == 1: 
        all_GL = all_GL.applymap(lambda x: (0,0) if pd.isna(x) else tuple(json.loads("[" + x + "]")))
    elif ploidy == 2: 
        all_GL = all_GL.applymap(lambda x: (0,0,0) if x == '.' or pd.isna(x) else tuple(json.loads("[" + x + "]")))
    #print("unphred finished")
    pos = (all_GL.applymap(len)).apply(any, axis = 1)
    assert all(pos) #check for len = 0 "truthy" means 0 is false and everything else is true
    assert np.all(all_GL.index == all_dosage.index) #check SNPs
    SNPs = all_dosage.index #output 1

    obsGLMerge=all_GL.applymap(unphred)
    
    
    print("obsGL:", len(gl), "left merged GL:", len(all_GL), "merged dosage:",len(all_dosage), "individual dosages", [len(df) for df in dosage_list])
    if outer: 
        assert len(all_dosage) == len(all_GL) and np.all([len(all_dosage) > len(df) for df in dosage_list])
    else:
        assert len(all_dosage) == len(all_GL) and np.all([len(all_dosage) < len(df) for df in dosage_list])

    #create allelic dosages
    M=all_dosage.shape[0] 
    N = len(dicto)
    sample_list = dicto.keys()
    final_dosages = [all_dosage.iloc[:, N*(k-1):N*k] for k in range(1, len(dosage_list) + 1)]
    
    allelic_dosages = np.zeros((len(dosage_list), 2, N, M)) #could change this to be smaller for ploidy 1 
    
    

#loop over samples
    for k, df in enumerate(final_dosages): 
        df.columns = list(dicto.keys())
        for sample in sample_list:
            if ploidy ==1: 
                allelic_dosages[k][0][sample_map(sample)] = df[sample].apply(lambda x: x) 
            elif ploidy == 2: 
                allelic_dosages[k][0][sample_map(sample)] = df[sample].apply(lambda x: x[0]) #(2,1)
                allelic_dosages[k][1][sample_map(sample)] = df[sample].apply(lambda x: x[1])#(2,2)

   #clip values to prevent underflow
    ad_clipped = np.clip(allelic_dosages, epsilon, 1 - epsilon)
    
    return SNPs, dicto, obsGLMerge, ad_clipped

# %% [markdown]
# def read_vcfs_genK(GL_path, *DS_paths, outer = True): 
#     #get PL --> unphred ---> pandas df #get DS ---> create allelic dosage structure numpy
#     ## K reference panel
#     
#  
#  #use bcftools to extract data
#     dicto = get_sample_names(DS_paths[0])
#     gl = query_bcftools(GL_path, dicto, "PL")
#     
#     
#     #get ploidy
#     #ploidy = len(json.loads("[" + gl.iloc[1,1] + "]")) - 1 #assuming gl all have the same ploidy
#     for entry in gl.iloc[:, 1]:  # Adjust the column index if necessary
#         if not pd.isna(entry) and entry != ".":
#             try:
#             # Attempt to parse the first suitable JSON string found
#                 ploidy = len(json.loads("[" + entry + "]")) - 1
#                 break  # Exit the loop after finding the first valid entry
#             except Exception as e:
#             # Optional: Print error if JSON parsing fails
#                 print(f"Error parsing JSON from entry '{entry}': {e}")
#               # Skip to the next entry if parsing fails
#                 continue  # Skip to the next entry if parsing fails
#     print("ploidy is ...", ploidy)
#     
#
#     
#     if ploidy ==1:
#         dosage_list = [query_bcftools(DS_paths[k], dicto, "DS") for k in range(len(DS_paths))]
#     else: 
#         dosage_list = [query_bcftools(DS_paths[k], dicto, "AP") for k in range(len(DS_paths))]
#     
#     #embedded function that knows what dicto is, otherwise its undefined
#     def sample_map(sampleID):
#         return(dicto[sampleID] - 1)     
#         
#     for df in dosage_list + [gl]: 
#         df.pop("drop")
#         
#     if outer:
#         all_dosage = reduce(lambda  left,right: pd.merge(left,right,on="ID",
#                                             how='outer'), dosage_list)  
#     else:
#         all_dosage = reduce(lambda  left,right: pd.merge(left,right,on="ID",
#                                             how='inner'), dosage_list) 
#     
#     
#     pos2 = np.array([int(st.split(":")[1])for st in all_dosage.ID])
#     all_dosage["pos"] = pos2
#     all_dosage["ref"] = np.array([st.split(":")[2] for st in all_dosage.ID])
#     all_dosage["alt"] = np.array([st.split(":")[3] for st in all_dosage.ID])
#     all_dosage = all_dosage.sort_values(["pos", "ref", "alt"]) #sort values for write to vcf (#NEED to sort of ref:ALT as well)
#     all_dosage.drop(["pos", "ref", "alt"], axis = 1, inplace = True)
#     #print(all_dosage)
# #OR 
#  #inner join is imputed and outer - inner is saved to a vcf file which is then (merged?) to the one generated by the output function
# #write a function for this and test it    
#     
#     
# #left merge GL 
#     all_GL = pd.merge(all_dosage["ID"], gl, how = 'left', on = "ID")
#     all_GL.index = all_GL.pop("ID") #should really be merge GL
#     all_dosage.index = all_dosage.pop("ID")
#     #print(all_GL.describe())
#
# #all dosage str --> tuple 
#     if ploidy==2: #for 0 dosage case only applies in outer merge case but still need to convert to tuple
#         all_dosage = all_dosage.applymap(lambda x: (0,0) if x == '.' or pd.isna(x) else tuple(json.loads("[" + x + "]")))
#     
#     
# #GL tuple, unphred
#     if ploidy == 1: 
#         all_GL = all_GL.applymap(lambda x: (0,0) if pd.isna(x) else tuple(json.loads("[" + x + "]")))
#     elif ploidy == 2: 
#         all_GL = all_GL.applymap(lambda x: (0,0,0) if x == '.' or pd.isna(x) else tuple(json.loads("[" + x + "]")))
#     #print("unphred finished")
#     pos = (all_GL.applymap(len)).apply(any, axis = 1)
#     assert all(pos) #check for len = 0 "truthy" means 0 is false and everything else is true
#     assert np.all(all_GL.index == all_dosage.index) #check SNPs
#     SNPs = all_dosage.index #output 1
#
#     obsGLMerge=all_GL.applymap(unphred)
#     
#     
#     print("obsGL:", len(gl), "left merged GL:", len(all_GL), "merged dosage:",len(all_dosage), "individual dosages", [len(df) for df in dosage_list])
#     if outer: 
#         assert len(all_dosage) == len(all_GL) and np.all([len(all_dosage) > len(df) for df in dosage_list])
#     else:
#         assert len(all_dosage) == len(all_GL) and np.all([len(all_dosage) < len(df) for df in dosage_list])
#
#     #create allelic dosages
#     M=all_dosage.shape[0] 
#     N = len(dicto)
#     sample_list = dicto.keys()
#     final_dosages = [all_dosage.iloc[:, N*(k-1):N*k] for k in range(1, len(dosage_list) + 1)]
#     
#     allelic_dosages = np.zeros((len(dosage_list), 2, N, M)) #could change this to be smaller for ploidy 1 
#     
#     
#
# #loop over samples
#     for k, df in enumerate(final_dosages): 
#         df.columns = list(dicto.keys())
#         for sample in sample_list:
#             if ploidy ==1: 
#                 allelic_dosages[k][0][sample_map(sample)] = df[sample].apply(lambda x: x) 
#             elif ploidy == 2: 
#                 allelic_dosages[k][0][sample_map(sample)] = df[sample].apply(lambda x: x[0]) #(2,1)
#                 allelic_dosages[k][1][sample_map(sample)] = df[sample].apply(lambda x: x[1])#(2,2)
#
#    #clip values to prevent underflow
#     ad_clipped = np.clip(allelic_dosages, epsilon, 1 - epsilon)
#     
#     return SNPs, dicto, obsGLMerge, ad_clipped

# %% [raw]
# gl_path = "/net/fantasia/home/kiranhk/1kg30xEUR/gl/bcftoolsgenogvcfs4x.vcf.gz"
# ds1 = "/net/fantasia/home/kiranhk/software/GLIMPSE2_for_kiran_kumar/GLIMPSE_ligate/EUREURBdiploid_1xchr20.vcf.gz"
# ds2 =  "/net/fantasia/home/kiranhk/software/GLIMPSE2_for_kiran_kumar/GLIMPSE_ligate/EUREURAdiploid_1xchr20.vcf.gz" 
#
# SNPs, dicto, gl, ad  = read_vcfs_genK(gl_path, ds1, ds2)
#
#
# #check equivalence
# SNPs_test, dicto_test, gl_test, ad_test = read_vcfs(gl_path, ds1, ds2)
#
# np.all(SNPs_test == SNPs), dicto_test == dicto
#
# np.where(SNPs != SNPs_test)
#
# SNPs[461775], SNPs_test[461775]
#
# SNPs[247], SNPs_test[246]
#
# np.allclose(ad, ad_test)


# %% [raw]
# dicto_test = pickle.load(open("dicto.p", "rb"))
# print(test[1] == dicto_test)
#
# SNPs_test = np.load("230913ASWSNPs.npy", allow_pickle = True)
#
# print(np.all(test[0] == SNPs_test))
#
# #test allelic dosages and obs gl 
# #ad_test = np.load("230918_haploid_ASWallelicdosages.npy")
# #np.allclose(test[3], ad_test)
#
# gl_test = pickle.load(open("ASWGL.csv", "rb"))
# #gl_test.index = test[2].index
#
#
# import pickle
# import numpy as np
# samples = pickle.load(open('ASWmetaimputed_allsamples.p', 'rb'))
# SNP_path = "230813ASWSNPs.npy"
# SNPs = np.load(SNP_path, allow_pickle = True)
#
# write_vcf(samples, SNPs, "example")

# %% [raw]
# def read_vcfs(GL_path, DS_path1, DS_path2, outer = True): 
#     #get PL --> unphred ---> pandas df #get DS ---> create allelic dosage structure numpy
#     
#
#  
#  #use bcftools to extract data
#     dicto = get_sample_names(DS_path1)
#     gl = query_bcftools(GL_path, dicto, "PL")
#     
#     
#     #get ploidy
#     #ploidy = len(json.loads("[" + gl.iloc[1,1] + "]")) - 1 #assuming gl all have the same ploidy
#     for entry in gl.iloc[:, 1]:  # Adjust the column index if necessary
#         if not pd.isna(entry) and entry != ".":
#             try:
#             # Attempt to parse the first suitable JSON string found
#                 ploidy = len(json.loads("[" + entry + "]")) - 1
#                 break  # Exit the loop after finding the first valid entry
#             except Exception as e:
#             # Optional: Print error if JSON parsing fails
#                 print(f"Error parsing JSON from entry '{entry}': {e}")
#               # Skip to the next entry if parsing fails
#                 continue  # Skip to the next entry if parsing fails
#     print("ploidy is ...", ploidy)
#     
#
#     
#     if ploidy ==1:
#         eur_dosage = query_bcftools(DS_path2, dicto, "DS")
#         afr_dosage = query_bcftools(DS_path1, dicto, "DS")
#     else: 
#         eur_dosage = query_bcftools(DS_path2, dicto, "AP")
#         afr_dosage = query_bcftools(DS_path1, dicto, "AP")
#         
#     #embedded function that knows what dicto is, otherwise its undefined
#     def sample_map(sampleID):
#         return(dicto[sampleID] - 1)     
#         
#
#
#     afr_dosage.pop("drop")
#     eur_dosage.pop("drop")
#     gl.pop("drop")
#
# #merge Dosage files
#     #all_dosage = pd.merge(eur_dosage, afr_dosage, on="ID")
# #Need to consider all variants now: 
#     if outer:
#         all_dosage = pd.merge(eur_dosage, afr_dosage, on="ID", how='outer')
#     
#     pos2 = np.array([int(st.split(":")[1])for st in all_dosage.ID])
#     all_dosage["pos"] = pos2
#     all_dosage["ref"] = np.array([st.split(":")[2] for st in all_dosage.ID])
#     all_dosage["alt"] = np.array([st.split(":")[3] for st in all_dosage.ID])
#     all_dosage = all_dosage.sort_values(["pos", "ref", "alt"]) #sort values for write to vcf (#NEED to sort of ref:ALT as well)
#     all_dosage.drop(["pos", "ref", "alt"], axis = 1, inplace = True)
# #OR 
#  #inner join is imputed and outer - inner is saved to a vcf file which is then (merged?) to the one generated by the output function
# #write a function for this and test it    
#     
#     
# #left merge GL 
#     all_GL = pd.merge(all_dosage["ID"], gl, how = 'left', on = "ID")
#     all_GL.index = all_GL.pop("ID") #should really be merge GL
#     all_dosage.index = all_dosage.pop("ID")
#     #print(all_GL.describe())
#
# #all dosage str --> tuple    
#     if ploidy==2: #for 0 dosage case
#         all_dosage = all_dosage.applymap(lambda x: (0,0) if pd.isna(x) else tuple(json.loads("[" + x + "]")))
#     
#     
# #GL tuple, unphred
#     if ploidy == 1: 
#         all_GL = all_GL.applymap(lambda x: (0,0) if pd.isna(x) else tuple(json.loads("[" + x + "]")))
#     elif ploidy == 2: 
#         all_GL = all_GL.applymap(lambda x: (0,0,0) if x == '.' or pd.isna(x) else tuple(json.loads("[" + x + "]")))
#     #print("unphred finished")
#     pos = (all_GL.applymap(len)).apply(any, axis = 1)
#     assert all(pos) #check for len = 0 "truthy" means 0 is false and everything else is true
#     assert np.all(all_GL.index == all_dosage.index) #check SNPs
#     SNPs = all_dosage.index #output 1
#
#     obsGLMerge=all_GL.applymap(unphred)
#
#     print("obsGL:", len(gl), "left merged GL:", len(all_GL), "merged dosage:",len(all_dosage), "african:", len(afr_dosage), "eur:", len(eur_dosage))
#     assert len(all_dosage) == len(all_GL) #and len(afr_dosage) > len(eur_dosage) and len(afr_dosage) < len(all_dosage) and len(eur_dosage) < len(all_dosage) 
#
#     #create allelic dosages
#     M=all_dosage.shape[0] 
#     N = len(dicto)
#
#     edosage = all_dosage.iloc[:,0:N] #2
#     adosage = all_dosage.iloc[:,N:N*2] #1
#
#     assert sum(edosage.index==adosage.index)==M
#     assert sum(edosage.index==all_GL.index)==M #check SNPs are the same in dosages and GLs
#
#     edosage.columns = list(dicto.keys())
#     adosage.columns = list(dicto.keys())
#
#
#     sample_list = dicto.keys()
#
#     allelic_dosages = np.zeros((2, 2, N, M)) #could change this to be smaller for ploidy 1 
#     
#
# #loop over samples
#     for sample in sample_list:
#         if ploidy ==1: 
#             allelic_dosages[1][0][sample_map(sample)] = edosage[sample].apply(lambda x: x) #2
#             allelic_dosages[0][0][sample_map(sample)] = adosage[sample].apply(lambda x: x) #1
#         elif ploidy == 2: 
#             allelic_dosages[1][0][sample_map(sample)] = edosage[sample].apply(lambda x: x[0]) #(2,1)
#             allelic_dosages[0][0][sample_map(sample)] = adosage[sample].apply(lambda x: x[0]) #(1,1)
#             allelic_dosages[1][1][sample_map(sample)] = edosage[sample].apply(lambda x: x[1])#(2,2)
#             allelic_dosages[0][1][sample_map(sample)] = adosage[sample].apply(lambda x: x[1])#(1,2)
#
#    #clip values to prevent underflow
#     ad_clipped = np.clip(allelic_dosages, epsilon, 1 - epsilon)
#             
#     return(SNPs, dicto, obsGLMerge, ad_clipped)
