# MetaGLIMPSE
Meta Imputation of Low Coverage Sequencing

***Overview***

This method takes the results from two or more single panel GLIMPSE2 imputations and combines the output using weights estimated via HMM tailored for each individual and marker. 

The output of this method is a vcf file at the union set of markers in each input file with the estimated genotype dosage for each sample and marker.

We have shown this method is better than the best single panel imputation in the cases where imputation with the mega panel is better than the other single panel imputations (which is most cases). See pre-print for more information: https://www.biorxiv.org/content/10.1101/2025.06.24.660721v1

***1. Installation***

git clone https://github.com/karinkumar/MetaGLIMPSE2.git

cd MetaGLIMPSE2/

Once you enter the MetaGLIMPSE2 folder, the executable is RunMetaGLIMPSE.py. The following options are required:


-- dosages:  paths of imputed genotypes vcf files (output from GLIMPSE2), 

-- gl:  vcf file with genotype likelihoods for each position in the _union_ set of markers of the dosage files, 

-- out:  prefix of outfiles. 

***2. Run Example***

See the example folder for African American input files derived from 1000 Genomes and downsampled to 1x and run the following code once you have installed the program and also have access to python. 

python3.8 RunMetaGLIMPSE.py --dosages ASWbcftoolsEURdiploid_1xchr20.vcf.gz ASWbcftoolsAFRdiploid_1xchr20.vcf.gz --gl bcftoolsgenogvcfs1x.vcf.gz --out ASWchr20

**2.1 Notes for Running GLIMPSE2**

To run GLIMPSE2 please check out the GLIMPSE tutorial https://odelaneau.github.io/GLIMPSE/
And please download the ap option branch (as this provides the necessary inputs for MetaGLIMPSE) (https://github.com/odelaneau/GLIMPSE/tree/ap-field)
Also use the following options --main 1 --burnin 19 --keep-monomorphic
These are *essential* in obtaining the correct results in both this branch of GLIMPSE2 and in MetaGLIMPSE

***3 Ligate*** 

MetaGLIMPSE2 produces meta-imputed chunks. In order to be turned into one vcf file for an entire chromosome, they need to be ligated used bcftools. The following code ligates the chunks in the example. 

ls -v ASWchr20*.vcf.gz > list.txt

bcftools concat -f list.txt -Oz -o fullASWchr20.vcf.gz

bcftools index fullASWchr20.vcf.gz



