# MetaGLIMPSE
Meta Imputation of Low Coverage Sequencing

***Overview***

This method takes the results from two or more single panel GLIMPSE2 imputations and combines the output using weights estimated via HMM for each individual and marker. 

The output of this method is a vcf file at the union set of markers in each input file with the estimated genotype dosage for each sample and marker.

See pre-print for more information: https://www.biorxiv.org/content/10.1101/2025.06.24.660721v1

***1. Installation***

git clone https://github.com/karinkumar/MetaGLIMPSE2.git

cd MetaGLIMPSE/

Once you enter the MetaGLIMPSE folder, the executable is RunMetaGLIMPSE.py. The following options are required:


-- dosages:  paths of imputed genotypes vcf files (output from GLIMPSE2), 

-- gl:  vcf file with genotype likelihoods for each position in the _union_ set of markers of the dosage files, 

-- out:  prefix of outfiles. 

***2. Run Example***

See the example folder for African American input files derived from 1000 Genomes and downsampled to 1x and run the following code once you have installed the program and also have access to python. 
```

python3.8 RunMetaGLIMPSE.py --dosages ASWbcftoolsEURdiploid_1xchr20.vcf.gz ASWbcftoolsAFRdiploid_1xchr20.vcf.gz --gl bcftoolsgenogvcfs1x.vcf.gz --out ASWchr20
```


***3 Ligate*** 

MetaGLIMPSE produces meta-imputed chunks. In order to be turned into one vcf file for an entire chromosome, they need to be ligated used bcftools. The following code ligates the chunks in the example. 
```
ls -v ASWchr20*.vcf.gz > list.txt

bcftools concat -f list.txt -Oz -o fullASWchr20.vcf.gz

bcftools index fullASWchr20.vcf.gz

```
**4 Running MetaGLIMPSE with your own data**
Once you've run the example, and it works. You will need to have a genotype likelihoods file for your data that calls the UNION set of bi-allelic SNPs between the reference panels you want to meta-impute with. You will also need to remove any SNPs that are in all reference panels are bi-allelic in those reference panels, but have different minor alleles. Also make sure to remove any singletons AC=1 or variants that have AC=0 from you reference panels.

Now you are ready to run GLIMPSE2 with each of your reference panels. Please download the ap option branch (as this provides the necessary inputs for MetaGLIMPSE) (https://github.com/odelaneau/GLIMPSE/tree/ap-field). Once you have chunked your code with GLIMPSE and binarized the reference panel (see tutorial here). You can run a particular chunk as follows: 

    LINE=$(sed -n "${{NUM}}p" YOUR_CHUNK_FILE.txt
    #stuff
    printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
    IRG=$(echo $LINE | cut -d" " -f3)
    ORG=$(echo $LINE | cut -d" " -f4)
    CHR=$(echo ${{LINE}} | cut -d" " -f2)
    REGS=$(echo ${{IRG}} | cut -d":" -f 2 | cut -d"-" -f1)
    REGE=$(echo ${{IRG}} | cut -d":" -f 2 | cut -d"-" -f2)

    # Execute GLIMPSE2_phase for the current chunk
    ./phase/bin/GLIMPSE2_phase --input-gl YOUR_GL_FILE.vcf.gz --reference YOUR_REF_PANEL.bin --output YOUR_OUTPUT_NAME_${{CHR}}_${{REGS}}_${\{REGE}}.bcf --out-ap-field --main 1 --burnin 19 --keep-monomorphic

These last four options are *essential* in obtaining the correct results in both this branch of GLIMPSE2 and in MetaGLIMPSE

Once you've run GLIMPSE2. Simply plug the imputed dosages and genotype likelihoods into MetaGLIMPSE as in Step 2. And then ligate the chunks as in Step 3. 
