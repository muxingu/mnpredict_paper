########## Important:
########## To run the code, you need to modify the "root" variable in "directory.txt" to the local directory that contains this file
########## Also you need to download the following files from the UKB return (information will be provided upon UKB's approval)
##########  - ids_geno.txt
##########  - ids_geno_vaf.txt
##########  - ids_pheno.txt
##########  - ids_pileup.txt
##########  - ids_cnv.txt
##########  - ids_dd.txt
########## These file must be put in the root directory
########## R and packages are also required

########## Execute following command in root directory
###### AML model
Rscript ./model/cox_stepwise-forward.r aml1 AML
Rscript ./model/final_param.r aml1
Rscript ./model/cox_final.r genopheno aml1

###### MDS model
Rscript ./model/cox_stepwise-forward.r mds1 MDS
Rscript ./model/final_param.r mds1
Rscript ./model/cox_final.r genopheno mds1

###### MPN model
Rscript ./model/cox_stepwise-forward.r mpn1 MPN
Rscript ./model/final_param.r mpn1
Rscript ./model/cox_final.r genopheno mpn1

