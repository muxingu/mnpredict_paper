########## Check README before executing this script

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

