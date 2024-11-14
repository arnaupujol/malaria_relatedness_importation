library(dcifer)

library(tidyverse)


#Specify the file name of the allele data table
sfile <- "/home/apujol/isglobal/projects/genmoz/data/Pipeline_Results/all_2022_allele_data_global_max_div_50cov_loci_pseudocigar.csv"

#Read allele data table
dsmp <- readDat(sfile, svar = "sampleID", lvar = "locus", avar = "new_allele")

#Specify Dcifer parameters
lrank <- 2

#Estimate complexity of infection and allele frequencies
coi   <- getCOI(dsmp, lrank = lrank)
afreq <- calcAfreq(dsmp, coi, tol = 1e-5) 
str(afreq, list.len = 2)

#Run Dcifer to calculate identity-by-descent
dres0 <- ibdDat(dsmp, coi, afreq, pval = TRUE, confint = TRUE, rnull = 0, 
                alpha = 0.05, nr = 1e3)  

#Save results in files
#IBD resimates
write_csv(data.frame(dres0[,,"estimate"]), file = "/home/apujol/isglobal/projects/genmoz/data/Pipeline_Results/ibd_dcifer/all22_div_global_max_50cov_v2_ibd_results_pseudocigar.csv")
#P-values
write_csv(data.frame(dres0[,,"p_value"]), file = "/home/apujol/isglobal/projects/genmoz/data/Pipeline_Results/ibd_dcifer/all22_div_global_max_50cov_v2_ibd_results_p_value_pseudocigar.csv")
#All results in RDS
saveRDS(dres0, file = "/home/apujol/isglobal/projects/genmoz/data/Pipeline_Results/ibd_dcifer/all22_div_global_max_50cov_v2_ibd_results_pseudocigar.rds")