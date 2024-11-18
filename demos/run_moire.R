library("moire")

#Specify the file name of the allele data table
sfile <- "/home/apujol/isglobal/manuscripts/importation_relatedness/data/data_filtered_for_study.csv"
#Load data
INPUT_DF <- read.csv(sfile)

#Rename the second column according to the needs of MOIRE
colnames(INPUT_DF)[2] <- "sample_id"

# set MOIRE parameters
dat_filter <- moire::load_long_form_data(INPUT_DF)
burnin <- 1e4
num_samples <- 1e4
pt_chains <- seq(1, .5, length.out = 30)

# run moire
mcmc_results <- moire::run_mcmc(
  dat_filter, is_missing = dat_filter$is_missing,
  verbose = TRUE, burnin = burnin, samples_per_chain = num_samples,
  pt_chains = pt_chains, pt_num_threads = length(pt_chains),
  thin = 10)

#Save results
saveRDS(mcmc_results, "/home/apujol/isglobal/manuscripts/importation_relatedness/data/moire_output.RDS") 

#Sample metrics can be extracted from these commands
eff_coi <- moire::summarize_effective_coi(mcmc_results)
naive_coi <- moire::summarize_coi(mcmc_results)
relatedness <- moire::summarize_relatedness(mcmc_results)
