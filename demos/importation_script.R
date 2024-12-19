# Load necessary libraries
library(readr)
library(dplyr)

simmetrize_matrix <- function(matrix, diagonal = FALSE) {
  # This function fills the empty values of a matrix from their symmetric values.
  
  # Set NA values to 0
  matrix[is.na(matrix)] <- 0
  
  # Add the matrix to its transpose
  matrix <- matrix + t(matrix)
  
  # If diagonal is FALSE, set diagonal elements to NA
  if (!diagonal) {
    diag(matrix) <- NA
  }
  
  return(matrix)
}

get_relatedness_to_population <- function(ibd_res_meta, ibd_pval_meta, ibd_threshold = 0.2, p_value = 0.05, 
                                           location = 'province', variable_name = 'rel_origin', pop_location = 'province') {
  # This function estimates the fraction of related pairs for each sample with specified populations.
  
  # Define filled IBD results
  ibd_matrix <- as.matrix(ibd_res_meta[, 1:nrow(ibd_res_meta)])
  class(ibd_matrix) <- "numeric"
  ibd_matrix <- simmetrize_matrix(ibd_matrix)
  
  pval_matrix <- as.matrix(ibd_pval_meta[, 1:nrow(ibd_pval_meta)])
  class(pval_matrix) <- "numeric"
  pval_matrix <- simmetrize_matrix(pval_matrix)
  
  # Defining samples to analyze
  sample_list <- unique(ibd_res_meta$sampleID)  # TO DO: make it for only people traveling
  ibd_res_meta[[variable_name]] <- NA
  
  # Loop over samples
  for (sample in sample_list) {
    # Location origin of sample
    target_location <- ibd_res_meta[ibd_res_meta$sampleID == sample, location][1]
    
    # Defining samples from origin location
    target_samples <- (ibd_res_meta$sampleID != sample) & (ibd_res_meta[[pop_location]] == target_location)
    
    # Fraction of pairs with IBD > ibd_threshold
    high_ibd <- ibd_matrix[target_samples, ibd_res_meta$sampleID == sample] >= ibd_threshold
    low_pval <- pval_matrix[target_samples, ibd_pval_meta$sampleID == sample] <= p_value
    rel_fraction <- mean(high_ibd & low_pval, na.rm = TRUE)
    
    # Saving relatedness with origin population
    ibd_res_meta[ibd_res_meta$sampleID == sample, variable_name] <- rel_fraction
    ibd_pval_meta[ibd_pval_meta$sampleID == sample, variable_name] <- rel_fraction
  }
  
  # Ensure the new variable is numeric
  ibd_res_meta[[variable_name]] <- as.numeric(ibd_res_meta[[variable_name]])
  ibd_pval_meta[[variable_name]] <- as.numeric(ibd_pval_meta[[variable_name]])
  
  return(list(ibd_res_meta = ibd_res_meta, ibd_pval_meta = ibd_pval_meta))
}

get_relatedness_origin_travels <- function(ibd_res_meta, ibd_pval_meta, 
                                           ibd_threshold = 0.2, p_value = 0.05, 
                                           travel2 = FALSE) {
  # This function computes the fraction of related pairs for each sample, with its origin population,
  # and destination populations for up to two travels reported.
  
  # Compute relatedness to origin population
  result <- get_relatedness_to_population(ibd_res_meta, ibd_pval_meta, 
                                          ibd_threshold = ibd_threshold, p_value = p_value, 
                                          location = 'province', variable_name = 'rel_origin', 
                                          pop_location = 'province')
  ibd_res_meta <- result$ibd_res_meta
  ibd_pval_meta <- result$ibd_pval_meta
  
  # Compute relatedness to the first destination population
  result <- get_relatedness_to_population(ibd_res_meta, ibd_pval_meta, 
                                          ibd_threshold = ibd_threshold, p_value = p_value, 
                                          location = 'travel_prov', variable_name = 'rel_dest1', 
                                          pop_location = 'province')
  ibd_res_meta <- result$ibd_res_meta
  ibd_pval_meta <- result$ibd_pval_meta
  
  # If second travel is included, compute relatedness to the second destination population
  if (travel2) {
    result <- get_relatedness_to_population(ibd_res_meta, ibd_pval_meta, 
                                            ibd_threshold = ibd_threshold, p_value = p_value, 
                                            location = 'travel_prov2', variable_name = 'rel_dest2', 
                                            pop_location = 'province')
    ibd_res_meta <- result$ibd_res_meta
    ibd_pval_meta <- result$ibd_pval_meta
  }
  
  return(list(ibd_res_meta = ibd_res_meta, ibd_pval_meta = ibd_pval_meta))
}

r_importation_prob <- function(ibd_res_meta, travel_1 = 'travel_prov', travel_2 = NULL, 
                               r_origin = 'rel_origin', r_des1 = 'rel_dest1', r_des2 = 'rel_dest2', 
                               class_name = 'prob_imported') {
  # This function estimates the probability of cases being imported at the individual level 
  # from the fraction of the genetic relatedness of the sample between the origin and the destination 
  # populations, as P(I) = sum_i (r_dest_i) / (sum_i (r_dest_i) + r_origin)
  
  # Initially setting all as unclassified
  ibd_res_meta[[class_name]] <- NA
  
  # If no travel, P(I) = 0
  no_travel_mask <- is.na(ibd_res_meta[[travel_1]])
  if (!is.null(travel_2)) {
    no_travel_mask <- no_travel_mask & is.na(ibd_res_meta[[travel_2]])
  }
  ibd_res_meta[no_travel_mask, class_name] <- 0
  
  # If only travel 1, use only r_des1 and r_origin
  trav1_mask <- !is.na(ibd_res_meta[[travel_1]])
  if (!is.null(travel_2)) {
    trav1_mask <- trav1_mask & is.na(ibd_res_meta[[travel_2]])
  }
  ibd_res_meta[trav1_mask, class_name] <- ibd_res_meta[trav1_mask, r_des1] / 
    (ibd_res_meta[trav1_mask, r_des1] + ibd_res_meta[trav1_mask, r_origin])
  
  if (!is.null(travel_2)) {
    # If only travel 2, use only r_des2 and r_origin
    trav2_mask <- is.na(ibd_res_meta[[travel_1]]) & !is.na(ibd_res_meta[[travel_2]])
    ibd_res_meta[trav2_mask, class_name] <- ibd_res_meta[trav2_mask, r_des2] / 
      (ibd_res_meta[trav2_mask, r_des2] + ibd_res_meta[trav2_mask, r_origin])
    
    # If both travels not null, use r_des1, r_des2, and r_origin
    trav12_mask <- !is.na(ibd_res_meta[[travel_1]]) & !is.na(ibd_res_meta[[travel_2]])
    ibd_res_meta[trav12_mask, class_name] <- (ibd_res_meta[trav12_mask, r_des1] + ibd_res_meta[trav12_mask, r_des2]) / 
      (ibd_res_meta[trav12_mask, r_des1] + ibd_res_meta[trav12_mask, r_des2] + ibd_res_meta[trav12_mask, r_origin])
  }
  
  return(ibd_res_meta)
}



# Metadata
pipeline_results_path <- "/home/apujol/isglobal/projects/genmoz/data/Pipeline_Results/"
all_data_div_covered_samples <- read_csv(paste0(pipeline_results_path, 'all_2022_allele_data_global_max_div_50cov_loci_v2.csv'))

# Select and drop duplicates for metadata
metadata <- all_data_div_covered_samples %>%
  select(sampleID, nida, province, travel, travel_prov, travel_days_728, source) %>%
  distinct()

# Specific pre-processing
metadata$province <- ifelse(metadata$province == 'Maputo Provincia', 'Maputo', metadata$province)
metadata$province <- ifelse(metadata$province == 'Maputo Cidade', 'Maputo City', metadata$province)

metadata$travel_prov <- ifelse(metadata$travel_prov == 'Maputo Provincia', 'Maputo', metadata$travel_prov)
metadata$travel_prov <- ifelse(metadata$travel_prov == 'Maputo Province', 'Maputo', metadata$travel_prov)
metadata$travel_prov <- ifelse(metadata$travel_prov == 'Maputo Cidade', 'Maputo City', metadata$travel_prov)

# Set the path for IBD data
ibd_data_path <- paste0(pipeline_results_path, "ibd_dcifer/")

# Read the IBD results and p-value files
ibd_res_filename <- paste0(ibd_data_path, "all22_div_global_max_50cov_v2_ibd_results_pseudocigar.csv")
ibd_res <- read_csv(ibd_res_filename)
sapply(ibd_res, as.numeric)
ibd_pval_filename <- paste0(ibd_data_path, "all22_div_global_max_50cov_v2_ibd_results_p_value_pseudocigar.csv")
ibd_pval <- read_csv(ibd_pval_filename)
sapply(ibd_pval, as.numeric)

# Set the row names (index) to match the columns in the data frames
rownames(ibd_res) <- colnames(ibd_res)
rownames(ibd_pval) <- colnames(ibd_pval)

# Define the prevalence of malaria in children per province from Inquérito Demográfico de Saúde 2022-2023
tdr_prevalence <- c('Maputo' = 0.003, 
                    'Maputo City' = 0.0004,  # Officially 0.0, but we assign the higher possible prevalence
                    'Inhambane' = 0.158, 
                    'Zambezia' = 0.349,
                    'Sofala' = 0.332, 
                    'Manica' = 0.102, 
                    'Tete' = 0.195, 
                    'Niassa' = 0.338, 
                    'Cabo Delgado' = 0.381, 
                    'Nampula' = 0.547, 
                    'Gaza' = 0.057)

# Merge ibd_res with metadata using 'sampleID' from metadata and rownames as the index for ibd_res
ibd_res$sampleID = rownames(ibd_res)
ibd_res_meta <- merge(ibd_res, metadata, by = "sampleID")
# Remove the Row.names column that was created during merge
#ibd_res_meta <- ibd_res_meta[, -1]
#ibd_res_meta$sampleID <- rownames(ibd_res)

# Merge ibd_pval with metadata using 'sampleID' from metadata and rownames as the index for ibd_pval
ibd_pval$sampleID = rownames(ibd_pval)
ibd_pval_meta <- merge(ibd_pval, metadata, by = "sampleID")

# Remove the Row.names column that was created during merge
#ibd_pval_meta <- ibd_pval_meta[, -1]
#ibd_pval_meta$sampleID <- rownames(ibd_res)

# Define thresholds to define related pairs
ibd_threshold <- 0.1
p_value <- 0.05

# Calculate all relatedness fractions with origin and destination populations
result <- get_relatedness_origin_travels(ibd_res_meta, ibd_pval_meta, 
                                         ibd_threshold = ibd_threshold, 
                                         p_value = p_value)
ibd_res_meta <- result[[1]]
ibd_pval_meta <- result[[2]]

# Filling missing data of relatedness at destination with relatedness at origin,
# so that both probabilities are equal
ibd_res_meta$rel_dest1[is.na(ibd_res_meta$rel_dest1)] <- ibd_res_meta$rel_origin[is.na(ibd_res_meta$rel_dest1)]
ibd_pval_meta$rel_dest1[is.na(ibd_pval_meta$rel_dest1)] <- ibd_pval_meta$rel_origin[is.na(ibd_pval_meta$rel_dest1)]

# Renormalising when relatedness is 0 for both origin and destination
mask2renorm <- (ibd_res_meta$rel_dest1 == 0) & (ibd_res_meta$rel_origin == 0)
ibd_res_meta$rel_origin[mask2renorm] <- 0.00001
ibd_res_meta$rel_dest1[mask2renorm] <- 0.00001

mask2renorm <- (ibd_pval_meta$rel_dest1 == 0) & (ibd_pval_meta$rel_origin == 0)
ibd_pval_meta$rel_origin[mask2renorm] <- 0.00001
ibd_pval_meta$rel_dest1[mask2renorm] <- 0.00001

# Define prevalences of origin and travel destinations
ibd_res_meta$origin_PR <- NA
ibd_res_meta$dest_PR <- NA

for (i in 1:nrow(ibd_res_meta)) {
    if (is.character(ibd_res_meta$travel_prov[i])) {
        ibd_res_meta$dest_PR[i] <- tdr_prevalence[ibd_res_meta$travel_prov[i]]
    }
    if (is.character(ibd_res_meta$province[i])) {
        ibd_res_meta$origin_PR[i] <- tdr_prevalence[ibd_res_meta$province[i]]
    }
}


# Imputing mean travel time when missing
travelled <- ibd_res_meta$travel == 1
mean_travel_time <- mean(ibd_res_meta$travel_days_728, na.rm = TRUE)
ibd_res_meta$travel_days_728[is.na(ibd_res_meta$travel_days_728) & travelled] <- mean_travel_time

# Obtain importation probabilities using imputed travel time
# Defining RxTxPR variables
ibd_res_meta$rel_dest1xTimp <- ibd_res_meta$rel_dest1 * ibd_res_meta$travel_days_728
ibd_res_meta$rel_origxTimp <- ibd_res_meta$rel_origin * (21 - ibd_res_meta$travel_days_728)
ibd_res_meta$rel_origxTimp[ibd_res_meta$rel_origxTimp < 0] <- 0
ibd_res_meta$rel_dest1xTimpxPR <- ibd_res_meta$rel_dest1xTimp * ibd_res_meta$dest_PR
ibd_res_meta$rel_origxTimpxPR <- ibd_res_meta$rel_origxTimp * ibd_res_meta$origin_PR

# Obtain importation probability using the function (assuming r_importation_prob is defined elsewhere)
ibd_res_meta <- r_importation_prob(ibd_res_meta, r_origin = 'rel_origxTimpxPR',
                                  r_des1 = 'rel_dest1xTimpxPR', class_name = 'prob_imp_timp_pr')


# Create mask for travelled cases
mask <- travelled
total_cases <- sum(!is.na(ibd_res_meta$travel))  # Count total cases where 'travel' is not NA
total_reported_travel <- sum(mask)  # Count total reported travel

# Calculating imported stats
prov_names <- unique(ibd_res_meta$travel_prov)
prov_names <- prov_names[!is.na(prov_names)]
col_names = c('N. cases', '% cases') 
imported_stats <- data.frame(matrix(ncol = 2, nrow = length(prov_names)))
colnames(imported_stats) <- col_names
rownames(imported_stats) <- prov_names

for (dest in prov_names) {
  print(dest)
  mask_dest <- ibd_res_meta$travel_prov == dest  # Filter for current destination
  imported <- sum(ibd_res_meta$prob_imp_timp_pr[mask_dest], na.rm = TRUE)  # Sum of imported probabilities for the destination
  
  # Store results in the data frame
  if(!is.na(imported)){
  imported_stats[dest, 'N. cases'] <- imported
  imported_stats[dest, '% cases'] <- (imported / total_cases) * 100}
}

# Calculate total stats
imported_stats['Total', 'N. cases'] <- sum(imported_stats[, 'N. cases'], na.rm = TRUE)
imported_stats['Total', '% cases'] <- sum(imported_stats[, '% cases'], na.rm = TRUE)


# Load necessary library
library(ggplot2)
library(tidyverse)

# Create a histogram for 'prob_imp_timp_pr' with specified breaks and range
plot_data <- data.frame(ibd_res_meta[ibd_res_meta$prob_imp_timp_pr >= 0.000001,])
ggplot(plot_data, aes(x = prob_imp_timp_pr)) +
  geom_histogram(bins = 20, 
                 fill = "blue", color = "black") +
  theme_minimal() +
  labs(x = "Probability of Imported Case", y = "Density") +
  ggtitle("Histogram of Imported Case Probabilities")


ggplot(ibd_res_meta, aes(x = prob_imp_timp_pr)) +
  geom_histogram(bins = 20, 
                 fill = "blue", color = "black") +
  xlim(0.000001, 1) +
  theme_minimal() +
  labs(x = "Probability of Imported Case", y = "Density") +
  ggtitle("Histogram of Imported Case Probabilities")


