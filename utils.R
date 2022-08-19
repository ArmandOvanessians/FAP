

### Make CNA data set using WES and WGS
merge_wes_wgs_cna <- function(wes, wgs){
  # get important features from raw wes data
  wes_gene_cna <- wes %>% dplyr::select(Tumor_Sample_Barcode, gene, CopyNumber)
  colnames(wes_gene_cna) <- c("sample", "feature","value")
  # get important features from raw wgs data
  wgs_gene_cna <- wgs %>% dplyr::select(Tumor_Sample_Barcode, gene, CopyNumber)
  colnames(wgs_gene_cna) <- c("sample", "feature","value")
  # wes data that contains samples that wgs data does not
  wes_gene_cna <-  subset(wes_gene_cna, !sample %in% unique(wgs_gene_cna$sample)) # Samples of WES where we don't have WGS information
  # make a combines wes and wgs dataset
  cna <- rbind(wgs_gene_cna, wes_gene_cna)
  cna$view <- "CNA"
  
  return(cna)
}


### SNV Preprocessing

accumulate_snv_mutations <- function(pan_can_genes, snvs){
  # get unique genes of snv dataset
  snvs_genes <- unique(snvs$feature)
  # combine unique genes with pancancer driver genes to make a gene dataset
  all_genes <-sort(unique(c(pan_can_genes, snvs_genes)))
  # Turn the sample and feature columns to factors and specify all the genes at interest
  snvs$feature <- factor(snvs$feature, levels = all_genes)
  snvs$sample <- factor(snvs$sample)
  #Makes a df where a frequency of each mutation is added in another column
  snv_df <- as.data.frame(table(snvs[,c("sample", "feature")]))
  colnames(snv_df) <- c("samples", "feature", "value")
  
  return(snv_df)
}

merge_wes_wgs <- function(wes_snv_maf, wgs_snv_maf){
  
  # Select the sample and gene coloumn of the whole exome sequenceing data for snv
  wes_snv <- wes_snv_maf %>% dplyr::select(Tumor_Sample_Barcode, Hugo_Symbol)
  colnames(wes_snv) <- c("sample", "feature")
  # Select the sample and gene coloumn of the whole genome sequenceing data for snv
  wgs_snv <- wgs_snv_maf %>% dplyr::select(Tumor_Sample_Barcode, Hugo_Symbol)
  colnames(wgs_snv) <- c("sample", "feature") 
  #Samples of WGS where we don't have WES information
  wgs_snv <-  subset(wgs_snv, !sample %in% unique(wes_snv$sample))
  #Compiled SNV data that will be used for the remainder of the pre processing
  snvs <- rbind(wes_snv, wgs_snv)
  
  return (snvs) 
}

# I replaced this function with SNV_thresh
#df_above_n <- function(snv_df, n){
# remove samples where number of mutations for the gene is equal to 0 
#non_zero_snv <- snv_df[snv_df$value != 0,] 
# get unique non zero genes across samples
#non_zero_snv_genes <- non_zero_snv %>% dplyr::select(feature) 
# make a df where each genes is mapped to the quantity of samples it has a mutation in
#genes_quant_df <- as.data.frame(table(non_zero_snv_genes[,c("feature")]))
#colnames(genes_quant_df) <- c("gene", "frequency")
# get all genes that are mutated in at least n samples
#genes_above_thresh <- genes_quant_df[genes_quant_df$frequency >= n ,]
# get a subset of the orginal snv_df where the genes are above the threshold
#df <- subset(snv_df, feature %in% genes_above_thresh$gene)
#colnames(df) <- c("sample", "feature", "value")

# return(df)
#}

make_snv <- function(wes, wgs, panCancer){
  snvs <- merge_wes_wgs(wes_snv_maf, wgs_snv_maf)
  snv_df <- accumulate_snv_mutations(genes, snvs)
  #snv <- df_above_n(snv_df, threshold)
  snv$view <- "SNV"
  return(snv)
}

# get genes that have mutations across at least thresh*100 % of samples
snv_thresh <- function(df, thresh){
  #convert long format to wide
  df_wide <- df %>% dplyr::select(sample, feature, value)
  df_wide <- reshape(df_wide, idvar = "feature", timevar = "sample", direction = "wide")
  #replace non zero mutations with 1 in order to calculate the percent of mutations present across the sampes
  df_wide[-1][df_wide[-1] != 0] <- 1
  df_wide$mut_across_samp <- rowMeans(df_wide[-1])
  # Sort df in descending order
  sorted_df_wide <- df_wide[order(-df_wide$mut_across_samp),]
  # plot the ordered df in terms of rank and variance
  sorted_df_wide$counter <- seq(0, by=1, to = nrow(sorted_df_wide) -1)
  plot(sorted_df_wide$counter, sorted_df_wide$mut_across_samp, xlab= "Gene Rank", ylab= "Percent")
  abline(h = thresh, col = "red")
  
  
  # get genes with a percent of at least n
  sorted_df_wide_n <- sorted_df_wide[sorted_df_wide$mut_across_samp >= thresh, ]
  # subset rows of original df that contain only the n most variable genes
  n_most <- subset(df, feature %in% sorted_df_wide_n$feature)
  
  return (n_most)
}


### RNA seq pre processing

make_rna_seq_df <- function(rna_seq){
  # turns rna seq from a wide to long df
  rna_seq_df <- melt(setDT(rna_seq), id.vars = "Name", variable.name = "sample")
  rna_seq_df <- rna_seq_df %>% dplyr::select(sample, Name, value) 
  colnames(rna_seq_df) <- c("sample", "feature", "value")
  # removes the the "." and characters after
  rna_seq_df$sample <- gsub("\\..*", "", rna_seq_df$sample)
  # removes - within name
  rna_seq_df$sample <- gsub("-", "", rna_seq_df$sample)
  # removes the _a present in some sample names
  rna_seq_df$sample <- gsub("_a", "", rna_seq_df$sample)
  # removes the _b present in some sample name 
  rna_seq_df$sample <- gsub("_b", "", rna_seq_df$sample)
  # replace samples of F with EP and samples of G with JP
  rna_seq_df$sample <- gsub("F", "EP", rna_seq_df$sample)
  rna_seq_df$sample <- gsub("G", "JP", rna_seq_df$sample)
  # Edge Cases
  rna_seq_df$sample <- gsub("EP_Tra_001", "EPTrans1", rna_seq_df$sample)
  rna_seq_df$sample <- gsub("JP_Dsc_003", "JP_Desc_3", rna_seq_df$sample)
  #remove
  # creates another column caled view for MOFA purposes
  rna_seq_df$view <- "rna-seq"
  return(rna_seq_df)
}


### Get the most varied genes in the given df
n_most_var_df<- function(df, n){
  df <- df %>% dplyr::select(sample, feature, value)
  # convert long df to wide df
  df_wide <- reshape(df, idvar = "feature", timevar = "sample", direction = "wide")
  #error check 
  df_wide
  #Calculate the variance of each feature across the samples
  df_wide$row_var = rowVars(as.matrix(df_wide[,c(-1)]))
  # Sort df in descending order
  sorted_df_wide <- df_wide[order(-df_wide$row_var),]
  # plot the ordered df in terms of rank and variance
  sorted_df_wide$counter <- seq(0, by=1, to = nrow(sorted_df_wide) -1)
  plot(sorted_df_wide$counter, sorted_df_wide$row_var, xlab= "Gene Rank", ylab= "Variance")
  abline(v = n, col = "red")
  # subset n most variable genes
  sorted_df_wide_n <- head(sorted_df_wide, n)
  # subset rows of original df that contain only the n most variable genes
  n_most <- subset(df, feature %in% sorted_df_wide_n$feature)
  return(n_most)
}


### Configure meta data




