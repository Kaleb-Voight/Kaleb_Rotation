library(dplyr)
library(purrr)
library(readr)
library(biomaRt)
library(EnsDb.Hsapiens.v86)
library(tidyr)
library(stringr)

setwd("/Users/kalebvoight/Desktop/Kaleb\ Rotation")

BreitmeyerC1R1 <- ("BreitmeyerC1R1.csv")
BreitmeyerC1R2 <- ("BreitmeyerC1R2.csv")
BreitmeyerC1R3 <- ("BreitmeyerC1R3.csv")
BreitmeyerC2R1 <- ("BreitmeyerC2R1.csv")
BreitmeyerC2R2 <- ("BreitmeyerC2R2.csv")
BreitmeyerC2R3 <- ("BreitmeyerC2R3.csv")
GuaquetaR1 <- ("GuaquetaR1.csv")
GuaquetaR2 <- ("GuaquetaR2.csv")
GuaquetaR3 <- ("GuaquetaR3.csv")
GuaquetaR4 <- ("GuaquetaR4.csv")
Rammonds <- ("Rammonds.csv")
Penney <- ("Penney.csv")

file_path <- c(BreitmeyerC1R1, BreitmeyerC1R2, BreitmeyerC1R3, BreitmeyerC2R1, BreitmeyerC2R2, BreitmeyerC2R3, GuaquetaR1, GuaquetaR2, GuaquetaR3, GuaquetaR4, Rammonds, Penney)

#Ensembl Conversion
is_gene_id <- function(id) {
  grepl("^ENSG", id)
}

clean_gene_name <- function(name) {
  gsub("\\..*$", "", name)
}

convert_ids_to_names <- function(gene_ids) {
  gene_ids_stripped <- clean_gene_name(gene_ids)
  
  edb <- EnsDb.Hsapiens.v86
  
  print(gene_ids_stripped)
  
  gene_info <- genes(edb, filter = GeneIdFilter(gene_ids_stripped))
  
  gene_df <- as.data.frame(gene_info)
  
  id_to_name <- setNames(gene_df$gene_name, clean_gene_name(gene_df$gene_id))
  return(id_to_name)
}

process_csv_file <- function(file_path) {
  tryCatch({
    df <- read_csv(file_path)
    message("Reading file: ", file_path)
    
    if (!"Gene_Name" %in% colnames(df)) {
      stop("Column 'Gene_Name' not found in file: ", file_path)
    }
    
    df$Gene_Name <- sapply(df$Gene_Name, clean_gene_name)
    message("Cleaned gene names")
    
    gene_names <- df$Gene_Name
    gene_ids <- gene_names[is_gene_id(gene_names)]
    
    if (length(gene_ids) > 0) {
      id_to_name <- convert_ids_to_names(gene_ids)
      message("Converted gene IDs to names")
      
      df$Gene_Name <- sapply(df$Gene_Name, function(x) {
        if (is_gene_id(x)) {
          return(id_to_name[x])
        } else {
          return(x)
        }
      })
    }
    
    processed_file_path <- paste0("P_", basename(file_path))
    write_csv(df, processed_file_path)
    message("P_", processed_file_path)
  }, error = function(e) {
    message("Error processing file ", file_path, ": ", e$message)
  })
}

map(file_path, process_csv_file)

#Combination of Files
MG <- ("Combined_MG_Dataset.csv")
CN <- ("Combined_CN_iPSC_Dataset.csv")

file_list <- c(MG, CN)

process_file <- function(file) {
  df <- read_csv(file)
  
  df$Gene_Name <- gsub("\\..*$", "", df$Gene_Name)
  df$Gene_Name <- str_trim(df$Gene_Name)
  df$Gene_Name <- toupper(df$Gene_Name)
  
  df <- df %>% distinct(Gene_Name, .keep_all = TRUE)
}

data_list <- lapply(file_list, process_file)

combined_data <- reduce(data_list, inner_join, by = "Gene_Name")

write_csv(combined_data, "Combined_MG_CN_iPSC_Dataset.csv")
