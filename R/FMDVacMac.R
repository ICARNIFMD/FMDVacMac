#' R Package for Prediction of Foot-and-Mouth Disease Virus Vaccine Matching Score (r1-value)
#'
#' @param field_isolate A `DNAStringSet` object containing one or more nucleotide sequences of field strains. 
#' A field Isolate represents a viral isolate collected from an infected animal in the field, and the comparison against the vaccine strain 
#' helps predict whether the vaccine will provide protection against the field isolate. The input should also be in standard FASTA format.
#'Example input: 
#'>ICFMD182/2022
#'ACAACCTCCACAGGTGAGTCGGCTAATCCCGTGACTGCCACCGTTGAAAACTACGGAGGCGAGACA....
#'>K154/12
#'ACCACTGCGACGGGAGAGTCAGCAGACCCTGTTACCACCACCGTTGAGAACTACGGCGGTGAAACA.....
#'>IND 168/2004
#'ACCACCACAACCGGTGAGTCGGCGGACCCGGTGACAACCACGGTTGAGAACTACGGAGGAACTACG.....
#' @param serotype the package trained for serotype "o".
#' @param vaccine_strain A `DNAStringSet` object containing the nucleotide sequence of the vaccine strain. 
#' The vaccine strain represents the reference strain used in the vaccine, and its sequence will be compared to the field strains
#' to predict the matching score. The input sequence should be in standard FASTA format, containing only nucleotide bases: A, T, G, and C.
#'Example input: 
#'>VP1_R2/1975
#'ACCACCTCCCCGGGTGAGTCAGCTGACCCCGTGACCGCCACTGTTGAAAACTACGGCGGTGAGACACAGG...
#' @param model_choice choice Specifies which machine learning model to use for prediction. Options are:
#'"xgboost"`: Uses the XGBoost model, a gradient boosting method known for its high performance in structured/tabular data.
#'"rf"`: Uses the Random Forest model, an ensemble method based on decision trees, often effective in handling complex datasets with multiple features.
#'"svm"`: Uses the Support Vector Machine (SVM) model, which is commonly applied for classification tasks and works well in high-dimensional spaces. 
#' "GBM" :ses Gradient Boosting Machine, an ensemble method that builds decision trees sequentially to improve predictive accuracy.
#'The choice of model will affect the prediction of the vaccine matching score based on the field and vaccine strain sequences.
#'
#' @return A `data.frame` with the following columns:
#'Sl. No.	Isolate_ID	r ̂_1-value	Remarks
#'1	ICFMD182/2022	0.45	    Protection
#'2	K154/12	        0.67	    Protection
#'3	IND 168/2004	0.37	    Protection
#'Sl. No.: The serial number of each field Isolate sequence.
#'`Isolate_ID`: The identifier for each field Isolate sequence (e.g., "ICFMD182/2022").
#'`r̂_1-value`: The predicted vaccine matching score for the field Isolate.
#'`Remarks`: A label indicating whether the match is "Protected" or "Non-Protected" based on the matching score. 
#'The output data.frame includes the following remarks:
#'"Protected": if the r̂_1-value ≥ 0.3, indicating cross-protection between the field and vaccine strains.
#'"Non-Protected": if the r̂_1-value < 0.3, indicating lack of cross-protection.
#' @export
#'
#' @examples FMDVacMac (field_isolate, vaccine_Strain, "rf")
#'field_isolate: 
#'>ICFMD182/2022
#'ACAACCTCCACAGGTGAGTCGGCTAATCCCGTGACTGCCACCGTTGAAAACTACGGAGGCGAGACA....
#'>K154/12
#'ACCACTGCGACGGGAGAGTCAGCAGACCCTGTTACCACCACCGTTGAGAACTACGGCGGTGAAACA.....
#' vaccine_strain: 
#'>VP1_R2/1975
#'ACCACCTCCCCGGGTGAGTCAGCTGACCCCGTGACCGCCACTGTTGAAAACTACGGCGGTGAGACACAGG...
#' model_choice : Options are "xgboost", "rf", or "svm"
#' @import Biostrings
#' @import e1071
#' @import openxlsx
#' @import randomForest
#' @import stringr
#' @import xgboost
#' @import caret

FMDVacMac <- function(field_isolate, serotype, vaccine_Strain, model_choice) {
  if (tolower(serotype) != "o") {
    stop("Invalid serotype. This function only supports serotype 'O'. Please input the correct serotype.")
  }
  
  # Read vaccine sequence
  if (vaccine_Strain == "R21975") {
    VPR<- system.file("extdata", "VP1_R21975_1.txt", package = "FMDVacMac")
    vaccine_Strain_seq <- readDNAStringSet(VPR, format = "fasta")
  } else if (vaccine_Strain == "Manisha") {
    Manisa<- system.file("extdata", "O1_Manisa_Turkey69.txt", package = "FMDVacMac")
    vaccine_Strain_seq <- readDNAStringSet(Manisa, format = "fasta")
  } else if (vaccine_Strain == "Ethopia") {
    ETH<- system.file("extdata", "ETH382005_VP1.txt", package = "FMDVacMac")
    vaccine_Strain_seq <- readDNAStringSet(ETH, format = "fasta")
  } else {
    stop("Invalid vaccine_Strain. Choose from: 'R21975', 'Manisha', 'Ethopia'.")
  }
  
  generate_vaccine_matching_features <- function(field_isolate, vaccine_Strain, k = 3) {
    field_features <- oligonucleotideFrequency(field_isolate, width = k)
    vaccine_features <- oligonucleotideFrequency(vaccine_Strain, width = k)
    
    calculate_distances <- function(field, vaccine) {
      apply(field, 1, function(row) sqrt(sum((vaccine - row)^2)))
    }
    distance_vector <- calculate_distances(field_features, vaccine_features[1, ])
    
    field_AA_seq <- translate(field_isolate, if.fuzzy.codon = "X")
    vaccine_AA_seq <- translate(vaccine_Strain, if.fuzzy.codon = "X")
    
    combined_sequences <- c(as.character(field_AA_seq), as.character(vaccine_AA_seq))
    seq_matrix <- do.call(rbind, strsplit(combined_sequences, ""))
    polymorphic_sites <- which(apply(seq_matrix, 2, function(col) length(unique(col)) > 1))
    
    polymorphic_matrix <- seq_matrix[, polymorphic_sites, drop = FALSE]
    coded_differences <- t(apply(polymorphic_matrix, 2, function(col) ifelse(col == col[1], 0, 1)))
    coded_differences_df <- as.data.frame(t(coded_differences))
    last_sequence <- coded_differences_df[nrow(coded_differences_df), ]
    
    poly_dist <- apply(coded_differences_df[-nrow(coded_differences_df), ], 1, function(sequence) {
      combined <- rbind(sequence, last_sequence)
      distance_matrix <- dist(combined, method = "binary")
      as.numeric(distance_matrix[1])
    })
    
    sequences <- as.character(c(field_AA_seq, vaccine_AA_seq))
    poisson_corrected_distance_variable_length <- function(seq1, seq2) {
      max_length <- max(nchar(seq1), nchar(seq2))
      seq1_padded <- str_pad(seq1, max_length, side = "right", pad = "N")
      seq2_padded <- str_pad(seq2, max_length, side = "right", pad = "N")
      num_diff <- sum(strsplit(seq1_padded, NULL)[[1]] != strsplit(seq2_padded, NULL)[[1]])
      p_hat <- num_diff / max_length
      if (p_hat <= 1) -log(1 - p_hat) else NA
    }
    
    vaccine_idx <- length(sequences)
    field_sequences <- sequences[1:(vaccine_idx - 1)]
    vaccine_sequence <- sequences[vaccine_idx]
    
    AA_Distance <- sapply(field_sequences, function(field_seq) {
      poisson_corrected_distance_variable_length(field_seq, vaccine_sequence)
    })
    
    find_nglycosylation_sites <- function(sequence) {
      pattern <- "N[^P][ST]"
      matches <- gregexpr(pattern, toupper(sequence), perl = TRUE)
      if (length(matches[[1]]) > 0 && matches[[1]][1] != -1) {
        substring(sequence, matches[[1]], matches[[1]] + 2)
      } else {
        character(0)
      }
    }
    
    ngly_sites <- lapply(sequences, find_nglycosylation_sites)
    field_isolates <- ngly_sites[1:(length(ngly_sites) - 1)]
    vaccine_sites <- ngly_sites[[length(ngly_sites)]]
    
    jaccard_index <- function(A, B) {
      intersection <- length(intersect(A, B))
      union <- length(union(A, B))
      if (union == 0) 0 else intersection / union
    }
    
    vaccine_field_distances <- sapply(field_isolates, function(field_sites) {
      jaccard_index(field_sites, vaccine_sites)
    })
    
    combined_df <- data.frame(
      "AMINO ACID DIST" = AA_Distance,
      "N-GLYC SITE DIS" = vaccine_field_distances,
      "AMINO ACID POLYM DIST" = poly_dist,
      "Kmer dist" = distance_vector
    )
    return(combined_df)
  }
  
  features_df <- generate_vaccine_matching_features(field_isolate, vaccine_Strain_seq)
  result_matrix <- as.matrix(features_df)
  rownames(result_matrix) <- names(field_isolate)
  
  vaccine_strain <- vaccine_Strain  # fix the case mismatch
  
  if (vaccine_strain == "R21975" && model_choice == "rf") {
    model_path <- system.file("extdata", "O_VP1_R21975_RF.RDS", package = "FMDVacMac")
    model <- readRDS(model_path)
    predictions <- predict(model, result_matrix, type = "response")
  } else if (vaccine_strain == "R21975" && model_choice == "svm") {
    model_path <- system.file("extdata", "O_VP1_R21975_svm.RDS", package = "FMDVacMac")
    model <- readRDS(model_path)
    predictions <- predict(model, result_matrix, type = "pred")
  } else if (vaccine_strain == "R21975" && model_choice == "xgboost") {
    model_path <- system.file("extdata", "O_VP1_R21975_xgboost.RDS", package = "FMDVacMac")
    model <- readRDS(model_path)
    predictions <- predict(model, result_matrix)
  } else if (vaccine_strain == "R21975" && model_choice == "GBM") {
    model_path <- system.file("extdata", "O_VP1_R21975_gbm.RDS", package = "FMDVacMac")
    model <- readRDS(model_path)
    predictions <- predict(model, result_matrix)
  } else if (vaccine_strain == "Manisha" && model_choice == "rf") {
    model_path <- system.file("extdata", "O_O1_Manisa_Turkey69_RF.RDS", package = "FMDVacMac")
    model <- readRDS(model_path)
    predictions <- predict(model, result_matrix, type = "response")
  } else if (vaccine_strain == "Manisha" && model_choice == "svm") {
    model_path <- system.file("extdata", "O_O1_Manisa_Turkey69_svm.RDS", package = "FMDVacMac")
    model <- readRDS(model_path)
    predictions <- predict(model, result_matrix, type = "pred")
  } else if (vaccine_strain == "Manisha" && model_choice == "xgboost") {
    model_path <- system.file("extdata", "O_O1_Manisa_Turkey69_xgboost.RDS", package = "FMDVacMac")
    model <- readRDS(model_path)
    predictions <- predict(model, result_matrix)
  } else if (vaccine_strain == "Manisha" && model_choice == "GBM") {
    model_path <- system.file("extdata", "O_O1_Manisa_Turkey69_gbm.RDS", package = "FMDVacMac")
    model <- readRDS(model_path)
    predictions <- predict(model, result_matrix)
  } else if (vaccine_strain == "Ethopia" && model_choice == "rf") {
    model_path <- system.file("extdata", "O_ETH382005_VP1_RF.RDS", package = "FMDVacMac")
    model <- readRDS(model_path)
    predictions <- predict(model, result_matrix, type = "response")
  } else if (vaccine_strain == "Ethopia" && model_choice == "svm") {
    model_path <- system.file("extdata", "O_ETH382005_VP1_svm.RDS", package = "FMDVacMac")
    model <- readRDS(model_path)
    predictions <- predict(model, result_matrix, type = "pred")
  } else if (vaccine_strain == "Ethopia" && model_choice == "xgboost") {
    model_path <- system.file("extdata", "O_ETH382005_VP1_xgboost.RDS", package = "FMDVacMac")
    model <- readRDS(model_path)
    predictions <- predict(model, result_matrix)
  } else if (vaccine_strain == "Ethopia" && model_choice == "GBM") {
    model_path <- system.file("extdata", "O_ETH382005_VP1_gbm.RDS", package = "FMDVacMac")
    model <- readRDS(model_path)
    predictions <- predict(model, result_matrix)
  } else {
    stop("Invalid vaccine_Strain or model_choice. Choose vaccine_Strain from: 'R21975', 'Manisha', 'Ethopia' AND model_choice from:'rf', 'svm', 'xgboost', 'GBM'.")
  }
  
  output <- data.frame(
    Accession = rownames(result_matrix),
    Prediction_Score = predictions,
    Remarks = ifelse(predictions > 0.3, "Protected", "Non-Protected")
  )
  print(output)
  return(output)
}
