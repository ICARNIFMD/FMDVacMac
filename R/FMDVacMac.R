#' R Package for Prediction of Foot-and-Mouth Disease Virus Vaccine Matching Score (r1-value)
#'
#' @param vaccine_strain A DNAStringSet object containing the vaccine strain sequence
#' @param field_isolate A DNAStringSet object containing the field strain sequences
#' @param model_choice Specifies which model to use for prediction. Options are "xgboost", "rf", or "svm"
#'
#' @return A `data.frame` with the following columns:
#' Field_Strain_ID: The identifier for each field strain sequence.
#' Matching_Score: The predicted vaccine matching score for the field strain.

#' @export
#'
#' @examples
#' vaccine_strain: R2/1992 TACTACTTCGCAGATTTAGAGGTGGCAGTGAAACACGAGGGGAACCTCACCTGGGTCCCGAACGGGGCGC
#' field_isolate: ACCACTGCGACGGGAGAGTCAGCAGACCCTGTTACCACCACCGTTGAGAACTACGGCGGTGAAACACAGG
#' model_choice : Options are "xgboost", "rf", or "svm"
#' FMDVacMac (vaccine_strain, field_isolate, "rf")
Fmdvacmac <- function(vaccine_strain, field_isolate, model_choice) {

  ############## Feature Engineering Function ##############
  generate_vaccine_matching_features <- function(field_Strain, vaccine_Strain, k = 3) {
    field_features <- oligonucleotideFrequency(field_Strain, width = k)
    vaccine_features <- oligonucleotideFrequency(vaccine_Strain, width = k)

    calculate_distances <- function(field, vaccine) {
      apply(field, 1, function(row) sqrt(sum((vaccine - row) ^ 2)))
    }
    distance_vector <- calculate_distances(field_features, vaccine_features[1, ])

    field_AA_seq <- translate(field_Strain, if.fuzzy.codon = "X")
    vaccine_AA_seq <- translate(vaccine_Strain, if.fuzzy.codon = "X")

    # Polymorphic Distance
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

    # Amino Acid Distance
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
    field_indices <- 1:(vaccine_idx - 1)
    vaccine_sequence <- sequences[vaccine_idx]
    field_sequences <- sequences[field_indices]

    AA_Distance <- sapply(field_sequences, function(field_seq) {
      poisson_corrected_distance_variable_length(field_seq, vaccine_sequence)
    })

    # N-glycosylation Sites
    find_nglycosylation_sites <- function(sequence) {
      pattern <- "N[^P][ST]"
      matches <- gregexpr(pattern, toupper(as.character(sequence)), perl = TRUE)
      if (length(matches[[1]]) > 0 && matches[[1]][1] != -1) {
        substring(sequence, matches[[1]], matches[[1]] + 2)
      } else {
        character(0)
      }
    }

    nucl_to_amino_acid_seq_final <- c(as.character(field_AA_seq), as.character(vaccine_AA_seq))
    ngly_sites <- lapply(nucl_to_amino_acid_seq_final, find_nglycosylation_sites)

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

  ############## Feature Extraction ##############
  features_df <- generate_vaccine_matching_features(field_Strain, vaccine_Strain)
  result_matrix <- as.matrix(features_df)
  rownames(result_matrix) <- names(field_Strain)

  ############## Prediction based on Model Choice ##############
  if (model_choice == "xgboost") {
    model <- readRDS("O_VP1_R21975_xgboost.RDS")
    predictions <- predict(model, result_matrix)
  } else if (model_choice == "rf") {
    model <- readRDS("O_VP1_R21975_RF.RDS")
    predictions <- predict(model, result_matrix, type = "response")
  } else if (model_choice == "svm") {
    model <- readRDS("O_VP1_R21975_svm.RDS")
    predictions <- predict(model, result_matrix, type = "pred")
  } else {
    stop("Invalid model_choice. Choose from 'xgboost', 'rf', or 'svm'.")
  }

  ############## Final Output ##############
  output <- data.frame(
    Accession = rownames(result_matrix),
    Prediction_Score = predictions,
    Remarks = ifelse(predictions > 0.3, "Protected", "Non-Protected")
  )

  print(output)
  return(output)
}
