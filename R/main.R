# ===========================================================================
#  m6APrediction -- R/main.R
#
#  This file contains ONLY the three prediction functions carried over from
#  Practical 5 and their roxygen2 documentation. No setwd(), library(),
#  read.csv() or devtools:: calls belong here: every command that is not a
#  function definition must live in a separate script outside the package.
# ===========================================================================


#' Encode DNA Sequences as Positional Nucleotide Factors
#'
#' Splits each equally long DNA sequence into single nucleotides and returns
#' one factor column per nucleotide position (\code{nt_pos1}, \code{nt_pos2},
#' ...). Every column carries the same four levels, \code{"A"}, \code{"T"},
#' \code{"C"} and \code{"G"}, in that fixed order.
#'
#' The level order must not be changed: \code{rf_fit} was trained on features
#' that were factorised with exactly this level set, so reordering them would
#' silently invalidate the model matrix. This function is an internal helper
#' used by \code{prediction_multiple()} and \code{prediction_single()}; it is
#' deliberately not exported.
#'
#' @param dna_strings A character vector of DNA sequences of equal length,
#'   for example \code{c("GAACA", "GGACC")}. For the m6A model each sequence
#'   is a 5-mer: two nucleotides upstream of the candidate adenine, the
#'   adenine itself, and two nucleotides downstream. The length of the first
#'   element determines the number of output columns.
#'
#' @return A \code{data.frame} with one row per input sequence and one factor
#'   column per nucleotide position, named \code{nt_pos1} to \code{nt_posn}.
#'   Each column is a factor on the levels \code{c("A", "T", "C", "G")}.
#'
#' @keywords internal
dna_encoding <- function(dna_strings){
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(unlist(strsplit(dna_strings, "")), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}


#' Predict m6A Status for Multiple Candidate Sites
#'
#' Scores a whole table of candidate sites with the random forest model
#' trained in Practical 5. Each \code{DNA_5mer} is expanded into five
#' positional nucleotide factors by \code{dna_encoding()}, the encoded columns
#' are appended to the supplied feature table, and the model is asked for the
#' probability of the positive (m6A) class. Probabilities are then converted
#' into a binary status using \code{positive_threshold}.
#'
#' @param ml_fit A fitted random forest model, normally the object shipped in
#'   \code{inst/extdata/rf_fit.rds} and loaded with
#'   \code{readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))}.
#'   It must be the model trained in Practical 5; it expects the seven
#'   features described under \code{feature_df} plus the five \code{nt_pos}
#'   columns produced by \code{dna_encoding()}.
#'
#' @param feature_df A \code{data.frame} with one row per candidate site and
#'   at least the seven required columns \code{gc_content}, \code{RNA_type},
#'   \code{RNA_region}, \code{exon_length}, \code{distance_to_junction},
#'   \code{evolutionary_conservation} and \code{DNA_5mer}. \code{DNA_5mer}
#'   must be a five-nucleotide character string with a central \code{A}.
#'   \code{exon_length} and \code{distance_to_junction} are supplied already
#'   log2-transformed and must not be transformed a second time.
#'
#' @param positive_threshold A single number, default \code{0.5}. A site
#'   whose predicted probability is strictly greater than this value is
#'   labelled \code{"Positive"}; a site exactly equal to it is labelled
#'   \code{"Negative"}.
#'
#' @return The input \code{data.frame} with two columns appended:
#'   \code{predicted_m6A_prob}, the predicted probability of the positive
#'   class (numeric), and \code{predicted_m6A_status}, a character column
#'   containing \code{"Positive"} or \code{"Negative"}. All original columns
#'   are preserved, in their original order.
#'
#' @examples
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds",
#'                               package = "m6APrediction"))
#' example_df <- read.csv(system.file("extdata", "m6A_input_example.csv",
#'                                    package = "m6APrediction"))
#' prediction_multiple(ml_fit, example_df[1:3, ], positive_threshold = 0.6)[,
#'   c("predicted_m6A_prob", "predicted_m6A_status")]
#'
#' @export
#' @import randomForest
#' @importFrom stats predict
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction",
                  "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df))) # Check errors if incorrect column names of input data.frame
  feature_df$RNA_type <- factor(feature_df$RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR"))
  infered_prob <- predict(ml_fit, newdata = cbind(feature_df, dna_encoding(feature_df$DNA_5mer)), type = "prob")
  feature_df$predicted_m6A_prob <- infered_prob[,"Positive"]
  feature_df$predicted_m6A_status <- "Negative"
  feature_df$predicted_m6A_status[feature_df$predicted_m6A_prob > positive_threshold] <- "Positive"
  return(feature_df) # Return a data.frame with predicted m6A prob and predicted m6A status
}


#' Predict m6A Status for a Single Candidate Site
#'
#' A convenience wrapper around \code{prediction_multiple()} that takes the
#' seven feature values of one candidate site as separate arguments, so that
#' a single site can be scored without building a \code{data.frame} by hand.
#' The arguments are assembled into a one-row table, passed to the batch
#' function, and the resulting probability and status are returned together as
#' a named vector.
#'
#' @param ml_fit A fitted random forest model, normally the object shipped in
#'   \code{inst/extdata/rf_fit.rds} and loaded with
#'   \code{readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))}.
#'
#' @param gc_content A number between 0 and 1, the GC content of the sequence
#'   window around the candidate site.
#'
#' @param RNA_type A character string, one of \code{"mRNA"}, \code{"lincRNA"},
#'   \code{"lncRNA"} or \code{"pseudogene"}.
#'
#' @param RNA_region A character string, one of \code{"CDS"},
#'   \code{"intron"}, \code{"3'UTR"} or \code{"5'UTR"}.
#'
#' @param exon_length A number, the length of the hosting exon, supplied
#'   already log2-transformed.
#'
#' @param distance_to_junction A number, the distance from the candidate site
#'   to the nearest exon junction, supplied already log2-transformed.
#'
#' @param evolutionary_conservation A number between 0 and 1, the
#'   evolutionary conservation score of the 5-mer context.
#'
#' @param DNA_5mer A five-nucleotide character string with a central
#'   \code{A}, for example \code{"GGACA"}.
#'
#' @param positive_threshold A single number, default \code{0.5}. A
#'   probability strictly greater than this value gives \code{"Positive"};
#'   a probability exactly equal to it gives \code{"Negative"}.
#'
#' @return A named character vector of length two:
#'   \code{predicted_m6A_prob}, the probability rendered as text (for example
#'   \code{"0.51"}), and \code{predicted_m6A_status}, either
#'   \code{"Positive"} or \code{"Negative"}. Both elements are character
#'   because a number and a label are combined in one vector; apply
#'   \code{as.numeric()} to the first element if a numeric probability is
#'   needed.
#'
#' @examples
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds",
#'                               package = "m6APrediction"))
#' prediction_single(ml_fit, gc_content = 0.6, RNA_type = "mRNA",
#'                   RNA_region = "CDS", exon_length = 12,
#'                   distance_to_junction = 5,
#'                   evolutionary_conservation = 0.8, DNA_5mer = "GGACA",
#'                   positive_threshold = 0.5)
#'
#' @export
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length,
                              distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){
  pred_outcome <- prediction_multiple(ml_fit = ml_fit,
                                      data.frame(gc_content = gc_content,
                                                 RNA_type = factor(RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")),
                                                 RNA_region = factor(RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR")),
                                                 exon_length = exon_length,
                                                 distance_to_junction = distance_to_junction,
                                                 evolutionary_conservation = evolutionary_conservation,
                                                 DNA_5mer = DNA_5mer),
                                      positive_threshold = positive_threshold)
  returned_vector <- c(pred_outcome$predicted_m6A_prob, pred_outcome$predicted_m6A_status)
  names(returned_vector) <- c("predicted_m6A_prob", "predicted_m6A_status")
  return(returned_vector) # Return a named vector with predicted m6A prob and predicted m6A status
}
