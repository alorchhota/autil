#' Functions to filter gene expression data
#' 
#' Each function filters gene expression data (gene x sample) based on a defined criteria.
#' @param expr.df data.frame. processed expression data (gene x sample). See details.
#' @param tpm.df data.frame. raw (unprocessed) TPM data (gene x sample). See details.
#' @param count.df data.frame. raw (uncrocessed) read-count data (gene x sample). See details.
#' @param raw.df data.frame. any type of raw (uncrocessed) data (gene x sample). See details.
#' @param annot.gene data.frame. gene annotation data frame with columns 'gene_id', 'gene_name', 'chr', 'start_pos', 'end_pos', etc.  See details.
#' @param annot.mappability data.frame. gene mappability data frame with at columns 'gene_id', 'mappability', etc. See details.
#' @param min.mappability numeric. minimum gene mappability.
#' @param chr.include character vector. chromosomoes to include.
#' @param chr.exclude character vector. chromosomoes to exclude.
#' @param chr.col character or numeric. column name or column number containing chromosomes of genes in the gene annotation data frame.
#' @param mappability.col character or numeric. column name or column number containing mappability of genes.
#' @param type.col character or numeric. column name or column number containing gene types in the gene annotation data frame.
#' @param type.values character vector. target gene types (e.g., protein_coding) to include.
#' @param min.tpm numeric. minimum TPM.
#' @param min.count numeric. minimum read count.
#' @param min.samples numeric. minimum number of samples.
#' @param min.var numeric. minimum variance of expression per gene across all samples.
#' @param min.mean numeric. minimum mean of expression per gene across all samples.
#' @param n maximum number of genes.
#' @details 
#' Above functions filter processed gene expression data (\code{expr.df}) based on several different criteria.
#' Typical filtering scenarios include:
#' \itemize{
#'   \item keeping only autosomal genes, excluding genes in sex and mitochondrial chromosomes using \code{filter_expr_on_chr},
#'   \item including genes with mappability >= 0.8 for an analysis using \code{filter_expr_on_mappability},
#'   \item inlcluding only protein-coding genes in an analysis using \code{filter_expr_on_gene_type},
#'   \item include genes with TPM > 0.1 in \eqn{\ge} 10 samples and read count > 6 in \eqn{\ge} 10 samples using \code{filter_expr_on_tpm_read},
#'   \item take only 500 most variant genes using \code{filter_expr_on_variance}, and
#'   \item take only 500 genes with highest coefficient of variation using \code{filter_expr_on_coeff_of_variation}.
#' }
#' 
#' Data formats are important for these filter functions.
#' 
#' Any expression-like data (\code{expr.df}, \code{raw.df}, \code{tpm.df}, or \code{count.df}) must be in a data.frame or a matrix-like object, 
#' where each row represents a gene and each column represents a sample (or an individual).
#' The data frame must have row-names with unique gene ids and column-names with unique sample ids.
#' 
#' Any gene annotation data (\code{annot.gene}, or \code{annot.mappability}) must be in a data.frame or a matrix-like object, 
#' where each row represents a gene, and each column represents a property of the gene.
#' The data frame must have row-names with unique gene ids matching to gene ids in expression data.
#' While necessary columns in an annotation data frame can be configured in each function,
#' the annotation data frame may contain more than necessary columns.
#' @describeIn filter_expr filteres to include or explore genes from given chromosomes (if not NULL). Chromosome names with or without "chr" prefix are gracefully handled.
#' @return Each function returns a data frame or matrix like object (similar to input type) with filtered expression data.
#' @export
filter_expr_on_chr <- function(expr.df, annot.gene, chr.include=as.character(1:22), chr.exclude=NULL, chr.col='chr'){
  stopifnot(length(chr.col) == 1 && ((is.character(chr.col) && chr.col %in% colnames(annot.gene)) || (is.numeric(chr.col) && chr.col <= ncol(annot.gene))))
  features.passed = rownames(expr.df)
  if(!is.null(chr.include)){
    chr.include = unique(c(extend_chr(chr.include), shorten_chr(chr.include)))
    features.passed = features.passed[annot.gene[features.passed, chr.col] %in% chr.include]
  }
  if(!is.null(chr.exclude)){
    chr.exclude = unique(c(extend_chr(chr.exclude), shorten_chr(chr.exclude)))
    features.passed = features.passed[!annot.gene[features.passed, chr.col] %in% chr.exclude]
  }
  expr.df = expr.df[features.passed,,drop=F]
  return(expr.df)
}

#' @describeIn filter_expr filters out genes with mappability less than a threshold.
#' @export
filter_expr_on_mappability <- function(expr.df, annot.mappability, min.mappability=0.97, mappability.col='mappability'){
  stopifnot(length(mappability.col) == 1 && ((is.character(mappability.col) && mappability.col %in% colnames(annot.mappability)) || (is.numeric(mappability.col) && mappability.col <= ncol(annot.mappability)) ))
  features = rownames(expr.df)
  mappabilities = annot.mappability[features, mappability.col]
  features.passed = features[mappabilities >= min.mappability]
  expr.df = expr.df[features.passed,,drop=F]
  return(expr.df)
}

#' @export
#' @describeIn filter_expr filters expression to include specific types of genes only.
filter_expr_on_gene_type <- function(expr.df, annot.gene, type.col = 'gene_type', type.values = 'protein_coding'){
  stopifnot(length(type.col) == 1 && ((is.character(type.col) && type.col %in% colnames(annot.gene)) || (is.numeric(type.col) && type.col <= ncol(annot.gene))))
  features = rownames(expr.df)
  gene_types = annot.gene[features, type.col]
  features.passed = features[!is.na(gene_types) & gene_types %in% type.values]
  expr.df = expr.df[features.passed,,drop=F]
  return(expr.df)
}

#' @export
#' @describeIn filter_expr filters expression data to include each gene with TPM > \code{min.tpm} in \eqn{\ge} \code{min.samples} samples and with read count > \code{min.count} in \eqn{\ge} \code{min.samples} samples.
filter_expr_on_tpm_read <- function(expr.df, tpm.df, count.df, min.tpm = 0.1, min.count = 6, min.samples = 10){
  ### make sure matrices have the same order as expr.df
  features = rownames(expr.df)
  samples = colnames(expr.df)
  
  stopifnot(all(features %in% rownames(tpm.df)))
  stopifnot(all(samples %in% colnames(tpm.df)))
  stopifnot(all(features %in% rownames(count.df)))
  stopifnot(all(samples %in% colnames(count.df)))
  
  tpm.df <- tpm.df[features, samples,drop=F]
  count.df <- count.df[features, samples,drop=F]
  
  ### filter
  n_samples_w_min_tpm <- apply(tpm.df>min.tpm, 1, sum)
  n_samples_w_min_count <- apply(count.df>min.count, 1, sum)
  has.min.samples <- (n_samples_w_min_tpm >= min.samples) & (n_samples_w_min_count >= min.samples)
  features.passed <- names(has.min.samples[has.min.samples])
  expr.df <- expr.df[features.passed,,drop=F]
  
  return(expr.df)
}

#' @export
#' @describeIn filter_expr filters expression data to include most variant \code{n} genes with mean \eqn{\ge} \code{min.mean} and variance \eqn{\ge} \code{min.var}
filter_expr_on_variance <- function(expr.df, raw.df, n, min.var=1e-6, min.mean=-Inf){
  ### make sure matrices have the same order as expr.df
  features = rownames(expr.df)
  samples = colnames(expr.df)
  
  stopifnot(all(features %in% rownames(raw.df)))
  stopifnot(all(samples %in% colnames(raw.df)))
  
  raw.df <- raw.df[features, samples,drop=F]
  
  sds = apply(raw.df, 1, sd)
  means = apply(raw.df, 1, mean)
  sds = sds[sds>=sqrt(min.var) & means >=min.mean]
  sds.sorted = sort(sds, decreasing = T)
  features.passed = names(sds.sorted[1:min(n, length(sds.sorted))])
  expr.df = expr.df[features.passed,,drop=F]
  return(expr.df)
}

#' @export
#' @describeIn filter_expr filters expression data to include genes with \code{n} highest coefficient of variation (variance/mean) and mean \eqn{\ge} \code{min.mean} and variance \eqn{\ge} \code{min.var}
filter_expr_on_coeff_of_variation <- function(expr.df, raw.df, n, min.var=1e-6, min.mean=1e-2){
  ### make sure matrices have the same order as expr.df
  features = rownames(expr.df)
  samples = colnames(expr.df)
  
  stopifnot(all(features %in% rownames(raw.df)))
  stopifnot(all(samples %in% colnames(raw.df)))
  
  raw.df <- raw.df[features, samples,drop=F]
  
  sds = apply(raw.df, 1, sd)
  means = apply(raw.df, 1, mean)
  sds1 = sds[sds>=sqrt(min.var) & means >=min.mean]
  means1 = means[sds>=sqrt(min.var) & means >=min.mean]
  coeff.var = sds1 / means1
  coeff.var.sorted = sort(coeff.var, decreasing = T)
  features.passed = names(coeff.var.sorted[1:min(n, length(coeff.var.sorted))])
  expr.df = expr.df[features.passed,,drop=F]
  return(expr.df)
}

