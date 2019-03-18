#' Filter cross-mappability data
#' 
#' Filter cross-mappability data to include gene pairs where 1) each gene is within a given set of genes, 2) the cross-mappability score is within a given range, or 3) the genes are not positionally overlapped.
#' 
#' @rdname filter_crossmap
#' @param crossmap_df Cross-mappability data frame. First two columns should contain genes (character) and the 3rd column should contain cross-mappability score (numeric) from the 1st gene to the 2nd gene.
#' @param overlap_df Position overlap data frame. First two columns should contain genes (character).
#' @param incl.genes character vector. Set of genes. Gene pairs within these genes are included in the filtered data.
#' @param min.crossmap numeric. Minimum cross-mappability score.
#' @param max.crossmap numeric. Maximum cross-mappability score.
#' @return A filttered cross-mappability data frame.
#' @export
#' @details 
#' \code{filter_crossmap_by_genes} filters cross-mappability data to include pairs where each gene is within a given set of genes.
#' \code{filter_crossmap_by_score} filters cross-mappability data to include pairs where the cross-mappability score is within a given range.
#' \code{filter_crossmap_by_overlap} filters cross-mappability data to include pairs where the genes are not positionally overlapped.
filter_crossmap_by_genes <- function(crossmap_df, incl.genes){
  # crossmap_df : cross-mappability data frame where first two columns are genes
  stopifnot(class(crossmap_df)=='data.frame' && ncol(crossmap_df)>=2)
  stopifnot(class(crossmap_df[,1])=='character' && class(crossmap_df[,2])=='character')
  stopifnot(class(incl.genes) == 'character')
  
  gene_fac = as.factor(crossmap_df[,1])
  gene_levels = levels(gene_fac)
  gene_levels_included = gene_levels %in% incl.genes
  gene1_included = gene_levels_included[as.integer(gene_fac)]
  
  gene_fac = as.factor(crossmap_df[,2])
  gene_levels = levels(gene_fac)
  gene_levels_included = gene_levels %in% incl.genes
  gene2_included = gene_levels_included[as.integer(gene_fac)]
  
  return(crossmap_df[gene1_included & gene2_included, , drop=F])
}

#' @rdname filter_crossmap
#' @export
filter_crossmap_by_score <- function(crossmap_df, min.crossmap=-Inf, max.crossmap=Inf){
  # 3rd column contains cross-mappability
  stopifnot(class(crossmap_df)=='data.frame' && ncol(crossmap_df)>=3)
  stopifnot(class(crossmap_df[,1])=='character' && class(crossmap_df[,2])=='character')
  stopifnot(is.numeric(crossmap_df[,3]))
  
  if (min.crossmap > -Inf)
    crossmap_df = crossmap_df[crossmap_df[,3]>= min.crossmap, ]
  if (max.crossmap < Inf)
    crossmap_df = crossmap_df[crossmap_df[,3] <= max.crossmap, ]
  
  return(crossmap_df)
}

#' @rdname filter_crossmap
#' @export
filter_crossmap_by_overlap <- function(crossmap_df, overlap_df){
  
  stopifnot(class(crossmap_df)=='data.frame' && ncol(crossmap_df)>=2)
  stopifnot(class(crossmap_df[,1])=='character' && class(crossmap_df[,2])=='character')
  stopifnot(class(overlap_df)=='data.frame' && ncol(overlap_df)>=2)
  stopifnot(class(overlap_df[,1])=='character' && class(overlap_df[,2])=='character')
  
  dataframe_diff <- function(df1, df2){
    return(df1[!duplicated(rbind(df2, df1))[-seq_len(nrow(df2))], ])
  }
  
  df1 = crossmap_df[,1:2]
  df2 = overlap_df[,1:2]
  colnames(df2) = colnames(df1)
  filtered_df = dataframe_diff(df1, df2)
  
  df2 = overlap_df[,2:1]
  colnames(df2) = colnames(df1)
  filtered_df = dataframe_diff(filtered_df, df2)
  
  filtered_df = merge(crossmap_df, filtered_df)  # add other columns including mappability
  return(filtered_df)
}
