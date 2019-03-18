#' Undirected cross-mappability
#' 
#' \code{directed_to_undirected_crossmap} generates undirected cross-mappability scores from directed cross-mappability scores.
#' @param crossmap_df Directed cross-mappability data frame. First two columns should contain genes (character) and the 3rd column should contain cross-mappability score (numeric) from the 1st gene to the 2nd gene.
#' @param FUN function to summarize mappability scores from a vector of scores.
#' @return A cross-mappability data frame.
#' @details In a directed cross-mappability data frame,  cross-mappability score from gene A to gene B can be different from cross-mappability from gene B to gene A. This function generates an undirected version of cross-mappability score for each gene pair using a function (min, mean, max etc.).
#' @export
directed_to_undirected_crossmap <- function(crossmap_df, FUN=min){
  require(data.table)
  
  stopifnot(class(crossmap_df)=='data.frame' && ncol(crossmap_df)>=3)
  stopifnot(class(crossmap_df[,1])=='character' && class(crossmap_df[,2])=='character')
  stopifnot(is.numeric(crossmap_df[,3]))
  
  genes = sort(unique(c(crossmap_df[,1], crossmap_df[,2])))
  genes1 = factor(crossmap_df[,1], levels = genes)
  genes2 = factor(crossmap_df[,2], levels = genes)
  min_idx = pmin(as.integer(genes1), as.integer(genes2))
  max_idx = pmax(as.integer(genes1), as.integer(genes2))
  
  crossmap_dt = data.table(gene1=genes[min_idx], gene2=genes[max_idx], crossmap = crossmap_df[,3])
  symmetric_crossmap_dt = crossmap_dt[ ,.(symmetric_crossmap=FUN(crossmap)), by=.(gene1, gene2)]
  original_colnames = colnames(crossmap_df)
  crossmap_df = as.data.frame(symmetric_crossmap_dt)
  colnames(crossmap_df) = original_colnames
  return(crossmap_df)
}
