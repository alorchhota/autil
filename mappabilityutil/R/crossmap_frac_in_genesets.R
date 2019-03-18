#' Function to compute fraction of cross-mappable gene pairs
#' 
#' \code{crossmap_frac_in_genesets} computes fraction of cross-mappable gene pairs among all possible pairs of the given set of genes.
#' @param genes character vector. Set of genes.
#' @param crossmap_df Cross-mappability data frame. First two columns should contain genes (character) and the 3rd column should contain cross-mappability score (numeric) from the 1st gene to the 2nd gene.
#' @param mappability_df Mappability data frame. First columns should contain genes (character) and the second columnd should contain corresponding mappability (numeric).
#' @param overlap_df Positional overlap data. First two columns should contain genes (character) with positional overlap.
#' @param crossmap_threshold numeric. A pair of genes are cross-mappable if cross-mappability score between them (in either direction) is greater or equal to this score.
#' @return The fraction of cross-mappable gene pairs in a given gene set.
#' @details 
#' A pair of genes are called "cross-mappable" if cross-mappability score between them (in either direction) is greater or equal to a threshold score. 
#' If there are \code{n} genes in a set, there could have \code{n(n-1)/2} gene pairs. 
#' This function computes what fractio of all these pairs are cross-mappable.
#' Genes with mappability of NA are excluded from the gene set before computing the fraction,
#' as cross-mappability with those genes are undefined.
#' If positionally overlapped genes are provided (if overlap_df is not NULL), 
#' overlapped gene pairs are excluded from consideration.
#' @export
crossmap_frac_in_genesets <- function(genes, crossmap_df, mappability_df, overlap_df=NULL, crossmap_threshold=1e-6){
  # genes: vector of ensembl ids
  # crossmap_df: cross-mappability, 1st two columns are genes, 3rd column cross-mappability
  # mappability_df: 1st column gene name, 2nd column mappability
  # overlap_df: first two columns are genes, 3rd column is # ofoverlap positions
  
  ### stop for different reasons
  stopifnot(class(genes)=='character')
  stopifnot(class(crossmap_df) == 'data.frame' && is.character(crossmap_df[,1]) && is.character(crossmap_df[,2]) && is.numeric(crossmap_df[,3]))
  stopifnot(class(mappability_df) == 'data.frame' && is.character(mappability_df[,1]) && is.numeric(mappability_df[,2]))
  stopifnot(is.null(overlap_df) || ( class(overlap_df) == 'data.frame' && is.character(overlap_df[,1]) && is.character(overlap_df[,2])))
  
  
  ### remove genes with NA or unmeasures mappability
  rownames(mappability_df) = mappability_df[,1]
  genes = unique(genes)
  unmeasured_genes = setdiff(genes, mappability_df[,1])
  genes = setdiff(genes, unmeasured_genes)
  na_mappable_genes = genes[is.na(mappability_df[genes,2])]
  genes = setdiff(genes, na_mappable_genes)
  
  ### filter crossmap for genes and cross-mappability
  filtered_crossmap_df = filter_crossmap_by_genes(crossmap_df, incl.genes = genes)
  filtered_crossmap_df = filter_crossmap_by_crossmappability(filtered_crossmap_df, min.crossmap = crossmap_threshold)
  filtered_crossmap_df = directed_to_undirected_crossmap(filtered_crossmap_df, FUN = max)
  
  ### filter crossmap and bg gene-pairs for gene-overlap
  n_gene = length(genes)
  total_pairs = (n_gene * (n_gene-1)/2)
  if(!is.null(overlap_df)){
    filtered_overlap_df = filter_crossmap_by_genes(overlap_df, incl.genes = genes)
    filtered_overlap_df[,3] = 1 # put a number on the 3rd column to use directed_to_undirected_crossmap()
    filtered_overlap_df = directed_to_undirected_crossmap(filtered_overlap_df[,1:3], FUN = max)
    filtered_crossmap_df = filter_crossmap_by_overlap(filtered_crossmap_df, filtered_overlap_df)
    total_pairs = total_pairs - nrow(filtered_overlap_df)  # exclude overlapped pairs
  }
  
  ### compute crossmap fraction
  crossmap_frac = nrow(filtered_crossmap_df) / total_pairs
  
  return(list(crossmap_frac = crossmap_frac,
              genes = genes,
              unmapped_genes = c(unmeasured_genes, na_mappable_genes) ))
}
