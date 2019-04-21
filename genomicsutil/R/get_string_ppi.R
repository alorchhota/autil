#' Protein-protein interactions (PPIs) within given set of genes.
#' 
#' This function returns PPIs available in STRING database within a given set of genes.
#' @param genes A character vector with gene symbols.
#' @param directed logical. Directed interactions?
#' @param string_version character. Version of STRING database.
#' @param species numeric. Species code for STRING (9606 for homo sapiens).
#' @param score_threshold numeric. Interactions in STRING with score >= score_threshold are considered true interactions.
#' @param string_dir character. Local path where STRING database will be stored locally so that it can be used offline. If \code{string_dir = ""}, a temporary directory is used.
#' @param max_homology_bitscore numeric. Maximum homology score. \code{0} to filter all interactions between homologous proteins, \code{Inf} not to filter any interaction.
#' @param benchmark_pathway character. Pathway type to benchmark. \code{NULL} to use all pathway types. Examples: \code{"KEGG", "REACTOME", "BIOCARTA", "Disease", "Pfam", "NCI", "ECOCYC".} See \code{pathwayType} parameter in \code{STRINGdb::benchmark_ppi()} function for details.
#' @details 
#' Note: all the genes are expected as HGNC gene symbols.
#' @return Returns a data frame with interactions between the given genes.
#' @export
#' @examples 
#' genes = c("TP53", "RBM3", "SF3", "LIM12", "MDM4", "TMEM160", "TP53BP2", "MDM2", "PDR", "MEG3", "EGFR")
#' interactions = get_string_ppi(genes)

get_string_ppi <- function(genes, directed = FALSE, string_version="10", species=9606, score_threshold=400, string_dir="", max_homology_bitscore=Inf, benchmark_pathway=NULL){
  
  require('STRINGdb')
  require('miscutil')
  string_db <- STRINGdb$new( version=string_version, species=species, score_threshold=score_threshold, input_directory=string_dir)
  
  # map gene names to string ids
  bg_df = data.frame(gene = genes)
  bg_string_df = string_db$map( bg_df, 'gene', removeUnmappedRows = TRUE, quiet = T)
  bg_string_ids = unique(bg_string_df$STRING_id)
  
  # get background interactions
  all_interactions_df = string_db$get_interactions(bg_string_ids)
  all_interactions_df = unique(all_interactions_df[,c('from','to', 'combined_score'), drop=F])
  colnames(all_interactions_df) = c('proteinA','proteinB','score')
  if(is.finite(max_homology_bitscore)){
    all_interactions_df = string_db$remove_homologous_interactions(interactions_dataframe = all_interactions_df, bitscore_threshold = max_homology_bitscore)
  }
  if(!is.null(benchmark_pathway)){
    all_interactions_df = string_db$benchmark_ppi(interactions_dataframe = all_interactions_df, pathwayType = benchmark_pathway, max_homology_bitscore = max_homology_bitscore)
    all_interactions_df = all_interactions_df[, c('proteinA', 'proteinB', 'score'), drop=F]
  }
  if(!directed)
    all_interactions_df = rminmax_df(all_interactions_df, col1 = 'proteinA', col2 = 'proteinB', unique = T)
  
  toret = merge(all_interactions_df, bg_string_df, by.x = 'proteinA', by.y = 'STRING_id')
  toret = merge(toret, bg_string_df, by.x = 'proteinB', by.y = 'STRING_id', suffixes = c('A','B'))
  toret = toret[,c('geneA','geneB','proteinA','proteinB','score'), drop = F]
  return(toret)
}
