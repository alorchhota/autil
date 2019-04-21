#' Map genes to STRING protein IDs.
#' 
#' This function maps HGNC gene symbols to protein ids in STRING.
#' @param genes A character vector with HGNC gene symbols.
#' @param string_version character. Version of STRING database.
#' @param species numeric. Species code for STRING (9606 for homo sapiens).
#' @param string_dir character. Local path where STRING database will be stored locally so that it can be used offline. If \code{string_dir = ""}, a temporary directory is used.
#' @return a character vector with mapped STRING ids.
#' @export
#' @examples 
#' genes = c('TP53', 'RBM3', 'SF3', 'LIM12', 'ATM', 'TMEM160', 'BCL2L1')
#' string_ids = map_genes_to_string_ids(genes)
map_genes_to_string_ids <- function(genes, string_version="10", species=9606, string_dir=""){
  require('STRINGdb')

  string_db <- STRINGdb$new(version=string_version, species=species, input_directory=string_dir)
  
  # map gene names to string ids
  genes_df = data.frame(gene = unique(genes), stringsAsFactors = F)
  genes_string_df = string_db$map(genes_df, 'gene', removeUnmappedRows = F, quiet = T)
  rownames(genes_string_df) = genes_string_df$gene
  string_ids = genes_string_df[genes, 'STRING_id']
  names(string_ids) = genes
  return(string_ids)
}
