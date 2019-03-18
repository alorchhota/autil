#' Enrichment of protein-protein interactions (PPIs)
#' 
#' This function computes enrichment of know PPIs in a given set of interactions.
#' @param interaction_df A data.frame with gene-gene interactions. It must have at least 2 columns containing gene symbols.
#' @param bg A character vector with background gene symbols. It should contain all genes in the interaction data frame and other genes analyzed to discover the interactions.
#' @param directed logical. Directed interactions? If TRUE, it is assumed that the g1_col gene regulates the g2_col gene.
#' @param g1_col character. Column name of interaction_df that contains one gene of an interaction.
#' @param g2_col character. Column name of interaction_df that contains another gene of an interaction.
#' @param n_iter numeric. Number of random iterations to compute the background distribution.
#' @param string_version character. Version of STRING database.
#' @param species numeric. Species code for STRING (9606 for homo sapiens).
#' @param score_threshold numeric. Interactions in STRING with score >= score_threshold are considered true interactions.
#' @param string_dir character. Local path where STRING database will be stored locally so that it can be used offline.
#' @details This function computes enrichment of known protein-protein interactions (PPIs) in the given set of interactions. 
#' It computes enrichment p-value as the probably that the number of known PPIs in the given set 
#' is greater than or equal to that in a random set of same size.
#' 
#' Note: all the genes (both in \code{interaction_df} and \code{bg}) are expected as HGNC gene symbols.
#' @return Returns a list object with the following items.
#' \item{pvalue}{numeric, enrichment p-value}
#' \item{random_known_ppi}{numeric vector, empirical number of known PPIs.}
#' \item{observed_known_ppi}{numeric, observed number of known PPIs in the given interaction set.}
#' @export
#' @examples 
#' gene1s = c('TP53', 'TP53', 'RBM3', 'SF3', 'LIM12', 'MDM4', 'TMEM160')
#' gene2s = g2=c('TP53BP2', 'MDM2', 'PDR', 'MEG3', 'TP53', 'MDM2', 'EGFR')
#' interaction_df = data.frame(g1 = gene1s, g2 = gene2s, stringsAsFactors = F)
#' bg = unique(c(gene1s, gene2s, c('CD96', 'KEAP1', 'SRSF1', 'TSEN2')))
#' enrich_res = string_ppi_enrichment(interaction_df, bg)

string_ppi_enrichment <- function(interaction_df, bg, directed = FALSE, g1_col=colnames(interaction_df)[1], g2_col=colnames(interaction_df)[2], n_iter=1000, string_version="10", species=9606, score_threshold=400, string_dir=""){
  require('STRINGdb')
  require('miscutil')
  
  stopifnot(g1_col %in% colnames(interaction_df) && g2_col %in% colnames(interaction_df))
  stopifnot(n_iter > 0)
  
  string_db <- STRINGdb$new( version=string_version, species=species, score_threshold=score_threshold, input_directory=string_dir)
  
  # map gene names to string ids
  bg_df = data.frame(gene = bg)
  bg_string_df = string_db$map( bg_df, 'gene', removeUnmappedRows = TRUE, quiet = T)
  bg_string_ids = unique(bg_string_df$STRING_id)
  
  mapped_interaction_df <- merge(interaction_df, bg_string_df, by.x = g1_col, by.y = 'gene')
  mapped_interaction_df <- merge(mapped_interaction_df, bg_string_df, by.x = g2_col, by.y = 'gene', suffixes = c('_1', '_2'))
  if(!directed){
    mapped_interaction_df = rminmax_df(mapped_interaction_df[,c('STRING_id_1', 'STRING_id_2')])
  } 
  mapped_interaction_df = unique(mapped_interaction_df[,c('STRING_id_1', 'STRING_id_2')])
  
  
  # get background interactions
  all_interactions_df = string_db$get_interactions(bg_string_ids)
  all_interactions_df = unique(all_interactions_df[,c('from','to'), drop=F])
  if(!directed)
    all_interactions_df = rminmax_df(all_interactions_df[,c('from', 'to'), drop=F], unique = T)
  all_interactions_srings = paste(all_interactions_df$from, all_interactions_df$to, sep = '|')
  
  # compute empirical p-value
  get_num_known_ppi <- function(cur_df, col1, col2, sorted_minmax=F){
    if(!sorted_minmax)
      cur_df = rminmax_df(cur_df[,c(col1, col2)])
    cur_interactions_string = paste(cur_df[,col1], cur_df[,col2], sep = '|')
    n_known_ppi = length(intersect(cur_interactions_string, all_interactions_srings))
    return(n_known_ppi)
  }
  
  total_interactions = nrow(mapped_interaction_df)
  
  generate_random_interactions <- function(bg_genes, n_interections, directed){
    bg_size <- length(bg_genes)
    
    if(directed){
      bg_int = matrix(1:(bg_size^2), nrow = bg_size)[!diag(bg_size)]
    } else {
      bg_int = matrix(1:(bg_size^2), nrow = bg_size)
      bg_int = bg_int[upper.tri(bg_int, diag = F)]
    }
    
    vals <- sample(bg_int, size = n_interections, replace = F)
    rand1 = bg_genes[vals %% bg_size]
    rand2 = bg_genes[vals %/% bg_size + 1]
    rand_df = data.frame(gene1 = rand1, gene2 = rand2, stringsAsFactors = F)
    return(rand_df)
  }
  
  random_known_ppi = sapply(1:n_iter, function(iter_idx){
    rand_df = generate_random_interactions(bg_string_ids, n_interections = total_interactions, directed = directed)
    n_rand_ppi = get_num_known_ppi(rand_df, col1 = 'gene1', col2 = 'gene2', sorted_minmax = F)
    return(n_rand_ppi)
  })
  
  observed_known_ppi = get_num_known_ppi(mapped_interaction_df, col1 = 'STRING_id_1', col2 = 'STRING_id_2', sorted_minmax = F)
  n_big_rand = sum(random_known_ppi >= observed_known_ppi)
  pvalue = (n_big_rand + 1) / (n_iter + 1)
  
  return(list(pvalue = pvalue, random_known_ppi = random_known_ppi, observed_known_ppi = observed_known_ppi))
  
}
