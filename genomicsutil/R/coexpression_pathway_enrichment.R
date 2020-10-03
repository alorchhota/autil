#' Enrichment of pathways (gene sets) in a co-expression network.
#' 
#' \code{coexpression_pathway_enrichment} computes enrichment p-values of pathways in a co-expression network.
#' @param net matrix. A gene x gene square matrix where each entry represents the weight of the edge between corresponding genes.
#' @param gene_sets list. List of pathways where each entry contains the genes in each pathway. \code{names(gene_sets)} represents the pathway name and each name must be unique.
#' @param min.gene integer. Minimum number of genes from a pathway available in the network.
#' @param max.gene integer. Maximum number of genes from a pathway available in the network.
#' @param iter integer. Number of random gene sets to compute the null distribution.
#' @details To compute the enrichment p-value of a pathway in a co-expression network,
#' we define the score of a pathway as the the sum of weights in the co-expression network
#' of all possible edges between the genes in the pathway. The enrichment p-value is 
#' then defined as the probability that the score of the pathway is at least as big as
#' a random gene set with same number of genes. 
#' 
#' To get the null distribution, we select a number of (\code{iter}) random gene sets where 
#' each gene set consists of the same number of randomly selected genes, compute
#' their scores, and fit a normal distribution.
#' 
#' Enrichment of a pathway is computed only if at least \code{min.gene} and 
#' at most \code{max.gene} genes from the pathway are available in the network.
#' This cirteria helps to avoid too small or to large pathways.
#' 
#' Note: higher values in \code{net} should correspond to higher confidence in the edge. 
#' If the sign of the values mean positive or negative association between genes,
#' you probably should provide absolute values.
#' 
#' @return Returns a \code{data.frame} including with the following columns.
#' \item{pathway}{Pathway name.}
#' \item{n.gene}{Number of genes available in the network.}
#' \item{p}{p-value for the pathway.}
#' \item{p.empirical}{Empirical p-value for the pathway.}
#' \item{fdr}{False discovery rate computed using Benjamini-Hochberg method.}
#' \item{fdr.empirical}{Empirical false discovery rate computed using Benjamini-Hochberg method.}
#' 
#' @export
#' @examples 
#' genes = c('TP53', 'RBM3', 'SF3', 'LIM12', 'ATM', 'TMEM160', 'BCL2L1', 'MDM2', 'PDR', 'MEG3', 'EGFR', 'CD96', 'KEAP1', 'SRSF1', 'TSEN2')
#' dummy_net = matrix(rnorm(length(genes)^2), nrow = length(genes), dimnames = list(genes, genes))
#' dummy_net = abs((dummy_net + t(dummy_net))/2) # symmetric undirected nework
#' dummy_gene_sets = list(pathway1=c('TP53', 'RBM3', 'SF3', 'SF5'),
#'                        pathway2=c('LIM12', 'MDM2', 'BCL2L1', 'TMEM160'),
#'                        pathway3=c('EGFR', 'TP53', 'CD96', 'SRSF1', 'RBM14'))
#' enrich_res = coexpression_pathway_enrichment(net = dummy_net, gene_sets = dummy_gene_sets)
#' head(enrich_res)
#' n_sig = sum(enrich_res$fdr <= 0.05)
#' print(sprintf('Number of significanly enriched pathways: %d', n_sig))

coexpression_pathway_enrichment <- function(net, gene_sets, min.gene=3, max.gene=100, iter=10000){
  require(MASS)
  
  stopifnot(is.list(gene_sets))
  stopifnot(setequal(rownames(net), colnames(net)))
  if(sum(!is.finite(net)) > 0)
    stop('every edge weight in the network must be finite.')
  if(length(unique(names(gene_sets))) != length(gene_sets))
    stop('each pathway (gene_sets) must have a unique name.')
  if(sum(net < 0))
    stop('the network has negative weight(s).')
  
  ### diagonal weights should be 0 (needed to compute sum-of-weights)
  diag(net) = 0
  
  ### process given gene sets to contain a valid number of background genes only
  bg_genes = unique(c(rownames(net), colnames(net)))
  gene_sets = lapply(gene_sets, intersect, y=bg_genes)
  gene_sets_lengths = sapply(gene_sets, length)
  gene_sets= gene_sets[gene_sets_lengths>=min.gene & gene_sets_lengths <= max.gene]
  
  
  ### compute null distribution of sum-of-weights
  gene_sets_lengths = sapply(gene_sets, length)
  all_pvalues = lapply(unique(gene_sets_lengths), function(gslen){
    ### compute null distribution for given geneset length
    null_scores = sapply(1:iter, function(it){
      iter_genes = sample(bg_genes, size = gslen)
      weight_sum = sum(net[iter_genes, iter_genes])
      return(weight_sum)
    })
    
    null_dist =  fitdistr(null_scores, densfun = 'normal')
    null_mean = null_dist$estimate['mean']
    null_sd = null_dist$estimate['sd']
    
    ### compute pvalue for each gene set given length
    gene_set_names = names(gene_sets_lengths[gene_sets_lengths==gslen])
    gene_set_pvalues = sapply(gene_set_names, function(gsname){
      gs_genes = gene_sets[[gsname]]
      observed_score = sum(net[gs_genes, gs_genes])
      empirical_p = (sum(null_scores >= observed_score)+1)/(length(null_scores)+1)
      normal_p = 1 - pnorm(observed_score, mean = null_mean, sd = null_sd)
      return(list(pathway = gsname,
                  n.gene = gslen,
                  p = normal_p,
                  p.empirical = empirical_p))
      
    })
    
    ret_df = as.data.frame(t(gene_set_pvalues))
    ret_df$ngene = gslen
    return(ret_df)
  })
  
  all_pvalues_df = as.data.frame(do.call(rbind, args = all_pvalues))
  all_pvalues_df$fdr = p.adjust(all_pvalues_df$p, method = 'BH')
  all_pvalues_df$fdr.empirical = p.adjust(all_pvalues_df$p.empirical, method = 'BH')
  return(all_pvalues_df)
}

