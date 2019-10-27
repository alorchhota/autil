#' AUC of a co-expression networks using STRING protein-protein interactions (PPIs).
#' 
#' This function computes area under the curve (AUC) of a co-expression network using protein-protein interactions (PPIs) from STRING.
#' @param net matrix. A gene x gene matrix where each entry represents the weight of the edge between corresponding genes. For a directed network, rows correspond to source genes and columns correspond to target genes. 
#' Row and column names of \code{net} should be HGNC gene symbols or gene ids. If they are not HGNC gene symbols, gene annotations (\code{gene_annot}) must be provided.
#' @param gene_annot data.frame. Gene annotations, where the row names correspond to gene ids and one of it columns (\code{symbol_col}) contains HGNC gene symbols.
#' @param symbol_col character or integer. Name or index of the column containg HGNC gene symbol in \code{gene_annot}. 
#' @param directed logical. Directed interactions? If TRUE, it is assumed that the g1_col gene regulates the g2_col gene.
#' @param string_version character. Version of STRING database.
#' @param species numeric. Species code for STRING (9606 for homo sapiens).
#' @param string_dir character. Local path where STRING database will be stored locally so that it can be used offline. If \code{string_dir = ""}, a temporary directory is used.
#' @param max_homology_bitscore numeric. Maximum homology score. \code{0} to filter all interactions between homologous proteins, \code{Inf} not to filter any interaction.
#' @param benchmark_pathway character. Pathway type to benchmark. \code{NULL} to use all pathway types. Examples: \code{"KEGG" / "REACTOME" / "BIOCARTA" / "Disease" / "Pfam" / "NCI" / "ECOCYC"}
#' @details This function computes the area under the ROC curve and the precision-recall curve. The area is computed using the \code{PRROC} package using weighted PPI scores. 
#' Note: higher values in \code{net} should correspond to higher confidence in the edge. If the sign of values mean positive or negative association between genes, you probably should provide absolute values. 
#' 
#' This function does not yet handle directed networks.
#' 
#' @return Returns a list object with the following items.
#' \item{pr}{Precision curve object.}
#' \item{roc}{Roc curve object.}
#' @export
#' @examples 
#' require('reshape2')
#' genes = c('TP53', 'RBM3', 'SF3', 'LIM12', 'ATM', 'TMEM160', 'BCL2L1', 'MDM2', 'PDR', 'MEG3', 'EGFR', 'CD96', 'KEAP1', 'SRSF1', 'TSEN2')
#' dummy_net = matrix(rnorm(length(genes)^2), nrow = length(genes), dimnames = list(genes, genes))
#' dummy_net = (dummy_net + t(dummy_net))/2 # symmetric undirected nework
#' auc_res = coexpression_ppi_auc(net = dummy_net)
#' print(sprintf('AUC under ROC curve: %s', auc_res$roc$auc))
#' print(sprintf('AUC under PR curve: %s', auc_res$pr$auc.integral))
#' plot(auc_res$roc)
#' plot(auc_res$pr)
coexpression_ppi_auc <- function(net, gene_annot = NULL, symbol_col = "gene_name", directed = FALSE, string_version="10", species=9606, string_dir="", max_homology_bitscore=Inf, benchmark_pathway=NULL){
  require('PRROC')
  require('reshape2')
  
  if(!is.null(gene_annot))
    stop('coexpression_ppi_auc() does not yet handle gene annotations. net should have HGNC gene symbols as row and column names.')
  
  if(directed == T)
    stop('coexpression_ppi_auc() does not yet handle directed networks.')
  if(directed == F && !isSymmetric(net))
    stop('undirected network must be symmetric.')
  if(sum(net < 0))
    warning('the network has negative weight(s).')
  
  # take only genes mappable to string
  bg = unique(c(rownames(net), colnames(net)))
  bg_string_ids = map_genes_to_string_ids(bg, string_version = string_version, species = species, string_dir = string_dir)
  bg = names(bg_string_ids[!is.na(bg_string_ids)])
  net = net[intersect(rownames(net), bg), intersect(colnames(net),bg)]
  
  # get ppi scores
  known_ppis_df = get_string_ppi(genes = bg, directed = directed, string_version = string_version, species = species, score_threshold = 0, max_homology_bitscore = max_homology_bitscore, benchmark_pathway = benchmark_pathway, string_dir = string_dir)
  known_ppi_score_mat = acast(data = known_ppis_df, formula = geneA ~ geneB, value.var = "score", fun.aggregate = max, na.rm = T)
  known_ppi_score_mat[is.infinite(known_ppi_score_mat)] = NA
  
  ppi_score_mat = matrix(NA, nrow = nrow(net), ncol = ncol(net), dimnames = list(rownames(net), colnames(net)))
  ppi_score_mat[rownames(known_ppi_score_mat), colnames(known_ppi_score_mat)] = known_ppi_score_mat
  if(!directed)
    ppi_score_mat = pmax(ppi_score_mat, t(ppi_score_mat), na.rm = T)  # made symmetric matrix
  ppi_score_mat[is.na(ppi_score_mat)] = 0
  ppi_score_mat = ppi_score_mat / 1000.0 # converted to probability
  
  probj = pr.curve(scores.class0 = net[lower.tri(net, diag = F)], weights.class0 = ppi_score_mat[lower.tri(net, diag = F)], curve = T)
  rocobj = roc.curve(scores.class0 = net[lower.tri(net, diag = F)], weights.class0 = ppi_score_mat[lower.tri(net, diag = F)], curve = T)
  
  return(list(pr = probj, roc = rocobj))
}
