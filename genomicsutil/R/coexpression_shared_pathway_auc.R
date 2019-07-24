#' AUC of a co-expression networks using known pathways (gene sets).
#' 
#' This function computes area under the curve (AUC) of a co-expression network using known pathways (gene sets).
#' @param net matrix. A gene x gene matrix where each entry represents the weight of the edge between corresponding genes. For a directed network, rows correspond to source genes and columns correspond to target genes. 
#' Row and column names of \code{net} should be HGNC gene symbols or gene ids. If they are not HGNC gene symbols, gene annotations (\code{gene_annot}) must be provided.
#' @param gene_sets list. List of pathways where each entry contains the genes (HGNC gene symbols) in each pathway.
#' @param gene_annot data.frame. Gene annotations, where the row names correspond to gene ids and one of it columns (\code{symbol_col}) contains HGNC gene symbols.
#' @param symbol_col character or integer. Name or index of the column containg HGNC gene symbol in \code{gene_annot}. 
#' @param directed logical. Directed interactions? If TRUE, it is assumed that the g1_col gene regulates the g2_col gene.
#' @details This function computes the area under the ROC curve and the precision-recall curve. 
#' The area is computed using the \code{PRROC} package considering an interaction between 
#' a pair of genes known if the genes have at least one pathway in common. 
#' Note: higher values in \code{net} should correspond to higher confidence in the edge. 
#' If the sign of values mean positive or negative association between genes, 
#' you probably should provide absolute values. 
#' 
#' This function does not yet handle directed networks.
#' 
#' @return Returns a list object with the following items.
#' \item{pr}{Precision curve object.}
#' \item{roc}{Roc curve object.}
#' \item{adjacency}{An adjacency matrix representing if a pair of genes have at least one common pathway.}
#' @export
#' @examples 
#' genes = c('TP53', 'RBM3', 'SF3', 'LIM12', 'ATM', 'TMEM160', 'BCL2L1', 'MDM2', 'PDR', 'MEG3', 'EGFR', 'CD96', 'KEAP1', 'SRSF1', 'TSEN2')
#' dummy_net = matrix(rnorm(length(genes)^2), nrow = length(genes), dimnames = list(genes, genes))
#' dummy_net = abs((dummy_net + t(dummy_net))/2) # symmetric undirected nework
#' gsea_fn = "your_gene_sets_file.gmt" # examples available at http://software.broadinstitute.org/gsea/msigdb/collections.jsp (download gene symbols GMT files)
#' gene_sets = read_gsea_db(gsea_fn)
#' auc_res = coexpression_shared_pathway_auc(net = dummy_net, gene_sets = gene_sets)
#' print(sprintf('AUC under ROC curve: %s', auc_res$roc$auc))
#' print(sprintf('AUC under PR curve: %s', auc_res$pr$auc.integral))
#' plot(auc_res$roc)
#' plot(auc_res$pr)
#' 

coexpression_shared_pathway_auc <- function(net, gene_sets, gene_annot = NULL, symbol_col = "gene_name", directed = FALSE){
  require('PRROC')
  
  if(!is.null(gene_annot))
    stop('coexpression_shared_pathway_auc() does not yet handle gene annotations. net should have HGNC gene symbols as row and column names.')
  
  if(directed == T)
    stop('coexpression_shared_pathway_auc() does not yet handle directed networks.')
  if(directed == F && !isSymmetric(net))
    stop('undirected network must be symmetric.')
  if(sum(net < 0))
    warning('the network has negative weight(s).')
  
  ### create shared pathway matrix
  shared_pathway_mat = matrix(0, nrow = nrow(net), ncol = ncol(net), dimnames = list(rownames(net), colnames(net)))
  bg = unique(c(rownames(net), colnames(net)))
  tmp <- lapply(gene_sets, function(genes){
    common_genes = intersect(genes, bg)
    if(length(common_genes) < 2)
      return(NULL)
    shared_pathway_mat[common_genes, common_genes] <<- 1
    return(NULL)
  })
  diag(shared_pathway_mat) <- 0
  
  ### compute AUC
  probj = pr.curve(scores.class0 = net[lower.tri(net, diag = F)], weights.class0 = shared_pathway_mat[lower.tri(net, diag = F)], curve = T)
  rocobj = roc.curve(scores.class0 = net[lower.tri(net, diag = F)], weights.class0 = shared_pathway_mat[lower.tri(net, diag = F)], curve = T)
  
  return(list(pr = probj, roc = rocobj, adjacency = shared_pathway_mat))
}
