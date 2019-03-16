#' Function to read GSEA gene sets
#' 
#' This function reads a file with gene sets in gmt format. Representive examples are found \href{http://software.broadinstitute.org/gsea/msigdb/collections.jsp}{here}.
#' @param fn either a character string or a connection object to read.
#' @return Returns a list of gene sets. Each item in the list is a character vector representing a gene set. The name of each gene set is available as the name of each item in the list. 
#' Note: the 2nd column in the file containing an url of the gene set is ignored.
#' @export
read_gsea_db <- function(fn){
  gsea_lines = readLines(fn)
  parssed_lines = strsplit(gsea_lines, split = '\t')
  gs_titles = sapply(parssed_lines, function(x) x[1])
  genesets = lapply(parssed_lines, function(x) x[-(1:2)])
  names(genesets) = gs_titles
  return(genesets)
}
