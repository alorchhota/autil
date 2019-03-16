#' TMM-normalization of RNA-seq data
#' 
#' This function normalizes a gene expression matrix using TMM.
#' @param expr_mat matrix. Data matrix where rows and columns represent genes and samples, respectively.
#' @return Returns a transformed matrix.
#' @examples 
#' y <- matrix( rpois(1000, lambda=5), nrow=200 )
#' tmm_y = tmm_normalize(y)
#' @export
#' @importFrom Rdpack reprompt
#' @references 
#' \insertRef{Robinson2010}{genomicsutil}
tmm_normalize <- function(expr_mat){
  require('edgeR')
  
  # expr_mat: gene x sample matrix
  stopifnot(class(expr_mat) == 'matrix')
  tmm_factors = calcNormFactors(expr_mat, method = "TMM")
  expr_mat = t(t(expr_mat)/tmm_factors)
  return(expr_mat)
}
