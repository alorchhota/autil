#' Transformation data to inverse normal distribution
#' 
#' This function transforms values in each row of a matrix to a rank based inverse normal transformation.
#' @param expr_mat matrix. Data matrix.
#' @return Returns a transformed matrix.
#' @examples 
#' x = matrix(1:100, nrow=5, byrow=T, dimnames = list(paste0('gene', 1:5), paste0('sample', 1:20)))
#' inv_norm_x = to_inv_normal(x)
#' @export
to_inv_normal <- function(expr_mat){
  stopifnot(is.matrix(expr_mat))
  # expr_mat : gene x sample matrix
  transformed_data <- apply(expr_mat, 1, function(x){
    qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
  })
  transformed_data = t(transformed_data)
  return(transformed_data)
}
