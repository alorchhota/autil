#' Correct expression data for covariates
#' 
#' This function removes effects of covariates from a gene expression matrix using a linear model.
#' @param expr_df data.frame or matrix. Data matrix where rows and columns represent genes and samples, respectively.
#' @param cov_df data.frame. Covariate data where rows and columns represent samples and covariates, respectively.
#' @return Returns a corrected data matrix.
#' @examples 
#' y <- matrix( rpois(5000, lambda=5), ncol=100 )
#' cov <- data.frame(C1 = rnorm(100), C2 = sample(LETTERS[1:5], size = 100, replace = T) )
#' corrected_y = rm_cov_from_expr(y, cov)
#' @export
rm_cov_from_expr <- function(expr_df, cov_df){
  require(limma)
  stopifnot(ncol(expr_df) == nrow(cov_df))
  
  if(setequal(colnames(expr_df), rownames(cov_df))){
    cov_df = cov_df[colnames(expr_df),,drop=F]
  } else {
    warning('samples names are not same in expr_df and cov_df.')
  }
  
  if(!is.matrix(expr_df))
    expr_df = as.matrix(expr_df)  # gene x sample
  if(!is.data.frame(cov_df))
    cov_df = as.data.frame(cov_df)  # sample x cov
  
  cov_mat = model.matrix( ~ ., data=cov_df)
  fit = lmFit(expr_df, cov_mat)
  residual_mat = expr_df - fit$coefficients %*% t(cov_mat)
  return(residual_mat)
}
