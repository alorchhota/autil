% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rm_cov_from_expr.R
\name{rm_cov_from_expr}
\alias{rm_cov_from_expr}
\title{Correct expression data for covariates}
\usage{
rm_cov_from_expr(expr_df, cov_df)
}
\arguments{
\item{expr_df}{data.frame or matrix. Data matrix where rows and columns represent genes and samples, respectively.}

\item{cov_df}{data.frame. Covariate data where rows and columns represent samples and covariates, respectively.}
}
\value{
Returns a corrected data matrix.
}
\description{
This function removes effects of covariates from a gene expression matrix using a linear model.
}
\examples{
y <- matrix( rpois(5000, lambda=5), ncol=100 )
cov <- data.frame(C1 = rnorm(100), C2 = sample(LETTERS[1:5], size = 100, replace = T) )
corrected_y = rm_cov_from_expr(y, cov)
}
