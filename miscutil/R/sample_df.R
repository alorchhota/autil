#' Random sampling of a data frame
#' 
#' \code{sample_df} takes a specified number of rows and/or columns from \code{x} either with or without replacement.
#' @param x data.frame or similar object.
#' @param size.row numeric. number of rows if \code{unit == 'total'}, otherwise percentage of rows.
#' @param size.col numeric. number of columns if \code{unit == 'total'}, otherwise percentage of columns.
#' @param unit either 'total' or 'percent'.
#' @param replace.row logical. should rows be sampled with replacement?
#' @param replace.col logical. should columns be sampled with replacement?
#' @return a data frame with sampled rows and/or columns.
#' @export
#' @examples
#' x = matrix(1:500, nrow = 50, dimnames = list(paste0('r',1:50),paste0('c',1:10)))
#' x1 = sample_df(x, size.row = 4)
#' x2 = sample_df(x, size.row = 4, unit = 'percent')
#' x3 = sample_df(x, size.row = 4, size.col = 5, unit = 'total')
sample_df <- function(x, size.row=ifelse(unit=='total', nrow(x), 100), size.col=ifelse(unit=='total', ncol(x), 100), unit='total', replace.row = F, replace.col = F){
  stopifnot(nrow(x) > 0 && ncol(x)>0)
  stopifnot(unit %in% c('total', 'percent'))
  stopifnot(size.row > 0 && size.col>0)
  stopifnot( !(replace.row == F && unit == 'total') || size.row <= nrow(x))
  stopifnot( !(replace.col == F && unit == 'total') || size.col <= ncol(x))
  stopifnot( !(replace.row == F && unit == 'percent') || size.row <= 100)
  stopifnot( !(replace.col == F && unit == 'percent') || size.col <= 100)
  
  nr = ifelse(unit == 'total', size.row, round(nrow(x) * size.row / 100))
  nc = ifelse(unit == 'total', size.col, round(ncol(x) * size.col / 100))

  if(nr < nrow(x)){
    rows = sample(x=1:nrow(x), size = nr, replace = replace.row)
  } else {
    rows = 1:nrow(x)
  }
  
  if(nc < ncol(x)){
    cols = sample(x=1:ncol(x), size = nc, replace = replace.col)
  } else {
    cols = 1:ncol(x)
  }
  
  sampled_df = x[rows, cols, drop=F]
  
  return(sampled_df)
}
