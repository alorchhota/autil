#' Return the first or the last part of a data.frame
#' 
#' \code{head_df} returns top-left (when \code{pos = 'left'}) or top-right (when \code{pos = 'right'}) part of a data.frame or a matrix-like object.
#' Similarly, \code{tail_df} returns bottom-left or bottom-right part of a data frame.
#' @param x data frame or matrix-like object.
#' @param pos "left" or "right"
#' @param nrow numeric. number of rows.
#' @param ncol numeric. number of columns.
#' @return returns a data-frame or a matrix-like object (depending on the input) with head or tail of the data-frame.
#' @examples 
#' x = matrix(1:100, nrow = 10, byrow = F, dimnames = list(paste0('r',1:10), paste0('c',1:10)))
#' head_df(x)
#' head_df(x, 'right')
#' tail_df(x, nrow = 2)
#' tail_df(x, 'right', nrow = 2, ncol = 5)


#' @export
#' @rdname head_df
head_df <- function(x, pos = 'left', nrow = 3, ncol=4){
  stopifnot(pos %in% c('left', 'right'))
  stopifnot(is.numeric(nrow(x)))
  stopifnot(is.numeric(ncol(x)))
  stopifnot(nrow>0 && ncol>0)
  
  if(nrow(x) ==0){
    rows = NULL
  } else {
    rows = 1:min(nrow, nrow(x))
  }
  
  if(ncol(x) == 0){
    rows = NULL
  } else if(pos == 'left'){
    cols = 1:min(ncol, ncol(x))
  } else if(pos == 'right'){
    cols = max(1, ncol(x)-ncol+1):ncol(x)
  }
  
  return(x[rows, cols, drop=F])
}


#' @export
#' @rdname head_df
tail_df <- function(x, pos = 'left', nrow = 3, ncol=4){
  stopifnot(pos %in% c('left', 'right'))
  stopifnot(is.numeric(nrow(x)))
  stopifnot(is.numeric(ncol(x)))
  stopifnot(nrow>0 && ncol>0)
  
  if(nrow(x) ==0){
    rows = NULL
  } else {
    rows = max(1, nrow(x)-nrow+1):nrow(x)
  }
  
  if(ncol(x) == 0){
    rows = NULL
  } else if(pos == 'left'){
    cols = 1:min(ncol, ncol(x))
  } else if(pos == 'right'){
    cols = max(1, ncol(x)-ncol+1):ncol(x)
  }
  
  return(x[rows, cols, drop=F])
}
