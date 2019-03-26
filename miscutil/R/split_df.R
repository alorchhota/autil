#' Function to split data.
#' 
#' \code{split_df} splits a data frame or similar object randomly with specified fractions of columns (or rows) in each split.
#' @param x data.frame or similar object.
#' @param frac numeric vector. fraction of data in different parts. if \code{sum(frac)!=1}, thne \code{frac} will be normalized.
#' @param by.col logical. If True, split columns, otherwise split rows.
#' @return list of data frames with splitted data. 
#' if \code{frac} vector has names, so does the returned list.
#' @export
#' @examples
#' x = matrix(1:50, nrow = 5, dimnames = list(paste0('r',1:5),paste0('c',1:10)))
#' splitted_x = split_df(x, frac = c(train=0.6, validation=0.2, test=0.2))
split_df <- function(x, frac, by.col = T){
  stopifnot(all(frac>=0))
  stopifnot(sum(frac)>0)
  stopifnot(nrow(x) > 0 && ncol(x)>0)
  
  sample_breaks = c(0,floor(cumsum(frac / sum(frac) * ifelse(by.col==T, ncol(x), nrow(x)))))
  stopifnot(length(unique(sample_breaks)) == length(sample_breaks))
  
  if(by.col==T){
    random_cols = sample(seq_len(ncol(x)))
    splits = lapply(seq_len(length(frac)), function(i) return(x[,random_cols[(sample_breaks[i]+1):sample_breaks[i+1]], drop=F]))
  } else {
    random_rows = sample(seq_len(nrow(x)))
    splits = lapply(seq_len(length(frac)), function(i) return(x[random_rows[(sample_breaks[i]+1):sample_breaks[i+1]], , drop=F]))
  }
  names(splits) = names(frac)
  
  return(splits)
}
