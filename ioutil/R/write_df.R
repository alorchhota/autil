#' Function to write a data frame into a file
#' 
#' This function is wrapper of the write.table() function with default values set to write a tab-delimited table with both row names and column names.
#' @param x data.frame
#' @param file either a character string or a connection open for writing.
#' @param sep character. The separator character.
#' @param quote logical. Should character or factor values be surrounded by double quotes?
#' @param row.names either a logical value indicating whether the row names of x are to be written along with x, or a character vector of row names to be written.
#' @param col.names either a logical value indicating whether the column names of x are to be written along with x, or a character vector of column names to be written. If col.names = NA and row.names = TRUE a blank column name is added.
#' @export
write_df <- function(x, file, sep = "\t", quote = F, row.names = T, col.names = NA){
  write.table(x = x, file = file, sep = sep, quote = quote, row.names = row.names, col.names = col.names)
}
