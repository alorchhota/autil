#' Functions to read or write a feather file.
#' 
#' \code{read_feather_df} reads a feather file in a data frame. \code{write_feather_df} writes a data.frame to a feather file.
#' @param fn character. File name.
#' @param rownames.col integer or character. Column name or index containing row names. Default: NULL to avoid row names.
#' @param columns integer or character vector. Columns (indexes or names) to read. Default: NULL to read all columns.
#' @param df data.frame to write.
#' @param rownames.title \code{NULL} or character. If \code{NULL}, row names are not written. 
#' Otherwise, row names are written as a column with title \code{rownames.title}.
#' @return \code{read_feather_df} returns a data frame with data in the file. \code{write_feateher_df} returns nothing.
#' @export
#' @examples
#' my_df = read_feather_df(feather::feather_example("mtcars.feather"))
#' rownames(my_df) = paste0('car', 1:nrow(my_df))
#' write_feather_df(my_df, fn = "my_df.feather", rownames.title = "car_id")
#' my_df2 = read_feather_df("my_df.feather", rownames.col = "car_id")
#' @rdname read_write_feather
read_feather_df <- function(fn, columns = NULL, rownames.col=NULL){
  require(feather)
  if(!file.exists(fn))
    stop(sprintf('file does not exist: %s', fn))
  if(!(is.null(rownames.col) || is.numeric(rownames.col) || is.character(rownames.col)))
    stop('rownames.col must be an integer between 1 and the total number of columns, or a valid column name, or NULL.')
  
  data_df = as.data.frame(feather::read_feather(fn, columns = columns))
  
  if(!is.null(rownames.col)){
    if(is.numeric(rownames.col)){
      if(rownames.col < 1 || rownames.col > ncol(data_df))
        stop('rownames.col must be an integer between 1 and the total number of columns, or a valid column name, or NULL.')
      rowname.idx = as.integer(rownames.col)
    } else if (is.character(rownames.col)){
      if(! rownames.col %in% colnames(data_df)){
        stop('rownames.col must be an integer between 1 and the total number of columns, or a valid column name, or NULL.')
      }
      rowname.idx = which(colnames(data_df) == rownames.col)
    } else {
      stop('rownames.col must be an integer between 1 and the total number of columns, or a valid column name, or NULL.')
    }
    
    rownames(data_df) = data_df[,rowname.idx]
    data_df = data_df[,-rowname.idx, drop = F]
  }
  
  return(data_df)
}


#' @export
#' @rdname  read_write_feather
write_feather_df <- function(df, fn, rownames.title = NULL){
  require(feather)
  
  if(!(is.null(rownames.title) || is.character(rownames.title)))
    stop('rownames.title must be a character or NULL.')
  if(is.null(rownames.title)){
    write_feather(df, path = fn)
  } else {
    existing_cols = colnames(df)
    if(rownames.title %in% existing_cols)
      stop('rownames.title must be different from any existing column name.')
    df[,rownames.title] = rownames(df)
    write_feather(df[,c(rownames.title, existing_cols)], path = fn)
  }
  invisible()
}
