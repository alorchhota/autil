#' Function to read a data frame
#' 
#' Reads a file in table format and creates a data frame. It is very much like read.table() function, but it is faster than read.table().
#' @param fn character. File name.
#' @param sep character. The filed seprator character. Default: a tab character.
#' @param quote character. The set of quote characters. Default: "".
#' @param row.names logical. Does each row has a name. Default: TRUE.
#' @param stringsAsFactors logical. Should characters should be converted to factors? Default: FALSE.
#' @param check.names logical. Should column names be checked to be valid names? Default: FALSE
#' @param lessColInHeader logial. Does the header line have one column less than the following lines? Default: FALSE.
#' @param skip integer. Number of lines to skip before beginning to read data. Default: 0.
#' @details The default values of this function are set to read a standard table with both row names and column names
#' Note: the function may not work properly if the parameters are not appropriate for the input file. So, please use the function cautiously.
#' @return A data frame with data in the file.
#' @export
#' @examples
#' my_df = read_df('path/to/my/file')

read_df <- function(fn, sep = '\t',  header = T, quote = "", row.names=T, stringsAsFactors = F, check.names = F, lessColInHeader=F, skip = 0){
  require(data.table)
  if(header==T && lessColInHeader==T){
    header_line = readLines(fn, n = skip+1)[skip+1]
    headers = strsplit(header_line, split = sep)[[1]]
    skip = skip + 1
  }
  
  data_df = fread(fn, 
                  sep = sep,
                  header = header & !lessColInHeader, 
                  skip = skip,
                  quote = quote,    # not available in old package
                  stringsAsFactors = stringsAsFactors, 
                  check.names = check.names, 
                  data.table = FALSE)
  if (row.names == TRUE){
    rownames(data_df) = data_df[,1]
    data_df = data_df[,-1, drop=F]
  }
  if(header==T && lessColInHeader==T){
    colnames(data_df) = headers
  }
  return(data_df)
}