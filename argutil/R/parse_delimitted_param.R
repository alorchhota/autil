#' Function to parse delimitted parameter
#' 
#' This function splits the given string by a delimiter. It can also remove empty strings.
#' @param param_str character. String to be parsed.
#' @param delim character. The delimiter.
#' @param rm.empty logical. If TRUE, removes every empty string after parsing.
#' @return A vector of parsed character strings.
#' @export
#' @examples 
#' x = "one,two,three,"
#' parsed_strings = parse_delimitted_param(x, delim = ',', rm.empty = T)
parse_delimitted_param <- function(param_str, delim=',', rm.empty=T){
  stopifnot(class(param_str)=='character' && length(param_str)==1)
  parts = strsplit(param_str, split = delim)[[1]]
  if(rm.empty==T){
    is_valid_parts = sapply(parts, function(s) nchar(s)>0)
    parts = parts[is_valid_parts]
  }
  return(parts)
}
