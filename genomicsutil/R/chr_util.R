#' Functions to handle chromosome names
#' 
#' \code{extend_chr} function adds a prefeix to a chromosome name if the prefix does not exist.  
#' \code{shorten_chr} function removes a prefix, if exists, from a chromosome names.  
#' 
#' @param chr character vector. chromosome names.
#' @param pfx character. prefix of a chromosome name.
#' @return character vector of (extended or shortened) chromosome names.
#' @export
#' @rdname chr_name
#' @examples 
#' chromosomes = c("1", "chr1", "5", "chr4", "chrX", "Y")
#' extend_chr(chromosomes)
#' shorten_chr(chromosomes)
#' 
extend_chr <- function(chr, pfx='chr'){
  if (!is.character(chr)) 
    chr <- as.character(chr)
  prefixes = substr(chr, 1, 3)
  chr[prefixes!=pfx] = paste0(pfx, chr[prefixes!=pfx])
  return(chr)
}

#' @rdname chr_name
#' @export
shorten_chr <- function(chr, pfx='chr'){
  if (!is.character(chr))
    chr <- as.character(chr)
  prefixes = substr(chr, 1, 3)
  chr[prefixes==pfx] = substr(chr[prefixes==pfx], 4, nchar(chr[prefixes==pfx]))
  return(chr)
}
