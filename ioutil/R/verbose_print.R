#' Verbose print
#' 
#' This function prints message along with current time.
#' @param msg character. message to print.
#' @param verbose numeric. print message if verbose > 0.
#' @details 
#' This function is useful for debugging. Whether the messages should be printed or not can be controlled by a single paramter verbose. If verbose = TRUE, the messages are printed, otherwise not. 
#' In addition, the function prints the current time which could be helpful to determine total runtime.
#' @export
#' @examples 
#' verbose_print('this message will be printed with current time.', verbose = 1)
#' verbose_print('this message will NOT be printed.', verbose = 0)
verbose_print <- function(msg, verbose=1){
  if(verbose > 0){
    print(sprintf("[%s] %s", format(Sys.time(), "%D %T"), msg))
  }
}
