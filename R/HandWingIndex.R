#' Hand-wing Index
#' 
#' Calculates hand-wing index from Kipp's distance and wing length values.
#'
#' @param k Kipp's distance (distance between tip of the first secondary feather 
#' and the tip of the longest primary feather). Must be in same units as `l`.
#' @param l wing length. Must be in same units as `k`.
#'
#' @return numeric value
#' @export
#'
#' @examples
HandWingIndex <- function(k = NULL,l = NULL){
  if(!is.null(k) & !(is.null(l))){
    hwi <- 100*k/l 
    return(hwi)
  }else{
    message("HandWingIndex function: please provide valid arguments")
    return(NULL)
  }
}