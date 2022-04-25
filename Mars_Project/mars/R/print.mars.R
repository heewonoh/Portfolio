#' Print a Mars Object
#'
#' @param x a `mars` object
#' @param ... additional arguments to print.mars, currently unused
#' @family methods
#'
#'
#' @export
#'
#' @examples mm <- mars(y~x1+x2,data=marstestdata)
#' print(mm)
print.mars<-function(x,...){
  print(x$call) #prints function called
  print(coefficients(x)) #prints coefficients and their values
  invisible(x)
}
