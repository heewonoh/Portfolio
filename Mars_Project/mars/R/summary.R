


#' Summary of Mars Object
#' @description summary method for class "mars". Prints a summary of hinges that make up each basis function in the optimal model, along with the coefficient of said basis functions
#' @param object an object of class "mars", usually a result of a call to [mars]
#' @param  ... additional arguments to print.mars, currently unused
#'
#' @return a summary of the mars object
#' @family methods
#' @export
#'
#' @author Gareth Bennett, Heewon Oh, Jusung Lee
#' @examples
#' mm<-mars(y ~.,data=mars::marstestdata)
#' summary(mm)
#'

summary.mars <- function(object,...){
  print(object$call)

  df<-as.data.frame(object$coefficients)
  rownames(df)[1]<-"(Intercept)"
  colnames(df)[1]<-"Coefficients"
  print(df) ##prints coefficients and intercept of function
  summary(object$residuals)
  for(i in 2:length(object$Bfuncs)){
    for(j in 1:nrow(object$Bfuncs[[i]])){
      cat("Hinge value:",object$x_names[object$Bfuncs[[i]][j,2]],"split at value t:",object$Bfuncs[[i]][j,3],"\n")
    } #prints hinge value
    cat("\n")
  }
}

