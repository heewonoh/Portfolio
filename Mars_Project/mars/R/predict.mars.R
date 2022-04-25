#' Predict Function
#'
#' @param object a mars object
#' @param newdata optional new data for which predictions are required
#' @param ... additional arguments to predict()--currently not used
#'
#' @return the corresponding matrix of basis functions
#' @family methods
#' @export
#'
#' @examples mm<- mars(y~x1+x2,dat=mars::marstestdata)
#' predict(mm,newdata=data.frame(x1=rnorm(100),x2=rnorm(100)))
predict.mars <- function(object,newdata,...) {
  if(missing(newdata) || is.null(newdata)) { #resolving null
    B <- as.matrix(object$B)
  }
  else {
    tt <- terms(object$formula,data = newdata)
    tt <- delete.response(tt)
    mf <- model.frame(tt,newdata)
    mt <- attr(mf, "terms")
    X <- model.matrix(mt, mf)[,-1] #removing intercept
    B <- make_B(X,object$Bfuncs)
  }
  beta <- object$coefficients
  drop(B %*% beta)
}

#B constructor for predictor
make_B <- function(X, Bfuncs) {
  len <- length(Bfuncs)
  ret <- init_B(nrow(X), len-1)
  for (i in 2:len) {
    temp <- 1
    for (j in 1:nrow(Bfuncs[[i]])) {
      temp <- temp*h(X[,Bfuncs[[i]][j,"v"]], Bfuncs[[i]][j,"s"], Bfuncs[[i]][j,"t"])
    }

    ret[,i] <- temp
  }
  return(as.matrix(ret))
}
