#' Multivariate Adaptive Regression Splines (MARS)
#' @description Fit Friedman's Multivariate Adaptive Regression Splines (MARS) model.
#'
#' @usage mars(formula, data, control = NULL)
#' @param formula an R formula
#' @param data a data frame containing your data
#' @param control an optional object of class 'mars.control'
#'
#' @details The MARS function attempts to fit non-linear data as accurately as possible. To do this the function defines one or more hinge-points that act as connectors between two or more linear functions.This process is completed through a multistage process including a forward pass and backwards pass.
#'
#' The forward pass algorithm takes a formula and data as input and attempts to output a matrix of basis functions of size Mmax. To do this, the algorithm considers each data point for each predictor and generates candidate basis function pairs. These candidates are selected based on reduction to the model's residual sum-of-squares.
#' The basis function matrix contains a row for each observations and a column for each basis function.
#'
#' The forward pass identifies every basis function that reduces the model's residual sum-of-squares, but this can produce an over-fit model. To create a model that generalized better with new data, the matrix is handed to the backwards pass algorithm. The algorithm 'prunes' any basis function that adds no significant predictive accuracy to the overall model.
#' The backwards algorithm compares subsets of basis function using generalized cross-validation scores in which 'd' controls the amount of penalization there is for extra model terms. By tracking GCV scores for each possible combination of basis function, the best possible model can be selected.
#' The backwards pass algorithm updates the basis matrix to contain only the function from the optimal model.

#' @return an object of class 'mars' which can be passed to plot,predict,summary and print methods.
#' @export
#'
#'
#' @examples
#' mm<-mars(y~.,dat=mars::marstestdata)
#' @import stats
#' @author Gareth Bennett, Heewon Oh, Jusung Lee
#' @references Jerome H. Friedman. Multivariate Adaptive Regression Splines (with discussion).Annals of Statistics 19/1, 1991. \url{https://statistics.stanford.edu/research/multivariate-adaptive-regression-splines.}
#' @seealso [mars.control] for constructing control objects
#' @seealso [plot.mars] for plotting results
#' @seealso [predict.mars] for predictions
#' @seealso [summary.mars] for summarizing mars objects
#' @seealso [print.mars] for printing mars objects

# MARS FUNCTION


mars <- function(formula,data,control=NULL){
  cc <- match.call() # save the call
  mf <- model.frame(formula,data)
  y <- model.response(mf)
  mt <- attr(mf,"terms")
  x <- model.matrix(mt,mf)[,-1,drop=FALSE]
  if(is.null(control))control<-mars.control() # validate control object
  control<-validate_mars.control(control)
  fwd <- fwd_stepwise(y,x,control)  # fwd
  bwd <- bwd_stepwise(fwd,control)  # bwd takes in fwd output
  mdl <- lm(y~.-1, data = data.frame(y = y, bwd$B)) # create model with bwd and drop intercept
  out <- c(list(call = cc, formula = formula, y = y, B = bwd$B, Bfuncs = bwd$Bfuncs, x_names=colnames(x)), mdl)
  class(out) <- c("mars", class(mdl))
  out
}


# FWD STEPWISE

fwd_stepwise <- function(y,x,control){
  # init constants
  Mmax <- control$Mmax
  N <- length(y)
  names(y) <- c(paste0(1:N))
  n <- ncol(x)
  B <- init_B(N,Mmax) # create B
  Bfuncs <- vector(mode="list", length=Mmax+1) # create Bfuncs
  splits <- data.frame(m=rep(NA,Mmax),v=rep(NA,Mmax),t=rep(NA,Mmax)) # create data frame for best splits
  for(i in 1:(Mmax/2)) {
    M <- 2*i - 1
    lof_best <- Inf   # to determine best split
    for(m in 1:M) {
      diff <- setdiff(1:n, Bfuncs[[m]][,"v"]) # remove best v
      for(v in diff){
        tt <- split_points(x[,v],B[,m]) # find split points
        for(t in tt) {
          Bnew <- data.frame(B[,1:M],
                             Btem1=B[,m]*h(x[,v], +1, t),
                             Btem2=B[,m]*h(x[,v], -1, t))
          gdat <- data.frame(y=y,Bnew)
          lof <- lof(y~.-1,gdat)  # remove intercept
          if(lof < lof_best) {
            lof_best <- lof;
            split_best <- c(m=m,v=v,t=t) # record current best split
          }
        }
      }

    }

    m <- split_best["m"]; v <- split_best["v"]; t <- split_best["t"];
    # use best split index to record B and Bfuncs
    B[,M+1] <- B[,m]*h(x[,v], -1,t)
    B[,M+2] <- B[,m]*h(x[,v], +1,t)

    Bfuncs[[M+1]] <- rbind(Bfuncs[[m]], c(s=-1,v,t))
    Bfuncs[[M+2]] <- rbind(Bfuncs[[m]], c(s=+1,v,t))
  }

  colnames(B) <- paste0("B",(0:(ncol(B)-1)))
  return(list(y=y,B=B, Bfuncs=Bfuncs))
}


# BWD STEPWISE
bwd_stepwise <- function(fwd, control) {
  Mmax <- ncol(fwd$B) - 1
  JStar <- 2:(Mmax+1) # skip intercept
  KStar <- JStar
  dat <- data.frame(y = fwd$y, fwd$B)
  LOFStar <- LOF(y~., dat, control) # obtain GCV from LOF
  for(M in (Mmax+1):2) {
    b <- Inf
    L <- KStar
    for(m in L) {
      K <- setdiff(L,m) # skip current m
      dat <- data.frame(y = fwd$y, fwd$B[,K])
      lof <- LOF(y~., dat, control) # obtain GCV for K in B
      if(lof < b) { # record current best lof
        b <- lof
        KStar <- K
      }
      if(lof < LOFStar) { # record best lof
        LOFStar <- lof
        JStar <- K
      }
    }
  }
  JStar <- c(1,JStar)
  return(list(y = fwd$y, B = fwd$B[,JStar], Bfuncs = fwd$Bfuncs[JStar]))
}

# Initialize B
init_B <- function(N,Mmax) {
  B <- data.frame(matrix(NA,nrow=N,ncol=(Mmax+1)))
  B[,1] <- 1
  names(B) <- c("B0",paste0("B",1:Mmax))
  return(B)
}

# FWD lof
lof <- function(form,data) {
  ff <- lm(form,data)
  return(sum(residuals(ff)^2)) # return sum of residuals squared
}

# BWD LOF
LOF <- function(form,data, marsControl) {
  # formula from tutorial 8 in Stat 361
  ff <- lm(form,data)
  RSS <- sum(residuals(ff)^2)
  N <- nrow(data)
  NC <- length(coefficients(ff))-1
  CM <- sum(diag(hatvalues(ff)))
  GCV <- RSS * N / (N - (CM + marsControl$d*NC))^2
  return (GCV)
}

# Hinge Function
h <- function(x,s,t) {
  # return 0 or positive
  return(pmax(0, s*(x-t)))
}


#Split points
split_points <- function(xv,Bm) {
  # return a filtered list excluding last element
  out <- sort(unique(xv[Bm>0]))
  return(out[-length(out)])
}



#' Mars Control Object
#' @description Constructor for `mars.control` objects
#'
#' This function constructs a `mars.control` object that specifies
#' parameters used in the model fitting procedure
#'
#' @param Mmax Maximum number of basis functions. Should be an even interger Default value is 1.
#' @param d The coefficient in the penalty term of the generalized cross validation measure. Default is 3.
#' @param trace Should we print information about the fitting? Default is `FALSE`
#'
#' @return a `mars.control` object
#' @export
#'
#' @examples mc<-mars.control(Mmax=10)


mars.control <- function(Mmax=2, d=3, trace=FALSE) {
  # use helper functions to validate or create control object
  Mmax<-as.integer(Mmax)
  control<-list(Mmax=Mmax,d=d,trace=trace)
  control<-validate_mars.control(control)
  control<-new_mars.control(control)
}

# Constructor
new_mars.control <- function(control) {
  structure(control,class="mars.control")
}

#Validator & Helper
validate_mars.control <- function(control) {
  # validate control object to be within bounds and fix if it is not
  stopifnot(is.integer(control$Mmax),is.numeric(control$d),is.logical(control$trace))
  if(control$Mmax<2){
    warning("Mmax must be >=2; setting to 2")
    control$Mmax<-2
  }
  if(control$Mmax%%2>0){
    control$Mmax<-2*ceiling(control$Mmax/2)
    warning("Mmax should be an even integer. Setting to ",control$Mmax)
  }
  control
}
