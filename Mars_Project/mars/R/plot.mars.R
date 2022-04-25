
#' Plot a mars object
#'
#' Plots the fitted basis functions that depend on
#' explanatory variable(main effects) or two explanatory variables
#' (two-way interactions). Use `predict.lm()` to see the residual plots
#'
#'
#' @param x a mars x
#'
#' @param ... additional arguments to pass to plot()
#'
#' @family methods
#' @export
#'
#' @examples mm<-mars(y~x1+x2,data=marstestdata,mars.control(Mmax=4))
#' plot(mm,col="red")
#' class(mm)<-"lm";plot(mm) #Standard diagnostic plots
#'
plot.mars<-function(x,...){
  data<-eval(x$call$data)
  tt<-terms(x$formula,data=data)
  tt<-delete.response(tt)
  mf<-model.frame(tt,data)
  mt<-attr(mf,"terms")
  X<-model.matrix(mt,mf)[,-1] #remove intercept
  Bf<-x$Bfuncs
  singleB<-which(sapply(Bf,function(x) NROW(x)==1))
  doubleB<-which(sapply(Bf,function(x) NROW(x)==2))
  nn<-ceiling(sqrt(length(singleB)+length(doubleB)))
  opar<-graphics::par(mfrow=c(nn,nn),mar=c(2,2,2,2))
  on.exit(graphics::par(opar))#use exit handler to reset pars
  for(i in singleB){
    vv<-Bf[[i]][1,"v"];varname<-x$x_names[[vv]]
    xx<-seq(from=min(X[,vv]),to=max(X[,vv]),length=100)
    bb<-h(xx,Bf[[i]][1,"s"],Bf[[i]][1,"t"])
    plot(xx,bb,type="l",xlab=varname,main=varname,...)#plot line graph
  }
}




