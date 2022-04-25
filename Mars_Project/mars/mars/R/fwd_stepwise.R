fwd_stepwise <-
function(formula,data){    ##3a changes
  
  fwd_list<- lm(formula,data)
 
  return(fwd_list)
  
}
