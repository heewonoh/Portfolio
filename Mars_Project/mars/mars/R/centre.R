centre <-
function(x,method) {
  switch(method,mean=mean(x),median=median(x),
         {message("\nmethod ",method,
                  " not implemented, using mean\n"); 
           mean(x)})
}
