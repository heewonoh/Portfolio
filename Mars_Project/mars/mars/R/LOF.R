LOF <-
function(form,data) {
  ff <- lm(form,data)
  return(sum(residuals(ff)^2))
  }
