#--- Getting empirical direct semivariogram
extractPointVgm <- function(g) {
  if(!is(g, "ataKrigVgm")) return(NULL)
  
  if(hasName(g, "pointVariogram")) {
    return(g$pointVariogram)
  } else {
    for(id in names(g)) {
      if(hasName(g[[id]], "pointVariogram"))
        g[[id]] <- g[[id]]$pointVariogram
    }
    return(g)
  }
}
