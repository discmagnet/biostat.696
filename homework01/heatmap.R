heatmap = function(X,Y,Z){
  # INPUTS: coordinates of each measurement and corresponding value
  library(akima)
  surface <- interp(X,Y,Z)
  Nx <- length(surface$x)
  Ny <- length(surface$y)
  temp <- as.matrix(surface$z)
  foo <- apply(t(temp),2,rev) # counter-clockwise rotation
  output <- data.frame(rep(surface$x, each = Ny),
                       rep(rev(surface$y), Nx),
                       as.vector(foo))
  colnames(output) <- c("lon","lat","obs")
  return(output)
  # OUTPUT: data frame ready to use in ggplot
}