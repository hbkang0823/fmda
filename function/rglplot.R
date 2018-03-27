## Inputs ##
# toPlot: the object to be plotted
# ref: reference face for the plot to be on (P x 3 matrix)
# qbreaks: break points at which 'toPlot' is separated and different colors are applied
# color: which colors to be applied

## Outputs ##
# plot

rglplot <- function(toPlot, ref, qbreaks, color) {

  require(akima)
  require(rgl)

  # plot

  intpp_z <- interp(x=ref[,1], y=ref[,2], z=ref[,3], nx=500,ny=500)
  intpp_p <- interp(x=ref[,1], y=ref[,2], z=toPlot, nx=500,ny=500)
  Pcol  = cut(intpp_p$z, breaks=qbreaks, include.lowest=TRUE)
  colorK <- color[match(Pcol, levels(Pcol))]

  open3d()
  surface3d(z=intpp_z$z, x=intpp_z$x, y=intpp_z$y, color=colorK, 
			xlab="",ylab="",zlab="",box=FALSE,axis=FALSE)
  viewdata <- data.frame(theta = c(3,25,45,90), phi=c(10,15,5,5))
  v=1; rgl.viewpoint(viewdata$theta[v], viewdata$phi[v])

}