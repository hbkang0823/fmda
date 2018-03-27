## Inputs ##
# face: P x 3 single face
# ifsave: T/F. if you want to save. default is FALSE
# plotname: if ifsave=TRUE, have plotname to save

## Outputs ##
# plot

plot_face <- function(face, ifsave=FALSE, plotname="tmp.png") {
  require(akima)
  require(rgl)

  viewdata <- data.frame(theta = c(3,25,45,90), phi=c(10,15,5,5))
  intpp <- interp(x=face[,1], y=face[,2], z=face[,3], nx=500,ny=500)
  open3d()
  surface3d(z=intpp$z,x=intpp$x,y=intpp$y,color=intpp$z,xlab="",ylab="",zlab="",box=FALSE,axis=FALSE)
  v=1; rgl.viewpoint(viewdata$theta[v], viewdata$phi[v])
  if (ifsave == TRUE) {
    rgl.snapshot(plotname)
  }
}
