## Inputs ##
# M0: reference in d-dimension after 'dimred'
# whichinner: if you want to define points that should go as 'inner' mesh that is finer
# obval: a control parameter to set outer boundary
# inbval: a control parameter to set inner boundary

## Outputs ##
# p: node
# t: edge

# requires 'INLA' package
# if not installed, use
# install.packages("INLA", repos="https://inla.r-inla-download.org/R/stable")
# Refernece: http://www.r-inla.org/

# also requires 'splancs' package
# if not installed, use
# install.packages("splancs", repos="http://R-Forge.R-project.org")

create_mesh <- function(M0, whichinner = c(), obval = -0.05, inbval = -0.05, ifplot=FALSE) {
  require(INLA)

  x.range <- max(M0[,1]) - min(M0[,1])
  maxedg <- c((x.range/32),(x.range/22))

  if (length(whichinner != 0)) {
    innerbound <- inla.nonconvex.hull(as.matrix(M0[whichinner,]), inbval, inbval, resolution=c(200,200))
    boundn <- inla.nonconvex.hull(as.matrix(M0),obval, obval, resolution=c(200,200))
    premesh <- inla.mesh.2d(boundary=list(innerbound, boundn), max.edge=maxedg, min.angle=5)
  } else {
    boundn <- inla.nonconvex.hull(as.matrix(M0), obval, obval, resolution=c(200,200))
    premesh <- inla.mesh.2d(boundary=boundn, max.edge=mean(maxedg), min.angle=5)
  }
  mesh_p <- premesh$loc[,1:2]
  mesh_t <- premesh$graph$tv
  
  if (ifplot) { plot(premesh,main=""); title("Mesh Created from M0") }
  return(list(p = mesh_p, t = mesh_t))
}