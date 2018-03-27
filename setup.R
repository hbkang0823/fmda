# Load necessary packages. If not installed, install the package.
if (!("tensorA" %in% rownames(installed.packages()))) { 
  install.packages("tensorA")
}
require(tensorA)

if (!("CompQuadForm" %in% rownames(installed.packages()))) { 
  install.packages("CompQuadForm")
}
require(CompQuadForm)

if (!("R.matlab" %in% rownames(installed.packages()))) { 
  install.packages("R.matlab")
}
require(R.matlab)

if (!("matlabr" %in% rownames(installed.packages()))) { 
  devtools::install_github("muschellij2/matlabr")
}
require(matlabr)

if (!("INLA" %in% rownames(installed.packages()))) { 
  install.packages("INLA", repos="https://inla.r-inla-download.org/R/stable")
  # refernece: http://www.r-inla.org/
}
require(INLA)

if (!("splancs" %in% rownames(installed.packages()))) { 
  install.packages("splancs", repos="http://R-Forge.R-project.org")
}
require(splancs)

if (!("akima" %in% rownames(installed.packages()))) { 
  install.packages("akima")
}
require(akima)

if (!("rgl" %in% rownames(installed.packages()))) { 
  install.packages("rgl")
}
require(rgl)

if (!("expm" %in% rownames(installed.packages()))) { 
  install.packages("expm")
}
require(expm)

if (!("colorRamps" %in% rownames(installed.packages()))) { 
  install.packages("colorRamps")
}
require(colorRamps)

# require(R.matlab)
# install.packages("R.matlab")

# Load functions
source("function/simulate_face.R")
source("function/create_mesh.R")
source("function/fpca2.R")
source("function/plot_face.R")
source("function/get_penalty.R")
source("function/MS_penreg.R")
source("function/get_orth_vec.R")
source("function/rglplot.R")
source("function/plot_PC.R")
source("function/plot_beta.R")
source("function/plot_ptsig.R")
source("function/plot_allsig.R")

# Load data for simulation
load("./rundata/for_simulation.Rdata")

