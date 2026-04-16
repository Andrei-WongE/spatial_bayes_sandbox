# INLA
install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/testing"),dep=TRUE,type="binary")

# OPENBUGS

install.packages(c("R2OpenBUGS", "CARBayes", "nimble"))

# In addition, it is useful to load the following R packages: akima, DCluster, shapefiles, mba, tmap, sf, sp and spdep. fillmap is available as R code within the participant files:
  
source("fillmap.R")

# DCluster can be obtained with command:
remotes::install_github("becarioprecario/DCluster")