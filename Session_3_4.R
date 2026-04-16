require(akima)
require(tmap)
# remotes::install_github('r-tmap/tmap')
require(sp)
# remotes::install_github("carrollrm/fillmap")
require(fillmap)
require(here)

#R commands: image, persp, contour, interp
asd<-read.table("surfXYZ.txt")
x<-asd$V1
y<-asd$V2
z<-asd$V3


asdP<-interp(x,y,z)
image(asdP)
contour(asdP,add=T)
persp(asdP)



##  MBA  plotting 

library(MBA)
data(LIDAR)
mba.int <- mba.surf(LIDAR, 300, 300, extend=TRUE)$xyz.est
image(mba.int, xaxs="r", yaxs="r")

contour(mba.int,add=T)

##  Polygons

source("fillmap.R")
# fillmap uses polygons objects to define maps areas

############################################################
##### fillmap code ########################################
## from GEoBUGS 

library(DCluster)
library(sf)
polySC<-readSplus("SC_geobugsSPlus.txt")

# polySC is a polygon object
# To use with fillmap you need to convert the format:

SCpoly<-st_as_sfc(polySC)
fillmap(SCpoly,"title",rnorm(46,2,2),n.col=5)

### from shapefiles 
library(sf)
library(sp)

polySC<-st_read("SC_county_alphasort.shp")
polySC<-st_geometry(polySC)
plot(polySC)

rand<-rnorm(46,1,1)
fillmap(polySC,"noise",rand,n.col=5)

# Using tmap

SCpoly<-st_read("SC_county_alphasort.shp")
SCpoly2<-as_Spatial(SCpoly)
plot(SCpoly2)

###SCcongen90 is dataframe including the vectors obs,expe,pov,inc######
#######################################################################
SCcongen90<-list(obs=c(0,7,1,5,1,1,5,16,0,17,4,0,0,1,1,7,1,3,0,0,8,2,13,7,0,8,0,3,2,4,1,11,0,1,2,3,3,8,6,14,3,11,6,0,1,5),
                 expe=c(1.0807,6.3775,0.622,6.6854,0.9142,1.0744,5.6518,8.1682,0.5749,18.0989,2.174,1.6619,1.9321,1.6148,
                        1.6713,3.0819,1.7562,4.9952,0.9362,1.2001,6.1293,2.5604,15.8589,2.9437,1.0399,7.276,0.9739,2.064,2.7206,2.8275,0.9425,
                        8.828,0.3644,1.775,1.5111,1.5111,2.5321,4.5836,3.9647,15.0264,0.732,10.8292,5.9848,1.4357,1.9949,6.9807),
                 pov=c(13.6,13.8,32.3,11,24.2,19.9,12.3,13.3,17,15.4,14.1,15.7,18,24.3,21.8,20.2,24.9,12,17.4,18.1,18.7,17.5,
                       10.6,13.8,22.8,13.7,21,12.6,14,14,26.9,9.4,17.8,23.1,23.1,14.4,10.8,22.1,10.1,13.6,15.9,11.2,18.3,14,26.4,10.6),
                 inc=c(36.786,38.534,20.727,37.205,24.3,27.607,45.822,40.161,32.247,38.458,33.232,31.715,29.505,25.896,28.919,
                       30.776,25.552,42.886,34.297,29.96,34.009,35.008,41.658,34.109,27.65,34.654,27.117,39.04,33.698,32.32,25.144,
                       45.14,29.805,25.008,25.993,32.231,36.912,28.624,37.054,39.587,31.324,37.092,31.948,30.801,23.748,44.619))
########################################################################
SCcon<-data.frame(SCcongen90)

#graphics....using data.frame and SCpoly2

areaID<-as.character(seq(1:46))
#area<-paste(c("area"),areaID,sep="")
attr<-data.frame(SCcon,row.names=areaID)
spg<-sp::SpatialPolygonsDataFrame(SCpoly2, attr, match.ID = TRUE)

qtm(spg,fill=c("obs","expe","pov","inc"),fill.palette="Blues",ncol=2)+tm_layout(legend.position=c("LEFT","BOTTOM"),legend.height=0.4)

####################spplot###############################################
#########################################################################

spplot(spg,zcol=c("obs","expe","pov","inc"),
       col.regions=grey(seq(0.9,0.1,length=40)),layout=c(4,1))


#########################################################################
#########################################################################

library(R2OpenBUGS)
Sys.setenv("OpenBUGS_PATH" = "C:/Users/Andre/Downloads/OpenBUGS323/OpenBUGS323")

parameters<-c("theta","a","b")
inits<-c(rep("gamma_poisson_BUGS_inits.txt",2))
data<-c("gamma_poisson_BUGS_data.txt")
model.file<-c("gamma_poisson_BUGS_model.txt")

GP_result <- bugs(data, inits, parameters, model.file,n.chains=2, n.iter=10000,n.burnin=8000)
print(GP_result)
plot(GP_result)

GP_result$mean
GP_result$summary
GP_result$DIC

GP_result$sims.matrix

######################## optional plotting #############################
library(sf)
library(sp)
source("fillmap.R")
polySC<-st_read("SC_county_alphasort.shp")
polySC<-as_Spatial(polySC)
fillmap(polySC,"",GP_result$mean$theta,n.col=5)
########################################################################

### lognormal model ###########################
parameters<-c("theta","v")
inits<-c(rep("log_normal_BUGS_inits.txt",2))
data<-c("log_normal_BUGS_data.txt")
model.file<-c("log_normal_BUGS_model.txt")

LN_result<-bugs(data, inits, parameters, model.file,n.chains=2, n.iter=10000,n.burnin=8000)
print(LN_result)
plot(LN_result)
LN_result$summary
LN_result$DIC

######## convergence checking #########################################
library(coda)
#length<-length(as.mcmc(LN_result$sims.matrix[,93]))  
#lengthmin <- ifelse(length>8000,length-8000, 0)             
geweke.plot(as.mcmc(LN_result$sims.matrix[1:2000,93]),auto.layout=FALSE)  
traceplot(as.mcmc(LN_result$sims.matrix[1:2000,93]))
######################################################################

### convolution model 
parameters<-c("theta","v","u")
inits<-c(rep("CONV_BUGS_inits.txt",2))
data<-c("CONV_BUGS_data.txt")
model.file<-c("CONV_BUGS_model.txt")


CONV_result<-bugs(data, inits, parameters, model.file,n.chains=2, n.iter=10000,n.burnin=8000)

print(CONV_result)
plot(CONV_result)
CONV_result$mean
CONV_result$summary
CONV_result$DIC

#################################################################################################
#################################################################################################

library(R2OpenBUGS)
parameters<-c("theta","v","a0","p","u","deviance")
inits<-c(rep("SMix_BUGS_inits.txt",2))
data<-c("SMix_BUGS_data.txt")
model.file<-c("SMix_BUGS_model.txt")

SMixRes<-bugs(data, inits, parameters, model.file,n.chains=2, n.iter=10000,n.burnin=8000)
print(SMixRes)
plot(SMixRes)
SMixRes$summary
SMixRes$DIC
theta<-SMixRes$mean$theta
u<-SMixRes$mean$u
v<-SMixRes$mean$v
P<-SMixRes$mean$p

library(coda)
Deviance<-SMixRes$sims.matrix[1:2000,"deviance"]
geweke.plot(as.mcmc(Deviance))

########## Leroux model  ##########################################
parameters<-c("theta","a0","s","rho","deviance")
model.file<-c("Leroux_BUGS_model_Q1RD.txt")
data<-c("Leroux_BUGS_data_Q1RDmodel.txt")
inits<-c(rep("Leroux_BUGS_inits_Q1RDmodel.txt",1))

LerRes<-bugs(data,inits,parameters,model.file,n.chains=1,n.iter=10000,n.burnin=8000)
print(LerRes)
plot(LerRes)
LerRes$summary
LerRes$DIC

#####################################################################
############ BYM2 model  ############################################
parameters<-c("theta","UCorr","Corr","rho","sigma","deviance")

model.file<-c("BYM2_BUGS_model.txt")
data<-c("BYM2_BUGS_data.txt")
inits<-c(rep("BYM2_BUGS_inits.txt",1))

BYM2Res<-bugs(data,inits,parameters,model.file,n.chains=1,n.iter=10000,n.burnin=8000)
print(BYM2Res)
plot(BYM2Res)
BYM2Res$summary
BYM2Res$DIC


Corr<-BYM2Res$mean$Corr
UCorr<-BYM2Res$mean$UCorr
rho<-BYM2Res$mean$rho
sigma<-BYM2Res$mean$sigma
theta<-BYM2Res$mean$theta


BYM2M<-cbind(Corr,UCorr,theta)

SCpoly<-st_read("SC_county_alphasort.shp")
SCpoly<-st_geometry(SCpoly)


fillmaps(SCpoly,"",BYM2M,n.col=4,leg.cex=0.9,
         lay.m=matrix(c(1,2,3),ncol=3,nrow=1))

#######################################################################

##### SMix maps ########################################################

SMixMat<-cbind(u,v,P,theta)
fillmaps(SCpoly,"",SMixMat,n.col=4,leg.loc="bottomleft",leg.cex=0.8,lay.m=matrix(c(1,2,3,4),ncol=2,nrow=2))











