require(nimble)

GPModel<-nimbleCode({
  for(i in 1:m){
    Y1998[i]~dpois(mu[i])
    mu[i]<-Exp98[i]*theta[i]
    theta[i]~dgamma(a,b)
  }
  a~dexp(10)
  b~dexp(10)
})

GPdata<-list(Y1998=c( 18,90,10,120,12,14,76,96,10,256,37,23,40,29,36,
                      55,21,63,15,19,129,47,260,60,10,184,22,45,43,44,10,171,11,
                      34,22,34,51,63,90,201,10,202,87,25,25,91))
GPconsts<-list(m=46,Exp98=c(19.334642001146, 105.221991510865, 8.9954123633133, 126.211287025262, 
                            12.9499400671852, 17.0850039703209, 85.5262771111914, 107.178846922884, 
                            11.0291918950188, 248.419380066852, 38.5954996425929, 27.0027208298727, 
                            32.2453350684913, 24.1871410613557, 29.3284980403873, 52.0933278275436, 
                            23.3496100847714, 69.1791167378613, 15.7011547559647, 17.5779462883105, 
                            98.0421453601469, 42.1724712080047, 277.747093167242, 49.9402374163248, 
                            15.0708479385354, 137.177683720537, 13.3400552455942, 38.1425892644401, 
                            46.222761591486, 49.646669857522, 16.011990994697, 161.116783742905, 
                            7.49225226944375, 27.1667732892036, 23.2255895652772, 27.0506021696774, 
                            50.282471254929, 68.9687528187193, 84.0568694371842, 241.020535657027, 
                            13.3636034454982, 194.239681727817, 84.0882670370562, 23.9367452023769, 
                            29.1377576211652, 121.126445726))
GPinits<-list(a=10,b=10,theta=c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
                                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
                                1, 1, 1, 1, 1))
######################################################################################
GPmod<-nimbleModel(code=GPModel,name="GPM",constants=GPconsts,data=GPdata,inits=GPinits)
################ compilation  #####################################################
CGP<-compileNimble(GPmod,showCompilerOutput = TRUE)


CGPcon<-configureMCMC(CGP,print=TRUE,enableWAIC=TRUE)
CGPcon$addMonitors(c("a","b","theta"))
GPMCMC<-buildMCMC(CGPcon)
CGPMCMC<-compileNimble(GPMCMC)

CGPMCMC$run(10000)

#or

niter=10000
samples<-runMCMC(CGPMCMC,niter=niter,nburnin=9000,nchains=1,summary=TRUE) 

samples$summary
theta<-samples$summary[3:48,1]

#or with wrapper
nimbleMCMC(code=GPModel,constants=GPconsts,data=GPdata,inits=GPinits,nchains=2,niter=10000,summary=TRUE)

#################################################################
#### graphics   ################################################
###############################################################

library(sf)
library(sp)
library(spdep)
SCpoly<-st_read("SC_county_alphasort.shp")
SCmap<-as_Spatial(SCpoly)
plot(SCmap)
GPM<-cbind(theta)
SCcon<-data.frame(GPM)
#graphics....using data.frame and SCmap
library(tmap)
areaID<-as.character(seq(1:46))
attr<-data.frame(SCcon,row.names=areaID)
spg<-SpatialPolygonsDataFrame(SCmap, attr, match.ID = TRUE)
labels<-c("theta")

qtm(spg,fill=c("theta"),fill.palette="Blues")+
  tm_layout(title="",panel.labels=labels,legend.position=c("LEFT","BOTTOM"),legend.height=0.5)

########################## alternate call ############################
mcmc.out<-nimbleMCMC(code=GPModel,constants=GPconsts,
                     data=GPdata,inits=GPinits,nchains=2,niter=10000,summary=TRUE)



##################################################################################
############################CV model###############################################

library(nimble)
CONVmodel<-nimbleCode({
  for(i in 1:m){
    Y1998[i]~dpois(mu[i])
    mu[i]<-Exp98[i]*theta[i]
    log(theta[i])<-a0+v[i]+u[i]
    prex[i]<-step(theta[i]-1)
    v[i]~dnorm(0,tauV)
    LDev[i]<--2*(Y1998[i]*log(mu[i]+0.001)-(mu[i]+0.001)-lfactorial(Y1998[i]))
  }
  u[1:m]~dcar_normal(adj[1:L],wei[1:L],num[1:m],tauU,zero_mean=1)
  for(k in 1:L) {wei[k] <- 1 }
  Dev<-sum(LDev[1:m])
  a0~dnorm(0,tau0)
  tauV~dgamma(2,0.5)
  tau0~dgamma(2,0.5)
  tauU~dgamma(2,0.5)   
})
#########################################################################################
CONVdata<-list(Y1998=c( 18,90,10,120,12,14,76,96,10,256,37,23,40,29,36,
                        55,21,63,15,19,129,47,260,60,10,184,22,45,43,44,10,171,11,
                        34,22,34,51,63,90,201,10,202,87,25,25,91))
CONVconsts<-list(m=46,Exp98=c(19.334642001146, 105.221991510865, 8.9954123633133, 126.211287025262, 
                              12.9499400671852, 17.0850039703209, 85.5262771111914, 107.178846922884, 
                              11.0291918950188, 248.419380066852, 38.5954996425929, 27.0027208298727, 
                              32.2453350684913, 24.1871410613557, 29.3284980403873, 52.0933278275436, 
                              23.3496100847714, 69.1791167378613, 15.7011547559647, 17.5779462883105, 
                              98.0421453601469, 42.1724712080047, 277.747093167242, 49.9402374163248, 
                              15.0708479385354, 137.177683720537, 13.3400552455942, 38.1425892644401, 
                              46.222761591486, 49.646669857522, 16.011990994697, 161.116783742905, 
                              7.49225226944375, 27.1667732892036, 23.2255895652772, 27.0506021696774, 
                              50.282471254929, 68.9687528187193, 84.0568694371842, 241.020535657027, 
                              13.3636034454982, 194.239681727817, 84.0882670370562, 23.9367452023769, 
                              29.1377576211652, 121.126445726),num = c(5, 5, 4, 5, 5, 4, 3, 6, 5, 4, 
                                                                       3, 4, 5, 6, 7, 5, 4, 4, 4, 6, 
                                                                       8, 5, 5, 6, 5, 3, 2, 7, 5, 7, 
                                                                       5, 6, 3, 5, 4, 7, 2, 9, 3, 6, 
                                                                       5, 4, 6, 7, 5, 4
                              ),adj = c(
                                33, 30, 24, 23, 4, 
                                41, 38, 32, 19, 6, 
                                25, 15, 6, 5, 
                                39, 37, 30, 23, 1, 
                                38, 25, 15, 6, 3, 
                                38, 5, 3, 2, 
                                27, 25, 15, 
                                45, 38, 22, 18, 14, 10, 
                                43, 40, 38, 32, 14, 
                                22, 18, 15, 8, 
                                46, 44, 42, 
                                46, 44, 29, 20, 
                                35, 31, 29, 28, 16, 
                                45, 43, 38, 21, 9, 8, 
                                38, 25, 18, 10, 7, 5, 3, 
                                35, 31, 28, 21, 13, 
                                35, 34, 26, 21, 
                                38, 15, 10, 8, 
                                41, 33, 24, 2, 
                                44, 40, 36, 29, 28, 12, 
                                45, 43, 35, 34, 31, 17, 16, 14, 
                                45, 34, 26, 10, 8, 
                                42, 39, 30, 4, 1, 
                                41, 36, 33, 30, 19, 1, 
                                27, 15, 7, 5, 3, 
                                34, 22, 17, 
                                25, 7, 
                                43, 40, 31, 29, 20, 16, 13, 
                                46, 28, 20, 13, 12, 
                                44, 42, 36, 24, 23, 4, 1, 
                                43, 28, 21, 16, 13, 
                                41, 40, 38, 36, 9, 2, 
                                24, 19, 1, 
                                45, 26, 22, 21, 17, 
                                21, 17, 16, 13, 
                                44, 41, 40, 32, 30, 24, 20, 
                                39, 4, 
                                32, 18, 15, 14, 9, 8, 6, 5, 2, 
                                37, 23, 4, 
                                43, 36, 32, 28, 20, 9, 
                                36, 32, 24, 19, 2, 
                                44, 30, 23, 11, 
                                40, 31, 28, 21, 14, 9, 
                                46, 42, 36, 30, 20, 12, 11, 
                                34, 22, 21, 14, 8, 
                                44, 29, 12, 11
                              ),L = 228)
CONVinits<-list(a0=-1.0,tau0=0.1,tauV=0.1,tauU=0.1,
                v=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                u=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))




CONVmod<-nimbleModel(code=CONVmodel,name="CONVM",constants=CONVconsts,data=CONVdata,inits=CONVinits)
CCONV<-compileNimble(CONVmod)

CCONVcon<-configureMCMC(CCONV,print=TRUE,enableWAIC=TRUE)
CCONVcon$addMonitors(c("a0","v","u","theta","Dev"))
CONVMCMC<-buildMCMC(CCONVcon)
CCONVMCMC<-compileNimble(CONVMCMC)
  
  niter=11000
  samples<-runMCMC(CCONVMCMC,niter=niter,nburnin=10000,nchains=1,summary=TRUE,WAIC=TRUE) 
  samples$summary
  
  Dev<-samples$summary[1,1]
  Dev
  DevSD<-samples$summary[1,3]
  pD=(DevSD**2)/2
  pD
  DIC=Dev+pD
  DIC


library(coda)
geweke.plot(as.mcmc(samples$samples[,"Dev"]),main="Deviance")
traceplot(as.mcmc(samples$samples[,"Dev"]),main="deviance")

theta<-samples$summary[6:51,1]
U<-samples$summary[52:97,1]
V<-samples$summary[98:143,1]
library(sf)
library(sp)
library(spdep)

SCpoly<-st_read("SC_county_alphasort.shp")
SCmap<-as_Spatial(SCpoly)
plot(SCmap)
CONVM<-cbind(V,U,theta)
SCcon<-data.frame(CONVM)

library(tmap)
areaID<-as.character(seq(1:46))
attr<-data.frame(SCcon,row.names=areaID)
spg<-SpatialPolygonsDataFrame(SCmap, attr, match.ID = TRUE)
labels<-c("V","U","theta")

qtm(spg,fill=c("V","U","theta"),fill.palette="Blues")+
  tm_layout(title="",panel.labels=labels,legend.position=c("LEFT","BOTTOM"),legend.height=0.5)

samples$WAIC

# To get the WAIC in nimble you must
# enableWAlC=TRUE in the configureMCMC step
# WAIC=TRUE in the runMCMC step

#####################################################################################
#####################################Adding regresion parameters###################################################
library(nimble)
CONVmodel<-nimbleCode({
  for(i in 1:m){
    Y1998[i]~dpois(mu[i])
    mu[i]<-Exp98[i]*theta[i]
    log(theta[i])<-a0+a1*Perc[i]+v[i]+u[i]
    prex[i]<-step(theta[i]-1)
    v[i]~dnorm(0,tauV)
    LDev[i]<--2*(Y1998[i]*log(mu[i]+0.001)-(mu[i]+0.001)-lfactorial(Y1998[i]))
  }
  u[1:m]~dcar_normal(adj[1:L],wei[1:L],num[1:m],tauU,zero_mean=1)
  for(k in 1:L) {wei[k] <- 1 }
  Dev<-sum(LDev[1:m])
  a0~dnorm(0,tau0)
  a1~dnorm(0,tau1)
  tauV~dgamma(2,0.5)
  tau0~dgamma(2,0.5)
  tauU~dgamma(2,0.5)
  tau1~dgamma(2,0.5)   
})
#########################################################################################
CONVdata<-list(Y1998=c( 18,90,10,120,12,14,76,96,10,256,37,23,40,29,36,
                        55,21,63,15,19,129,47,260,60,10,184,22,45,43,44,10,171,11,
                        34,22,34,51,63,90,201,10,202,87,25,25,91))
CONVconsts<-list(m=46,Exp98=c(19.334642001146, 105.221991510865, 8.9954123633133, 126.211287025262, 
                              12.9499400671852, 17.0850039703209, 85.5262771111914, 107.178846922884, 
                              11.0291918950188, 248.419380066852, 38.5954996425929, 27.0027208298727, 
                              32.2453350684913, 24.1871410613557, 29.3284980403873, 52.0933278275436, 
                              23.3496100847714, 69.1791167378613, 15.7011547559647, 17.5779462883105, 
                              98.0421453601469, 42.1724712080047, 277.747093167242, 49.9402374163248, 
                              15.0708479385354, 137.177683720537, 13.3400552455942, 38.1425892644401, 
                              46.222761591486, 49.646669857522, 16.011990994697, 161.116783742905, 
                              7.49225226944375, 27.1667732892036, 23.2255895652772, 27.0506021696774, 
                              50.282471254929, 68.9687528187193, 84.0568694371842, 241.020535657027, 
                              13.3636034454982, 194.239681727817, 84.0882670370562, 23.9367452023769, 
                              29.1377576211652, 121.126445726),num = c(5, 5, 4, 5, 5, 4, 3, 6, 5, 4, 
                                                                       3, 4, 5, 6, 7, 5, 4, 4, 4, 6, 
                                                                       8, 5, 5, 6, 5, 3, 2, 7, 5, 7, 
                                                                       5, 6, 3, 5, 4, 7, 2, 9, 3, 6, 
                                                                       5, 4, 6, 7, 5, 4
                              ),adj = c(
                                33, 30, 24, 23, 4, 
                                41, 38, 32, 19, 6, 
                                25, 15, 6, 5, 
                                39, 37, 30, 23, 1, 
                                38, 25, 15, 6, 3, 
                                38, 5, 3, 2, 
                                27, 25, 15, 
                                45, 38, 22, 18, 14, 10, 
                                43, 40, 38, 32, 14, 
                                22, 18, 15, 8, 
                                46, 44, 42, 
                                46, 44, 29, 20, 
                                35, 31, 29, 28, 16, 
                                45, 43, 38, 21, 9, 8, 
                                38, 25, 18, 10, 7, 5, 3, 
                                35, 31, 28, 21, 13, 
                                35, 34, 26, 21, 
                                38, 15, 10, 8, 
                                41, 33, 24, 2, 
                                44, 40, 36, 29, 28, 12, 
                                45, 43, 35, 34, 31, 17, 16, 14, 
                                45, 34, 26, 10, 8, 
                                42, 39, 30, 4, 1, 
                                41, 36, 33, 30, 19, 1, 
                                27, 15, 7, 5, 3, 
                                34, 22, 17, 
                                25, 7, 
                                43, 40, 31, 29, 20, 16, 13, 
                                46, 28, 20, 13, 12, 
                                44, 42, 36, 24, 23, 4, 1, 
                                43, 28, 21, 16, 13, 
                                41, 40, 38, 36, 9, 2, 
                                24, 19, 1, 
                                45, 26, 22, 21, 17, 
                                21, 17, 16, 13, 
                                44, 41, 40, 32, 30, 24, 20, 
                                39, 4, 
                                32, 18, 15, 14, 9, 8, 6, 5, 2, 
                                37, 23, 4, 
                                43, 36, 32, 28, 20, 9, 
                                36, 32, 24, 19, 2, 
                                44, 30, 23, 11, 
                                40, 31, 28, 21, 14, 9, 
                                46, 42, 36, 30, 20, 12, 11, 
                                34, 22, 21, 14, 8, 
                                44, 29, 12, 11
                              ),L = 228,
                 Perc= c( 19.1,17.7,40.4,18.7,27.6,30.4,13,15,18.5,18.7,20.6,20.8,23.7,27.8,22.6,23.2,26.2,11.7,20.6,
                          23.4,21.7,20.3,15.4,22.1,23.7,19.5,25.9,18.3,19.7,19.5,27.1,12.9,21,29.5,33.1,19,13.6,26.7,17,17,19.3,17.2,20.9,18.7,32.2,13.1 ))
CONVinits<-list(a0=-1.0,a1=0.1,tau0=0.1,tau1=0.1,tauV=0.1,tauU=0.1,
                v=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                u=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))




CONVmod<-nimbleModel(code=CONVmodel,name="CONVM",constants=CONVconsts,data=CONVdata,inits=CONVinits)
CCONV<-compileNimble(CONVmod)

CCONVcon<-configureMCMC(CCONV,print=TRUE,enableWAIC=TRUE)
CCONVcon$addMonitors(c("a0","a1","v","u","theta","Dev"))
CONVMCMC<-buildMCMC(CCONVcon)
CCONVMCMC<-compileNimble(CONVMCMC)

niter=21000
samples<-runMCMC(CCONVMCMC,niter=niter,nburnin=20000,nchains=1,summary=TRUE,WAIC=TRUE) 

Dev<-samples$summary[1,1]
Dev
DevSD<-samples$summary[1,3]
pD=(DevSD**2)/2
pD
DIC=Dev+pD
DIC

library(coda)
geweke.plot(as.mcmc(samples$samples[,"Dev"]),main="Deviance")
traceplot(as.mcmc(samples$samples[,"Dev"]),main="deviance")

output2<-samples$summary
theta.list = output2[ grep( "^theta\\[", rownames(output2) ), ]
theta<-theta.list[,1]
U.list = output2[ grep( "^u\\[", rownames(output2) ), ]
U<-U.list[,1]
V.list= output2[ grep( "^v\\[", rownames(output2) ), ]
V<-V.list[,1]

a1<-mean(samples$samples[,"a1"])
plot(density(samples$samples[,"a1"]),main="a1")
a0<-mean(samples$samples[,"a0"]);X11()
plot(density(samples$samples[,"a0"]),main="a0")
# Is a1 important? coefficient of poverty, but lower limit of CI is negative, suggesting the poverty is not significant related to respiratory cancer at county level


library(sf)
library(sp)
library(spdep)
SCpoly<-st_read("SC_county_alphasort.shp")
SCmap<-as_Spatial(SCpoly)
plot(SCmap)
CONVM<-cbind(V,U,theta)
SCcon<-data.frame(CONVM)
library(tmap)

areaID<-as.character(seq(1:46))
attr<-data.frame(SCcon,row.names=areaID)
spg<-SpatialPolygonsDataFrame(SCmap, attr, match.ID = TRUE)
labels<-c("V","U","theta")

qtm(spg,fill=c("V","U","theta"),fill.palette="Blues")+
  tm_layout(title="",panel.labels=labels,legend.position=c("LEFT","BOTTOM"),legend.height=0.5)


#####################################################################################
#########################Varible selection###############################################
library(nimble)
# Standarized variables, more efficient for sampling

VARSelect<-nimbleCode({
  for(i in 1:159){
    x1c[i]<-(x1[i]-mean(x1[]))/sd(x1[])
    x2c[i]<-(x2[i]-mean(x2[]))/sd(x2[])
    x3c[i]<-(x3[i]-mean(x3[]))/sd(x3[])
    x4c[i]<-(x4[i]-mean(x4[]))/sd(x4[])
    x5c[i]<-(x5[i]-mean(x5[]))/sd(x5[])
    y1[i]~dbin(p1[i],n[i])
    #y2[i]~dbin(p2[i],n[i])
    logit(p1[i])<-b0+b[1]*psi[1]*x1c[i]+b[2]*psi[2]*x2c[i]+b[3]*psi[3]*x3c[i]+b[4]*psi[4]*x4c[i]+b[5]*psi[5]*x5c[i]+v[i]
    v[i]~dnorm(0,tauV)
  }
  for( j in 1: 5){
    psi[j]~dbern(p[j])
    p[j]~dbeta(0.5,0.5)}
  
  tauV<-pow(sdV,-2)
  sdV~dunif(0,2)
  b0~dnorm(0,taub0)
  for(j in 1:5){b[j]~dnorm(0,taub[j])
    taub[j]<-pow(sdb[j],-2)
    sdb[j]~dunif(0,2)
    #taub[j]~dgamma(2,1)
  }
  #taub0~dgamma(2,1)
  taub0<-pow(sdb0,-2)
  sdb0~dunif(0,2)
  
})
VARSelectdata<-list(LBW2007=c( 35,8,18,7,65,17,93,115,55,26,
                               336,25,13,24,34,93,36,34,10,79,18,132,84,13,395,11,28,
                               225,136,8,578,17,952,72,93,99,31,110,17,54,20,8,59,1286,
                               31,22,231,192,26,5,58,24,47,16,21,71,163,182,26,1558,33,
                               4,100,62,56,22,1108,43,199,14,32,47,20,9,221,189,22,75,
                               20,27,23,14,15,39,29,13,80,27,147,6,11,204,23,19,28,9,38,
                               10,41,7,41,23,10,9,37,372,154,23,12,164,36,31,27,17,69,
                               11,22,6,17,21,376,121,8,16,11,156,34,6,60,10,2,35,12,13,
                               33,86,79,43,7,16,114,23,11,20,44,97,90,62,9,26,43,0,10,
                               27,148,10,13,19,31 ), VLBW2007=c( 6,1,2,1,11,0,11,17,
                                                                 6,1,58,9,1,5,8,18,5,8,0,19,4,35,20,2,66,2,4,32,30,0,
                                                                 140,5,185,15,12,22,6,17,3,6,2,3,9,276,5,7,43,37,7,0,11,
                                                                 5,10,2,3,7,35,33,4,311,7,0,18,14,9,4,212,4,31,3,7,7,5,
                                                                 3,47,37,3,15,6,5,5,3,4,9,4,3,23,4,24,2,6,42,3,2,8,2,6,
                                                                 1,5,1,9,4,3,4,4,79,35,5,3,26,10,7,3,3,13,0,2,5,0,2,90,
                                                                 27,2,4,1,30,7,0,9,3,0,8,2,4,6,11,16,12,0,1,21,1,2,5,8,
                                                                 14,20,10,3,7,4,0,1,5,22,1,0,3,6 ))
names(VARSelectdata)=c("y1","y2")
VARSelectconsts<-list(
  Birth2007=c( 293,146,185,33,576,203,1237,1529,352,
               274,2449,174,207,261,461,939,350,302,98,894,187,1720,825
               ,141,4012,215,351,3297,1590,42,5270,119,11274,727,827,
               1444,287,1722,155,390,178,246,442,12028,287,156,1613,
               2019,187,87,764,261,412,226,239,966,1473,2513,305,13802,
               420,36,1170,857,430,180,14092,612,3226,108,414,338,337,
               131,2749,2115,164,869,189,241,275,133,96,381,227,159,773,
               382,1724,94,137,1793,358,168,353,108,318,155,370,88,377,
               270,116,243,700,3376,1624,401,159,2161,330,361,258,196,
               738,122,271,35,191,98,3384,1209,63,193,116,979,375,60,
               501,79,18,358,127,166,150,681,701,452,100,85,1033,158,
               119,216,378,888,1225,556,80,270,534,19,92,302,1905,100,
               124,156,293 ),popden2007=c( 35.29,
                                           24.33,36.87,11.02,178.2,70.85,414,202.06,70.1,36.96,
                                           618.94,56.61,34.74,33.1,68.22,97.02,27.4,127.32,21.76,
                                           77.3,42.73,224.39,383.66,13.59,567.14,37.91,85.52,482.35,
                                           944.31,16.43,1908.69,8.64,2034.12,66.93,81.15,376.19,71.75,
                                           268.71,38.41,80.8,92.53,101.8,47.83,2748.19,40.06,29.51,
                                           290.33,624.66,23.15,10.13,105.81,55.66,32.76,62.22,58.54,
                                           538.67,186.34,703.78,82.77,1876.7,66.53,19.22,177.41,146.38,
                                           54.66,40.34,1794.14,151.96,457.69,20.22,101.8,62.7,104.39,
                                           38.47,576.48,347.75,27.84,173.08,36.88,39.87,31.18,24.57,
                                           31.33,69.15,91.77,42.54,58.5,92.9,116.56,38.36,28.19,201.88,
                                           93.35,33.54,98.68,19.14,82.96,26.35,45.2,21.77,47.15,63.56,
                                           36.93,51.95,118.07,864.91,347.35,168.91,31.65,408.08,169.95,
                                           131.34,52.09,78.78,133.25,39.78,58.78,17.59,44.52,16.99,
                                           609.1,628.13,24.6,23.19,38.15,317.38,140.96,10.13,67.04,
                                           16.8,9.64,47.92,23.15,30.3,30.59,82.5,156.99,75.88,65.37,
                                           34.58,153.5,32.41,28.53,65.01,84.68,144.55,252.58,39.71,
                                           20.69,30.77,45.06,10.71,22.94,103.57,322.01,22.65,21.77,
                                           22.54,37.36 ),black2007=c( 19.19,18.46,17.46,49.33,42.82,
                                                                      3.88,11.54,9.65,33.37,11.88,50.35,26.52,4.77,36.53,15.26,
                                                                      28.99,50.35,26.53,60.09,20.04,24.82,17.49,2.7,28.51,40.37,
                                                                      27.83,10.77,5.87,25.36,61.71,62.22,30.01,23.17,26.51,23.18,
                                                                      15.08,27.98,17.5,22.86,43.61,1.65,1.64,40.03,54.31,29.86,
                                                                      49.24,64.07,34,49.52,7.65,13.95,29.71,32.54,31.96,0.92,
                                                                      18.91,13.34,3.66,8.92,42.88,0.9,9.71,25.62,3.66,28.52,
                                                                      39.11,21.13,5.04,6.97,75.71,5.92,18.79,19.25,10.77,
                                                                      31.16,27.36,26.03,7.83,22.99,14.97,55.46,40.48,39.45,
                                                                      23.43,28.79,25.18,34.96,18.62,42.05,32.37,24.04,34.25,
                                                                      2.17,58.73,8.85,33.44,37.68,33.44,40.03,29.73,46.86,
                                                                      26.39,26.42,25.13,1.72,46.56,35.23,6.74,18.31,15.07,
                                                                      43.3,2.35,10.39,13.73,12.87,33.17,27.17,46.92,1.87,
                                                                      60.6,52.35,38.3,27.6,44.21,34.62,31.98,11.98,58.96,
                                                                      49.59,56.83,58.76,28.65,40.48,42.41,60.66,37.08,27.53,
                                                                      24.51,1.48,32.06,33.11,41.65,40.35,1.56,27.79,4.32,
                                                                      15.39,27.66,56.11,52.73,19.96,44.23,35.2,2.75,3.66,
                                                                      36.67,41.74,40.67,29.38 ),Inc2007=c( 33275,30431,31424,
                                                                                                           41795,39180,43813,51283,48366,31781,34460,36954,38711,
                                                                                                           38316,32302,60879,34861,33142,46610,28348,46583,31239,
                                                                                                           43196,47022,39270,45124,45308,33356,63518,36158,26051,
                                                                                                           43674,30119,64665,33666,33024,68224,30607,58627,42626,
                                                                                                           30385,40654,56964,32650,51753,32393,30266,33073,55659,
                                                                                                           28703,32430,54132,33344,26223,34526,35710,79166,42685,
                                                                                                           85318,36616,58052,39902,36994,46260,42769,33060,38835,
                                                                                                           64005,41696,51764,27426,38722,61783,36667,42076,63395,
                                                                                                           52911,33902,49820,44695,32530,30967,28152,28310,50652,
                                                                                                           40923,34041,36092,
                                                                                                           55446,40993,36260,37334,38666,46582,27742,42627,34453,
                                                                                                           36645,36026,35093,31895,30732,50423,34961,47245,35817,
                                                                                                           41095,49995,70872,43514,59828,39506,50617,37173,51043,
                                                                                                           37814,39076,44128,28877,39188,26487,36944,55642,37613,
                                                                                                           32638,31094,39310,33965,27687,34882,31515,27021,33003,
                                                                                                           29657,28950,30394,39665,37335,32480,40383,29442,40867,
                                                                                                           30582,36096,42302,35046,37568,53022,34286,31964,34760,
                                                                                                           38773,33287,31707,41325,43408,31697,32461,35367,37643 ),
  Pov2007=c( 29.7,31.7,29,39.7,24.2,19.1,17,16.8,32.6,29.1,
             35.6,26.4,25.7,35.1,12.9,24.5,34.5,17.9,36.3,15.1,32.4,
             20.2,15.2,28.3,23.7,18.5,24.8,7.5,28.4,49.4,20.6,32.3,
             12.3,30,31.2,8.4,31.5,13.1,26.3,40.3,16.8,14.2,32.9,21.3,
             29.7,35.3,35,14.8,41.9,36.4,12.7,24.3,40.2,32.4,24,5.9,23.3,
             6.1,23,20.8,21.9,18.5,21.9,19.9,31.1,34.5,11.5,20.9,16.7,
             34.8,23.5,11.6,23.9,23.5,9.8,16.4,29.1,16.2,24.7,28.9,32,
             37.7,39.3,15.5,22.4,29.6,30.6,12.8,24.6,25.4,31.5,29.6,
             18.9,37.2,18.7,32.3,27.8,29.5,33.3,31.8,34.1,15.9,28,17.9,
             21.1,27.6,16,8,18.9,9.1,29.7,18.1,27.6,15.7,24.4,25.3,24.9,
             33.4,22,35.5,34.3,16.5,23.5,29.5,35.7,27.2,24.7,40.9,35.5,35,
             38.1,32.6,33,33.7,43.6,25.3,40.7,32.5,20.6,35.7,23.6,37.9,
             28.5,20.7,25.8,22.7,17.6,29.7,33.5,29,25.7,31.6,34.2,19.9,
             20.1,34.4,26.9,24,27.7 ),UER2007=c( 5.3,5.7,4.8,5,5.6,3,
                                                 4.2,4.9,6.9,4.4,5.2,5.3,4.5,4.3,3.5,4.4,7.8,5.2,6.6,4,4.5,
                                                 4.8,3.6,4.8,4,9.5,5.9,3.6,3.9,5.6,5.7,5.3,4,6,4.5,3.8,5.7,
                                                 4,4.7,5.8,4.3,3.7,6.1,4.8,5,5.3,5.6,4.7,5.6,3.1,3.3,5.7,5.2,
                                                 4.3,3.8,3.9,4.7,3.3,4.9,4.8,3.6,5.1,3.6,4.7,4,5.4,4,4,3.6,
                                                 8.4,4.7,3.6,6.5,5,4.4,3.9,6.1,4.2,4.8,6.8,7.2,8.9,5.9,4.3,
                                                 5.1,3.7,5.1,3.7,5.2,5.8,3.7,3.9,4.2,7.3,3.8,4.7,6.4,4.1,5.9,
                                                 4.4,5.2,4.1,4.7,4.4,4.6,5.2,5.3,3,3.9,4.1,5.5,3.8,4,4.3,4.5,
                                                 4,4.4,5.4,5.9,6.5,6.1,5,6,4.9,5.5,5.7,4.9,6.3,7.1,6.4,7.9,
                                                 4.8,5.9,7.8,6,4.1,5.3,4.8,3.7,5.8,5.7,6.7,5.7,3.6,6.4,4.5,
                                                 4.6,4.9,7.6,5.7,5.4,4.5,5.2,3.6,4.6,5.9,6.7,5.4,5.4 ))
names(VARSelectconsts)<-c("n","x1","x2","x3","x4","x5")

##############################################################
VARSelectinits<-list(b0=-100,b=c(-0.1,-0.1,-0.1,-0.1,-0.1),sdb0=0.1,sdV=0.5,sdb=c(0.1,0.1,0.1,0.1,0.1),
                     v=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0),psi=c(1,1,1,1,1),p=c(0.5,0.5,0.5,0.5,0.5))

################################################################
VARCON<-nimbleModel(code=VARSelect,data=VARSelectdata,constants=VARSelectconsts,inits=VARSelectinits)
CLN<-compileNimble(VARCON)

CLNcon<-configureMCMC(CLN,print=TRUE,enableWAIC=TRUE)
CLNcon$addMonitors(c("b","psi","p"))

LNMCMC<-buildMCMC(CLNcon)
CLNMCMC<-compileNimble(LNMCMC)

niter=100000
samples<-runMCMC(CLNMCMC,niter=niter,nburnin=90000,nchains=1,summary=TRUE,WAIC=TRUE) 
#library(coda)
# geweke.plot(as.mcmc(samples$samples[,"psi"]),main="psi")

samples$summary
# Only those estimates >0.5 are to be considered

##############################################################################################
###########################################CARBayes and INLA##################################
#SC respiratory cancers
library(sp)
library(spdep)
library(CARBayes)
library(sf)

#CARBayes UH model.....have to assume leroux with rho=0

SClist<-list(Y1998=c( 18,90,10,120,12,14,76,96,10,256,37,23,40,29,36,
                      55,21,63,15,19,129,47,260,60,10,184,22,45,43,44,10,171,11,
                      34,22,34,51,63,90,201,10,202,87,25,25,91 ),
             
             Exp98=c(19.334642001146, 105.221991510865, 8.9954123633133, 126.211287025262, 
                     12.9499400671852, 17.0850039703209, 85.5262771111914, 107.178846922884, 
                     11.0291918950188, 248.419380066852, 38.5954996425929, 27.0027208298727, 
                     32.2453350684913, 24.1871410613557, 29.3284980403873, 52.0933278275436, 
                     23.3496100847714, 69.1791167378613, 15.7011547559647, 17.5779462883105, 
                     98.0421453601469, 42.1724712080047, 277.747093167242, 49.9402374163248, 
                     15.0708479385354, 137.177683720537, 13.3400552455942, 38.1425892644401, 
                     46.222761591486, 49.646669857522, 16.011990994697, 161.116783742905, 
                     7.49225226944375, 27.1667732892036, 23.2255895652772, 27.0506021696774, 
                     50.282471254929, 68.9687528187193, 84.0568694371842, 241.020535657027, 
                     13.3636034454982, 194.239681727817, 84.0882670370562, 23.9367452023769, 
                     29.1377576211652, 121.126445726))

SCresp<-data.frame(SClist)



SCpoly<-st_read("SC_county_alphasort.shp")
SCmap<-st_geometry(SCpoly)
W.nb <- poly2nb(SCmap)
W.mat <- nb2mat(W.nb, style="B")

## UH model
form<-Y1998~1+offset(log(Exp98))
modelUH<-S.CARleroux(form,family="poisson",data=SCresp, W=W.mat,rho=0, burnin=10000, 
                     n.sample=11000)
modelUH$summary
print(modelUH)
modelUH$samples
modelUH$modelfit

sample.fitted<-modelUH$samples$fitted

m=46;L=1000
ypred<-matrix(0,ncol=m,nrow=L);ppl<-matrix(0,ncol=m,nrow=L)

for (j in 1:m){
  for (i in 1:L){
    ypred[i,j]<-rpois(1,sample.fitted[i,j])}}
ycopy<-matrix(0,ncol=m,nrow=L)
for(i in 1:m){
  for(j in 1:L){
    ycopy[j,i]<-SCresp$Y1998[i]}}
for(i in 1:m){
  for( j in 1:L){
    a<-(ycopy[j,i]-ypred[j,i])**2
    ppl[j,i]<-a}}

mspeI<-rep(0,m)
mspeI<-colMeans(ppl)
mspe<-mean(mspeI)
mspe
sd(mspeI)


###### sample theta  ###############################
L =1000
sample.theta<-matrix(0,nrow=L,ncol=m)
for (i in 1:L){
  for(j in 1:m){
    sample.theta[i,j]<-sample.fitted[i,j]/SCresp$Exp98[j]
  }}

theta.est<-colMeans(sample.theta[,])

source("fillmap.R")

fillmap(SCmap,"",theta.est,n.col=4)

UHeffect<-log(theta.est)-0.0005;x11()
fillmap(SCmap,"",UHeffect,n.col=4)


#### ICAR model #############################################
form<-Y1998~1+offset(log(Exp98))
modelCH<-S.CARleroux(form,family="poisson",data=SCresp, W=W.mat,rho=1, burnin=10000, 
                     n.sample=11000)
# Rho equal to 1 is the ICAR model, rho=0 is the UH model
modelCH$summary
print(modelCH)
modelCH$samples
modelCH$modelfit

sample.theta<-matrix(0,nrow=L,ncol=m)
for (i in 1:L){
  for(j in 1:m){
    sample.theta[i,j]<-sample.fitted[i,j]/SCresp$Exp98[j]
  }}
theta.est<-colMeans(sample.theta[,])
CHeffect<-log(theta.est)-0.001
fillmap(SCmap,"",CHeffect,n.col=5)


## Conv model

form<-Y1998~1+offset(log(Exp98))
modelCONV<-S.CARbym(form,family="poisson",data=SCresp, W=W.mat, burnin=40000, 
                    n.sample=45000)
modelCONV$summary
modelCONV$modelfit
sample.theta<-matrix(0,nrow=L,ncol=m)
for (i in 1:L){
  for(j in 1:m){
    sample.theta[i,j]<-sample.fitted[i,j]/SCresp$Exp98[j]
  }}
theta.est<-colMeans(sample.theta[,])
fillmap(SCmap,"",theta.est,n.col=5)



sample.fitted<-modelCONV$samples$fitted
ypred<-matrix(0,ncol=m,nrow=L);ppl<-matrix(0,ncol=m,nrow=L)
for (j in 1:m){
  for (i in 1:L){
    ypred[i,j]<-rpois(1,sample.fitted[i,j])}}
ycopy<-matrix(0,ncol=m,nrow=L)
for(i in 1:m){
  for(j in 1:L){
    ycopy[j,i]<-SCresp$Y1998[i]}}
for(i in 1:m){
  for( j in 1:L){
    a<-(ycopy[j,i]-ypred[j,i])**2
    ppl[j,i]<-a}}
mspeI<-rep(0,m)
mspeI<-colMeans(ppl)
mspe<-mean(mspeI)
mspe
sd(mspeI)

## Leroux model, balances UH and UC     
form<-Y1998~1+offset(log(Exp98))
modelLER<-S.CARleroux(form,family="poisson",data=SCresp, W=W.mat, burnin=50000, 
                      n.sample=55000)

modelLER$summary
modelLER$modelfit

sample.fitted<-modelLER$samples$fitted
sample.theta<-matrix(0,nrow=L,ncol=m)
for (i in 1:L){
  for(j in 1:m){
    sample.theta[i,j]<-sample.fitted[i,j]/SCresp$Exp98[j]
  }}
theta.est<-colMeans(sample.theta[,])

fillmap(SCmap,"",theta.est,n.col=4,leg.loc="bottomleft",leg.cex=1.0)

ypred<-matrix(0,ncol=m,nrow=L);ppl<-matrix(0,ncol=m,nrow=L)
for (j in 1:m){
  for (i in 1:L){
    ypred[i,j]<-rpois(1,sample.fitted[i,j])}}
ycopy<-matrix(0,ncol=m,nrow=L)
for(i in 1:m){
  for(j in 1:L){
    ycopy[j,i]<-SCresp$Y1998[i]}}
for(i in 1:m){
  for( j in 1:L){
    a<-(ycopy[j,i]-ypred[j,i])**2
    ppl[j,i]<-a}}
mspeI<-rep(0,m)
mspeI<-colMeans(ppl)
mspe<-mean(mspeI)
mspe
sd(mspeI)


##################################################################################
############################INLA##################################################

library(INLA)
#SC respiratory cancers 
#INLA UH model.....
SClist<-list(Y1998=c( 18,90,10,120,12,14,76,96,10,256,37,23,40,29,36,
                      55,21,63,15,19,129,47,260,60,10,184,22,45,43,44,10,171,11,
                      34,22,34,51,63,90,201,10,202,87,25,25,91 ),
             Exp98=c(19.334642001146, 105.221991510865, 8.9954123633133, 126.211287025262, 
                     12.9499400671852, 17.0850039703209, 85.5262771111914, 107.178846922884, 
                     11.0291918950188, 248.419380066852, 38.5954996425929, 27.0027208298727, 
                     32.2453350684913, 24.1871410613557, 29.3284980403873, 52.0933278275436, 
                     23.3496100847714, 69.1791167378613, 15.7011547559647, 17.5779462883105, 
                     98.0421453601469, 42.1724712080047, 277.747093167242, 49.9402374163248, 
                     15.0708479385354, 137.177683720537, 13.3400552455942, 38.1425892644401, 
                     46.222761591486, 49.646669857522, 16.011990994697, 161.116783742905, 
                     7.49225226944375, 27.1667732892036, 23.2255895652772, 27.0506021696774, 
                     50.282471254929, 68.9687528187193, 84.0568694371842, 241.020535657027, 
                     13.3636034454982, 194.239681727817, 84.0882670370562, 23.9367452023769, 
                     29.1377576211652, 121.126445726))

SCresp<-data.frame(SClist)

m=46
region<-seq(1:m);region2<-region;idx<-region2

#UH model
formulaUH = Y1998~ f(region, model = "iid",param=c(2,0.5))
resUH = inla(formulaUH,family="poisson",data=SCresp,control.compute=list(dic=TRUE,cpo=TRUE,waic=TRUE),E=Exp98)
summary(resUH)

cpo<-resUH$cpo$cpo
locdic<-resUH$dic$local.dic
LPML<-sum(log(cpo))
LPML

resUH$dic$dic
resUH$dic$p.eff

UHeffect<-resUH$summary.random$region[,2]

####### graphics #############################################
source("fillmap.R")

SCpoly<-st_read("SC_county_alphasort.shp")   # SC counties 46#
SCpoly<-as_Spatial(SCpoly)
plot(SCpoly)

fillmap(SCpoly,"UH",UHeffect,n.col=5);x11()

par(mfrow=c(1,2))
fillmap(SCpoly,"cpo",cpo,n.col=4)
fillmap(SCpoly,"local dic",locdic,n.col=4)


#############################################################################

###### adj and num from SC counties here #####################
#  if you have the adj and num vectors already then you can create the graph file 
#inla.geobugs2inla(adj, num, graph.file=".........txt")
###############################################################

### ICAR model   ###############################################
region<-seq(1:m);region2<-region;idx<-region
formICAR<-Y1998~1+f(region,model="besag",param=c(2,0.5),graph="SCgraph.txt") 
resICAR = inla(formICAR,family="poisson",data=SCresp,control.compute=list(dic=TRUE,cpo=TRUE,waic=TRUE),E=SClist$Exp98)
cpo<-resICAR$cpo$cpo
LPML<-sum(log(cpo))
LPML
dic<-resICAR$dic$dic
pD<-resICAR$dic$p.eff
CHeffect<-resICAR$summary.random$region[,2]
par(mfrow=c(1,2))
fillmap(SCpoly,"ICAR",CHeffect,n.col=5)
fillmap(SCpoly,"CPO",cpo,n.col=5)

#### # Conv model ############################################################################

formCONV = Y1998~1+f(region, model = "iid",param=c(2,0.5))+f(region2, model = "besag",
                                                             param=c(2,0.5),graph="SCgraph.txt")
resCONV = inla(formCONV,family="poisson",data=SCresp,
               control.compute=list(dic=TRUE,cpo=TRUE,waic=TRUE),E=Exp98)
summary(resCONV)

cpo<-resICAR$cpo$cpo
LPML<-sum(log(cpo))
LPML
resCONV$dic$dic
resCONV$dic$p.eff

UHeffect<-resCONV$summary.random$region[,2]
CHeffect<-resCONV$summary.random$region2[,2]
par(mfrow=c(2,2))
fillmap(SCpoly,"UH",UHeffect,n.col=5)
fillmap(SCpoly,"CH",CHeffect,n.col=5)
fillmap(SCpoly,"CPO",cpo,n.col=5)

####################################################################################################

#### Leroux Model################################################################
#library(INLABMA)
#library(spdep)
#W.nb <- poly2nb(SCpoly)
#W.mat <- nb2mat(W.nb, style="B")

#rlambda <- seq(0.03, 0.8, length.out = 20)
#errorhyper <- list(prec = list(prior = "loggamma",
#param = c(1, 0.01), initial = log(1), fixed = FALSE))
#form2 <- Y1998 ~ 1+offset(log(Exp98))
#lerouxmodels <- mclapply(rlambda, function(lambda) {
#leroux.inla(form2, d =SCresp, W = W.mat,
#lambda = lambda, improve = TRUE,
#family = "poisson",
#control.predictor = list(compute = TRUE),
#control.compute = list(dic = TRUE, cpo = TRUE,waic=TRUE))
#})
#resLER<- INLABMA(lerouxmodels, rlambda, 0, impacts = FALSE)

#cpo<-resLER$cpo
#LPML<-sum(log(cpo$cpo))
#LPML
#resLER$dic$dic
#resLER$dic$p.eff
#resLER$rho

####################################################################################
#################### BYM2  #########################################################
region2<-region

formBYM2 = Y1998~1+f(region2, model = "bym2",scale.model=TRUE,param=c(2,0.5),graph="SCgraph.txt")
resBYM2 = inla(formBYM2,family="poisson",data=SCresp,control.compute=list(dic=TRUE,cpo=TRUE,waic=TRUE),E=Exp98)
summary(resBYM2)

cpo<-resBYM2$cpo
LMPL<-sum(log(cpo$cpo))
LMPL

resBYM2$dic$dic
resBYM2$dic$p.eff
rand<-resBYM2$summary.random$region2[,2]
UH<-rand[1:46]
CH<-rand[47:92]
################################# using fillmap #######################

fillmap(SCpoly,"",UH,n.col=5)
fillmap(SCpoly,"",CH,n.col=5)

############## using tmap ############################################ #####################################################################
library(sp)
library(spdep)
#setwd(".........................................................")
SCpoly<-st_read("SC_county_alphasort.shp")
SCmap<-as_Spatial(SCpoly)
LERM<-cbind(UH,CH)
SCcon<-data.frame(LERM)
#graphics....using data.frame and SCmap
library(tmap)
areaID<-as.character(seq(1:46))
attr<-data.frame(SCcon,row.names=areaID)
spg<-SpatialPolygonsDataFrame(SCmap, attr, match.ID = TRUE)
labels<-c("UH","CH")
qtm(spg,fill=c("UH","CH"),fill.palette="Blues")+
  tm_layout(title="",panel.labels=labels,legend.position=c("LEFT","BOTTOM"),legend.height=0.5)





