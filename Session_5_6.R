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



