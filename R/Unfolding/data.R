rm(list=ls(all=TRUE))
#
setwd("~/Dropbox/encuesta_cep73_abril2015")
library(smacof)
library(coda)
library(sm)
library(Rmpfr)
library(foreign)
library(rgl)
library(plyr)
dyn.load("bayesianunfoldingmac.so")

set_globals <- function(nslice,nburn,nrowX,ncolX,NS,N,NDIM,UNFOLD,NMISSING,X,CONSTRAINTS) {
  res <-.C("copyFromR",
           as.integer(nslice),
           as.integer(nburn),
           as.integer(nrowX),
           as.integer(ncolX),
           as.integer(NS),
           as.integer(N),
           as.integer(NDIM),
           as.integer(UNFOLD),
           as.integer(NMISSING),
           as.double(X),
           as.double(CONSTRAINTS))
}
do_lbfgs <- function(kpnp,kpnq,yrotate,rmatrix){
  .C("mainlbfgs",
     as.integer(kpnp),
     as.integer(kpnq),
     as.double(yrotate),
     as.double(rmatrix))
}
do_sliceu <- function(theta,thetanow2,theta1000,ssenow,XTRUE,thetaLeft,thetaRight,WW,PP,XCOORDS,SIGMAPRIOR){
  .C("sliceunfolding",
     as.double(theta),
     as.double(thetanow2),
     as.double(theta1000),
     as.double(ssenow),
     as.double(XTRUE),
     as.double(thetaLeft),
     as.double(thetaRight),
     as.double(WW),
     as.integer(PP),
     as.double(XCOORDS),
     as.double(SIGMAPRIOR))
}
#
#


#cep73 <- read.spss("Encuesta CEP 73 Abril 2015.sav", use.value.labels = TRUE, to.data.frame = TRUE)
#CARGAR DATOS
cep73 <- read.csv("~/Dropbox/encuesta_cep73_abril2015/Encuesta CEP 73 Abril 2015.csv")

#CAMBIAR NOMBRES DE VARIABLES TERMOMETRO POR POLITICOS
library(plyr)
cep73 <- rename(cep73, 
       c("MB_P19_A"="A Allamand", 
         "MB_P19_B"="I Allende", 
         "MB_P19_C"="O Andrade", 
         "MB_P19_D"="A Arenas",
         "MB_P19_E"="M Bachelet",
         "MB_P19_F"="M Enriquez-O", 
         "MB_P19_G"="C Escalina", 
         "MB_P19_H"="A Espina", 
         "MB_P19_I"="N Eyzaguirre", 
         "MB_P19_J"="G Girardi", 
         "MB_P19_K"="R LagosE", 
         "MB_P19_L"="R LagosW", 
         "MB_P19_M"="E Matthei", 
         "MB_P19_N"="C Mockenber", 
         "MB_P19_O"="C Montes", 
         "MB_P19_P"="MJ Ossandon", 
         "MB_P19_Q"="R Penailillo", 
         "MB_P19_R"="L Perez", 
         "MB_P19_S"="S Pinera", 
         "MB_P19_T"="J Pizarro", 
         "MB_P19_U"="J Quintana", 
         "MB_P19_V"="X Rincon", 
         "MB_P19_W"="G Tellier", 
         "MB_P19_X"="C Toha", 
         "MB_P19_Y"="C Vallejo", 
         "MB_P19_Z"="A Velasco", 
         "MB_P19_AA"="Von Baer", 
         "MB_P19_BB"="I Walker", 
         "MB_P19_CC"="P Walker", 
         "MB_P19_DD"="Ma Nunez", 
         "MB_P19_EE"="F Kast"))
colnames(cep73[35:65])


#MATRIZ DE TERMOMETRO
T <- cep73[,35:65]
T[T > 5] <- NA
colSums(is.na(T))

#Politicos con mas notas
a <- as.data.frame(colSums(!is.na(T)))
a[order(a[,1]),] 
cutoff <- 850
T <- T[,colSums(!is.na(T))>=cutoff]
#Agregar a Kast
kast <- subset(cep73,select=c("F Kast"))
T <- cbind(T,kast)
T[T > 5] <- NA
#fin de agregar
T <- as.matrix(T)
a <- colnames(TT)
b<- c("AA","IA","MB","MEO","RLE","RLW","EM","SP","CT","CV","AV","FK")

#sacar personas sin notas
cutoff <- 12
T <- T[rowSums(!is.na(T))>=cutoff,]
T <- T*20
T <- (100-T)/50

########PARAMETROS
nrowX <- nrow(T)
ncolX <- ncol(T)
nburn <- 500
nslice <- 1500
NS <- 2
N <- NS*(nrowX+ncolX) - ((NS*(NS+1))/2)
NDIM <- NS*(nrowX+ncolX) - (NS-1)
UNFOLD <- 1
#  ROWS MUST HAVE LESS THAN 7 MISSING ENTRIES -  COLUMNS MUST HAVE AT LEAST 1/4 NON-MISSING ENTRIES (THIS IS HARD WIRED IN THE C CODE)
NMISSING <- 7


#
TT <- T
TT[is.na(TT)] <- -999.0 # No est?? funcionando
X <- as.vector(t(TT))


#  CONSTRAINTS
#
CONSTRAINTS <- rep(1,NS*(nrowX+ncolX))
#
if (NS==1){
  CONSTRAINTS[NS*(nrowX+ncolX)] <- 0
}
if (NS==2){
  CONSTRAINTS[(NS*(nrowX+ncolX)-NS):(NS*(nrowX+ncolX))] <- 0
}
if (NS==3){
  CONSTRAINTS[(NS*(nrowX+ncolX)-4):(NS*(nrowX+ncolX))] <- 0
  CONSTRAINTS[(NS*(nrowX+ncolX)-6)] <- 0
}
#


#  SMACOF
#
weightmat <- rep(1,nrowX*ncolX)
dim(weightmat) <- c(nrowX,ncolX)
weightmat[is.na(T)] <- 0
#smacofRect(delta, ndim = 2, circle = c("none","row","column"), weightmat = NULL, init = NULL, verbose = FALSE, itmax = 1000, reg = 1e-6, eps = 1e-6)
SMACOF.result <- smacofRect(delta=TT, ndim=NS,weightmat=weightmat, itmax=10000)
#
zmetric <- as.numeric(t(SMACOF.result$conf.col))
xmetric <- as.numeric(t(SMACOF.result$conf.row))
rmatrix <- c(zmetric,xmetric)
rmatrix[(NS*(nrowX+ncolX)-NS):(NS*(nrowX+ncolX))] <- 0
yrotate <- rep(0,(NS*(nrowX+ncolX)))

zz <- SMACOF.result$conf.col
xx <- SMACOF.result$conf.row


#GRAFICO
# plot(xx[,1],xx[,2],type="n",asp=1,
#      main="",
#      xlab="",
#      ylab="",
#      #xlim=c(-1,1),ylim=c(-1.5,1.5),
#      cex=1.2,font=2)
# mtext("Derecha - Izquierda",side=1,line=2.75,cex=1.2)
# # y-axis title
# #mtext("Social/Estilo de vida",side=2,line=2.5,cex=1.2)
# points(xx[,1],xx[,2],col="gray50")
# #points(xx[,1][presidential.vote==2],xx[,2][presidential.vote==2],pch='N',col="gray50",font=2)
# #points(xx[,1][presidential.vote==1],xx[,2][presidential.vote==1],pch='H',col="gray65",font=2)
# #points(xx[,1][presidential.vote==3],xx[,2][presidential.vote==3],pch='W',col="gray80",font=2)
# points(zz[,1],zz[,2],col="red")
# text(zz[,1],zz[,2],a,col="black",cex=1,font=1)


######################
######################
#  L-BFGS
#
set_globals(nslice,nburn,nrowX,ncolX,NS,N,NDIM,UNFOLD,NMISSING,X,CONSTRAINTS)
#
lbfgs.result <- do_lbfgs(nrowX,ncolX,yrotate,rmatrix)
lbfgs.coords <- lbfgs.result[[3]]
dim(lbfgs.coords) <- c(NS,(nrowX+ncolX))
X3 <- t(lbfgs.coords)
lbfgs.stimuli <- X3[1:ncolX,]
lbfgs.individuals <- X3[(ncolX+1):(nrowX+ncolX),]

# plot(lbfgs.individuals[,1],lbfgs.individuals[,2],type="n",asp=1,
#      main="",
#      xlab="",
#      ylab="",
#      xlim=c(-1.5,1.5),ylim=c(-1.5,1.5),font=2)
# points(lbfgs.individuals[,1],lbfgs.individuals[,2],col="gray50")
# a <- colnames(TT)
# text(lbfgs.stimuli[,1],lbfgs.stimuli[,2],a, col="black",cex=1,font=1)



#######################
#######################
#  BAYESIAN UNFOLDING
#
theta <- c(runif(NDIM,min=-.5,max=.5)) # Random starts for Slice Sampler
#theta <- lbfgs.result[[3]] # Non-random starts for Slice Sampler
theta2 <- theta
theta1000 <- rep(0,nslice*NDIM)
dim(theta1000) <-c(nslice*NDIM,1)
ssenow <-rep(0,(2*(nslice+nburn)))
dim(ssenow) <-c((2*(nslice+nburn)),1)
XTRUE <- lbfgs.result[[3]]
thetaL <- rep(-99.0, NDIM)
thetaR <- rep(99.0, NDIM)
dim(thetaL) <- dim(thetaR) <- c(NDIM, 1)
thetaL[NDIM] <- 0.10
thetaR[NDIM] <- 0.50
WW <- 1.0
PP <- 3.0
XCOORDS <- rep(0,(nrowX+ncolX)*NS)
SIGMAPRIOR <- 100.0
#
result4 <- do_sliceu(theta,theta2,theta1000,ssenow,XTRUE,thetaL,thetaR,WW,PP,XCOORDS,SIGMAPRIOR)
#
#  length(result4[[3]]) = NDIM*nslice
#
#  CALCULATE THE MEANS FROM THE SLICE SAMPLER --  mean(result4[[4]][2501:4000]) is the mean of the variance term
#
sigma_squared_hat <- mean(result4[[4]][2501:4000])
sigma_squared_hat_sd <- sd(result4[[4]][2501:4000])
#
#  SAMPLES
#
samples <- matrix(result4[[3]], ncol=NDIM, byrow=TRUE)
#
stimuli <- vector("list",NS)
for(i in 1:NS){
  stimuli[[i]] <- samples[,(seq(i, ncolX*NS, by=NS))]
  stimuli[[i]] <- as.mcmc(stimuli[[i]])
}
#
individuals <- vector("list",NS)
for(i in 1:NS){
  individuals[[i]] <- samples[,(seq((ncolX*NS+i), NDIM, by = NS))]
  individuals[[i]] <- as.mcmc(individuals[[i]])
}
#
#  BAYESIAN UNFOLDING PLOT OF INDIVIDUALS AND STIMULI
#
z.onedim <- summary(stimuli[[1]])[[1]][,1]
z.twodim <- summary(stimuli[[2]])[[1]][,1]
x.onedim <- summary(individuals[[1]])[[1]][,1]
x.twodim <- summary(individuals[[2]])[[1]][,1]
x.twodim[length(x.onedim)] <- 0
#
#quartz()
plot(x.onedim, x.twodim, type="n", asp=1,
     main="",
     xlab="",
     ylab="",
 xlim=c(5.3,6.7), ylim=c(0.7,1.3), 
      font=2)
mtext("Izquierda - Derecha",side=1,line=2.75,cex=1.2)
mtext("Social/Estilo de vida",side=2,line=2.5,cex=1.2)
points(x.onedim,x.twodim,col="gray50",cex=0.5)
points(z.onedim,z.twodim,col="red",cex=0.5)#pch=c(1:12)
text(jitter(z.onedim),jitter(z.twodim),b,col="black",cex=0.5,font=2)
#text(z.onedim,z.twodim,b,col="black",cex=0.8,font=2)
#legend("topleft", colnames(T), col="black",text.font=1, pt.cex=0.1, inset=.01, bty="n")

#  BAYESIAN UNFOLDING PLOT OF STIMULI WITH 90% CREDIBLE INTERVALS
#
z.lower.onedim <- apply(stimuli[[1]], 2, function(x){quantile(x,0.05)})
z.upper.onedim <- apply(stimuli[[1]], 2, function(x){quantile(x,0.95)})
z.lower.twodim <- apply(stimuli[[2]], 2, function(x){quantile(x,0.05)})
z.upper.twodim <- apply(stimuli[[2]], 2, function(x){quantile(x,0.95)})

#  BAYESIAN UNFOLDING VARIANCE (FIRST DIMENSION) PLOT
x.lower.onedim <- apply(individuals[[1]], 2, function(x){quantile(x,0.05)})
x.upper.onedim <- apply(individuals[[1]], 2, function(x){quantile(x,0.95)})
x.lower.twodim <- apply(individuals[[2]], 2, function(x){quantile(x,0.05)})
x.upper.twodim <- apply(individuals[[2]], 2, function(x){quantile(x,0.95)})


# plot(density(x.upper.onedim-x.lower.onedim),
#      xlab="",
#      ylab="",
#      main="")
# # Main title
# mtext("First Dimension Uncertainty Estimates\nIndividuals on Top, Stimuli on Bottom",side=3,line=1.25,cex=1.2,font=2)
# # x-axis title
# mtext("Difference between Upper and Lower 90% Credible Intervals",side=1,line=2.75,cex=1.2)
# # y-axis title
# mtext("Density",side=2,line=2.5,cex=1.2)
# rug(z.upper.onedim-z.lower.onedim)
# #
# #
# #  BAYESIAN UNFOLDING VARIANCE (SECOND DIMENSION) PLOT
# #
# quartz()
# plot(density(x.upper.twodim-x.lower.twodim),
#      xlab="",
#      ylab="",
#      main="")
# # Main title
# mtext("Second Dimension Uncertainty Estimates\nIndividuals on Top, Stimuli on Bottom",side=3,line=1.25,cex=1.2,font=2)
# # x-axis title
# mtext("Difference between Upper and Lower 90% Credible Intervals",side=1,line=2.75,cex=1.2)
# # y-axis title
# mtext("Density",side=2,line=2.5,cex=1.2)
# rug(z.upper.twodim-z.lower.twodim)



#
#  BAYESIAN UNFOLDING PLOT OF RESPONDENT 6'S POSTERIOR DENSITY
#
nresp <- 300
resp.2.samples <- cbind(individuals[[1]][,nresp],
                        individuals[[2]][,nresp])
#par(mar = c(2.5,2,1,0.5))
sm.density(resp.2.samples,col="lightgray",xlab="Izquierda - Derecha",ylab="Social",phi=30,theta=60)#,xlim=c(-0.5,2),ylim=c(-2.25,0.25))
#
