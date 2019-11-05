####################################################
## Code to Fit ST Unc Model to Bronchiolitis Data ##
####################################################

###################
## Preliminaries ##
###################

##Load Libraries
library(LatticeKrig)
library(data.table)
library(parallel)

## Source in Helper Functions
source("../Source/ModelFitFunctions.R")
source("../Source/AMCMCUpdate.R")

## Set jittering length and grid size
t.jit.len <- 7 #Not true value (masked for confidentiality)
grid.box.len <- 0.0018/2 #1/2 side length of spatial discretized grid box
s.jit.len <- 0.02 #Not true value (masked for confidentiality)

## Which model to run?
include.X.unc <- TRUE
include.ST.unc <- TRUE

## Set Number of cores to use for parallel processing
ncores.unc.update <- 1
ncores.param.update <- 1
one.year <- 24  #Number of time periods per year

## Load Census information
load("../Data/CensusData.RData") #Water locations are NA

############################################################
## Load City Grid  and keep only grid cells that are land ##
############################################################
load('../Data/Norfolk200mgrid.Rdata')
water <- which(is.na(rowSums(NorfolkCensus)))
lats <- c(36.82, 36.98)
lon <- c(-76.346, -76.133)
kp.grid <- which(apply(grid.cent,1,function(x){
  length(which(x[1]==grid.cent.land[,1] & x[2]==grid.cent.land[,2]))>0
}))
lon.locs <- unique(grid.cent[,1])
lat.locs <- unique(grid.cent[,2])

####################################
## Format PM25/SO2 estimates ##
####################################
load('../Data/NorfolkFullPredictionsSO2.Rdata')
load('../Data/NorfolkFullPredictionsPM25.Rdata')

D <- rdist(grid.cent.land,coords) #assign to closest of grid coordinates
grp <- apply(D,1,which.min)

tm <- c(which(as.Date(dates)==as.Date("2003-09-01")):248)

so2 <- so2[grp,tm]
so2.var <- so2.var[grp,tm]
pm25 <- pm25[grp,tm]
pm25.var <- pm25.var[grp,tm]
so2 <- so2[,-(1:one.year)]
so2.var <- so2.var[,-(1:one.year)]
pm25 <- pm25[,-(1:one.year)]
pm25.var <- pm25.var[,-(1:one.year)]

#####################################
## Remove the first year of cases  ##
## Because we won't have counts of ##
## number of susceptibles          ## 
#####################################

## Temporal Windows
strt.dates <- c("01/01/","01/15/","02/01/","02/15/","03/01/",
                "03/15/","04/01/","04/15/","05/01/","05/15/","06/01/","06/15/","07/01/",
                "07/15/","08/01/","08/15/","09/01/","09/15/","10/01/",
                "10/15/","11/01/","11/15/","12/01/","12/15/")
fin.dates <- c("01/14/","01/31/","02/14/","02/28/","03/14/","03/31/",
               "04/14/","04/30/","05/14/","05/31/","06/14/","06/30/","07/14/","07/31/",
               "08/14/","08/31/","09/14/","09/30/","10/14/","10/31/",
               "11/14/","11/30/","12/14/","12/31/")
win.strt <- c()
win.end <- c()
for(y in 2003:2013){
  win.strt <- c(win.strt,paste(strt.dates,y,sep=""))
  win.end <- c(win.end,paste(fin.dates,y,sep=""))
}
date.cutoff <- which(win.strt=="05/01/2013")
date.strt <- which(win.strt=="08/15/2003")
win.strt <- win.strt[-(date.cutoff:length(win.strt))]
win.end <- win.end[-(date.cutoff:length(win.end))]
win.strt <- win.strt[-(1:date.strt)]
win.end <- win.end[-(1:date.strt)]
win.strt.cases <- win.strt[-(1:one.year)]
win.end.cases <- win.end[-(1:one.year)]

########################################
## Load in Data on Cases and Controls ##
########################################
cases <- read.table("../Data/SimulatedCasesInfo.txt", header=TRUE)
cases$JitteredBDay <- as.Date(cases$JitteredBDay)
cases$JitteredInfection <- as.Date(cases$JitteredInfection)
controls <- read.table("../Data/SimulatedControlInfo.txt", header=TRUE)
controls$JitteredBDay <- as.Date(controls$JitteredBDay)
N <- nrow(cases) + nrow(controls)

###############################################
## Calculate Overlap Probabilities for cases ##
###############################################
cases.sbt <- mclapply(split(cases, 1:nrow(cases)), 
                      get.case.overlap, 
                      mc.cores=ncores.param.update)

##################################################
## Calculate Overlap Probabilities for controls ##
##################################################
controls.sb <- mclapply(split(controls, 1:nrow(controls)), 
                        get.control.overlap, 
                        mc.cores=ncores.param.update)

##################################
## Useful Numbers for Reference ##
##################################
G <- nrow(grid.cent.land)
NT <- length(win.strt)
NT.case <- NT-one.year

####################
## Initiate Delta ##
####################
delta <- length(controls.sb)+length(cases.sbt)

#####################################
## Function for Quicker Tabulating ##
#####################################
b.case <- sapply(cases.sbt, function(x){x$bwin})
b.cntrl <- sapply(controls.sb, function(x){x$bwin})
s.case <- sapply(cases.sbt, function(x){x$gcell})
s.cntrl <- sapply(controls.sb, function(x){x$gcell})
cntrl.cnts <- quick.tab(c(b.case,b.cntrl),c(s.case,s.cntrl))
N0 <- matrix(0,nrow=NT.case,ncol=G)
for(i in 1:nrow(N0)){
  N0[i,] <- colSums(cntrl.cnts[(i+1):(i+one.year),])
}
t.case <- sapply(cases.sbt, function(x){x$twin})
case.cnts <- quick.tab(t.case,s.case)
case.cnts <- case.cnts[-(1:one.year),]
kp.spat <- which(colSums(N0)>0)

###########################
## Define Basis matrices ##
###########################

## M matrix for Spatial Piece
source("../Source/GetAdjMatrix.R")

## Moran Basis Functions
A <- AdjMat(length(lat.locs),length(lon.locs))
unitnum <- A$unitlabel
A <- A$A[kp.grid,kp.grid]
ones <- matrix(1,nrow=nrow(A),ncol=1)
P.orth <- diag(nrow(A))-(ones%*%t(ones))/sum(ones)
M <- eigen(P.orth%*%A%*%P.orth)
M$vectors <- M$vectors[,order(M$values,decreasing=TRUE)]
M$values <- sort(M$values,decreasing=TRUE)
which.pos <- which(M$values>0)
Xg <-M$vectors[,1:100]
priprec.g <- t(Xg)%*%(diag(rowSums(A))-A)%*%Xg
tau.g <- 0.01

## M matrix for Temporal Effects
Aw <- matrix(0,nrow=NT-one.year,ncol=NT-one.year)
for(r in 1:(nrow(Aw)-1)){
  Aw[r,r+1] <- 1
}
Aw <- Aw+t(Aw)
ones <- matrix(1,nrow=nrow(Aw),ncol=1)
P.orth <- diag(nrow(Aw))-(ones%*%t(ones))/sum(ones)
M <- eigen(P.orth%*%Aw%*%P.orth)
M$vectors <- M$vectors[,order(M$values,decreasing=TRUE)]
M$values <- sort(M$values,decreasing=TRUE)
Xw <-M$vectors[,1:75]
priprec.w <- t(Xw)%*%(diag(rowSums(Aw))-Aw)%*%Xw
tau.w <- 0.001

#####################
## Starting Values ##
#####################
beta.vec <- matrix(c(-5,rep(0,2+ncol(NorfolkCensus))),ncol=1)
tm.log.odds <- rowSums(case.cnts)/rowSums(N0)+.001 
tm.log.odds <- log(tm.log.odds/(1-tm.log.odds))-beta.vec[1]
eta.vec <- matrix(solve(t(Xw)%*%Xw+tau.w*priprec.w)%*%t(Xw)%*%tm.log.odds,nrow=ncol(Xw),ncol=1)
spat.log.odds <- colSums(case.cnts)/colSums(N0)+0.001
spat.log.odds <- log(spat.log.odds/(1-spat.log.odds))-beta.vec[1]
psi.vec <- matrix(solve(t(Xg[kp.spat,])%*%Xg[kp.spat,]+tau.g*priprec.g)%*%t(Xg[kp.spat,])%*%spat.log.odds[kp.spat],nrow=ncol(Xg),ncol=1)

###################
## Update Lambda ##
###################
system.time(Lambda <- update.lam(cntrl.cnts))

###################
## MCMC Settings ##
###################
source("../Source/AMCMCUpdate.R")
burn <- 1000000
num.it <- 1000000
thin <- 1
beta.ind <- 1:length(beta.vec)
eta.ind <- length(beta.vec)+(1:length(eta.vec))
psi.ind <- max(eta.ind)+(1:length(psi.vec))
amcmc <- list(mn=matrix(0,nrow=max(psi.ind),ncol=1),
              var=matrix(0,nrow=max(psi.ind),ncol=max(psi.ind)))
amcmc.it <- 250
kp <- 0
kpseq <- seq(burn+thin,burn+thin*num.it,by=thin)
save.seq <- round(seq(.02,.98,by=.02)*(burn))
printseq <- kpseq[round(seq(1,length(kpseq),length=20))]

############################
## Matrices to Hold Draws ##
############################
beta.draws <- matrix(0,nrow=num.it,ncol=length(beta.vec))
eta.draws <- matrix(0,nrow=num.it,ncol=nrow(eta.vec))
psi.draws <- matrix(0,nrow=num.it,ncol=nrow(psi.vec))

#####################
## Begin MCMC Loop ##
#####################
system.time({
  for(it in 1:(burn+thin*num.it)){
    
    if(include.X.unc){
      ## Draw new pollution - X matrix
      cur.so2 <- matrix(so2 + sqrt(so2.var)*rnorm(length(so2)),
                        nrow=nrow(so2),ncol=ncol(so2))
      cur.pm25 <- matrix(pm25 + sqrt(pm25.var)*rnorm(length(pm25)),
                         nrow=nrow(so2),ncol=ncol(so2))
    }
    
    ## Update lambda
    s.case <- sapply(cases.sbt,function(x) x$gcell)
    b.case <- sapply(cases.sbt,function(x) x$bwin)
    s.cntrl <- sapply(controls.sb,function(x) x$gcell)
    b.cntrl <- sapply(controls.sb,function(x) x$bwin) 
    cntrl.cnts <- quick.tab(c(b.case,b.cntrl),c(s.case,s.cntrl))
    Lambda <- update.lam(cntrl.cnts)
    
    if(include.ST.unc){
      ## Update Control Spatial Locations and Bdays
      controls.sb <- mclapply(controls.sb,update.both.cntrl,mc.cores=ncores.unc.update)
    
      ## Update Cases Spatial Locations, BDays, and RSV Dates jointly
      cases.sbt <- mclapply(cases.sbt,update.both.case,mc.cores=ncores.unc.update)
    }
    
    ## Update beta, eta and psi jointly
    prop.var <- (0.00001^2)*diag(nrow(amcmc$mn))
    if(it>amcmc.it){
      prop.var <- (2.4^2/nrow(amcmc$mn))*(amcmc$var+prop.var)
    }
    all.prop <- c(beta.vec,eta.vec,psi.vec)+t(chol(prop.var))%*%rnorm(nrow(amcmc$mn))
    prop.beta <- matrix(all.prop[beta.ind],ncol=1)
    prop.eta <- matrix(all.prop[eta.ind],ncol=1)
    prop.psi <- matrix(all.prop[psi.ind],ncol=1)
    bvecs <- list(prop=prop.beta,cur=beta.vec)
    etavecs <- list(prop=prop.eta,cur=eta.vec)
    psivecs <- list(prop=prop.psi,cur=psi.vec)
    cntrl.llike.diff <- sum(as.numeric(mclapply(controls.sb,control.like.diff,bvec=bvecs,eta=etavecs,psi=psivecs,mc.cores=ncores.param.update)))
    case.llike.diff <- sum(as.numeric(mclapply(cases.sbt,case.like.diff,bvec=bvecs,eta=etavecs,psi=psivecs,mc.cores=ncores.param.update)))
    prior.prop <- sum(dnorm(prop.beta,0,5,log=TRUE))-(tau.w/2)*t(prop.eta)%*%priprec.w%*%prop.eta-
      (tau.g/2)*t(prop.psi)%*%priprec.g%*%prop.psi
    prior.cur <- sum(dnorm(beta.vec,0,5,log=TRUE))-(tau.w/2)*t(eta.vec)%*%priprec.w%*%eta.vec-
      (tau.g/2)*t(psi.vec)%*%priprec.g%*%psi.vec
    MH.ratio <- cntrl.llike.diff+case.llike.diff+prior.prop-prior.cur
    if(log(runif(1,0,1))<MH.ratio){
      beta.vec <- prop.beta
      eta.vec <- prop.eta
      psi.vec <- prop.psi
    }
    amcmc <- AMCMC.update(rbind(beta.vec,eta.vec,psi.vec),amcmc$mn,amcmc$var,it)
    
    
    ## Update Delta
    delta <- rgamma(1,N+1/2,1)
    
    ## Retain Draw if Necessary
    if(it%in%kpseq){
      kp <- kp+1
      if(it%in%printseq){
        cat(paste(round(100*kp/num.it,2),"of Draws Obtained (Post Burn)\n"))
      }
      beta.draws[kp,] <- beta.vec
      eta.draws[kp,] <- eta.vec
      psi.draws[kp,] <- psi.vec
      save(file="./RSVResults.RData",list=c("beta.draws","eta.draws","psi.draws","Xg","Xw"))
    }
    
    if(it%in%save.seq){
      ## Save current draw so can start where left off if necessary
      save(file="./CurrentDraw.RData",list=c("beta.vec","psi.vec","eta.vec"))
    }
    
    
  }
})









