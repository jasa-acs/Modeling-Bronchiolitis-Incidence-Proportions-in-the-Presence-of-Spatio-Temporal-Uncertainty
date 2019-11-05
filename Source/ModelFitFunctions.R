###################################################
## Calculate Overlap Probabilities for cases     ##
## get.case.overlap(obj) takes a list of cases   ##
## and calculates the percentage overlap of each ##
## case with the spatial grid as defined by      ##
## Equation (6) in the manuscript                ##
###################################################
get.case.overlap <- function(obj){
  
  ## Spatial Overlap
  glocs <- which(rdist(obj[,c("JitteredLon","JitteredLat")],
                       grid.cent.land) < (sqrt(2)*(grid.box.len+s.jit.len)))
  glocs <- glocs[!glocs%in%water]
  gwgts <- rep(0,length(glocs))
  for(g in glocs){
    unif.pts <- cbind(runif(1000,grid.cent.land[g,1]-grid.box.len,grid.cent.land[g,1]+grid.box.len),
                      runif(1000,grid.cent.land[g,2]-grid.box.len,grid.cent.land[g,2]+grid.box.len))
    gwgts[g==glocs] <- mean(rdist(obj[,c("JitteredLon","JitteredLat")],unif.pts)<=s.jit.len)
  }
  gcell <- glocs[which.max(gwgts)]
  # plot(grid.cent.land[,1], grid.cent.land[,2])
  # points(grid.cent.land[glocs,1], grid.cent.land[glocs,2], pch=19, col="blue")
  # points(obj[,1], obj[,2], pch=19, col="red")
  
  ## Temporal Overlap
  bday.jit.set <- seq.Date(obj[,"JitteredBDay"]-t.jit.len, obj[,"JitteredBDay"]+t.jit.len, by="day")
  inf.jit.set <- seq.Date(obj[,"JitteredInfection"]-t.jit.len, obj[,"JitteredInfection"]+t.jit.len, by="day")
  wlocs <- which(abs(obj[,"JitteredBDay"]-as.Date(win.strt, format="%m/%d/%Y")) <= 8 |
                   abs(obj[,"JitteredBDay"]-as.Date(win.end, format="%m/%d/%Y")) <= 8)
  bt.wgt <- matrix(0, nrow=length(wlocs), ncol=one.year)
  for(w in wlocs){
    dseq <- seq.Date(as.Date(win.strt[w], format="%m/%d/%Y"),
                     as.Date(win.end[w], format="%m/%d/%Y"), by="day")
    doverlap <- length(base::intersect(dseq,bday.jit.set))/length(dseq)
    for(w2 in 0:(one.year-1)){
      if((w+w2)<=length(win.strt)){
        dseq <- seq.Date(as.Date(win.strt[w+w2], format="%m/%d/%Y"),
                         as.Date(win.end[w+w2], format="%m/%d/%Y"), by="day")
        bt.wgt[w==wlocs,w2+1] <- doverlap*length(base::intersect(dseq,inf.jit.set))/length(dseq)
      }
    }
  }
  
  ## Choose birthday and infection window
  bwin <- arrayInd(which.max(bt.wgt), .dim=dim(bt.wgt))
  twin <- wlocs[bwin[1]] + bwin[2] - 1
  bwin <- wlocs[bwin[1]]
  
  ## Return all the info
  return(list(glocs=glocs, gwgts=gwgts, bt.wgt=bt.wgt, gcell=gcell, bwin=bwin, twin=twin))
}

#######################################################
## Calculate Overlap Probabilities for controls      ##
## get.control.overlap(obj) takes a list of controls ##
## and calculates the percentage overlap of each     ##
## case with the spatial grid as defined by          ##
## Equation (6) in the manuscript                    ##
#######################################################
get.control.overlap <- function(obj){
  
  ## Spatial Overlap
  glocs <- which(rdist(obj[,c("JitteredLon","JitteredLat")],
                       grid.cent.land) < (sqrt(2)*(grid.box.len+s.jit.len)))
  glocs <- glocs[!glocs%in%water]
  gwgts <- rep(0,length(glocs))
  for(g in glocs){
    unif.pts <- cbind(runif(1000,grid.cent.land[g,1]-grid.box.len,grid.cent.land[g,1]+grid.box.len),
                      runif(1000,grid.cent.land[g,2]-grid.box.len,grid.cent.land[g,2]+grid.box.len))
    gwgts[g==glocs] <- mean(rdist(obj[,c("JitteredLon","JitteredLat")],unif.pts)<=s.jit.len)
  }
  gcell <- glocs[which.max(gwgts)]
  
  ## Temporal Overlap
  bday.jit.set <- seq.Date(obj[,"JitteredBDay"]-t.jit.len, obj[,"JitteredBDay"]+t.jit.len, by="day")
  wlocs <- which(abs(obj[,"JitteredBDay"]-as.Date(win.strt, format="%m/%d/%Y")) <= 8 |
                   abs(obj[,"JitteredBDay"]-as.Date(win.end, format="%m/%d/%Y")) <= 8)
  wwgts <- rep(0, length(wlocs))
  for(w in wlocs){
    dseq <- seq.Date(as.Date(win.strt[w], format="%m/%d/%Y"),
                     as.Date(win.end[w], format="%m/%d/%Y"), by="day")
    wwgts[w==wlocs] <- length(base::intersect(dseq,bday.jit.set))/length(dseq)
  }
  bwin <- wlocs[which.max(wwgts)]
  
  ## Return all the info
  return(list(glocs=glocs, gwgts=gwgts, wwgts=wwgts, gcell=gcell, bwin=bwin))
  
}

################################################
## Function for Quicker Tabulating            ##
## takes spatial windows and temporal windows ##
## and counts the number of obs within        ##
## each time x spatial window                 ##
################################################
quick.tab <- function(bdays,cells){
  tab <- matrix(0,nrow=NT,ncol=G)
  DT <- data.table(s=cells,b=bdays)
  agg.cnts <- DT[,.N,by=names(DT)]
  tab[(agg.cnts$s-1)*NT+agg.cnts$b] <- agg.cnts$N
  return(tab)
}

##################################################
## Functions to Evaluate Like within each       ##
## for a list of controls/cases, calculates the ##
## difference in likelihood given a proposed    ##
## beta, eta and psi vector                     ##
## as defined by Equation (10) in the           ##
## manuscript                                   ##
##################################################
control.like.diff <- function(x,bvec=beta.vec,eta=eta.vec,psi=psi.vec){
  ## Find valid times
  times <- which(win.strt.cases%in%win.strt[x$bwin+(0:(one.year-1))])
  
  ## Calculate logit(prob) from proposal
  logit.prob.prop <- bvec$prop[1]+bvec$prop[2]*cur.so2[x$gcell,times]+bvec$prop[3]*cur.pm25[x$gcell,times] + 
    c(bvec$prop[4:(3+ncol(NorfolkCensus))]%*%t(NorfolkCensus[x$gcell,])) +
    Xw[times,]%*%eta$prop+sum(Xg[x$gcell,]*psi$prop)
  
  ## Calculate logit(prob) from current
  logit.prob.cur <- bvec$cur[1]+bvec$cur[2]*cur.so2[x$gcell,times]+bvec$cur[3]*cur.pm25[x$gcell,times] + 
    c(bvec$cur[4:(3+ncol(NorfolkCensus))]%*%t(NorfolkCensus[x$gcell,])) +
    Xw[times,]%*%eta$cur+sum(Xg[x$gcell,]*psi$cur)
  
  ## Back transform to probability
  prob.prop <- exp(logit.prob.prop)/(1+exp(logit.prob.prop))
  prob.cur <- exp(logit.prob.cur)/(1+exp(logit.prob.cur))
  
  ## Sum the 1-prob because we didn't see a case here
  return(sum(log(1-prob.prop))-sum(log(1-prob.cur)))
}

case.like.diff <- function(x,bvec=beta.vec,eta=eta.vec,psi=psi.vec){
  ## Find valid times
  times <- which(win.strt.cases%in%win.strt[x$bwin+(0:(one.year-1))])
  if(length(times)==0){
    return(0) 
  } else {
    case.time <- which(win.strt.cases%in%win.strt[x$twin]) #-24
    if(length(case.time)>0){
      times <- times[times!=case.time]
    }
    
    ## Calculate logit(prob) under proposed values
    logit.prob.prop <- bvec$prop[1]+bvec$prop[2]*cur.so2[x$gcell,times]+bvec$prop[3]*cur.pm25[x$gcell,times]+
      c(bvec$prop[4:(3+ncol(NorfolkCensus))]%*%t(NorfolkCensus[x$gcell,]))+
      Xw[times,]%*%eta$prop+sum(Xg[x$gcell,]*psi$prop)
    if(length(case.time)!=0){
      c.logit.prob.prop <- bvec$prop[1]+bvec$prop[2]*cur.so2[x$gcell,case.time]+bvec$prop[3]*cur.pm25[x$gcell,case.time]+
        c(bvec$prop[4:(3+ncol(NorfolkCensus))]%*%t(NorfolkCensus[x$gcell,])) +
        sum(Xw[case.time,]*eta$prop)+sum(Xg[x$gcell,]*psi$prop)
      c.prob.prop <- exp(c.logit.prob.prop)/(1+exp(c.logit.prob.prop))
    } else {
      c.prob.prop <- 1
    }
    
    ## Calculate logit(prob) under current values
    logit.prob.cur <- bvec$cur[1]+bvec$cur[2]*cur.so2[x$gcell,times]+bvec$cur[3]*cur.pm25[x$gcell,times]+
      c(bvec$cur[4:(3+ncol(NorfolkCensus))]%*%t(NorfolkCensus[x$gcell,]))+
      Xw[times,]%*%eta$cur+sum(Xg[x$gcell,]*psi$cur)
    if(length(case.time)!=0){
      c.logit.prob.cur <- bvec$cur[1]+bvec$cur[2]*cur.so2[x$gcell,case.time]+bvec$cur[3]*cur.pm25[x$gcell,case.time]+
        c(bvec$cur[4:(3+ncol(NorfolkCensus))]%*%t(NorfolkCensus[x$gcell,])) +
        sum(Xw[case.time,]*eta$cur)+sum(Xg[x$gcell,]*psi$cur)
      c.prob.cur <- exp(c.logit.prob.cur)/(1+exp(c.logit.prob.cur))
    } else {
      c.prob.cur <- 1
    }
    
    ## Back transform to probability
    prob.prop <- exp(logit.prob.prop)/(1+exp(logit.prob.prop))
    prob.cur <- exp(logit.prob.cur)/(1+exp(logit.prob.cur))
    
    ## Sum the log(prob) @ time + log(1-prob) because we didn't see a case on other time periods
    return(log(c.prob.prop)+sum(log(1-prob.prop))-(log(c.prob.cur)+sum(log(1-prob.cur))))
  }
}

################################
## Update Lambda              ##
## Dirichlet update of Lambda ##
################################
update.lam <- function(cnts){
  gamvars <- matrix(rgamma(length(cnts),c(cnts)+1,1),nrow=nrow(cnts),ncol=ncol(cnts))
  return(gamvars/sum(gamvars))
}

############################################
## Draw time windows/regions for cases    ##
## a list of cases (x), draws new spatial ## 
## and temporal windows for true location ##
## as described in Section 4.1            ##
############################################
update.both.case <- function(x,bvec=beta.vec,eta=eta.vec,psi=psi.vec) { 
  ## Update Location
  case.time <- x$twin-24
  if (case.time <= 0) return(x) #if the case before 9/1/04
  logit.probs <- bvec[1]+bvec[2]*cur.so2[x$glocs,case.time]+bvec[3]*cur.pm25[x$glocs,case.time]+
    t(bvec[4:(3+ncol(NorfolkCensus))]%*%t(NorfolkCensus[x$glocs,]))+
    (Xw[case.time,]%*%eta)[1]+(Xg[x$glocs,]%*%psi)
  if(length(Xw[case.time,]%*%eta)>1) cat("Error: too long")
  probs <- (exp(logit.probs)/(1+exp(logit.probs)))*Lambda[x$bwin,x$glocs]*x$gwgts
  x$gcell <- sample(x$glocs,size=1,prob=probs)
  
  ##Update b/t
  if(length(x$wlocs)>1) {
    z <- c(x$wlocs[1]+0:23,x$wlocs[2]+0:23)
  } else { 
    z <- x$wlocs+0:23
  }
  
  ## get valid case/control times
  z.case <- z-24
  stf <- which(z.case>0 & z.case < 209) #only sample cases after 9/1/04
  z.case <- z.case[stf]
  if (length(z.case) < 2) return(x) #if only 1 time or no possible times
  
  stf2 <- which(z < 233) #only births before 5/1/13
  z <- z[stf2]
  
  ## get thetas/lambdas
  logit.theta <- bvec[1]+bvec[2]*cur.so2[x$gcell,z.case]+bvec[3]*cur.pm25[x$gcell,z.case]+
    bvec[4:(3+ncol(NorfolkCensus))]%*%t(NorfolkCensus[x$gcell,])+
    (Xw[z.case,]%*%eta)[1]+Xg[x$gcell,]%*%psi
  theta <- exp(logit.theta)/(1+exp(logit.theta))
  lambda <- Lambda[z,x$gcell]
  
  thta <- rep(0,48) #the impossible cases (to early/late) get probability of 0
  thta[stf] <- theta
  
  lmbda <- rep(0,48)
  lmbda[stf2] <- lambda 
  
  probs <- lmbda*thta*as.vector(t(x$bt.wgt))
  draw <- sample(1:length(probs),size=1,prob=probs) 
  ind <- arrayInd(draw,c(24,length(x$wlocs)))
  x$bwin <- x$wlocs[ind[2]]
  x$twin <- x$bwin+ind[1]-1
  
  return(x)
}

###############################################
## Draw time windows/regions for controls    ##
## a list of controls (x), draws new spatial ## 
## and temporal windows for true location    ##
## as described in Section 4.1               ##
###############################################
update.both.cntrl <- function(x) {
  probs <- x$wwgts*Lambda[x$wlocs,x$gcell]
  if (length(x$wlocs) > 1) x$bwin <- sample(x$wlocs,size=1,prob=probs) #times lambda?
  
  probs <- x$gwgts*Lambda[x$bwin,x$glocs]
  if (length(x$glocs) > 1) x$gcell <- sample(x$glocs,size=1,prob=probs)
  return(x)
}
