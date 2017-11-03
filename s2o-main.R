rm(list = ls())

library("ggplot2")
library("reshape2")
library("deSolve")
library("rootSolve")

setwd("C:/Users/rosherbal/Dropbox/Kings/results/s2o")

path <- "C:/Users/rosherbal/Dropbox/Kings/scripts/R"
spath <- "C:/Users/rosherbal/Dropbox/Kings/results/s2o/"


source(paste(path, "fanalysis.R", sep = "/"))
source(paste(path, "S2O-models.R", sep = "/"))

#####
#variables
#####

#global needs for analysis...
mod = s2o.mod[["f5"]]#function
modPar = s2o.par[["di_nd2"]]#parameters
pp <- modPar$pp #initial parameters
ic <- modPar$ic

## time course...
tt <- c(0.0, 200,0.001)#t0, tfinal, step

#bifurcation analysis...
mol <- "PP" #selected molecule to track
vn1 <- "nt" #names of parameter 1 to change
vv1 <- seq(0,20,0.05) #define parameter points to check 
jac <- T#calculate the jacobian matrix?
st <- seq(0.0,4,0.05) #set of initial condidions of state variables to check

#####
#exploring different inicial conditions and parameter sets
######

#matrix of initial conditions from st
ics <- matrix(data = NA, nrow = length(st), ncol = length(ic), dimnames = list(NULL,names(ic)))
ics[,"OO"] <- st
ics[,"PP"] <- rev(st)
ics[,"PO"] <- as.numeric(pp["tot"]) - (st + rev(st))
ics[which(ics[,2]<0),2] <- 0
ics<-ics[-1,]
ics<-ics[-dim(ics)[1],]

ic <- ics[38,] #fixed inicial conditions

#create all combinations of parameter set
pps <- expand.grid(0:1,0:1,0:1,0:1,0:1,0:1,0:1,0:1,c(0,0.05),c(0,0.05),c(0,0.05),c(0,0.05),c(0,0.05),c(0,0.05),c(0,0.05),c(0,0.05))
colnames(pps) <-c("p0", "p1", "p2", "p3", "d0", "d1", "d2", "d3", "bp0", "bp1", "bp2", "bp3", "bd0", "bd1", "bd2", "bd3")
#select set of parameters 
pps <- pps[which(pps[,"p1"]==1 & pps[,"p0"]==1 & pps[,"p2"]==0 & pps[,"p3"]==1 & 
                 pps[,"d1"]==1 & pps[,"d0"]==0 & pps[,"d2"]==1 & pps[,"d3"]==1 &
                 pps[,"bp1"]>0 & pps[,"bp2"]==0 & pps[,"bp0"]>0 &   pps[,"bp3"]>0 & 
                 pps[,"bd2"]>0 & pps[,"bd0"]==0 & pps[,"bd1"]==0 & pps[,"bd3"]>0),]#
pps <- do.call(rbind,apply(pps,1,function(r){sr <- length(r)-length(r[which(r!=0)]);if(sr==6){return(r)}}))

extpp <-cbind(rep(pp["nt"],dim(pps)[1]),rep(pp["nd"],dim(pps)[1]))
extpp <- cbind(extpp,rep(pp["tot"],dim(pps)[1]))
colnames(extpp) <- c("nt","nd", "tot")
pps <- as.data.frame(cbind(pps, extpp), stringsAsFactors = F)

########
#Analysis
########

#####
#### time course diagrams
##############
pp["nt"]<- 4
tc <- tcf(ic=ic,tt=tt,mod=mod,pp=pp)#funciton in fanalysis.R
tc$OP <- pp["tot"] - apply(tc[,-1],1, sum)
tc$OP[tc$OP<0.0] <- 0.0
gtc <- melt(as.data.frame(tc), id.vars="time")

##############
#### bifurcation analysis
##############

#calculate the matrix of points 
system.time(
  lsps <- ssa1p(initCon=ics,parSet=pps,mxInCon=st,val1par=vv1,nam1par=vn1, fmod=mod)
)

#calculate bistable region for each paramSet
bsregion <- lapply(lsps,BiOsreg, mol=mol,nam1par=vn1, cutoffval=0.5)

#know if they are switches or oscillators
type <- do.call(c,lapply(lsps,function(df){cls <- unique(df$st); print(cls);
if(all(2%in%cls & 0%in%cls)){if(all(4%in%cls & 3%in%cls)){cond = "so"} else if(4%in%cls){cond = "ds"}else if(3%in%cls){cond = "hs"}else{cond = "s"}}
else if(3%in%cls){if(4%in%cls){cond = "do"}else{cond = "o"}}
else{cond=NA}; return(cond)}))

#save raw data
ba1p <- list(lsps,bsregion,type,pps)
names(ba1p) <- c("lvalues","lbsregion","type","pps")
save(ba1p, file = "C:/Users/rosherbal/Dropbox/Kings/results/s2o/ba1p_pps8_bd1_up2one.RData")

#remove miscalculations for bistable models
lsps <- lapply(1:length(ba1p[["lvalues"]]), rmmiscal, lsps=ba1p[["lvalues"]], bsre=ba1p[["lbsregion"]], nam1par=vn1)
lsps <- lapply(1:length(lsps),function(n,l){cbind(l[[n]],gt=rep(n,dim(lsps[[n]])[1]))},lsps)

##############
#### amplitude in oscillatory systems
##############

#select oscillatory systems
oscmod <- ba1p[["lvalues"]][which(ba1p[["type"]]%in%c("so","o","do"))]#grep("bd1", dict[,2])
names(oscmod)<-which(ba1p[["type"]]%in%c("so","o","do"))
"
if(length(oscmod) > 0){
  #calculate A and T for each set 
  AT[i] <- lapply(i,function(i,pps,infr,tms,inc,thrs=5,period=F){#AT names(oscmod)
    i<- as.numeric(i)
    inf_up <- max(ba1p[["lvalues"]][[i]][which(ba1p[["lvalues"]][[i]]$st==3),"nt"]) + 0.1
    inf_dw <- min(ba1p[["lvalues"]][[i]][which(ba1p[["lvalues"]][[i]]$st==3),"nt"])- 0.1
    print("#############")
    print(i)
    print(c(inf_dw,inf_up))
    print("#############")

    if((inf_up-inf_dw)<=2){infr<- seq(inf_dw,inf_up,0.1)}
    else{infr <- c(seq(inf_dw,(inf_dw+2), 0.1),seq((inf_dw+2.1), inf_up,0.5))}
    infr <- seq(inf_dw,inf_up,0.5) 
    m <- ATfun(inframe=infr,parset=pps[i,],timeseq=tms,incon=inc,thrs,period=F)
    m[,"gt"] <- i
    return(m)},ba1p[["pps"]],ba1p[["bsregion"]],tt,ics[41,])
  
  
 save(AT,file=paste(spath,"ba1p_bd1_up2one_AT.Rdata", sep = "/")) 
  
  names(AT)<- names(oscmod)
}

##########
# ploting
#########

dfAT <- do.call(rbind,AT)
dfsps <- do.call(rbind,lsps)

cols2 <-c("OO" = "#ae0033", "PP" = "#006d8f", "OP" = "#9d9d9d", "PO" = "#5e5e5e")

#for oscillators
ggplot()+
  theme_bw(base_size = 20) +
  xlab("Phosphate donor (AU)") + ylab("Steady State Concentration (AU)") + ggtitle("BP system")
  geom_line(data=gSPS[which(gSPS$st%in%c(3) & gSPS$variable%in%c("OO","PP")),], aes(nt,value, colour=factor(variable)), size = 1, linetype = "dashed" ) +
  geom_line(data=gSPS[which(gSPS$st%in%c(2,4)&gSPS$variable%in%c("OO","PP")),], aes(nt,value, colour=factor(variable)),size = 1.5) +
  geom_line(data=dfAT[which(dfAT$n%in%c("OO","PP")),], aes(as.numeric(p), as.numeric(An), colour=factor(n)), size = 0.5) +
  geom_line(data=dfAT[which(dfAT$n%in%c("OO","PP")),], aes(as.numeric(p),as.numeric(Ax), colour=factor(n)), size = 0.5) + 
  facet_wrap(~gt, nrow = 4) +
  scale_y_continuous(breaks=seq(min(st),max(st),1), limit=c(0, as.numeric(pp["tot"]))) +
  scale_x_continuous(breaks=seq(0,max(gSPS$nt),1)) +
  scale_colour_manual(values = cols2, name="Form")


#for switches
ggplot()+
  theme_bw(base_size = 20) +
  xlab("Phosphate donor (AU)") + ylab("Steady State Concentration (AU)") + ggtitle("BP system")
  geom_line(data=gSPS[which(gSPS$st%in%c(0,3) & gSPS$variable%in%c("OO","PP")),], aes(nt,value, colour=factor(variable)), size = 1, linetype = "dashed" ) +
  geom_line(data=gSPS[which(gSPS$st%in%c(2,4)&gSPS$variable%in%c("OO","PP")&gSPS$value<=0.5),], aes(nt,value, colour=factor(variable)),size = 1.5) +
  geom_line(data=gSPS[which(gSPS$st%in%c(2,4)&gSPS$variable%in%c("OO","PP")&gSPS$value>=0.5),], aes(nt,value, colour=factor(variable)),size = 1.5) +
  facet_wrap(~gt, nrow = 4) +
  scale_y_continuous(breaks=seq(min(st),max(st),1), limit=c(0, as.numeric(pp["tot"]))) +
  scale_x_continuous(breaks=seq(0,max(gSPS$nt),1)) +
  scale_colour_manual(values = cols2, name="Form")





