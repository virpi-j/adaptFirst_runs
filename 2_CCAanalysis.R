rm(list=ls())
gc()
library(data.table)
setwd("/scratch/project_2000994/PREBASruns/finRuns/Rsrc/virpiSbatch/")
.libPaths(c("/projappl/project_2000994/project_rpackages", .libPaths()))
libpath <- .libPaths()[1]

parPath <- "/scratch/project_2000994/PREBASruns/metadata/paramUnc/"
load(paste0("/scratch/project_2000994/PREBASruns/finRuns/uncRuns/regRuns/XY",2015,".rdata"))


load(paste0(parPath,"pCROB_unc.rdata"))
parindCrob <- parind
pCrobdim <- nrow(pCROBbirch)
load(paste0(parPath,"pPREL_unc.rdata"))
parindPrel <- parind
pPreldim <- nrow(pPREL_unc)
data <- read.delim(paste0(parPath,"Yasso15.dat"), header = TRUE, sep="\t") 
pYasdim <- length(data[[1]])
pYas_unc <- matrix(0,pYasdim,35)
library(stringr)
for(ind in 1:pYasdim){ # read parameter values to matrix
    pYas_unc[ind,] <-as.numeric(unlist(str_split(data[[1]][[ind]], pattern = "  ")))[2:36]
}
load(paste0("/scratch/project_2000994/PREBASruns/finRuns/uncRuns/regRuns/loadParids.rdata"))

# soilCt, organic soil carbon stock
load(paste0("/scratch/project_2000994/PREBASruns/finRuns/uncRuns/regRuns/soilCt_organic"))
EC<-rbind(as.matrix(EC1),as.matrix(EC2))

pCROBASr<-data.table()
pPRELr<-data.table()
pYASr<-data.table()
for(ij in 1:nsim){ 
  pCROBASr<-rbind(pCROBASr, t(c(pCROBpine[parids1[ij],],
                                              pCROBspruce[parids1[ij],],
                                              pCROBbirch[parids1[ij],])
          ))
  pPRELr <- rbind(pPRELr,t(pPREL_unc[parids2[ij],]))
  pYASr <- rbind(pYASr,t(pYas_unc[parids3[ij],]))
}
  
NoHarv<-F # include NoHarv-scenario?
resfolder<-"/scratch/project_2000994/PREBASruns/finRuns/uncRuns/regRuns/"
#pdf(file=paste0("/scratch/project_2000994/PREBASruns/finRuns/uncRuns/regRuns/RdIn_NoHarv",NoHarv,".pdf"))#,  width = 880, height = 680)
CCAyear<-2015
for(CCAyear in c(2015, 2035, 2050)){
  load(paste0("/scratch/project_2000994/PREBASruns/finRuns/uncRuns/regRuns/XY",CCAyear,".rdata"))
  if(!NoHarv) X<-X[X$indHarv!=1,]  
  Y<-Y[,c(1:4,6)]
  reg_info<-c(reg_info[c(1,3:19)],"Whole Country")
  devtools::source_url("https://raw.githubusercontent.com/ForModLabUHel/utilStuff/master/RedundancyIndexCCA/RdInd_calc.r")
  xx<-X$indRCP
  nx1<-which(xx==1)
  nx2<-which(xx==2)
  nx3<-which(xx==3)
  nrep<-length(nx2)/nsim
  xx[nx1]<-rep(climModids[1:nsim],times=nrep)
  xx[nx2]<-rep(climModids[1:nsim],times=nrep)
  xx[nx3]<-rep(climModids[1:nsim],times=nrep)
  X<-cbind(X[,c(-1,-7,-8)],indGCM=xx,X[,c(7,8)])
  Ppar<-T
  if(Ppar){
    nrep<-length(xx)/nsim
    indCrob<-pCROBASr#rep(x=parids1[1:nsim],times=nrep)
    nCrob <- ncol(indCrob)
    indPrel<-pPRELr#rep(x=parids2[1:nsim],times=nrep)
    nPrel <- ncol(pPRELr)
    indPYas<-pYASr#rep(x=parids3[1:nsim],times=nrep)
    nYas <- ncol(indPYas)
    indECorg <- t(EC[,1:nsim])
    nECorg <- ncol(indECorg)+2 #last two columns from X are EC for N20 and CH4
    indCorg <- as.matrix(soilCt[1:nsim],nsim,1)
    nCorg <- ncol(indCorg)
    rHarv <- as.matrix(harvestLimsr[1:nsim,])
    nrHarv <- ncol(rHarv)
    X<-cbind(X,indECorg=indECorg,indCrob=indCrob,indPrel=indPrel,indPYas=indPYas,indCorg=indCorg,rHarv=rHarv)
  #  X<-cbind(X,indCrob=indCrob,indPrel=indPrel,indPYas=indPYas)
  }
  RdIn<-list()
  for(r_no in unique(X$indReg)){
  nn<-which(X$indReg==r_no)
  Xi<-X[nn,-5] # exclude constant region index
  Yi<-Y[nn,]
  RdIn[[r_no]]<-RdInd_calc(Xi,Yi)
}
#png(file=paste0("RdIn",CCAyear,".png"),  width = 880, height = 680)
#par(mfrow=c(2,3))
#for(ij in 1:ncol(Y)){
#  barplot(RdIn[[r_no]][ij,],names.arg = names(X)[c(1:4,6)], main = names(Y)[ij])
#}
#dev.off()

RdInY<-list()
for(ij in 1:ncol(Y)){
  Rdrr<-data.table()
  for(r_no in unique(X$indReg)){
    Rdrr<-rbind(Rdrr, data.table(r_no=r_no,t(RdIn[[r_no]][ij,])))
  }
  names(Rdrr)[2:6]<-names(X)[c(1:4,6)]
  names(Rdrr)[2:ncol(Rdrr)]<-names(X)[-5]
  regorder<-match(c(1:2,4:20),as.vector(Rdrr$r_no))
  RdInY[[ij]] <- Rdrr[regorder,]
}
names(RdInY)<-names(Y)

ffont<-1.8
library("colorspace")
palette1<-turbo(10) #(11)
#palette1<-palette1[c(1,3,5,7,9,11,2,4,6,8,10)]
#palette1 <- rainbow(15) 
#palette1 <- viridis(8)

for(ij in 1:ncol(Y)){
  filee<-paste0(resfolder,"CCA",colnames(Y)[ij],CCAyear,".pdf")
  if(CCAyear==2015){
    #png(file = filee, width = 600, height = 1000)
    pdf(file = filee, width = 6, height = 10)
  } else {
    pdf(file = filee, width = 4, height = 10)
  }
  xplot<-matrix(0,nrow=nrow(RdInY[[ij]]), ncol=10)
  xplot[,1:5]<-as.matrix(RdInY[[ij]][,2:6])
  xtmp <- as.matrix(RdInY[[ij]][,c(-6:-1)])
  xtmpECorg <- xtmp[,which(grepl("EC",colnames(xtmp), fixed = TRUE))]#xtmp[,1:nECorg]
  xtmpCrob <- xtmp[,which(grepl("Crob",colnames(xtmp), fixed = TRUE))]#xtmp[,(nECorg+1):(nECorg+nCrob)]
  xtmpPrel <- xtmp[,which(grepl("Prel",colnames(xtmp), fixed = TRUE))]#xtmp[,(nECorg+nCrob+1):(nECorg+nCrob+nPrel)]
  xtmpYas <- xtmp[,which(grepl("Yas",colnames(xtmp), fixed = TRUE))]#xtmp[,(nECorg+nCrob+nPrel+1):(nECorg+nCrob+nPrel+nYas)]
  xtmpCorg <- xtmp[,which(grepl("indCorg",colnames(xtmp), fixed = TRUE))]#xtmp[,(nECorg+nCrob+nPrel+nYas+1):(nECorg+nCrob+nPrel+nYas+nCorg)]
  xtmprHarv <- xtmp[,which(grepl("rHarv",colnames(xtmp), fixed = TRUE))]#xtmp[,(nECorg+nCrob+nPrel+nYas+nCorg):(nECorg+nCrob+nPrel+nYas+nCorg+nrHarv)]
  for(r_no in 1:length(RdInY[[ij]]$r_no)){
    #xplot[r_no,6:11] <-c(max(xtmp2[r_no,]),max(xtmp3[r_no,]),max(xtmp4[r_no,]),max(xtmp1[r_no,]),max(xtmp5[r_no]),max(xtmp6[r_no,]))
    xplot[r_no,6:10] <-c(max(xtmpCrob[r_no,]),max(xtmpPrel[r_no,]),max(xtmpYas[r_no,]),
                         max(xtmpECorg[r_no,]),max(xtmprHarv[r_no,]))
    #    xplot[r_no,6:8] <-c(sum(xtmp1[r_no,]),sum(xtmp2[r_no,]),sum(xtmp3[r_no,]))
  }
  par(mai=c(1.02,0.82*3,0.82,0.82*2))
  xplot<-data.table(xplot)
  names(xplot)<-c(names(X)[c(1:2)],"RCP","Harv","GCM","pCrob","pPrel","pYas","pECorg","pHarv")
  par(mai=c(1.02,0.82,0.02,0.42))
  xlabs <- NA*c(1:length(reg_info))
  ll<-F
  if(CCAyear==2015){ 
    par(mai=c(1.02,0.82*ffont+2,0.02,0.42))
    xlabs<-reg_info[regorder]
  }
  if(CCAyear==2050) ll<-T

  barplot(t(xplot), beside=T,
          xlab="RdInd", #yaxt="no",
          #fig.width = 10, fig.height = 5,
          cex.lab=ffont, cex.axis=ffont, 
          cex.main=ffont, cex.sub=ffont,
         # xlim=c(0,0.9),
          cex.names=ffont, width = 1.5,
          fig.width = 8, fig.height = 5,
          #main = paste(names(Y)[ij], CCAyear),
          names.arg = xlabs,
          las=1, #col = rainbow(nrow(new_mat)))
          horiz=T,legend.text = ll,
          args.legend = list(#legend=c(names(X[,2:5]),"par"),
                              x = "topright",
#                              fill=palette1,
#                              border=palette1,
                             inset = c(0.05, 0.05), 
                            cex=ffont,box.lty=0),
          col=palette1[1:(ncol(xplot))],
          border=palette1[1:(ncol(xplot))]
  )  
  dev.off()

  }
}
#dev.off()
