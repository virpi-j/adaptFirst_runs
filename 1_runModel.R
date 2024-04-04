#station_id <- 1 # This to localsettings! Station id 1 to 6
set.seed(1)
print(paste("station",station_id))
stations <- data.frame(name=c("Vantaa_Helsinki","Jokioinen","Jyvaskyla","Kajaani",
                              "Sodankyla","Utsjoki"),
              location=c("Helsinki_Vantaa_lentoasema", # Uusimaa 1
                              "Jokioinen_Ilmala", # Kanta-Hame 9
                              "Jyvaskyla_lentoasema", # Keski-Suomi 6 
                              "Kajaani_lentoasema", # Kainuu 16
                              "Sodankyla_Tahtela", # Lappi 8
                              "Utsjoki_Kevo"), # Lappi 8
              ID = c(1:6), 
              x = c(24.96, 23.48, 25.67, 27.67, 26.63, 27.01),
              y = c(60.33, 60.80, 62.4, 64.28, 67.37, 69.76))

r_nos_stations <- list()
r_nos_stations[[1]] <- c(1)#,9,11,13,15)
r_nos_stations[[2]] <- c(11,1,9)#,1,4,9,13,15)
r_nos_stations[[3]] <- c(6)#,4,13,17,7,3,12)
r_nos_stations[[4]] <- c(16,19)#,7,18,19)
r_nos_stations[[5]] <- c(8)
r_nos_stations[[6]] <- c(8)

r_no = region = r_nos_stations[[station_id]][1] # region ID

stat_name <- stations[station_id,"name"]
nYears <- 2050-2015

###### Calibrated PREBAS #################################
if(calibratedPREBAS){
  load("../modelParameters/pCROBAS_newVcalP_CN2703.rdata")
  load("../modelParameters/pPRELES_newVcalP_CN2703.rdata")
  pCROB_new <- pCROBAS_newVcalP_CN
  pCROB_cc <- pCROBAS_newVcalP_CN
  
  # Annikki's corrections
  pCROB_new[54,2] <- 1.5033  # tissue N
  pCROB_new[55,2] <- 0.0113  # tissue N
  pCROB_new[54,3] <- 3.00    # tissue N birch
  pCROB_new[55,3] <- 0.037   # tissue N
  pCROB_new[57,3] <- 0.2     # tissue N
  pCROB_new[58,3] <- 0.05     # tissue N
  
  pCROB_new[63,2] <- -4.866   # restricted N uptake, intercept
  pCROB_new[64,2] <- 0.993    # restricted N uptake, coeff
  pCROB_new[63,1] <- -4.596   # restricted N uptake, intercept
  pCROB_new[64,1] <- 0.887    # restricted N uptake, coeff
  pCROB_new[63,3] <- -4.696   # restricted N uptake, intercept
  pCROB_new[64,3] <- 1.282       # restricted N uptake, coeff
  pCROB_new[41,1:3] <- restrictionSwitch        # set restriction switch active, 1 = default. 
  #   - This should be automatic in the current New Version that can be downloaded
  
  pCROB <- pCrobasX <- pCROB_new
  
  pPREL_new <- pPRELES_newVcalP_CN
  pPREL <- pPRELES <- pPREL_new
}  

devtools::source_url("https://raw.githubusercontent.com/virpi-j/adaptFirst_runs/master/settings.R")
parsCN_new_alfar <- parsCN_alfar

if(calibratedPREBAS){
  parsCN_cc_alfar <- parsCN_alfar
  parsCN_new_alfar <- parsCN_alfar
  
  parsCN_cc_alfar <- parsCN_alfar
  parsCN_cc_alfar[1,3]<- 2*parsCN_alfar[1,1]
  parsCN_new_alfar <- parsCN_cc_alfar
} 

xy_UTM <- data.frame()
for(ij in 1:nrow(stations)){
  xy <- stations[ij,c("ID","x","y")]
  coordinates(xy) <- c("x","y")
  proj4string(xy) <- CRS("+proj=longlat +datum=WGS84") 
  res <- spTransform(xy, crsX)#CRS(paste("+proj=utm +zone=",35," ellps=WGS84",sep='')))
  xy_UTM <- rbind(xy_UTM, station_coords<-as.data.frame(res)[2:3])
}
colnames(xy_UTM)<-c("x_UTM","y_UTM")
stations <- cbind(stations,xy_UTM)
print(stations)
if(toRaster & station_id==1) write.csv(stations[,c("ID","name","x_UTM","y_UTM")],file="weather_stations.csv")

##source_url("https://raw.githubusercontent.com/virpi-j/adaptFirst_runs/master/05_create_CO2cols.R")
CO2_RCPyears <- read.csv2(file=paste0(climatepath,"co2_concentrations_PREBAS.csv"),header=T,sep = ";")
co2Names<-names(CO2_RCPyears)[2:ncol(CO2_RCPyears)]
names(CO2_RCPyears)[1]<-"year"

#load(paste0("/scratch/project_2000994/PREBASruns/finRuns/input/maakunta/maakunta_",r_no,"_IDsTab.rdata"))
#data.all <- cbind(data.all,data.IDs[match(data.all$segID, data.IDs$maakuntaID),4:5])
print(paste("PREBAS version",vPREBAS))

## Weather station coordinates:

xy <- stations[station_id,c("ID","x","y")]
coordinates(xy) <- c("x","y")
proj4string(xy) <- CRS("+proj=longlat +datum=WGS84") 
res <- spTransform(xy, crsX)#CRS(paste("+proj=utm +zone=",35," ellps=WGS84",sep='')))
station_coords<-as.data.frame(res)[2:3]
#station_coords <- LongLatToUTM(xy) # dd to UTM
print(station_coords)

d<- sqrt((data.all$x - station_coords$x)^2+(data.all$y - station_coords$y)^2)
nn.d <- order(d, decreasing=F)[1:nSitesRun]

print(range(d))

ops <- list(data.all[nn.d,])
if(toRaster){
  ndat <- sample(1:nrow(data.all),10000)
  plot(data.all$x[ndat],data.all$y[ndat])
  points(station_coords$x,station_coords$y,cex=2,col="red")
  ndat <- sample(1:nrow(ops[[1]]),1000)
  points(ops[[1]]$x[ndat],ops[[1]]$y[ndat],col="green")
}
print(paste("Weather station",stations[station_id,"name"]))
print(paste("Forest data: maximum distance from weather station",
            round(max(d[nn.d])/1000,2),"km"))
print(paste("Area of forest within the closest",nSitesRun,"segments is", round(sum(sum(ops[[1]]$area)),2),"hectares"))
if(toRaster){
  source_url("https://raw.githubusercontent.com/virpi-j/adaptFirst_runs/master/2_createRaster.R")
}
if(toRaster) dev.off()

#### Climate selection #############################################
rcps = "CurrClim" 

#############
#pathtoken = "/scratch/project_2000994/PREBASruns/finRuns/"
if(climScen > 0){
  climatepath = "/scratch/project_2000994/RCP/"
  climMod <- c("CanESM2.","CNRM.","GFDL.","HadGEM2.","MIROC.")
  rcpx <- c("rcp26","rcp45","rcp85")
  rcps <- rcpsFile <-paste0(climMod[ClimModid],rcpx[climScen])
  rcpsName <- rcps
  CO2fixed <- 0
} else if(climScen<0){
  if(CO2fixed==0){
    rcpsFile <- paste0(stat_name,"_1991_2100_constant_change_v3.csv")
    Co2Col<-which(co2Names=="X2005.fixed")
    rcpsName <- "constant"
  } else {
    rcpsFile <- paste0(stat_name,"_1991_2100_seasonally_perturbed_v1.csv")
    rcpsName <- "perturbed"
    Co2Col<-CO2fixed
  }
  weatherData<-read.csv2(file=paste0(climatepath,rcpsFile),sep = ",")
  print(paste("CO2scenario", names(CO2_RCPyears)[Co2Col+1]))
}
print(paste("Climate scenario",rcpsName))

if(climScen<0){
  deltaP <- unique(weatherData$Pchange)
  deltaT <- unique(weatherData$deltaT)
  if(outType=="testRun"){
    deltaT<-deltaT[c(1,2,length(deltaT))]
    deltaP<-deltaP[c(1,2,length(deltaP))]
  }
  
  
  deltaTP <- matrix(0,2,length(deltaP)*length(deltaT))
  index <- 1
  for(iT in 1:length(deltaT)){
    for(iP in 1:length(deltaP)){
      deltaTP[1,index] <- deltaT[iT]
      deltaTP[2,index] <- deltaP[iP]
      index <- index+1
    }
  }
  deltaTP0 <- which(deltaTP[1,]==0 & deltaTP[2,]==0)
  deltaTP <- deltaTP[,c(deltaTP0,setdiff(1:ncol(deltaTP),which(deltaTP[1,]==0 & deltaTP[2,]==0)))]
  print(paste("Run",ncol(deltaTP),"iterations for IRS"))
  toMem <- ls()
  
  deltaIDs <- 1:ncol(deltaTP)
} else if(climScen>0){
  deltaIDs <- 1
}

sampleID <- 1
if(outType=="testRun"){
  deltaID<-deltaIDs[1]
  easyInit=FALSE
  forceSaveInitSoil=F 
  cons10run = F
  initilizeSoil <- F
  procDrPeat=F
  coeffPeat1=-240 
  coeffPeat2=70
  coefCH4 = 0.34#g m-2 y-1
  coefN20_1 = 0.23
  coefN20_2 = 0.077#g m-2 y-1
  landClassUnman=NULL
  compHarvX = 0
  sampleX=NULL
  funPreb = regionPrebas
  initSoilCreStart=NULL
  outModReStart=NULL
  reStartYear=1
}
source_url("https://raw.githubusercontent.com/virpi-j/adaptFirst_runs/master/functions.R")
source_url("https://raw.githubusercontent.com/ForModLabUHel/IBCcarbon_runs/master/general/functions.r")


if(outType=="testRun"){
  # CurrClim scenario using the IBC-carbon settings to get soilC initialization
  sampleXs0 <- list()
  if(climScen>=0 | (CO2fixed==0 & harvscen=="Base" & harvinten=="Base")){
    outType<-"testRun"
    nYears<-2050-2015
    endingYear <- nYears + startingYear
    print(paste("Simulate soilC for",nYears,"years"))
    sampleXs0 <- runModelAdapt(1,
                               outType=outType,  
                               climScen=0,
                               rcps = "CurrClim",
                               CO2fixed=CO2fixed,
                               harvScen="Base",
                               harvInten="Base")
  }
  # IRS runs
  outType<-"dTabs"
  nYears <- 2100-1991#2015
  if(climScen > 0) nYears <- 2100-2015
  endingYear <- nYears + startingYear
  rcps <- rcpsFile 
  #source_url("https://raw.githubusercontent.com/virpi-j/adaptFirst_runs/master/settings.R")
  print(paste("Simulate for",nYears,"years"))
  source_url("https://raw.githubusercontent.com/virpi-j/adaptFirst_runs/master/functions.R")
  source("~/adaptFirst_runs/functions.R")
  source_url("https://raw.githubusercontent.com/ForModLabUHel/IBCcarbon_runs/master/general/functions.r")
  sampleXs <- lapply(deltaIDs, function(jx) { 
    runModelAdapt(jx,
                  outType=outType, climScen=climScen,
                  rcps = rcpsFile,
                  CO2fixed=CO2fixed,
                  harvScen=harvscen,#"Base" or #BaseTapio
           harvInten=harvinten)})

  if(climScen<10) sampleXs <- list(sampleXs0, sampleXs)
  
} else {
  # Baseline for soil and deadWood initialization
  sampleXs0 <- list()
  if(CO2fixed==0 & harvscen=="Base" & harvinten=="Base"){
    outType<-"testRun"
    nYears<-2050-2015
    if(climScen > 0) nYears <- 2100-2015
    endingYear <- nYears + startingYear
    print(paste("Simulate soilC for",nYears,"years"))
    sampleXs0 <- runModelAdapt(1,
                               outType=outType,  
                               rcps = "CurrClim",
                               CO2fixed=CO2fixed,
                               harvScen="Base",
                               harvInten="Base")
  }
  # IRS runs
  outType<-"dTabs"
  nYears <- 2100-1991#2015
  if(climScen > 0) nYears <- 2100-2015
  endingYear <- nYears + startingYear
  rcps <- rcpsFile 
  source_url("https://raw.githubusercontent.com/virpi-j/adaptFirst_runs/master/settings.R")
  print(paste("Simulate for",nYears,"years"))
  source_url("https://raw.githubusercontent.com/virpi-j/adaptFirst_runs/master/functions.R")
  source_url("https://raw.githubusercontent.com/ForModLabUHel/IBCcarbon_runs/master/general/functions.r")
  #  source("~/adaptFirst_runs/functions.R")
#  sampleXs <- runModelAdapt(1,
#                outType=outType,  
#                rcps = rcpsFile,
#                CO2fixed=CO2fixed,
#                climScen=climScen,
#                harvScen="Base",
#                harvInten="Base")
  sampleXs <- mclapply(deltaIDs, function(jx) {
    runModelAdapt(jx,
             outType=outType,  
             rcps = rcpsFile, #paste0(stat_name,"_1991_2100_constant_change_v1.csv"),
             CO2fixed=CO2fixed, climScen=climScen,
             #harvScen="baseTapio",#"Base" or baseTapio
             harvScen=harvscen,
             harvInten=harvinten)
    }, mc.cores = nCores,mc.silent=FALSE)      
  if(climScen<10){
    sampleXs <- list(sampleXs0, sampleXs)
  }
}

if(climScen<0){
  save(sampleXs,deltaTP,file = paste0("Results/outputs",station_id,".rdata"))
  print("Results saved as lists")
} else {
  output <- sampleXs[[2]][[1]]
  print(output[,1:4])
  save(output,file = paste0("Results/outputs",station_id,"_",rcps,
                              "_",harvscen,"_",harvinten,
                              "_Nrestrct",restrictionSwitch,".rdata"))
}

if(climScen<0){
  load(paste0("Results/outputs",station_id,".rdata"))
  output <- list()
  ndeltaTP <- ncol(deltaTP)
  m <- nrow(sampleXs[[2]][[1]])
  for(k in 1:m){
    tmp <- data.frame()
    index <- 0
    for(ij in 1:ndeltaTP){
      tmp <- rbind(tmp,cbind(t(deltaTP[,ij+index]),sampleXs[[2]][[ij]][k,2:ncol(sampleXs[[2]][[1]])]))
    }  
    names(tmp)[1:2] <- c("deltaT","deltaP")
    
    output[[k]] <- tmp 
    names(output)[k] <- sampleXs[[2]][[1]][k,1]
    
  }
  save(output,file = paste0("Results/outputs_",stat_name,"_",harvscen,"_",harvinten,"_",rcpsName,"_",co2Names[Co2Col],".rdata"))
  
  plotFigs <- TRUE
  if(plotFigs){
    pdf(file=paste0("Results/results_",stat_name,"_",harvscen,"_",harvinten,"_",rcpsName,"_",co2Names[Co2Col],".pdf"))
    for(k in 1:m){
      contourPlot <- TRUE
      if(contourPlot){
        zz1<-matrix(as.numeric(output[[k]][,3]),nrow=length(deltaT),ncol=length(deltaP),byrow = TRUE)
        zz7<-matrix(as.numeric(output[[k]][,10]),nrow=length(deltaT),ncol=length(deltaP),byrow = TRUE)
        par(mfrow=c(1,1))   
        nlev <- 10
        zrange <- range(cbind(zz1,zz7), finite = TRUE)
        filled.contour(deltaT, deltaP, zz1, 
                       zlim = zrange,
                       levels = pretty(zrange, nlev), nlevels = nlev,
                       col =  hcl.colors(20, "Spectral"),
                       xlab = "deltaT", ylab = "deltaP",  
                       main = paste0(names(output)[k],"/",rcpsName,"/",co2Names[Co2Col]," ",perStarts[1],"-",perEnds[1])
        )
        filled.contour(deltaT, deltaP, zz7,       
                       zlim = zrange,
                       levels = pretty(zrange, nlev), nlevels = nlev,
                       col =  hcl.colors(20, "Spectral"),
                       xlab = "deltaT", ylab = "deltaP",  
                       main = paste0(names(output)[k],"/",rcpsName,"/",co2Names[Co2Col]," ",
                                     perStarts[8],"-",perEnds[8])
        )
      } else {
        par(mfrow=c(1,2))   
        x<- as.numeric(output[[k]][,1])
        y<- as.numeric(output[[k]][,2])
        z1<- as.numeric(output[[k]][,3])
        z7<- as.numeric(output[[k]][,10])
        #scatter3D(x, y, z, colvar = z, col = NULL, add = FALSE)
        scatter3D(x, y, z1, pch = 18, cex = 2, 
                  theta = 20, phi = 20, ticktype = "detailed",
                  xlab = "deltaT", ylab = "deltaP", zlab = names(output)[k],  
                  #surf = list(x = x.pred, y = y.pred, z = z.pred,  
                  #            facets = NA, fit = fitpoints), 
                  main = paste0("per1:",perStarts[1],"-",perEnds[1])
        )
        scatter3D(x, y, z7, pch = 18, cex = 2, 
                  theta = 20, phi = 20, ticktype = "detailed",
                  xlab = "deltaT", ylab = "deltaP", zlab = names(output)[k],  
                  #surf = list(x = x.pred, y = y.pred, z = z.pred,  
                  #            facets = NA, fit = fitpoints), 
                  main = paste0("per7:",perStarts[8],"-",perEnds[8])
        )
      }
    }
    dev.off()
  }
}


# models outputs to NAs, outputDT, initSoilC and plots
#Sys.chmod(list.dirs("NAs"), "0777",use_umask=FALSE)
#f <- list.files("NAs", all.files = TRUE, full.names = TRUE, recursive = TRUE)
#Sys.chmod(f, (file.info(f)$mode | "0777"),use_umask=FALSE)

#Sys.chmod(list.dirs("outputDT"), "0777",use_umask=FALSE)
#f <- list.files("outputDT", all.files = TRUE, full.names = TRUE, recursive = TRUE)
#Sys.chmod(f, (file.info(f)$mode | "0777"),use_umask=FALSE)

#Sys.chmod(list.dirs("initSoilC"), "0777",use_umask=FALSE)
#f <- list.files("initSoilC", all.files = TRUE, full.names = TRUE, recursive = TRUE)
#Sys.chmod(f, (file.info(f)$mode | "0777"),use_umask=FALSE)

#Sys.chmod(list.dirs("plots"), "0777",use_umask=FALSE)
#f <- list.files("plots", all.files = TRUE, full.names = TRUE, recursive = TRUE)
#Sys.chmod(f, (file.info(f)$mode | "0777"),use_umask=FALSE)