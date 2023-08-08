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
              x = c(24.96, 23.50, 25.67, 27.67, 26.63, 27.01),
              y = c(60.33, 60.81, 62.4, 64.28, 67.37, 69.76))

r_nos_stations <- list()
r_nos_stations[[1]] <- c(1,9,11,13,15)
r_nos_stations[[2]] <- c(1,4,9,11,13,15)
r_nos_stations[[3]] <- c(6,4,13,17,7,3,12)
r_nos_stations[[4]] <- c(16,7,18,19)
r_nos_stations[[5]] <- c(8)
r_nos_stations[[6]] <- c(8)

r_no = region = r_nos_stations[[station_id]][1] # region ID

stat_name <- stations[station_id,"name"]
nYears <- 2050-2015
devtools::source_url("https://raw.githubusercontent.com/virpi-j/adaptFirst_runs/master/settings.R")

xy_UTM <- data.frame()
for(ij in 1:nrow(stations)){
  xy <- stations[ij,c("ID","x","y")]
  xy_UTM <- rbind(xy_UTM, LongLatToUTM(xy)[2:3])
}
colnames(xy_UTM)<-c("x_UTM","y_UTM")
stations <- cbind(stations,xy_UTM)
print(stations)

##source_url("https://raw.githubusercontent.com/virpi-j/adaptFirst_runs/master/05_create_CO2cols.R")
CO2_RCPyears <- read.csv2(file=paste0(climatepath,"co2_concentrations_PREBAS.csv"),header=T,sep = ";")
co2Names<-names(CO2_RCPyears)[2:ncol(CO2_RCPyears)]
names(CO2_RCPyears)[1]<-"year"

#load(paste0("/scratch/project_2000994/PREBASruns/finRuns/input/maakunta/maakunta_",r_no,"_IDsTab.rdata"))
#data.all <- cbind(data.all,data.IDs[match(data.all$segID, data.IDs$maakuntaID),4:5])
print(paste("PREBAS version",vPREBAS))

## Weather station coordinates:

xy <- stations[station_id,c("ID","x","y")]
station_coords <- LongLatToUTM(xy) # dd to UTM

d<- sqrt((data.all$x - station_coords$x)^2+(data.all$y - station_coords$y)^2)
nn.d <- order(d[d<100], decreasing=F)[1:nSitesRun]
print(range(d))

ops <- list(data.all[nn.d,])

print(paste("Weather station",stations[station_id,"name"]))
print(paste("Forest data: maximum distance from weather station",
            round(max(d[nn.d])/1000,2),"km"))
print(paste("Area of forest within the closest",nSitesRun,"segments is", round(sum(sum(ops[[1]]$area)),2),"hectares"))
if(toRaster){
  source_url("https://raw.githubusercontent.com/virpi-j/adaptFirst_runs/master/2_createRaster.R")
}
rcps = "CurrClim" 
#rcps <- paste0(stat_name,"_1991_2100_constant_change_v1.csv")
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

print(paste("Climate scenario",rcpsName))
print(paste("CO2scenario", names(CO2_RCPyears)[Co2Col+1]))

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
  if(CO2fixed==0 & harvscen=="Base" & harvinten=="Base"){
    outType<-"testRun"
    nYears<-2050-2015
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
  nYears <- 2100-2015
  endingYear <- nYears + startingYear
  rcps <- rcpsFile 
  source_url("https://raw.githubusercontent.com/virpi-j/adaptFirst_runs/master/settings.R")
  print(paste("Simulate for",nYears,"years"))
  source_url("https://raw.githubusercontent.com/virpi-j/adaptFirst_runs/master/functions.R")
  sampleXs <- lapply(deltaIDs, function(jx) { 
    runModelAdapt(jx,
                  outType=outType,  
                  rcps = rcpsFile,#"paste0(stat_name,"_1991_2100_constant_change_v1.csv"),
                  CO2fixed=CO2fixed,
                  harvScen=harvscen,#"Base" or #BaseTapio
           harvInten=harvinten)})
  sampleXs <- list(sampleXs0, sampleXs)
  
} else {
  # Baseline for soil and deadWood initialization
  sampleXs0 <- list()
  if(CO2fixed==0 & harvscen=="Base" & harvinten=="Base"){
    outType<-"testRun"
    nYears<-2050-2015
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
  nYears <- 2100-2015
  endingYear <- nYears + startingYear
  rcps <- rcpsFile 
  source_url("https://raw.githubusercontent.com/virpi-j/adaptFirst_runs/master/settings.R")
  print(paste("Simulate for",nYears,"years"))
  source_url("https://raw.githubusercontent.com/virpi-j/adaptFirst_runs/master/functions.R")
  sampleXs <- mclapply(deltaIDs, function(jx) {
    runModelAdapt(jx,
             outType=outType,  
             rcps = rcpsFile, #paste0(stat_name,"_1991_2100_constant_change_v1.csv"),
             CO2fixed=CO2fixed,
             #harvScen="baseTapio",#"Base" or baseTapio
             harvScen=harvscen,
             harvInten=harvinten)
    }, mc.cores = nCores,mc.silent=FALSE)      
  sampleXs <- list(sampleXs0, sampleXs)
}

save(sampleXs,deltaTP,file = paste0("outputs",station_id,".rdata"))
print("Results saved as lists")

load(paste0("outputs",station_id,".rdata"))
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
save(output,file = paste0("outputs_",stat_name,"_",harvscen,"_",harvinten,"_",rcpsName,"_",co2Names[Co2Col],".rdata"))

plotFigs <- TRUE
if(plotFigs){
  pdf(file=paste0("results_",stat_name,"_",harvscen,"_",harvinten,"_",rcpsName,"_",co2Names[Co2Col],".pdf"))
  for(k in 1:m){
    contourPlot <- TRUE
    if(contourPlot){
    zz1<-matrix(as.numeric(output[[k]][,3]),nrow=length(deltaT),ncol=length(deltaP),byrow = TRUE)
    zz7<-matrix(as.numeric(output[[k]][,9]),nrow=length(deltaT),ncol=length(deltaP),byrow = TRUE)
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
                   main = paste0(names(output)[k],"/",rcpsName,"/",co2Names[Co2Col]," ",perStarts[7],"-",perEnds[7])
    )
  } else {
    par(mfrow=c(1,2))   
    x<- as.numeric(output[[k]][,1])
    y<- as.numeric(output[[k]][,2])
    z1<- as.numeric(output[[k]][,3])
    z7<- as.numeric(output[[k]][,9])
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
            main = paste0("per7:",perStarts[7],"-",perEnds[7])
    )
  }
  }
  dev.off()
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