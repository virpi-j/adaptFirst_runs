#station_id <- 1 # This to localsettings! Station id 1 to 6
set.seed(1)
print(paste("station",station_id))
stations <- data.frame(name=c("Helsinki","Jokioinen","Jyvaskyla","Kajaani",
                              "Sodankyla","Utsjoki"),
              location=c("Helsinki_Vantaa_lentoasema", # Uusimaa 1
                              "Jokioinen_Ilmala", # Kanta-Hame 9
                              "Jyvaskyla_lentoasema", # Keski-Suomi 6 
                              "Kajaani_lentoasema", # Kainuu 16
                              "Sodankyla_Tahtela", # Lappi 8
                              "Utsjoki_Kevo"), # Lappi 8
              ID = c(1:6), 
              x = c(24.96, 23.5, 25.67, 27.67, 26.63, 27.01),
              y = c(60.33, 60.81, 62.4, 64.28, 67.37, 69.76))

r_nos_stations <- list()
r_nos_stations[[1]] <- c(1,9,11,13,15)
r_nos_stations[[2]] <- c(9,4,11,13)
r_nos_stations[[3]] <- c(6,4,13,17,7,3,12)
r_nos_stations[[4]] <- c(16,7,18,19)
r_nos_stations[[5]] <- c(8)
r_nos_stations[[6]] <- c(8)

r_no = region = r_nos_stations[[station_id]][1] # region ID
xy <- stations[station_id,c("ID","x","y")]
stat_name <- stations[station_id,"name"]
devtools::source_url("https://raw.githubusercontent.com/virpi-j/adaptFirst_runs/master/settings.R")
source_url("https://raw.githubusercontent.com/virpi-j/adaptFirst_runs/master/05_create_CO2cols.R")

#load(paste0("/scratch/project_2000994/PREBASruns/finRuns/input/maakunta/maakunta_",r_no,"_IDsTab.rdata"))
#data.all <- cbind(data.all,data.IDs[match(data.all$segID, data.IDs$maakuntaID),4:5])
print(paste("PREBAS version",vPREBAS))

## Weather station coordinates:

station_coords <- LongLatToUTM(xy) # dd to UTM

d<- sqrt((data.all$x - station_coords$x)^2+(data.all$y - station_coords$y)^2)
nn.d <- order(d, decreasing=F)[1:nSitesRun]

ops <- list(data.all[nn.d,])

print(paste("Weather station",stations[station_id,"name"]))
print(paste("Forest data: maximum distance from weather station",
            round(max(d[nn.d])/1000,2),"km"))
print(paste("Area of forest within the closest",nSitesRun,"segments is", round(sum(sum(ops[[1]]$area)),2),"hectares"))

rcps = "CurrClim" 
rcps <- paste0(stat_name,"_1991_2100_constant_change_v1.csv")

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
}
#source_url("https://raw.githubusercontent.com/virpi-j/adaptFirst_runs/master/functions.R")
if(outType=="testRun"){
  outType<-"dTabs"
  rcps = "CurrClim" 
  source_url("https://raw.githubusercontent.com/virpi-j/adaptFirst_runs/master/functions.R")
  sampleXs0 <- runModel(deltaIDs[1],
                        outType=outType, 
                        harvScen="Base",
                        harvInten="Base")
  
  rcps <- paste0(stat_name,"_1991_2100_constant_change_v1.csv")
  source_url("https://raw.githubusercontent.com/virpi-j/adaptFirst_runs/master/functions.R")
  sampleXs <- lapply(deltaIDs, function(jx) { 
    runModel(jx, 
           outType=outType, 
           CO2fixed=CO2fixed,
           harvScen="Base",
           harvInten="Base")})
  sampleXs <- list(sampleXs0, sampleXs)
  
} else {
  # Baseline for soil and deadWood initialization
  rcps = "CurrClim" 
  source_url("https://raw.githubusercontent.com/virpi-j/adaptFirst_runs/master/functions.R")
  sampleXs0 <- runModel(deltaIDs[1],
                        outType=outType, 
             harvScen="Base",
             harvInten="Base")
  # deltaT & deltaP climate runs  
  rcps <- paste0(stat_name,"_1991_2100_constant_change_v1.csv")
  source_url("https://raw.githubusercontent.com/virpi-j/adaptFirst_runs/master/functions.R")
  sampleXs <- mclapply(deltaIDs, function(jx) {
    runModel(jx,
             outType=outType,
             CO2fixed=CO2fixed,
             harvScen="Base",
             harvInten="Base")
    }, mc.cores = nCores,mc.silent=FALSE)      
  sampleXs <- list(sampleXs0, sampleXs)
}

save(sampleXs,deltaTP,file = "outputs.rdata")


# load("outputs.rdata")
output <- list()
ndeltaTP <- ncol(deltaTP)
m <- nrow(sampleXs[[1]])
for(k in 1:m){
  tmp <- data.frame()
  index <- 0
  for(ij in 1:ndeltaTP){
    if(ij==deltaTP0) {
      tmp <- rbind(tmp,cbind(t(deltaTP[,1]),sampleXs[[1]][k,2:ncol(sampleXs[[1]])]))
      index <- -1
    } else {
      tmp <- rbind(tmp,cbind(t(deltaTP[,ij+1+index]),sampleXs[[2]][[ij+index]][k,2:ncol(sampleXs[[1]])]))
    }
  }  
  names(tmp)[1:2] <- c("deltaT","deltaP")

  output[[k]] <- tmp 
  names(output)[k] <- sampleXs[[1]][k,1]
  
}
save(output,file = paste0("outputs_",stat_name,".rdata"))

plotFigs <- TRUE
if(plotFigs){
  pdf(file=paste0("results_",stat_name,".pdf"))
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
      main = paste0(names(output)[k]," ",perStarts[1],"-",perEnds[1])
    )
    filled.contour(deltaT, deltaP, zz7,       
                   zlim = zrange,
                   levels = pretty(zrange, nlev), nlevels = nlev,
                   col =  hcl.colors(20, "Spectral"),
                   xlab = "deltaT", ylab = "deltaP",  
            main = paste0(names(output)[k]," ",perStarts[7],"-",perEnds[7])
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