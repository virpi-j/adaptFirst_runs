#devtools::source_url("https://raw.githubusercontent.com/ForModLabUHel/IBCcarbon_runs/master/finRuns/Rsrc/settings.r")
#source_url("https://raw.githubusercontent.com/ForModLabUHel/IBCcarbon_runs/master/general/functions.r")
#setwd("/scratch/project_2000994/PREBASruns/finRuns/")
#rcpfile <- "CurrClim"
#pathFiles <- paste0("outputDT/forCent",r_no,"/")

data.IDs_rnos <- data.frame()
for(r_noi in r_nos_stations[[station_id]]){  
  load(paste0("/scratch/project_2000994/PREBASruns/finRuns/input/maakunta/maakunta_",r_noi,"_IDsTab.rdata"))
  #data.IDs$segID <- data.IDs$maakuntaID
  data.IDs_rnos <- rbind(data.IDs_rnos, data.IDs)
}
data.IDs<-data.IDs_rnos
rm("data.IDs_rnos")
gc()

data.IDs <- data.IDs[segID!=0]
setkey(data.IDs,segID)

if(!exists("varXs")) varXs <- c("age")
varX <- varXs[1]

for(varX in varXs){

  x<-load(paste0("../outputDT/weatherStation",station_id,"/",varX,"_harscen",harvScen,"_harInten",harvInten,"_",rcps,"_Nswitch",restrictionSwitch,".rdata"))
  outX <- get(x)
  #outX <- ops[[1]]
  #outX$segID <- outX$maakuntaID
  setkey(outX,segID)
  setkey(data.IDs,segID)
  
  tabX <- merge(data.IDs,outX)
  colnames(tabX)[colnames(tabX)=="x.x"]<-"x"
  colnames(tabX)[colnames(tabX)=="y.x"]<-"y"
  #dev.off()
  ndat <- sample(1:nrow(tabX),1000)
  #plot(tabX$x[ndat],tabX$y[ndat],col="blue")
  
  #tabX <- merge(outX,data.IDs)
  #colnames(tabX)[colnames(tabX)=="x.y"]<-"x"
  #colnames(tabX)[colnames(tabX)=="y.y"]<-"y"
  rm(outX);gc()
  
  filee <- paste0("../rasters/",stations$name[station_id],"_",varX,"_harscen",harvScen,"_harInten",harvInten,"_",rcps,"_Nswitch",restrictionSwitch)

  rastX <- rasterFromXYZ(tabX[,.(x,y,per1)])
  crs(rastX) <- crsX
  plot(rastX)
  writeRaster(rastX,filename = paste0(filee,"_per1.tiff"),overwrite=T)

  rastX <- rasterFromXYZ(tabX[,.(x,y,per7)])
  crs(rastX) <- crsX
  plot(rastX)
  writeRaster(rastX,filename = paste0(filee"_per7.tiff"),overwrite=T)
  
  #hist(rastX, main = paste(varX))
  
  rm(tabX);
  rm(rastX);gc()
  
  # if(varX!="DeadWoodVolume")  file.remove(paste0(pathFiles,fileXs))
  print(varX)
}
#dev.off()


###soilType raster

print("all rasters created")

# createRast outputs to rasters and plots
Sys.chmod(list.dirs("rasters"), "0777",use_umask=FALSE)
f <- list.files("rasters", all.files = TRUE, full.names = TRUE, recursive = TRUE)
Sys.chmod(f, (file.info(f)$mode | "0777"),use_umask=FALSE)

Sys.chmod(list.dirs("plots"), "0777",use_umask=FALSE)
f <- list.files("plots", all.files = TRUE, full.names = TRUE, recursive = TRUE)
Sys.chmod(f, (file.info(f)$mode | "0777"),use_umask=FALSE)
