devtools::source_url("https://raw.githubusercontent.com/ForModLabUHel/IBCcarbon_runs/master/finRuns/Rsrc/settings.r")
source_url("https://raw.githubusercontent.com/ForModLabUHel/IBCcarbon_runs/master/general/functions.r")
rcpfile <- rcps
#pathFiles <- paste0("outputDT/forCent",r_no,"/")

data.IDs_rnos <- data.frame()
for(r_noi in r_nos_stations[[station_id]]){  
  load(paste0("input/maakunta/maakunta_",r_noi,"_IDsTab.rdata"))
 # data.IDs$segID <- data.IDs$maakuntaID
  data.IDs_rnos <- rbind(data.IDs_rnos, data.IDs)
}
data.IDs<-data.IDs_rnos
rm("data.IDs_rnos")
gc()

data.IDs <- data.IDs[segID!=0]
setkey(data.IDs,segID)

varXs <- c("age")
varX <- varXs[1]
#nSamples <- ceiling(dim(data.all)[1]/nSitesRun)


#pdf(paste0("plots/histRast_",r_no,
#           "_harscen",harvScen,
#           "_harInten",harvInten,"_",
#           rcpfile,".pdf"))

#if(!exists("varXs")) varXs <- c(varNames[varSel], specialVars)

for(varX in varXs){
  # varX <- varXs[1]
  #fileXs <- list.files(path = paste0(pathtoken,pathFiles), pattern = paste0(varX,"_harscen",harvScen,"_harInten",harvInten,"_",rcps))
  #if(length(fileXs) != nSamples) stop(paste0(nSamples-length(fileXs)," files missing"))
  
  #outX <- data.table()
  #for(i in 1:length(fileXs)){
  #  load(paste0(pathFiles,fileXs[i]))
  #  outX <- rbind(outX,get(varX))
  #}
  
  outX <- ops[[1]]
  setkey(outX,segID)
  setkey(data.IDs,segID)
  
  tabX <- merge(outX,data.IDs)
  rm(outX);gc()
  colnames(tabX)[colnames(tabX)=="x.x"]<-"x"
  colnames(tabX)[colnames(tabX)=="y.x"]<-"y"
  # can make a loop 
  rastX <- rasterFromXYZ(tabX[,.(x,y,age)])
  crs(rastX) <- crsX
  writeRaster(rastX,filename = paste0("../adaptFirst/rasters/station",station_id,"/",
                                      varX,".tiff"),overwrite=T)
  hist(rastX, main = paste(varX))
  
  rm(tabX);gc()
  rm(rastX);gc()
  
  # if(varX!="DeadWoodVolume")  file.remove(paste0(pathFiles,fileXs))
  print(varX)
}
dev.off()


###soilType raster

print("all rasters created")

# createRast outputs to rasters and plots
Sys.chmod(list.dirs("rasters"), "0777",use_umask=FALSE)
f <- list.files("rasters", all.files = TRUE, full.names = TRUE, recursive = TRUE)
Sys.chmod(f, (file.info(f)$mode | "0777"),use_umask=FALSE)

Sys.chmod(list.dirs("plots"), "0777",use_umask=FALSE)
f <- list.files("plots", all.files = TRUE, full.names = TRUE, recursive = TRUE)
Sys.chmod(f, (file.info(f)$mode | "0777"),use_umask=FALSE)
