rm(list=ls())
gc()

station_id <- 1 # This to localsettings! Station id 1 to 6
stations <- data.frame(name=c("Helsinki","Jokioinen","Jyväskylä","Kajaani",
                              "Sodankylä","Utsjoki"),
              location=c("Helsinki_Vantaa_lentoasema",
                              "Jokioinen_Ilmala",
                              "Jyvaskyla_lentoasema",
                              "Kajaani_lentoasema",
                              "Sodankyla_Tahtela",
                              "Utsjoki_Kevo"),
              ID = c(1:6), x = c(24.96, 23.5, 25.67, 27.67, 26.63, 27.01),
               y = c(60.33, 60.81, 62.4, 64.28, 67.37, 69.76))

r_nos_stations <- c(1,5,3,6,5,4) # check these! the region of the weather station

r_no = region = r_nos_stations[station_id] # region ID
xy <- stations[station_id,c("ID","x","y")]
stat_name <- stations[station_id,"name"]
devtools::source_url("https://raw.githubusercontent.com/virpi-j/adaptFirst_runs/master/settings.R")
source_url("https://raw.githubusercontent.com/virpi-j/adaptFirst_runs/master/functions.R")
source_url("https://raw.githubusercontent.com/virpi-j/adaptFirst_runs/master/05_create_CO2cols.R")

load(paste0("/scratch/project_2000994/PREBASruns/finRuns/input/maakunta/maakunta_",r_no,"_IDsTab.rdata"))
data.all <- cbind(data.all,data.IDs[match(data.all$segID, data.IDs$maakuntaID),4:5])

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
toMem <- ls()

outType="testRun"

sampleID <- 1
if(outType=="testRun"){
easyInit=FALSE
forceSaveInitSoil=F 
cons10run = F
procDrPeat=F
coeffPeat1=-240 
coeffPeat2=70
coefCH4 = 0.34#g m-2 y-1
coefN20_1 = 0.23
coefN20_2 = 0.077#g m-2 y-1
landClassUnman=NULL
compHarvX = 0
}
sampleXs <- lapply(sampleID, function(jx) { 
  runModel(jx, 
           outType=outType, 
           harvScen="Base",
           harvInten="Base")})

#mclapply(sampleIDs, function(jx) {
#    runModel(jx,outType="testRun",harvScen=harvScen,harvInten=harvInten)
#}, mc.cores = nCores,mc.silent=FALSE)      

# models outputs to NAs, outputDT, initSoilC and plots
Sys.chmod(list.dirs("NAs"), "0777",use_umask=FALSE)
f <- list.files("NAs", all.files = TRUE, full.names = TRUE, recursive = TRUE)
Sys.chmod(f, (file.info(f)$mode | "0777"),use_umask=FALSE)

Sys.chmod(list.dirs("outputDT"), "0777",use_umask=FALSE)
f <- list.files("outputDT", all.files = TRUE, full.names = TRUE, recursive = TRUE)
Sys.chmod(f, (file.info(f)$mode | "0777"),use_umask=FALSE)

Sys.chmod(list.dirs("initSoilC"), "0777",use_umask=FALSE)
f <- list.files("initSoilC", all.files = TRUE, full.names = TRUE, recursive = TRUE)
Sys.chmod(f, (file.info(f)$mode | "0777"),use_umask=FALSE)

Sys.chmod(list.dirs("plots"), "0777",use_umask=FALSE)
f <- list.files("plots", all.files = TRUE, full.names = TRUE, recursive = TRUE)
Sys.chmod(f, (file.info(f)$mode | "0777"),use_umask=FALSE)