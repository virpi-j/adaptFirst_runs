library(data.table)

CO2_obs <- data.frame(year = 1991:2020, CO2 = c(355.68, 356.42, 357.13, 358.61, 360.67, 362.58, 363.48, 366.27, 368.38,
  369.64, 371.15, 373.15, 375.64, 377.44, 379.46, 381.59, 383.37, 385.46,
  386.95, 389.21, 391.85, 394.06, 396.74, 398.87, 401.01, 404.41, 406.76, 
  408.72, 411.66, 414.24))
RCPyears <- c(2000,2005,seq(from=2010,to=2100,by=10))
years <- 2000:2100

CO2_RCP <- data.frame(year = RCPyears,
                      CO2_RCP26 = c(368.9, 378.8, 389.3, 412.1, 430.8, 
                                440.2, 442.7, 441.7, 437.5, 431.6,
                                426.0, 420.9),
                      CO2_RCP45 = c(368.9, 378.8, 389.1, 411.1, 435.0, 460.8,
                                486.5, 508.9, 524.3, 531.1, 533.7, 538.4),
                      CO2_RCP85 = c(368.9, 378.8, 389.3, 415.8, 448.8, 489.4,
                                540.5, 603.5, 677.1, 758.2, 844.8, 935.9))

CO2_RCPyears <- data.table(year=1991:1999, matrix(CO2_obs[which(CO2_obs$year<2000),2],2000-CO2_obs[1,1],3))
colnames(CO2_RCPyears) <- colnames(CO2_RCP)
CO2_RCPyears <- rbind(CO2_RCPyears, CO2_RCP[1,])
for(id in 2:nrow(CO2_RCP)){
  nyears <- which(years> RCPyears[id-1] & years < RCPyears[id])
  tmpCO2 <- data.frame(year=years[nyears])
  for(k in 1:3){
    tmp <- CO2_RCP[id-1,1+k]+(CO2_RCP[id,1+k] - CO2_RCP[id-1,1+k])/(RCPyears[id]-RCPyears[id-1])*(years[nyears]-RCPyears[id-1])
    tmpCO2 <- cbind(tmpCO2, data.table(tmp))
  }
  colnames(tmpCO2) <- colnames(CO2_RCP)
  CO2_RCPyears <- rbind(CO2_RCPyears, tmpCO2)
  CO2_RCPyears <- rbind(CO2_RCPyears, CO2_RCP[id,])
}


