rm(list=ls())

library(terra)
library(raster)
library(dplyr)
library(sf)
library(spData)

stERS <- us_states[,"NAME"] %>%
  summarize(geometry = st_union(geometry)) #
# merge states (can remove if want states separate)

stGrid <- st_make_grid(stERS, 
                       n=c(500,500), 
                       what="polygons", 
                       square = F, 
                       flat_topped = T) %>%
  st_as_sf()

stGrid <- st_intersection(stGrid, stERS) %>%
  mutate(Hex = rownames(.))

# Land ----

rat <- raster("/Users/gabriela.quinlan/Desktop/camelinaLoc/createdData/landUseRasters/camMosaic19_23.tif", RAT=T)
#camCoord <- xyFromCell(rat, which.max(rat))
#write.csv(camCoord, "camelinaPts.csv") # 303135 locations /14846299153 cells 

rat2 <- raster("/Users/gabriela.quinlan/Desktop/camelinaLoc/createdData/landUseRasters/canolaMosaic19_23.tif", RAT=T)
#canCoord <- xyFromCell(rat2, which.max(rat2)) 
#write.csv(canCoord, "canolaPts.csv") # 35607992 locations /14846299153 cells 


camCoord <- read.csv("/Users/gabriela.quinlan/Desktop/camelinaLoc/allPointValues_30m/camelinaPts.csv", row.names = 1)
canCoord <- read.csv("/Users/gabriela.quinlan/Desktop/camelinaLoc/allPointValues_30m/canolaPts.csv", row.names = 1)


# Soil ----
di <- raster("/Users/gabriela.quinlan/Desktop/camelinaLoc/rawData/L48_2020_DI_240m.tif")
pi <- raster("/Users/gabriela.quinlan/Desktop/camelinaLoc/rawData/L48_2020_PI_240m.tif")

pi[pi == 255] <- NA # take out 255 (NA)
di[di == 255] <- NA

piCam <- cbind(terra::extract(pi, camCoord), camCoord)
diCam <- cbind(terra::extract(di, camCoord), camCoord)

piCan <- cbind(terra::extract(pi, canCoord), canCoord)
diCan <- cbind(terra::extract(di, canCoord), canCoord)

canolaSoil <- full_join(piCan, diCan) %>%
  rename(pi = `terra::extract(pi, canCoord)`, 
         di = `terra::extract(di, canCoord)`)

camelinaSoil <- full_join(piCam, diCam) %>%
  rename(pi = `terra::extract(pi, camCoord)`, 
         di = `terra::extract(di, camCoord)`)

#write.csv(canolaSoil, "canolaSoil.csv")
#write.csv(camelinaSoil, "camelinaSoil.csv")

## Models & metadata ----
modelList <- list()
sdEnvData <- vector()
medianEnvData <- vector() 
meanEnvData <- vector()
lengthEnvData <- vector()

modelList[[1]] <- ecdf(canolaSoil[,"di"])
modelList[[2]] <- ecdf(canolaSoil[,"pi"])
modelList[[3]] <- ecdf(camelinaSoil[,"di"])
modelList[[4]] <- ecdf(camelinaSoil[,"pi"])

sdEnvData[1] <- sd(canolaSoil[,"di"], na.rm=T)
sdEnvData[2] <- sd(canolaSoil[,"pi"], na.rm=T)
sdEnvData[3] <- sd(camelinaSoil[,"di"], na.rm=T)
sdEnvData[4] <- sd(camelinaSoil[,"pi"], na.rm=T)

medianEnvData[1] <- median(canolaSoil[,"di"], na.rm=T)
medianEnvData[2] <- median(canolaSoil[,"pi"], na.rm=T)
medianEnvData[3] <- median(camelinaSoil[,"di"], na.rm=T)
medianEnvData[4] <- median(camelinaSoil[,"pi"], na.rm=T)

meanEnvData[1] <- mean(canolaSoil[,"di"], na.rm=T)
meanEnvData[2] <- mean(canolaSoil[,"pi"], na.rm=T)
meanEnvData[3] <- mean(camelinaSoil[,"di"], na.rm=T)
meanEnvData[4] <- mean(camelinaSoil[,"pi"], na.rm=T)

lengthEnvData[1] <- length(canolaSoil[,"di"])
lengthEnvData[2] <- length(canolaSoil[,"pi"])
lengthEnvData[3] <- length(camelinaSoil[,"di"])
lengthEnvData[4] <- length(camelinaSoil[,"pi"])

metaDatSoil <- as.data.frame(cbind(medianEnvData, sdEnvData, meanEnvData, lengthEnvData,
                                   c("canola", "canola", "camelina", "camelina"), 
                                   c("di", "pi", "di", "pi")))

colnames(metaDatSoil) <- c("medianEnvData", "sdEnvData", "meanEnvData", "lengthEnvData", "crop", "envVar")

#saveRDS(modelList, "soilMods.rds")
#write.csv(metaDatSoil, "soilParams.csv")

# Need to format like everything else so can combine later. 
# did this by hand and went back and wrote code after. 
canolaModelParamsSoil <- metaDatSoil[which(metaDatSoil$crop == "canola"),] %>%
  dplyr::select(!c(crop, envVar)) # di then pi

camelinaModelParamsSoil <- metaDatSoil[which(metaDatSoil$crop == "camelina"),] %>%
  dplyr::select(!c(crop, envVar)) # di then pi
  
#write.csv(canolaModelParamsSoil, "soilCanolaModelParams.csv")
#write.csv(camelinaModelParamsSoil, "soilCamelinaModelParams.csv")

# Weather ----
years <- c(2019:2023, 2069:2073)
bioVar <- c("01", "03", "08", "09", "10", "11", "12", "16", "17", "18", "19")
scenerio <- c("lev", "wre")

out <- do.call(expand.grid, list(years, bioVar, scenerio)) %>%
  rename("year" = 1, "bioVar" = 2, "scenerio" = 3,)

weatherCam <- list()
weatherCan <- list()

projNew <- crs(rat2)
i=1
for(j in 1:4){ #length(bioVar)
  for(k in 1:length(years)){
    
    ratPath <- paste("/Users/gabriela.quinlan/Desktop/camelinaLoc/rawData/", scenerio[i], "/bio",
                     bioVar[j], "_equiv_", years[k], ".tiff", sep="")
    ratTemp <- raster(ratPath, RAT=T)
    ratTemp2 <- projectRaster(ratTemp, crs=projNew)
    
    index <- as.numeric(rownames(out[which(out$year == years[k] &
                                             out$bioVar == bioVar[j] &
                                             out$scenerio == scenerio[i]),]))
    print(index)
    weatherCam[[index]] <- terra::extract(ratTemp2, camCoord)
    
  }
}

#write.csv(do.call(cbind, weatherCam), "weatherCam.csv")


weatherCam <- list()

for(j in 5:9){ #length(bioVar)
  for(k in 1:length(years)){
    
    ratPath <- paste("/Users/gabriela.quinlan/Desktop/camelinaLoc/rawData/", scenerio[i], "/bio",
                     bioVar[j], "_equiv_", years[k], ".tiff", sep="")
    ratTemp <- raster(ratPath, RAT=T)
    ratTemp2 <- projectRaster(ratTemp, crs=projNew)
    
    index <- as.numeric(rownames(out[which(out$year == years[k] &
                                             out$bioVar == bioVar[j] &
                                             out$scenerio == scenerio[i]),]))
    print(index)
    weatherCam[[index]] <- terra::extract(ratTemp2, camCoord)
    
    
  }
}

#write.csv(do.call(cbind, weatherCam), "weatherCam1.csv")


weatherCam <- list()

for(j in 10:length(bioVar)){ #length(bioVar)
  for(k in 1:length(years)){
    
    ratPath <- paste("/Users/gabriela.quinlan/Desktop/camelinaLoc/rawData/", scenerio[i], "/bio",
                     bioVar[j], "_equiv_", years[k], ".tiff", sep="")
    ratTemp <- raster(ratPath, RAT=T)
    ratTemp2 <- projectRaster(ratTemp, crs=projNew)
    
    index <- as.numeric(rownames(out[which(out$year == years[k] &
                                             out$bioVar == bioVar[j] &
                                             out$scenerio == scenerio[i]),]))
    print(index)
    weatherCam[[index]] <- terra::extract(ratTemp2, camCoord)
    
    
  }
}

#write.csv(do.call(cbind, weatherCam), "weatherCam2.csv")

weatherCam <- list()

i=2
for(j in 1:4){ #length(bioVar)
  for(k in 1:length(years)){
    
    ratPath <- paste("/Users/gabriela.quinlan/Desktop/camelinaLoc/rawData/", scenerio[i], "/bio",
                     bioVar[j], "_equiv_", years[k], ".tiff", sep="")
    ratTemp <- raster(ratPath, RAT=T)
    ratTemp2 <- projectRaster(ratTemp, crs=projNew)
    
    index <- as.numeric(rownames(out[which(out$year == years[k] &
                                             out$bioVar == bioVar[j] &
                                             out$scenerio == scenerio[i]),]))
    print(index)
    weatherCam[[index]] <- terra::extract(ratTemp2, camCoord)
    
    
  }
}

#write.csv(do.call(cbind, weatherCam), "weatherCam3.csv")

weatherCam <- list()

for(j in 5:9){ #length(bioVar)
  for(k in 1:length(years)){
    
    ratPath <- paste("/Users/gabriela.quinlan/Desktop/camelinaLoc/rawData/", scenerio[i], "/bio",
                     bioVar[j], "_equiv_", years[k], ".tiff", sep="")
    ratTemp <- raster(ratPath, RAT=T)
    ratTemp2 <- projectRaster(ratTemp, crs=projNew)
    
    index <- as.numeric(rownames(out[which(out$year == years[k] &
                                             out$bioVar == bioVar[j] &
                                             out$scenerio == scenerio[i]),]))
    print(index)
    weatherCam[[index]] <- terra::extract(ratTemp2, camCoord)
    
    
  }
}

#write.csv(do.call(cbind, weatherCam), "weatherCam4.csv")


weatherCam <- list()

for(j in 10:length(bioVar)){ #length(bioVar)
  for(k in 1:length(years)){
    
    ratPath <- paste("/Users/gabriela.quinlan/Desktop/camelinaLoc/rawData/", scenerio[i], "/bio",
                     bioVar[j], "_equiv_", years[k], ".tiff", sep="")
    ratTemp <- raster(ratPath, RAT=T)
    ratTemp2 <- projectRaster(ratTemp, crs=projNew)
    
    index <- as.numeric(rownames(out[which(out$year == years[k] &
                                             out$bioVar == bioVar[j] &
                                             out$scenerio == scenerio[i]),]))
    print(index)
    weatherCam[[index]] <- terra::extract(ratTemp2, camCoord)
    
    
  }
}

#write.csv(do.call(cbind, weatherCam), "weatherCam5.csv")

i=1
for(j in 1:3){ #length(bioVar)
  for(k in 1:length(years)){
    
    ratPath <- paste("/Users/gabriela.quinlan/Desktop/camelinaLoc/rawData/", scenerio[i], "/bio",
                     bioVar[j], "_equiv_", years[k], ".tiff", sep="")
    ratTemp <- raster(ratPath, RAT=T)
    ratTemp2 <- projectRaster(ratTemp, crs=projNew)
    
    index <- as.numeric(rownames(out[which(out$year == years[k] &
                                             out$bioVar == bioVar[j] &
                                             out$scenerio == scenerio[i]),]))
    print(index)
    weatherCan[[index]] <- terra::extract(ratTemp2, canCoord)
    
  }
}

#write.csv(do.call(cbind, weatherCan), "weatherCan.csv")


weatherCan <- list()

for(j in 4:6){ #length(bioVar)
  for(k in 1:length(years)){
    
    ratPath <- paste("/Users/gabriela.quinlan/Desktop/camelinaLoc/rawData/", scenerio[i], "/bio",
                     bioVar[j], "_equiv_", years[k], ".tiff", sep="")
    ratTemp <- raster(ratPath, RAT=T)
    ratTemp2 <- projectRaster(ratTemp, crs=projNew)
    
    index <- as.numeric(rownames(out[which(out$year == years[k] &
                                             out$bioVar == bioVar[j] &
                                             out$scenerio == scenerio[i]),]))
    print(index)
    weatherCan[[index]] <- terra::extract(ratTemp2, canCoord)
    
    
  }
}

#write.csv(do.call(cbind, weatherCan), "weatherCan1.csv")


weatherCan <- list()

for(j in 7:9){ #length(bioVar)
  for(k in 1:length(years)){
    
    ratPath <- paste("/Users/gabriela.quinlan/Desktop/camelinaLoc/rawData/", scenerio[i], "/bio",
                     bioVar[j], "_equiv_", years[k], ".tiff", sep="")
    ratTemp <- raster(ratPath, RAT=T)
    ratTemp2 <- projectRaster(ratTemp, crs=projNew)
    
    index <- as.numeric(rownames(out[which(out$year == years[k] &
                                             out$bioVar == bioVar[j] &
                                             out$scenerio == scenerio[i]),]))
    print(index)
    weatherCan[[index]] <- terra::extract(ratTemp2, canCoord)
    
    
  }
}

#write.csv(do.call(cbind, weatherCan), "weatherCan2.csv")

weatherCan <- list()

for(j in 10:length(bioVar)){ #length(bioVar)
  for(k in 1:length(years)){
    
    ratPath <- paste("/Users/gabriela.quinlan/Desktop/camelinaLoc/rawData/", scenerio[i], "/bio",
                     bioVar[j], "_equiv_", years[k], ".tiff", sep="")
    ratTemp <- raster(ratPath, RAT=T)
    ratTemp2 <- projectRaster(ratTemp, crs=projNew)
    
    index <- as.numeric(rownames(out[which(out$year == years[k] &
                                             out$bioVar == bioVar[j] &
                                             out$scenerio == scenerio[i]),]))
    print(index)
    weatherCan[[index]] <- terra::extract(ratTemp2, canCoord)
    
    
  }
}

#write.csv(do.call(cbind, weatherCan), "weatherCan3.csv")

weatherCan <- list()

i=2
for(j in 1:3){ #length(bioVar)
  for(k in 1:length(years)){
    
    ratPath <- paste("/Users/gabriela.quinlan/Desktop/camelinaLoc/rawData/", scenerio[i], "/bio",
                     bioVar[j], "_equiv_", years[k], ".tiff", sep="")
    ratTemp <- raster(ratPath, RAT=T)
    ratTemp2 <- projectRaster(ratTemp, crs=projNew)
    
    index <- as.numeric(rownames(out[which(out$year == years[k] &
                                             out$bioVar == bioVar[j] &
                                             out$scenerio == scenerio[i]),]))
    print(index)
    weatherCan[[index]] <- terra::extract(ratTemp2, canCoord)
    
    
  }
}

#write.csv(do.call(cbind, weatherCan), "weatherCan4.csv")

weatherCan <- list()

for(j in 4:6){ #length(bioVar)
  for(k in 1:length(years)){
    
    ratPath <- paste("/Users/gabriela.quinlan/Desktop/camelinaLoc/rawData/", scenerio[i], "/bio",
                     bioVar[j], "_equiv_", years[k], ".tiff", sep="")
    ratTemp <- raster(ratPath, RAT=T)
    ratTemp2 <- projectRaster(ratTemp, crs=projNew)
    
    index <- as.numeric(rownames(out[which(out$year == years[k] &
                                             out$bioVar == bioVar[j] &
                                             out$scenerio == scenerio[i]),]))
    print(index)
    weatherCan[[index]] <- terra::extract(ratTemp2, canCoord)
    
    
  }
}

#write.csv(do.call(cbind, weatherCan), "weatherCan5.csv")


weatherCan <- list()

for(j in 7:9){ #length(bioVar)
  for(k in 1:length(years)){
    
    ratPath <- paste("/Users/gabriela.quinlan/Desktop/camelinaLoc/rawData/", scenerio[i], "/bio",
                     bioVar[j], "_equiv_", years[k], ".tiff", sep="")
    ratTemp <- raster(ratPath, RAT=T)
    ratTemp2 <- projectRaster(ratTemp, crs=projNew)
    
    index <- as.numeric(rownames(out[which(out$year == years[k] &
                                             out$bioVar == bioVar[j] &
                                             out$scenerio == scenerio[i]),]))
    print(index)
    weatherCan[[index]] <- terra::extract(ratTemp2, canCoord)
    
    
  }
}

#write.csv(do.call(cbind, weatherCan), "weatherCan6.csv")

weatherCan <- list()


for(j in 10:length(bioVar)){ #length(bioVar)
  for(k in 1:length(years)){
    
    ratPath <- paste("/Users/gabriela.quinlan/Desktop/camelinaLoc/rawData/", scenerio[i], "/bio",
                     bioVar[j], "_equiv_", years[k], ".tiff", sep="")
    ratTemp <- raster(ratPath, RAT=T)
    ratTemp2 <- projectRaster(ratTemp, crs=projNew)
    
    index <- as.numeric(rownames(out[which(out$year == years[k] &
                                             out$bioVar == bioVar[j] &
                                             out$scenerio == scenerio[i]),]))
    print(index)
    weatherCan[[index]] <- terra::extract(ratTemp2, canCoord)
    
    
  }
}

#write.csv(do.call(cbind, weatherCan), "weatherCan7.csv")

## Data to models ----
rm(list=ls())
# Make our functions and mean/sd 
filesAll <- list.files("/Users/gabriela.quinlan/Desktop/weatherConditions/", pattern = "Can")

modelList <- list()
sdEnvData <- vector()
medianEnvData <- vector() 
meanEnvData <- vector()
lengthEnvData <- vector()
start <- 0

# Take average of 5 years 
#for (j in 1: length(filesAll)) { #length(filesAll)
j=1
print(paste("on", j, "of", length(filesAll)))

dat <- read.csv(filesAll[[j]], row.names = 1) 

print(ncol(dat))

if(ncol(dat) == 10) { dat <- dat %>%
  mutate(V1 = rowMeans(select(dat, V1:V5), na.rm = TRUE), 
         V2 = rowMeans(select(dat, V6:V10), na.rm = TRUE)) %>%
  dplyr::select(V1:V2)}
if(ncol(dat) == 20) {dat <- dat %>%
  mutate(V1 = rowMeans(select(dat, V1:V5), na.rm = TRUE), 
         V2 = rowMeans(select(dat, V6:V10), na.rm = TRUE), 
         V3 = rowMeans(select(dat, V11:V15), na.rm = TRUE),
         V4 = rowMeans(select(dat, V16:V20), na.rm = TRUE)) %>%
  dplyr::select(V1:V4)}
if(ncol(dat) == 30) {dat <- dat %>%
  mutate(V1 = rowMeans(select(dat, V1:V5), na.rm = TRUE), 
         V2 = rowMeans(select(dat, V6:V10), na.rm = TRUE), 
         V3 = rowMeans(select(dat, V11:V15), na.rm = TRUE),
         V4 = rowMeans(select(dat, V16:V20), na.rm = TRUE),
         V5 = rowMeans(select(dat, V21:V25), na.rm = TRUE),
         V6 = rowMeans(select(dat, V26:V30), na.rm = TRUE)) %>%
  dplyr::select(V1:V6)}
if(ncol(dat) == 40) {dat <- dat %>%
  mutate(V1 = rowMeans(select(dat, V1:V5), na.rm = TRUE), 
         V2 = rowMeans(select(dat, V6:V10), na.rm = TRUE), 
         V3 = rowMeans(select(dat, V11:V15), na.rm = TRUE),
         V4 = rowMeans(select(dat, V16:V20), na.rm = TRUE),
         V5 = rowMeans(select(dat, V21:V25), na.rm = TRUE),
         V6 = rowMeans(select(dat, V26:V30), na.rm = TRUE),
         V7 = rowMeans(select(dat, V31:V35), na.rm = TRUE),
         V8 = rowMeans(select(dat, V36:V40), na.rm = TRUE)) %>%
  dplyr::select(V1:V8)}
if(ncol(dat) == 50) {dat <- dat %>%
  mutate(V1 = rowMeans(select(dat, V1:V5), na.rm = TRUE), 
         V2 = rowMeans(select(dat, V6:V10), na.rm = TRUE), 
         V3 = rowMeans(select(dat, V11:V15), na.rm = TRUE),
         V4 = rowMeans(select(dat, V16:V20), na.rm = TRUE),
         V5 = rowMeans(select(dat, V21:V25), na.rm = TRUE),
         V6 = rowMeans(select(dat, V26:V30), na.rm = TRUE),
         V7 = rowMeans(select(dat, V31:V35), na.rm = TRUE),
         V8 = rowMeans(select(dat, V36:V40), na.rm = TRUE), 
         V9 = rowMeans(select(dat, V41:V45), na.rm = TRUE), 
         V10 = rowMeans(select(dat, V46:V50), na.rm = TRUE)) %>%
  dplyr::select(V1:V10)}


for (i in 1: ncol(dat)) {
  start <- start + 1
  modelList[[start]] <- ecdf(dat[,i])
  sdEnvData[start] <- sd(dat[,i])
  medianEnvData[start] <- median(dat[,i])
  meanEnvData[start] <- mean(dat[,i])
  lengthEnvData[start] <- length(dat[,i])
} # i 
#} # j 


#saveRDS(modelList, paste( "canolaModel1.rds", sep = ""))
#write.csv(cbind(sdEnvData, medianEnvData, meanEnvData, lengthEnvData), "canolaModelParams1.csv")


## Create Metadata to match ----

out <- matrix(nrow = 44, ncol = 3)
colnames(out) <- c("scenerio", "envVar", "years")
start <- 0
years <- c("2019_2023", "2069_2073")

# Take same loop as above so get the index's 
for (i in 1:length(scenerio)) {
  for(j in 1:length(bioVar)) { 
    for(k in 1:length(years)) { 
      start <- start + 1
      
      out[start, "scenerio"] <- scenerio[i]
      out[start, "envVar"] <- bioVar[j]
      out[start, "years"] <- years[k]
    }
  }
}

out2 <- as.data.frame(out) %>%
  mutate(crop = "canola")

# index for chunks I saved them in 
metaData <- as.data.frame(out) %>%
  mutate(crop = "camelina") %>%
  rbind(out2) %>%
  mutate(subset = ifelse(scenerio == "lev" & envVar %in% bioVar[1:3], 1, NA)) %>%
  mutate(subset = ifelse(scenerio == "lev" & envVar %in% bioVar[4:6], 2, subset)) %>%
  mutate(subset = ifelse(scenerio == "lev" & envVar %in% bioVar[7:9], 3, subset)) %>%
  mutate(subset = ifelse(scenerio == "lev" & envVar %in% bioVar[10:11], 4, subset)) %>%
  mutate(subset = ifelse(scenerio == "wre" & envVar %in% bioVar[1:3], 5, subset)) %>%
  mutate(subset = ifelse(scenerio == "wre" & envVar %in% bioVar[4:6], 6, subset)) %>%
  mutate(subset = ifelse(scenerio == "wre" & envVar %in% bioVar[7:9], 7, subset)) %>%
  mutate(subset = ifelse(scenerio == "wre" & envVar %in% bioVar[10:11], 8, subset)) %>%
  mutate(subset = ifelse(crop == "camelina", "", subset)) 

#write.csv(metaData, 'metaData.csv')

soilMeta <- full_join(as.data.frame(cbind(scenerio = c(NA, NA), envVar = c("di", "pi"), years = c(NA, NA), crop = c("camelina", "camelina"))), 
          as.data.frame(cbind(scenerio = c(NA, NA), envVar = c("di", "pi"), years = c(NA, NA), crop = c("canola", "canola"))))

metaDatFull <- full_join(metaData, soilMeta)          

# Now combine summary files and add metadata 
filesAll <- list.files(path="/Users/gabriela.quinlan/Desktop/camelinaLoc/createdData/modelSummaries/", full.names = TRUE, pattern = ".csv")
result <- rbindlist(sapply(filesAll, fread,simplify = FALSE), idcol = 'filename', use.names=TRUE) %>%
  separate(filename, into= c("1", "2", "3", "4", "5", "6", "7", "8", "file", "10", "11")) %>% 
  dplyr::select(!c(1:8, 10:11)) %>%
  cbind(metaDatFull)

result <- rbind(result[1:88,], result[91:92,], result[89:90,]) # for some reason flipping soil 

#write.csv(result, "/Users/gabriela.quinlan/Desktop/camelinaLoc/createdData/modelSummaries/modelParamsAll.csv")


# Density at new hex size ----
library(sf)
camCoord2 <- st_as_sf(camCoord, coords = c("x","y"), crs = crs(pi))
camCoord2 <- st_transform(camCoord2, crs = crs(stGrid))

canCoord2 <- st_as_sf(canCoord, coords = c("x","y"), crs = crs(pi))
canCoord2 <- st_transform(canCoord2, crs = crs(stGrid))

data_sf_summary <- stGrid %>% 
  mutate(countsCam = lengths(st_intersects(., camCoord2)),
         countsCan = lengths(st_intersects(., canCoord2)))

df <- nc %>% st_set_geometry(NULL)

canCamArea$area <- as.numeric(st_area(stGrid))

canCamArea <- canCamArea %>%
  mutate(propCam = (countsCam*900) / area ,
         propCan = (countsCan*900) / area)

write.csv(canCamArea, "~/Desktop/camelinaLoc/createdData/landUse/canCamArea.csv")

