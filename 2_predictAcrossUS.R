# Make predictions ----

rm(list=ls())

library(tidyr); library(tibble); library(dplyr); library(ggplot2); 
library(ranger); library(sf); library(mgcv); library(spData);
library(patchwork); library(purrr)


# Make Grid to predict on ----
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


climAll <- read.csv("createdData/predictAcrossUS/climProjectAll.csv")
soilAll <- dplyr::select(read.csv("createdData/predictAcrossUS/soilAll.csv"), c("Hex", "di", "pi"))

metaDat <- read.csv("createdData/modelSummaries/modelParamsAll.csv", row.names = 1)
#soilMetaDat <- read.csv("createdData/modelSummaries/soilParams.csv", row.names = 1)

# Soils ----
soilModels <- readRDS("createdData/modelSummaries/soilMods.rds")

soilMetaDat <- metaDat %>% dplyr::filter(grepl("soil", file))

soilMetaDat$colName <- paste(soilMetaDat$crop, soilMetaDat$envVar)

denseOutPred <- matrix(nrow = 72406, ncol= length(soilModels), NA) 
highLowPred <- matrix(nrow = 72406, ncol= length(soilModels), NA)
probPNPred <- matrix(nrow = 72406, ncol= length(soilModels), NA)


for(i in 1:length(soilModels)) {
  varSD <- soilMetaDat[i,"sdEnvData"]
  medianTemp <- soilMetaDat[i,"medianEnvData"]
  predDense <- soilModels[[i]]
  var <- soilMetaDat[i,"envVar"]
  crop <- soilMetaDat[i,"crop"]
  
  soilTemp <- soilAll[, var]
  
  denseOutPred[,i] <- replace_na(predDense(soilTemp + varSD) -
                                   predDense(soilTemp - varSD) , 0)
  
  probPNPred[, i] <- ifelse(soilTemp > medianTemp, 1-denseOutPred[,i], denseOutPred[,i]-1)
}

denseOutPred <- as.data.frame(denseOutPred)
probPNPred <- as.data.frame(probPNPred)


colnames(denseOutPred) <- soilMetaDat$colName
colnames(probPNPred) <- soilMetaDat$colName

write.csv(denseOutPred, "soilDenseOutPred.csv")
write.csv(probPNPred, "soilprobPNPred.csv")

# Camelina Weather ----
camModels <- readRDS("createdData/modelSummaries/camelinaModel.rds")

camModels2 <- list(camModels[[1]], camModels[[1]], 
                     camModels[[3]], camModels[[3]],
                     camModels[[5]], camModels[[5]],
                     camModels[[7]], camModels[[7]],
                     camModels[[9]], camModels[[9]],
                     camModels[[11]], camModels[[11]],
                     camModels[[13]], camModels[[13]],
                     camModels[[15]], camModels[[15]],
                     camModels[[17]], camModels[[17]],
                     camModels[[19]], camModels[[19]],
                     camModels[[21]], camModels[[21]],
                   camModels[[1]], camModels[[1]], 
                   camModels[[3]], camModels[[3]],
                   camModels[[5]], camModels[[5]],
                   camModels[[7]], camModels[[7]],
                   camModels[[9]], camModels[[9]],
                   camModels[[11]], camModels[[11]],
                   camModels[[13]], camModels[[13]],
                   camModels[[15]], camModels[[15]],
                   camModels[[17]], camModels[[17]],
                   camModels[[19]], camModels[[19]],
                   camModels[[21]], camModels[[21]])

metaDatCam <- metaDat %>%
  mutate(years = recode(years, "2019_2023" = "19_23", "2069_2073" = "69_73")) %>%
  dplyr::filter(crop == "camelina") %>% 
  mutate(envVar = ifelse(nchar(envVar) == 1, 
                         paste("X0", envVar, sep=""), paste("X", envVar, sep=""))) %>%
  mutate(colName = ifelse(years == "19_23", "pres", "pred")) %>%
  mutate(colName = paste("cam", scenerio, colName, envVar, sep="_")) 

datTemp <- metaDatCam %>%
  filter(scenerio == "lev") %>%
  mutate(sdEnvData = ifelse(years == "69_73", NA, sdEnvData), 
         medianEnvData = ifelse(years == "69_73", NA, medianEnvData), 
         meanEnvData = ifelse(years == "69_73", NA, meanEnvData)) %>%
  fill(sdEnvData, medianEnvData, meanEnvData)

metaDatCam[1:22,3:6] <- datTemp[1:22,3:6]
metaDatCam[23:44,3:6] <- datTemp[1:22,3:6]
metaDatCam <- metaDatCam[1:44,]


denseOutPred <- matrix(nrow = 72406, ncol= length(camModels2), NA) # 44 = 11 var * 2 scenerios * 2 years
probPNPred <- matrix(nrow = 72406, ncol= length(camModels2), NA)

for(i in 1:length(camModels2)) {
  varSD <- metaDatCam[i,"sdEnvData"]
  medianTemp <- metaDatCam[i,"medianEnvData"]
  predDense <- camModels2[[i]]
  var <- metaDatCam[i,"envVar"]
  
  climTemp <- climAll[which(climAll$scenerio == metaDatCam[i,"scenerio"] &
                              climAll$years == metaDatCam[i, "years"]), var]
  
  denseOutPred[,i] <- replace_na(predDense(climTemp + varSD) -
                                   predDense(climTemp - varSD) , 0)
  
  probPNPred[, i] <- ifelse(climTemp > medianTemp, 1-denseOutPred[,i], denseOutPred[,i]-1)
}

denseOutPred <- as.data.frame(denseOutPred)
probPNPred <- as.data.frame(probPNPred)


colnames(denseOutPred) <- metaDatCam$colName[1:44]
colnames(probPNPred) <- metaDatCam$colName[1:44]

write.csv(denseOutPred, "camDenseOutPred.csv")
write.csv(probPNPred, "camprobPNPred.csv")

# Canola Weather ----

metaDat <- read.csv("createdData/modelSummaries/modelParamsAll.csv", row.names = 1)

metaDatCan2 <- metaDat %>%
  mutate(years = recode(years, "2019_2023" = "19_23", "2069_2073" = "69_73")) %>%
  dplyr::filter(crop == "canola") %>% 
  mutate(envVar = ifelse(nchar(envVar) == 1, 
                         paste("X0", envVar, sep=""), paste("X", envVar, sep=""))) %>%
  mutate(colName = ifelse(years == "19_23", "pres", "pred")) %>%
  mutate(colName = paste("can", scenerio, colName, envVar, sep="_"))


canFiles <- list.files("createdData/modelSummaries/", pattern = glob2rx("canola*.rds"))


canModelTemp <- readRDS(paste("createdData/modelSummaries/", canFiles[1], sep=""))

canModels <- c(canModelTemp[1], canModelTemp[1], 
                  canModelTemp[3], canModelTemp[3], 
                  canModelTemp[5], canModelTemp[5])

canModelTemp <- readRDS(paste("createdData/modelSummaries/", canFiles[2], sep=""))

canModels <- c(canModels, canModelTemp[1], canModelTemp[1], 
               canModelTemp[3], canModelTemp[3], 
               canModelTemp[5], canModelTemp[5])

canModelTemp <- readRDS(paste("createdData/modelSummaries/", canFiles[3], sep=""))

canModels <- c(canModels, canModelTemp[1], canModelTemp[1], 
               canModelTemp[3], canModelTemp[3], 
               canModelTemp[5], canModelTemp[5])

canModelTemp <- readRDS(paste("createdData/modelSummaries/", canFiles[4], sep=""))

canModels <- c(canModels, canModelTemp[1], canModelTemp[1], 
               canModelTemp[3], canModelTemp[3])

canModels <- c(canModels, canModels)

datTemp <- metaDatCan2 %>%
  filter(scenerio == "lev") %>%
  mutate(sdEnvData = ifelse(V1 %in% c(2, 4, 6), NA, sdEnvData), 
         medianEnvData = ifelse(V1 %in% c(2, 4, 6), NA, medianEnvData), 
         meanEnvData = ifelse(V1 %in% c(2, 4, 6), NA, meanEnvData)) %>%
  fill(sdEnvData, medianEnvData, meanEnvData)

metaDatCan2[1:22,3:6] <- datTemp[1:22,3:6]
metaDatCan2[23:44,3:6] <- datTemp[1:22,3:6]
metaDatCan2 <- metaDatCan2[1:44,]

  
  denseOutPred <- matrix(nrow = 72406, ncol= length(canModels), NA) # 11 var * 2 scenerios * 2 years
  probPNPred <- matrix(nrow = 72406, ncol= length(canModels), NA)
  
  for(i in 1:length(canModels)) {
    varSD <- metaDatCan2[i,"sdEnvData"]
    medianTemp <- metaDatCan2[i,"medianEnvData"]
    predDense <- canModels[[i]]
    var <- metaDatCan2[i,"envVar"] 
    #print(var)
    
    climTemp <- climAll[which(climAll$scenerio == metaDatCan2[i,"scenerio"] &
                                climAll$years == metaDatCan2[i, "years"]), var]
    
    denseOutPred[,i] <- replace_na(predDense(climTemp + varSD) -
                                     predDense(climTemp - varSD) , 0)
    
    probPNPred[, i] <- ifelse(climTemp > medianTemp, 1-denseOutPred[,i], denseOutPred[,i]-1)
  }


colnames(denseOutPred) <- metaDatCan2[1:44,]$colName
colnames(probPNPred) <- metaDatCan2[1:44,]$colName

write.csv(denseOutPred, "canDenseOutPred.csv")
write.csv(probPNPred, "canprobPNPred.csv")
