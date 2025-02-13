# Analyze from predictions 

rm(list=ls())

library(tidyr); library(tibble); library(dplyr); library(ggplot2); 
library(ranger); library(sf); library(mgcv); library(spData);
library(patchwork); library(purrr); library(forcats); library(scico)

# Equation for t.test from summary statistics 
t.test2 <- function(m1,m2,s1,s2,n1,n2,m0=0,equal.variance=FALSE)
{
  if( equal.variance==FALSE ) 
  {
    se <- sqrt( (s1^2/n1) + (s2^2/n2) )
    # welch-satterthwaite df
    df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
  } else
  {
    # pooled standard deviation, scaled by the sample sizes
    se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) ) 
    df <- n1+n2-2
  }      
  t <- (m1-m2-m0)/se 
  dat <- c(m1-m2, se, t, 2*pt(-abs(t),df))    
  names(dat) <- c("Difference of means", "Std Error", "t", "p-value")
  return(dat) 
}

# Model summaries
modelParams <- read.csv("createdData/modelSummaries/modelParamsAll.csv")

# Read in predictions across US ----
soilDenseOutPred <- read.csv("createdData/predictAcrossUS/soilDenseOutPred.csv")
soilprobPNPred <- read.csv("createdData/predictAcrossUS/soilprobPNPred.csv")

camDenseOutPred <- read.csv("createdData/predictAcrossUS/camDenseOutPred.csv")
camprobPNPred <- read.csv("createdData/predictAcrossUS/camprobPNPred.csv")

canDenseOutPred <- read.csv("createdData/predictAcrossUS/canDenseOutPred.csv")
canprobPNPred <- read.csv("createdData/predictAcrossUS/canprobPNPred.csv")

# Crop area
camBio <- read.csv("createdData/landUse/canCamArea.csv", row.names = 1)



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

# Format Data ----

modelParamsComp <- modelParams %>%
  mutate(years = ifelse(is.na(years), "2019_2023", years),
         scenerio = ifelse(is.na(scenerio), "lev", scenerio)) %>% 
  filter(years == "2019_2023" & scenerio == "lev") %>% 
  dplyr::select(sdEnvData, meanEnvData, lengthEnvData, envVar, crop)

## Predict - Limiting Factor ----

vars <- c("lev_pres", "lev_pred", "wre_pres", "wre_pred")
camLimit <- matrix(nrow = 72406, ncol = length(vars))
canLimit <- matrix(nrow = 72406, ncol = length(vars))

canLimitYN <- list()
camLimitYN <- list()


for (i in 1: length(vars)) {
denseOutPred <- cbind(dplyr::select(camDenseOutPred, contains(vars[i])),
                      dplyr::select(soilDenseOutPred, contains("camelina")))
#camLimit[,i] <- apply(denseOutPred, 1, min)
#camLimitYN[[i]] <- as.data.frame(denseOutPred == camLimit[,i])

tempDF  <- denseOutPred %>%
  mutate(across(everything(), ~ ifelse(. < 0.2, 0, .))) %>%
  rowwise() %>%
  mutate(mean = mean(c_across(starts_with("cam")), na.rm = TRUE))

camLimit[,i] <- tempDF$mean
camLimitYN[[i]] <- as.data.frame(tempDF[,-14] == 0)

}

colnames(camLimit) <- paste("cam", vars, sep="_")

for (i in 1: length(vars)) {
  denseOutPred <- cbind(dplyr::select(canDenseOutPred, contains(vars[i])),
                        dplyr::select(soilDenseOutPred, contains("canola")))
  #canLimit[,i] <- apply(denseOutPred, 1, min)
  #canLimitYN[[i]] <- as.data.frame(denseOutPred == canLimit[,i])
  
  tempDF  <- denseOutPred %>%
    mutate(across(everything(), ~ ifelse(. < 0.2, 0, .))) %>%
    rowwise() %>%
    mutate(mean = mean(c_across(starts_with("can")), na.rm = TRUE))
  
  canLimit[,i] <- tempDF$mean
  canLimitYN[[i]] <- as.data.frame(tempDF[,-14] == 0)
}

colnames(canLimit) <- paste("can", vars, sep="_")

predictDenseAll <- as.data.frame(cbind(camLimit, canLimit)) %>%
  mutate(Hex = row_number()) %>%
  mutate(canolaPctChgWRE = (can_wre_pred - can_wre_pres) / can_wre_pres, 
         canolaPctChgLEV = (can_lev_pred - can_lev_pres) / can_lev_pres,
         camelinaPctChgWRE = (cam_wre_pred - cam_wre_pres) / cam_wre_pres,
         camelinaPctChgLEV = (cam_lev_pred - cam_lev_pres) / cam_lev_pres)
  mutate(canolaChgWRE = can_wre_pred - can_wre_pres, 
        canolaChgLEV = can_lev_pred - can_lev_pres,
       camelinaChgWRE = cam_wre_pred - cam_wre_pres,
      camelinaChgLEV = cam_lev_pred - cam_lev_pres)


temp1 <- as.data.frame(colSums(canLimitYN[[1]])) %>%
  rownames_to_column(var = "var") %>%
  mutate(var = recode(var, "canola.di" = "can_lev_pres_di", "canola.pi" = "can_lev_pres_pi")) %>%
  rename("sum" = `colSums(canLimitYN[[1]])`)
temp2 <- as.data.frame(colSums(canLimitYN[[2]])) %>%
  rownames_to_column(var = "var") %>%
  mutate(var = recode(var, "canola.di" = "can_lev_pred_di", "canola.pi" = "can_lev_pred_pi")) %>%
  rename("sum" = `colSums(canLimitYN[[2]])`)
temp3 <- as.data.frame(colSums(canLimitYN[[3]])) %>%
  rownames_to_column(var = "var") %>%
  mutate(var = recode(var, "canola.di" = "can_wre_pres_di", "canola.pi" = "can_wre_pres_pi")) %>%
  rename("sum" = `colSums(canLimitYN[[3]])`)
temp4 <- as.data.frame(colSums(canLimitYN[[4]])) %>%
  rownames_to_column(var = "var") %>%
  mutate(var = recode(var, "canola.di" = "can_wre_pred_di", "canola.pi" = "can_wre_pred_pi")) %>%
  rename("sum" = `colSums(canLimitYN[[4]])`)
temp5 <- as.data.frame(colSums(camLimitYN[[1]])) %>%
  rownames_to_column(var = "var") %>%
  mutate(var = recode(var, "camelina.di" = "cam_lev_pres_di", "camelina.pi" = "cam_lev_pres_pi")) %>%
  rename("sum" = `colSums(camLimitYN[[1]])`)
temp6 <- as.data.frame(colSums(camLimitYN[[2]])) %>%
  rownames_to_column(var = "var") %>%
  mutate(var = recode(var, "camelina.di" = "cam_lev_pred_di", "camelina.pi" = "cam_lev_pred_pi")) %>%
  rename("sum" = `colSums(camLimitYN[[2]])`)
temp7 <- as.data.frame(colSums(camLimitYN[[3]])) %>%
  rownames_to_column(var = "var") %>%
  mutate(var = recode(var, "camelina.di" = "cam_wre_pres_di", "camelina.pi" = "cam_wre_pres_pi")) %>%
  rename("sum" = `colSums(camLimitYN[[3]])`)
temp8 <- as.data.frame(colSums(camLimitYN[[4]])) %>%
  rownames_to_column(var = "var") %>%
  mutate(var = recode(var, "camelina.di" = "cam_wre_pred_di", "camelina.pi" = "cam_wre_pred_pi")) %>%
  rename("sum" = `colSums(camLimitYN[[4]])`)

limDat <- rbind(temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8) %>%
  separate(var, into= c("Crop", "Scenario", "Years", "var")) %>%
  mutate(Years = recode(Years, "pres" = "2019-2023", "pred" = "2069-2073"), 
         Crop = recode(Crop, "can" = "Canola", "cam" = "Camelina")) %>%
  mutate(var=recode(var, 'X01'='Annual Mean Temp', 'X03'='Isothermality','X12'='Annual Precipitation',
                    'X08'='Mean Temp Wet Q', 'X09'='Mean Temp Dry Q', 
                    'X10'='Mean Temp Warm Q', 'X11'='Mean Temp Cold Q',
                    'X16'='Ppt Wet Q', 'X17'='Ppt Dry Q',
                    'X18'='Ppt Warm Q', 'X19'='Ppt Cold Q', 
                    "di"="Soil DI", "pi" = "Soil PI"), 
         Scenario = recode(Scenario, "lev" = "Stabilizing Emissions Scenario", 
                           "wre" = "High Emissions Scenario")) %>%
  mutate(Scenario = relevel(as.factor(Scenario), ref= "Stabilizing Emissions Scenario")) %>%
  mutate(title = paste(Crop, Scenario, Years))


canLimSum <- as.data.frame(cbind(rowSums(canLimitYN[[1]]),
                              rowSums(canLimitYN[[2]]),
                              rowSums(canLimitYN[[3]]),
                              rowSums(canLimitYN[[4]]))) %>%
  rownames_to_column("Hex") %>%
  rename("Canola Stabilizing Emissions Scenario 2019-2023" = "V1",
         "Canola Stabilizing Emissions Scenario 2069-2073" = "V2",
         "Canola High Emissions Scenario 2019-2023" = "V3", 
         "Canola High Emissions Scenario 2069-2073" = "V4") %>%
  pivot_longer(cols= `Canola Stabilizing Emissions Scenario 2019-2023`:`Canola High Emissions Scenario 2069-2073`,
               names_to = "title", values_to = "sum") 

limMapAllSum <- as.data.frame(cbind(rowSums(camLimitYN[[1]]),
                                 rowSums(camLimitYN[[2]]),
                                 rowSums(camLimitYN[[3]]),
                                 rowSums(camLimitYN[[4]]))) %>%
  rownames_to_column("Hex") %>%
  rename("Camelina Stabilizing Emissions Scenario 2019-2023" = "V1",
         "Camelina Stabilizing Emissions Scenario 2069-2073" = "V2",
         "Camelina High Emissions Scenario 2019-2023" = "V3", 
         "Camelina High Emissions Scenario 2069-2073" = "V4") %>%
  pivot_longer(cols= `Camelina Stabilizing Emissions Scenario 2019-2023`:`Camelina High Emissions Scenario 2069-2073`,
               names_to = "title", values_to = "sum") %>%
  full_join(canLimSum)



predictSp <- stGrid %>%
  mutate(Hex = as.numeric(Hex)) %>%
  left_join(predictDenseAll) 

# limits (pos/ neg)
soilprobPNPred2 <- soilprobPNPred %>%
  mutate(can_wre_pres_di = canola.di, 
         can_wre_pred_di = canola.di,
         can_lev_pres_di = canola.di,
         can_lev_pred_di = canola.di,
         can_wre_pres_pi = canola.pi, 
         can_wre_pred_pi = canola.pi,
         can_lev_pres_pi = canola.pi,
         can_lev_pred_pi = canola.pi,
         cam_wre_pres_di = camelina.di, 
         cam_wre_pred_di = camelina.di,
         cam_lev_pres_di = camelina.di,
         cam_lev_pred_di = camelina.di,
         cam_wre_pres_pi = camelina.pi, 
         cam_wre_pred_pi = camelina.pi,
         cam_lev_pres_pi = camelina.pi,
         cam_lev_pred_pi = camelina.pi) %>%
  dplyr::select(!c(camelina.pi, camelina.di, canola.pi, canola.di))

pnPredMap <- full_join(camprobPNPred, canprobPNPred) %>% 
  full_join(soilprobPNPred2) %>%
  rename("Hex" = "X") %>% 
  pivot_longer(-Hex) %>%
  separate(name, into = c("Crop", "Scenario", "Years", "var")) %>%
  mutate(Years = recode(Years, "pres" = "2019-2023", "pred" = "2069-2073"), 
         Crop = recode(Crop, "can" = "Canola", "cam" = "Camelina")) %>%
  mutate(var=recode(var, 'X01'='Annual Mean Temp', 'X03'='Isothermality','X12'='Annual Precipitation',
                    'X08'='Mean Temp Wet Q', 'X09'='Mean Temp Dry Q', 
                    'X10'='Mean Temp Warm Q', 'X11'='Mean Temp Cold Q',
                    'X16'='Ppt Wet Q', 'X17'='Ppt Dry Q',
                    'X18'='Ppt Warm Q', 'X19'='Ppt Cold Q', 
                    "di"="Soil DI", "pi" = "Soil PI"), 
         Scenario = recode(Scenario, "lev" = "Stabilizing Emissions Scenario", 
                           "wre" = "High Emissions Scenario")) %>%
  pivot_wider(values_from= value, names_from = var) %>%
  mutate(Scenario = relevel(as.factor(Scenario), ref= "Stabilizing Emissions Scenario")) %>%
  mutate(title = paste(Crop, Scenario, Years))

## Plot ----

## Get boundry ----

camelinaArea <- stGrid %>%
  mutate(Hex = as.numeric(Hex)) %>%
  left_join(camBio[,c("Hex", "propCam")]) %>%
  filter(propCam > quantile(na.omit(camBio$propCam), 0.95)) %>%
  st_union()
canolaArea <- stGrid %>%
  mutate(Hex = as.numeric(Hex)) %>%
  left_join(camBio[,c("Hex", "propCan")]) %>%
  filter(propCan > quantile(na.omit(camBio$propCan), 0.95)) %>%
  st_union()

### Fig. 5 ----

orderName <- limDat %>% group_by(var) %>% summarise(top = mean(sum)) %>% arrange(top) %>% dplyr::select(var)

limDat$var <- factor(limDat$var, levels=as.vector(orderName)$var)


ggplot(limDat, aes(x=var, y=sum, fill=(Crop))) +
  geom_col(position=position_dodge2(reverse = TRUE)) + 
  facet_grid(vars(Scenario), vars(Years)) +
  xlab("Environmental Variable") + 
  ylab("Limiting Factor (n hexbins)") + 
  theme_bw(base_size=15) + 
  scale_fill_manual(name = "Year", values = c("yellowgreen", 'goldenrod')) + 
  theme(axis.text.x = element_text(color="black"), 
        axis.text.y = element_text(color="black")) +
  coord_flip(ylim = c(28000, 64000)) 
#ggsave("fig5.tif", width=8, height=6.5)

### Fig. 4 ----
a <- ggplot() + ggtitle("E. Canola Change High Emission Scenario") + 
  theme_grey(base_size = 6) + 
  geom_sf(data = predictSp, aes(fill= canolaChgWRE), color = NA) + 
  geom_sf(data= us_states, fill = NA, color = 'grey', linewidth =0.07) + 
  geom_sf(data= canolaArea, fill = NA, linewidth=0.2, color = 'navy') + 
  scale_fill_gradient2(  low = "red",
                         mid = "white",
                         high = "blue", 
                         limits = c(-0.365, 0.221)) +
  theme(legend.title=element_blank(), 
        legend.position="bottom",legend.key.width = unit(0.5, "cm"),
        axis.text.x = element_text(color="black"), 
        axis.text.y = element_text(color="black"))
b <- ggplot() + ggtitle("D. Camelina Change High Emission Scenario") + 
  theme_grey(base_size = 6) + 
  geom_sf(data = predictSp, aes(fill= camelinaChgWRE), color = NA) + 
  geom_sf(data= us_states, fill = NA, color = 'grey', linewidth =0.07) + 
  geom_sf(data= camelinaArea, fill = NA, linewidth=0.2, color = 'navy') + 
  scale_fill_gradient2(  low = "red",
                         mid = "white",
                         high = "blue", 
                         limits = c(-0.365, 0.221)) +
  theme(legend.title=element_blank(), 
        legend.position="bottom", legend.key.width = unit(0.5, "cm"),
        axis.text.x = element_text(color="black"), 
        axis.text.y = element_text(color="black"))
c <- ggplot() + ggtitle("B. Canola Change Stabilizing Emission Scenario") + 
  theme_grey(base_size = 6) + 
  geom_sf(data = predictSp, aes(fill= canolaChgLEV), color = NA) + 
  geom_sf(data= us_states, fill = NA, color = 'grey', linewidth =0.07) + 
  geom_sf(data= canolaArea, fill = NA, linewidth=0.2, color = 'navy') + 
  scale_fill_gradient2(  low = "red",
                         mid = "white",
                         high = "blue", 
                         limits = c(-0.365, 0.221)) +
  theme(legend.title=element_blank(), 
        legend.position="bottom",legend.key.width = unit(0.5, "cm"),
        axis.text.x = element_text(color="black"), 
        axis.text.y = element_text(color="black"))
d <- ggplot() + ggtitle("A. Camelina Change Stabilizing Emission Scenario") + 
  theme_grey(base_size = 6) + 
  geom_sf(data = predictSp, aes(fill= camelinaChgLEV), color = NA) + 
  geom_sf(data= us_states, fill = NA, color = 'grey', linewidth =0.07) + 
  geom_sf(data= camelinaArea, fill = NA, linewidth=0.2, color = 'navy') + 
  scale_fill_gradient2(  low = "red",
                         mid = "white",
                         high = "blue", 
                         limits = c(-0.365, 0.221))  + #-0.005, 0.004
  theme(legend.title=element_blank(), 
        legend.position="bottom", legend.key.width = unit(0.5, "cm"), 
        axis.text.x = element_text(color="black"), 
        axis.text.y = element_text(color="black"))

e <- ggplot(predictDenseAll) +
  theme_bw(base_size = 6) + 
  theme(axis.text.x = element_text(color="black"), 
        axis.text.y = element_text(color="black"), 
        legend.position="bottom") + 
  #xlim(-0.05, 0.05) + 
  xlab("Change") + 
  ylab("Density") + 
  ggtitle("F. Change Under High\n Emission Scenario") + 
  geom_density(aes(camelinaChgWRE, color="cam"), bw=0.05,  size= 0.2) + #bw = .0005
  geom_density(aes(canolaChgWRE, color="can"),bw=0.05,  size=0.2) + #bw = .0005
  geom_vline(aes(xintercept=mean(camelinaChgWRE)),color= "yellowgreen", linetype= "dashed", size=0.5) + 
  geom_vline(aes(xintercept=mean(canolaChgWRE)), color= "goldenrod", linetype= "dashed", size=0.5) + 
  scale_color_manual(values = c(can="goldenrod", cam= "yellowgreen"), labels = c("Camelina", "Canola"), name="Crop") 

f <- ggplot(predictDenseAll) +
  theme_bw(base_size = 6) + 
  theme(axis.text.x = element_text(color="black"), 
        axis.text.y = element_text(color="black"), 
        legend.position="bottom") + 
  #xlim(-0.05, 0.05) + 
  xlab("Change") + 
  ylab("Density") + 
  ggtitle("C. Change Under Stabilizing\n Emission Scenario") + 
  geom_density(aes(camelinaChgLEV, color="cam"),bw=0.05, size= 0.2) + # bw = .0005
  geom_density(aes(canolaChgLEV, color="can"), bw=0.05,size=0.2) + # bw = .0005
  geom_vline(aes(xintercept=mean(camelinaChgLEV)),color= "yellowgreen", linetype= "dashed", size=0.5) + 
  geom_vline(aes(xintercept=mean(canolaChgLEV)), color= "goldenrod", linetype= "dashed", size=0.5) + 
  scale_color_manual(values = c(can="goldenrod", cam= "yellowgreen"), labels = c("Camelina", "Canola"), name="Crop") 


d + c + f + b + a +  e +  
  plot_layout(guides = "collect", widths = c(10,10, 5)) & theme(legend.position = 'bottom')
#ggsave("fig4.tif", width=9, height=5)

### Fig. 3 ----

a <- ggplot() + ggtitle("D. Canola 2069-2073 High Emission Scenario") + 
  geom_sf(data = predictSp, aes(fill= can_wre_pred), color = NA) + 
  geom_sf(data= us_states, fill = NA, color = 'grey', linewidth =0.07) + 
  geom_sf(data= canolaArea, fill = NA, linewidth=0.2, color = 'white') + 
  #scale_fill_gradient2(limits = c(0, 0.851), guide="none") +
  theme(legend.title=element_blank(), 
        axis.text.x = element_text(color="black"), 
        axis.text.y = element_text(color="black")) +
  scale_fill_viridis(option="magma", limits = c(0, 0.851))
b <- ggplot() + ggtitle("C. Camelina 2069-2073 High Emission Scenario") + 
  geom_sf(data = predictSp, aes(fill= cam_wre_pred), color = NA) + 
  geom_sf(data= us_states, fill = NA, color = 'grey', linewidth =0.07) + 
  geom_sf(data= camelinaArea, fill = NA,linewidth=0.2, color = 'white') + 
  #scale_fill_gradient2(limits = c(0, 0.851), guide="none") +
  theme(legend.title=element_blank(), 
        axis.text.x = element_text(color="black"), 
        axis.text.y = element_text(color="black"))  +
  scale_fill_viridis(option="magma", limits = c(0, 0.851))
c <- ggplot() + ggtitle("B. Canola 2069-2073 Stabilizing Emission Scenario") + 
  geom_sf(data = predictSp, aes(fill= can_lev_pred), color = NA) + 
  geom_sf(data= us_states, fill = NA, color = 'grey', linewidth =0.07) + 
  geom_sf(data= canolaArea, fill = NA, linewidth=0.2, color = 'white') + 
  #scale_fill_gradient2(limits = c(0, 0.851), guide="none") +
  theme(legend.title=element_blank(), 
        axis.text.x = element_text(color="black"), 
        axis.text.y = element_text(color="black"))  +
  scale_fill_viridis(option="magma", limits = c(0, 0.851))
d <- ggplot() + ggtitle("A. Camelina 2069-2073 Stabilizing Emission Scenario") + 
  geom_sf(data = predictSp, aes(fill= cam_lev_pred), color = NA) + 
  geom_sf(data= us_states, fill = NA, color = 'grey', linewidth =0.07) + 
  geom_sf(data= camelinaArea, fill = NA, linewidth=0.2, color = 'white') + 
  #scale_fill_gradient2(limits = c(0, 0.851), guide="none") +
  theme(legend.title=element_blank(), 
        axis.text.x = element_text(color="black"), 
        axis.text.y = element_text(color="black"))  +
  scale_fill_viridis(option="magma", limits = c(0, 0.851))

d + c + b + a + plot_layout(guides = "collect")
#ggsave("fig3Sup.tif", width=10.5, height=8)

d / c + plot_layout(guides = "collect")
#ggsave("fig3.tif", width=7, height=7)

### Fig. 2 ----
a <- ggplot() + ggtitle("Canola 2019-2023 High Emission Scenario") + 
  geom_sf(data = predictSp, aes(fill= can_wre_pres), color = NA) + 
  geom_sf(data= us_states, fill = NA, color = 'grey', linewidth =0.07) + 
  geom_sf(data= canolaArea, fill = NA,  color = 'white', linewidth = 0.2) + 
  #scale_fill_gradient2(limits = c(0, 0.851), guide="none") +
  theme(legend.title=element_blank(), 
        axis.text.x = element_text(color="black"), 
        axis.text.y = element_text(color="black")) +
  scale_fill_viridis(option="magma", limits = c(0, 0.851))
b <- ggplot() + ggtitle("Camelina 2019-2023 High Emission Scenario") + 
  geom_sf(data = predictSp, aes(fill= cam_wre_pres), color = NA) + 
  geom_sf(data= us_states, fill = NA, color = 'grey', linewidth =0.07) + 
  geom_sf(data= camelinaArea, fill = NA,  color = 'white', linewidth =0.2) + 
  #scale_fill_gradient2(limits = c(0, 0.851), guide="none") +
  theme(legend.title=element_blank(), 
        axis.text.x = element_text(color="black"), 
        axis.text.y = element_text(color="black")) +
  scale_fill_viridis(option="magma", limits = c(0, 0.851))
c <- ggplot() + ggtitle("C. Canola 2019-2023 Stabilizing Emission Scenario") + 
  geom_sf(data = predictSp, aes(fill= can_lev_pres), color = NA) + 
  geom_sf(data= us_states, fill = NA, color = 'grey', linewidth =0.07) + 
  geom_sf(data= canolaArea, fill = NA,  color = 'white', linewidth =0.2) + 
  #scale_fill_gradient2(limits = c(0, 0.851), guide="none") +
  theme(legend.title=element_blank(), 
        axis.text.x = element_text(color="black"), 
        axis.text.y = element_text(color="black")) +
  scale_fill_viridis(option="magma", limits = c(0, 0.851))
d <- ggplot() + ggtitle("A. Camelina 2019-2023 Stabilizing Emission Scenario") + 
  geom_sf(data = predictSp, aes(fill= cam_lev_pres), color = NA) + 
  geom_sf(data= us_states, fill = NA, color = 'grey', linewidth =0.07) + 
  geom_sf(data= camelinaArea, fill = NA,  color = 'white', linewidth =0.2) + 
  #scale_fill_gradient2(limits = c(0, 0.851), guide="none")  + # 0.0051
  theme(legend.title=element_blank(), 
        axis.text.x = element_text(color="black"), 
        axis.text.y = element_text(color="black")) +
  scale_fill_viridis(option="magma", limits = c(0, 0.851))

corEnv <- unique(camBio[,c("Hex", "propCam", "propCan")]) %>%
  mutate(Hex = as.numeric(Hex)) %>%
  full_join(predictDenseAll)

e <- ggplot(data= corEnv[which(corEnv$propCan > 0), ], aes(y= log(propCan), x= can_lev_pres)) + 
  geom_point(fill = "goldenrod", pch=21) + 
  annotate("text", x = 0.07, y = -0.05, label = bquote("R"^"2"~"= 0.34")) + 
  geom_smooth(method="lm", color="black") + 
  theme_bw() + 
  ylim(c(-12, 0)) + 
  ylab("Canola Proportion  (log-scale)") + 
  xlab("Probability of Canola") + 
  ggtitle("D. Canola 2019-2023") + 
  theme_bw() +
  theme(axis.text.x = element_text(color="black"), 
        axis.text.y = element_text(color="black"))

f <- ggplot(data= corEnv[which(corEnv$propCam > 0), ], aes(y= log(propCam), x= cam_lev_pres)) + 
  geom_point(fill = "yellowgreen", pch=21) + 
  annotate("text", x = 0.1, y = -2.45, label = bquote("R"^"2"~"= 0.13")) + 
  geom_smooth(method="lm", color="black") + 
  ylab("Camelina Proportion (log-scale)") + 
  xlab("Probability of Camelina") + 
  ggtitle("B. Camelina 2019-2023") + 
  theme_bw() +
  theme(axis.text.x = element_text(color="black"), 
        axis.text.y = element_text(color="black"))
  

d + f + c + e + plot_layout(guides = "collect", widths = c(10, 5)) & theme(legend.position = 'left')
#ggsave("fig2.tif", width=10.5, height=8)

### Fig. 1 ----
fig1Dat <- modelParamsComp %>%
  mutate(envVar=recode(envVar, '01'='Annual Mean Temp', '03'='Isothermality','12'='Annual Precipitation',
                       '08'='Mean Temp Wet Q', '09'='Mean Temp Dry Q', 
                       '10'='Mean Temp Warm Q', '11'='Mean Temp Cold Q',
                       '16'='Ppt Wet Q', '17'='Ppt Dry Q',
                       '18'='Ppt Warm Q', '19'='Ppt Cold Q', 
                       "di"="Soil DI", "pi" = "Soil PI"), 
         crop = recode(crop, "camelina" = "Camelina", "canola" = "Canola")) %>%
  mutate(envVar = as.factor(envVar)) %>%
  mutate(envVar=fct_relevel(envVar,c("Annual Mean Temp", "Mean Temp Wet Q", "Mean Temp Dry Q", "Mean Temp Warm Q", "Mean Temp Cold Q", 
                                     "Annual Precipitation", "Ppt Wet Q", "Ppt Dry Q", "Ppt Warm Q", "Ppt Cold Q", 
                                     "Isothermality", "Soil DI", "Soil PI")))

dat_text <- data.frame(label = rep("p < 0.01", 13),
                       crop = rep("Canola", 13),
                       envVar  = unique(fig1Dat$envVar))

p <- ggplot(fig1Dat, aes(x= crop, y= meanEnvData, color= crop)) + 
  geom_point(size=3) + 
  geom_errorbar(aes(ymin=meanEnvData-sdEnvData, ymax=meanEnvData+sdEnvData), size= 2, width=.2) + 
  facet_wrap(~factor(envVar), ncol=5, scales = "free") +
  scale_x_discrete(labels= c("Camelina", "Canola")) + 
  scale_color_manual(values=c("yellowgreen", "goldenrod"), guide="none") + 
  ylab("Value") + xlab("Crop") + 
  theme_bw(base_size = 15) + 
  theme(axis.text.x = element_text(color="black"), 
        axis.text.y = element_text(color="black"))
 p +  geom_text(data= dat_text, aes(x = crop, y = Inf, label = label), hjust = 1.2, vjust = 2, color= "black")  


#ggsave("fig1.tif", width = 11, height=7)

## Analayis ----
mod <- list()

envIndex <- unique(modelParamsComp$envVar)

for (i in 1:length(envIndex)) {
  m1 <- modelParamsComp[which(modelParamsComp$envVar == envIndex[i] & modelParamsComp$crop == "camelina"), ]$meanEnvData
  m2 <- modelParamsComp[which(modelParamsComp$envVar == envIndex[i] & modelParamsComp$crop == "canola"), ]$meanEnvData
  
  s1 <- modelParamsComp[which(modelParamsComp$envVar == envIndex[i] & modelParamsComp$crop == "camelina"), ]$sdEnvData
  s2 <- modelParamsComp[which(modelParamsComp$envVar == envIndex[i] & modelParamsComp$crop == "canola"), ]$sdEnvData
  
  n1 <- modelParamsComp[which(modelParamsComp$envVar == envIndex[i] & modelParamsComp$crop == "camelina"), ]$lengthEnvData
  n2 <- modelParamsComp[which(modelParamsComp$envVar == envIndex[i] & modelParamsComp$crop == "canola"), ]$lengthEnvData
  
  mod[[i]] <- t.test2(m1,m2,s1,s2,n1,n2) # get same results as t.test (Welch Two Sample t-test)
}

names(mod) <- recode(envIndex, '01'='Annual Mean Temp', '03'='Isothermality','12'='Annual Precipitation',
                     '08'='Mean Temp Wet Q', '09'='Mean Temp Dry Q', 
                     '10'='Mean Temp Warm Q', '11'='Mean Temp Cold Q',
                     '16'='Ppt Wet Q', '17'='Ppt Dry Q',
                     '18'='Ppt Warm Q', '19'='Ppt Cold Q', 
                     "di"="Soil DI", "pi" = "Soil PI")

mod

summary(lm(predictDenseAll$cam_lev_pres ~ predictDenseAll$cam_wre_pres))
summary(lm(predictDenseAll$can_lev_pres ~ predictDenseAll$can_wre_pres))

t.test(predictDenseAll$canolaChgLEV)
t.test(predictDenseAll$canolaChgWRE) # lower than 0
t.test(predictDenseAll$camelinaChgLEV)
t.test(predictDenseAll$camelinaChgWRE)

t.test(predictDenseAll$canolaChgLEV, predictDenseAll$camelinaChgLEV) # ns
t.test(predictDenseAll$canolaChgWRE, predictDenseAll$camelinaChgWRE) # sig 


summary(lm(log(propCan) ~ can_lev_pres, corEnv[which(corEnv$propCan >0),]))
summary(lm(log(propCam) ~ cam_lev_pres, corEnv[which(corEnv$propCam >0),]))

# Canola -- greater suitability, but difference is smaller in the future 
t.test(predictSp$cam_lev_pres,  predictSp$can_lev_pres) # 1.51%
t.test(predictSp$cam_wre_pres, predictSp$can_wre_pres) # 2.31%

t.test(predictSp$cam_lev_pred,  predictSp$can_lev_pred) # 1.07%
t.test(predictSp$cam_wre_pred, predictSp$can_wre_pred) # 1.25%

#### Fig S3 ----

limChg <- limDat %>%
  dplyr::select(Crop, Scenario, Years, var, sum) %>%
  pivot_wider( names_from = Years, values_from = sum) %>%
  mutate(chg = `2069-2073` - `2019-2023`)

orderName2 <- limChg %>% group_by(var) %>% summarise(top = mean(chg)) %>% arrange(top) %>% dplyr::select(var)
limChg$var <- factor(limChg$var, levels=as.vector(orderName2)$var)

ggplot(limChg, aes(x= var, y=chg, fill= Crop)) + 
  geom_col(position=position_dodge2(reverse = TRUE)) + 
  facet_grid(vars(Scenario)) +
  coord_flip() + 
  xlab("Environmental Variable") + 
  ylab("Change in Limiting Factor (n hexbins)") + 
  theme_bw(base_size=15) + 
  scale_fill_manual(name = "Crop", values = c("yellowgreen", 'goldenrod')) + 
  theme(axis.text.x = element_text(color="black"), 
        axis.text.y = element_text(color="black"))
#ggsave("figS3.tif", width=7, height=6.5)


#### Fig S2  ----

e <- ggplot(data= predictDenseAll, aes(y= can_wre_pres, x= can_lev_pres)) + 
  geom_point(fill = "goldenrod", pch=21) + 
  annotate("text", x = 0.15, y = 0.75, label = bquote("R"^"2"~"= 0.97")) + 
  geom_abline(slope= 1, intercept = 0) + 
  ylab("High Emissions Scenario") + 
  xlab("Stabilizing Emissions Scenario") + 
  ggtitle("B. Canola 2019-2023") + 
  theme_bw() +
  theme(axis.text.x = element_text(color="black"), 
        axis.text.y = element_text(color="black"))

f <- ggplot(data= predictDenseAll, aes(y= cam_wre_pres, x= cam_lev_pres)) + 
  geom_point(fill = "yellowgreen", pch=21) + 
  annotate("text", x = 0.2, y = 0.7, label = bquote("R"^"2"~"= 0.96")) + 
  geom_abline(slope= 1, intercept = 0) + 
  ylab("High Emissions Scenario") + 
  xlab("Stabilizing Emissions Scenario") + 
  ggtitle("A. Camelina 2019-2023") + 
  theme_bw() +
  theme(axis.text.x = element_text(color="black"), 
        axis.text.y = element_text(color="black"))

f + e 
#ggsave("figS2-2.tif", width=6, height=4)

### Fig S1 ----
# Distribution of crop
# all the Scenarios/ years are the same distribution 
a <- ggplot(camBio) +
  theme_bw(base_size = 20) + 
  theme(axis.text.x = element_text(color="black"), 
        axis.text.y = element_text(color="black"), 
        legend.position="bottom") + 
  xlab("Camelina Proportional Area") + 
  ylab("Frequency") + 
  ggtitle("A") + 
  #geom_density(aes(Camelina),color="black",  bw = .0005, size= 1) + 
  geom_histogram(aes(x = propCam), bins = 30, fill = "yellowgreen", color = "black") +
  geom_vline(xintercept=quantile(na.omit(camBio$propCam), 0.95)[[1]], 
             color= "black", linetype= "dashed", size=0.5) 

b <- ggplot(camBio) +
  theme_bw(base_size = 20) + 
  theme(axis.text.x = element_text(color="black"), 
        axis.text.y = element_text(color="black"), 
        legend.position="bottom") + 
  xlab("Canola Proportional Area") + 
  ylab("Frequency") + 
  ggtitle("B") + 
  #geom_density(aes(Canola), color="black",  bw = .0005, size= 1) + 
  geom_histogram(aes(x = propCan), bins = 30, fill = "goldenrod", color = "black") +
  geom_vline(xintercept=quantile(na.omit(camBio$propCan), 0.95)[[1]], 
             color= "black", linetype= "dashed", size=0.5) 
a + b 
#ggsave("figs1-2.tif", width=10, height=5)

### Fig. S4 ----
limSp <- stGrid %>%
  mutate(Hex = as.integer(Hex)) %>%
  left_join(pnPredMap) 

#cols <- c("#8A215D", "#503D8B", "#589C98")
#cols <- c("#B5C65C", "#E9C868", "#BE4F45")
#cols <- c("white", "lightblue", "black")
#cols <- c("white", "blue", "black")
#cols <- c("violet", "orange", "cyan")
#cols <- c("red", "yellow", "blue")
#colsTemp <- c("blue", "yellow", "red")


a <- ggplot() + ggtitle("A. Limiting Factor: Soil PI") + 
  theme_grey(base_size = 6) + 
  #scale_fill_gradientn(colours = cols, na.value=NA,
   #                    values = scales::rescale(c(-1, -0.75, 0, 0.75, 1))) +
  scale_fill_scico(palette = 'managua', limits = c(-1,1)) + 
  geom_sf(data = na.omit(limSp), aes(fill= `Soil PI`), color = NA) + 
  geom_sf(data= us_states, fill = NA, color = 'grey') + 
  #guides(fill="none") + 
  theme(legend.title=element_blank(), 
        legend.position="bottom",legend.key.width = unit(0.5, "cm"),
        axis.text.x = element_text(color="black"), 
        axis.text.y = element_text(color="black")) + 
  facet_wrap(~title, ncol= 2)

b <- ggplot() + ggtitle("B. Limiting Factor: Soil DI") + 
  theme_grey(base_size = 6) + 
  #scale_fill_gradientn(colours = cols, na.value=NA,
   #                    values = scales::rescale(c(-1, -0.75, 0, 0.75, 1))) +
  scale_fill_scico(palette = 'managua', limits = c(-1,1)) +
  geom_sf(data = na.omit(limSp), aes(fill= `Soil DI`), color = NA) + 
  geom_sf(data= us_states, fill = NA, color = 'grey') + 
  #guides(fill="none") + 
  theme(legend.title=element_blank(), 
        legend.position="bottom",legend.key.width = unit(0.5, "cm"),
        axis.text.x = element_text(color="black"), 
        axis.text.y = element_text(color="black")) + 
  facet_wrap(~title, ncol= 2)

c <- ggplot() + ggtitle("C. Limiting Factor: Annual Mean Temp") + 
  theme_grey(base_size = 6) + 
  #scale_fill_gradientn(colours = colsTemp, na.value=NA,
   #                    values = scales::rescale(c(-1, -0.75, 0, 0.75, 1))) +
  scale_fill_scico(palette = 'managua', limits = c(-1,1), direction=-1) +
  geom_sf(data = na.omit(limSp), aes(fill= `Annual Mean Temp`), color = NA) + 
  geom_sf(data= us_states, fill = NA, color = 'grey') + 
  #guides(fill="none") + 
  theme(legend.title=element_blank(), 
        legend.position="bottom",legend.key.width = unit(0.5, "cm"),
        axis.text.x = element_text(color="black"), 
        axis.text.y = element_text(color="black")) + 
  facet_wrap(~title, ncol= 2)

d <- ggplot() + ggtitle("D. Limiting Factor: Isothermality") + 
  theme_grey(base_size = 6) + 
  #scale_fill_gradientn(colours = cols, na.value=NA,
   #                    values = scales::rescale(c(-1, -0.75, 0, 0.75, 1))) +
  scale_fill_scico(palette = 'managua', limits = c(-1,1))+
  geom_sf(data = na.omit(limSp), aes(fill= `Isothermality`), color = NA) + 
  geom_sf(data= us_states, fill = NA, color = 'grey') + 
  #guides(fill="none") + 
  theme(legend.title=element_blank(), 
        legend.position="bottom",legend.key.width = unit(0.5, "cm"),
        axis.text.x = element_text(color="black"), 
        axis.text.y = element_text(color="black")) + 
  facet_wrap(~title, ncol= 2)

e <- ggplot() + ggtitle("E. Limiting Factor: Mean Temp Wet Q") + 
  theme_grey(base_size = 6) + 
  #scale_fill_gradientn(colours = colsTemp, na.value=NA,
   #                    values = scales::rescale(c(-1, -0.75, 0, 0.75, 1))) +
  scale_fill_scico(palette = 'managua', limits = c(-1,1), direction=-1) +
  geom_sf(data = na.omit(limSp), aes(fill= `Mean Temp Wet Q`), color = NA) + 
  geom_sf(data= us_states, fill = NA, color = 'grey') + 
  #guides(fill="none") + 
  theme(legend.title=element_blank(), 
        legend.position="bottom",legend.key.width = unit(0.5, "cm"),
        axis.text.x = element_text(color="black"), 
        axis.text.y = element_text(color="black")) + 
  facet_wrap(~title, ncol= 2)

f <- ggplot() + ggtitle("F. Limiting Factor: Mean Temp Dry Q") + 
  theme_grey(base_size = 6) + 
  #scale_fill_gradientn(colours = colsTemp, na.value=NA,
   #                    values = scales::rescale(c(-1, -0.75, 0, 0.75, 1))) +
  scale_fill_scico(palette = 'managua', limits = c(-1,1), direction=-1) +
  geom_sf(data = na.omit(limSp), aes(fill= `Mean Temp Dry Q`), color = NA) + 
  geom_sf(data= us_states, fill = NA, color = 'grey') + 
  #guides(fill="none") + 
  theme(legend.title=element_blank(), 
        legend.position="bottom",legend.key.width = unit(0.5, "cm"),
        axis.text.x = element_text(color="black"), 
        axis.text.y = element_text(color="black")) + 
  facet_wrap(~title, ncol= 2)

g <- ggplot() + ggtitle("G. Limiting Factor: Mean Temp Warm Q") + 
  theme_grey(base_size = 6) + 
  #scale_fill_gradientn(colours = colsTemp, na.value=NA,
   #                    values = scales::rescale(c(-1, -0.75, 0, 0.75, 1))) +
  scale_fill_scico(palette = 'managua', limits = c(-1,1), direction=-1) +
  geom_sf(data = na.omit(limSp), aes(fill= `Mean Temp Warm Q`), color = NA) + 
  geom_sf(data= us_states, fill = NA, color = 'grey') + 
  #guides(fill="none") + 
  theme(legend.title=element_blank(), 
        legend.position="bottom",legend.key.width = unit(0.5, "cm"),
        axis.text.x = element_text(color="black"), 
        axis.text.y = element_text(color="black")) + 
  facet_wrap(~title, ncol= 2)

h <- ggplot() + ggtitle("H. Limiting Factor: Mean Temp Cold Q") + 
  theme_grey(base_size = 6) + 
  #scale_fill_gradientn(colours = colsTemp, na.value=NA,
   #                    values = scales::rescale(c(-1, -0.75, 0, 0.75, 1))) +
  scale_fill_scico(palette = 'managua', limits = c(-1,1), direction=-1) +
  geom_sf(data = na.omit(limSp), aes(fill= `Mean Temp Cold Q`), color = NA) + 
  geom_sf(data= us_states, fill = NA, color = 'grey') + 
  #guides(fill="none") + 
  theme(legend.title=element_blank(), 
        legend.position="bottom",legend.key.width = unit(0.5, "cm"),
        axis.text.x = element_text(color="black"), 
        axis.text.y = element_text(color="black")) + 
  facet_wrap(~title, ncol= 2)

i <- ggplot() + ggtitle("I. Limiting Factor: Annual Precipitation") + 
  theme_grey(base_size = 6) + 
  #scale_fill_gradientn(colours = cols, na.value=NA,
   #                    values = scales::rescale(c(-1, -0.75, 0, 0.75, 1))) +
  scale_fill_scico(palette = 'managua', limits = c(-1,1)) +
  geom_sf(data = na.omit(limSp), aes(fill= `Annual Precipitation`), color = NA) + 
  geom_sf(data= us_states, fill = NA, color = 'grey') + 
  #guides(fill="none") + 
  theme(legend.title=element_blank(), 
        legend.position="bottom",legend.key.width = unit(0.5, "cm"),
        axis.text.x = element_text(color="black"), 
        axis.text.y = element_text(color="black")) + 
  facet_wrap(~title, ncol= 2)

j <- ggplot() + ggtitle("J. Limiting Factor: Ppt Wet Q") + 
  theme_grey(base_size = 6) + 
  #scale_fill_gradientn(colours = cols, na.value=NA,
   #                    values = scales::rescale(c(-1, -0.75, 0, 0.75, 1))) +
  scale_fill_scico(palette = 'managua', limits = c(-1,1)) +
  geom_sf(data = na.omit(limSp), aes(fill= `Ppt Wet Q`), color = NA) + 
  geom_sf(data= us_states, fill = NA, color = 'grey') + 
  #guides(fill="none") + 
  theme(legend.title=element_blank(), 
        legend.position="bottom",legend.key.width = unit(0.5, "cm"),
        axis.text.x = element_text(color="black"), 
        axis.text.y = element_text(color="black")) + 
  facet_wrap(~title, ncol= 2)

k <- ggplot() + ggtitle("K. Limiting Factor: Ppt Dry Q") + 
  theme_grey(base_size = 6) + 
  #scale_fill_gradientn(colours = cols, na.value=NA,
   #                    values = scales::rescale(c(-1, -0.75, 0, 0.75, 1))) +
  scale_fill_scico(palette = 'managua', limits = c(-1,1)) +
  geom_sf(data = na.omit(limSp), aes(fill= `Ppt Dry Q`), color = NA) + 
  geom_sf(data= us_states, fill = NA, color = 'grey') + 
  #guides(fill="none") + 
  theme(legend.title=element_blank(), 
        legend.position="bottom",legend.key.width = unit(0.5, "cm"),
        axis.text.x = element_text(color="black"), 
        axis.text.y = element_text(color="black")) + 
  facet_wrap(~title, ncol= 2)

l <- ggplot() + ggtitle("L. Limiting Factor: Ppt Warm Q") + 
  theme_grey(base_size = 6) + 
  #scale_fill_gradientn(colours = cols, na.value=NA,
   #                    values = scales::rescale(c(-1, -0.75, 0, 0.75, 1))) +
  scale_fill_scico(palette = 'managua', limits = c(-1,1)) +
  geom_sf(data = na.omit(limSp), aes(fill= `Ppt Warm Q`), color = NA) + 
  geom_sf(data= us_states, fill = NA, color = 'grey') + 
  #guides(fill="none") + 
  theme(legend.title=element_blank(), 
        legend.position="bottom",legend.key.width = unit(0.5, "cm"),
        axis.text.x = element_text(color="black"), 
        axis.text.y = element_text(color="black")) + 
  facet_wrap(~title, ncol= 2)

m <- ggplot() + ggtitle("M. Limiting Factor: Ppt Cold Q") + 
  theme_grey(base_size = 6) + 
  #scale_fill_gradientn(colours = cols, na.value=NA,
   #                    values = scales::rescale(c(-1, -0.75, 0, 0.75, 1))) +
  scale_fill_scico(palette = 'managua', limits = c(-1,1)) +
  geom_sf(data = na.omit(limSp), aes(fill= `Ppt Cold Q`), color = NA) + 
  geom_sf(data= us_states, fill = NA, color = 'grey') + 
  #guides(fill="none") + 
  theme(legend.title=element_blank(), 
        legend.position="bottom",legend.key.width = unit(0.5, "cm"),
        axis.text.x = element_text(color="black"), 
        axis.text.y = element_text(color="black")) + 
  facet_wrap(~title, ncol= 2)

a + b + c 
#ggsave("figS4-1.tif", width= 12, height=5)

d + e + f 
#ggsave("figS4-2.tif", width= 12, height=5)

g + h + i 
#ggsave("figS4-3.tif", width= 12, height=5)

j + k +l 
#ggsave("figS4-4.tif", width= 12, height=5)

m
#ggsave("figS4-5.tif", width= 4, height=5)


### Fig. S.5 ----
limSp <- stGrid %>%
  left_join(limMapAllSum) 

ggplot() + ggtitle("Number of Limiting Factors") + 
  theme_grey(base_size = 6) + 
  #scale_fill_gradientn(colours = c("white", "darkblue"), na.value=NA,
   #                    values = scales::rescale(c(1, 10))) +
  geom_sf(data = na.omit(limSp), aes(fill= sum), color = NA) + 
  geom_sf(data= us_states, fill = NA, color = 'grey') + 
  #guides(fill="none") + 
  theme(legend.title=element_blank(), 
        legend.position="bottom",legend.key.width = unit(0.5, "cm"),
        axis.text.x = element_text(color="black"), 
        axis.text.y = element_text(color="black")) + 
  facet_wrap(~title, ncol= 2) + 
  scale_fill_viridis(option="magma", direction = -1)

#ggsave("figS5.tif", width= 4, height=5)

# Precip Truncated? ----
dat <- read.csv("allPointValues_30m/weatherConditions/weatherCam1.csv", row.names = 1)
dat <- read.csv("allPointValues_30m/weatherConditions/weatherCam2.csv", row.names = 1)

if(ncol(dat) == 50) {dat <- dat %>%
  mutate(V1 = rowMeans(dplyr::select(dat, V1:V5), na.rm = TRUE), 
         V2 = rowMeans(dplyr::select(dat, V6:V10), na.rm = TRUE), 
         V3 = rowMeans(dplyr::select(dat, V11:V15), na.rm = TRUE),
         V4 = rowMeans(dplyr::select(dat, V16:V20), na.rm = TRUE),
         V5 = rowMeans(dplyr::select(dat, V21:V25), na.rm = TRUE),
         V6 = rowMeans(dplyr::select(dat, V26:V30), na.rm = TRUE),
         V7 = rowMeans(dplyr::select(dat, V31:V35), na.rm = TRUE),
         V8 = rowMeans(dplyr::select(dat, V36:V40), na.rm = TRUE), 
         V9 = rowMeans(dplyr::select(dat, V41:V45), na.rm = TRUE), 
         V10 = rowMeans(dplyr::select(dat, V46:V50), na.rm = TRUE)) %>%
  dplyr::select(V1:V10)}

if(ncol(dat) == 20) {dat <- dat %>%
  mutate(V1 = rowMeans(dplyr::select(dat, V1:V5), na.rm = TRUE), 
         V2 = rowMeans(dplyr::select(dat, V6:V10), na.rm = TRUE), 
         V3 = rowMeans(dplyr::select(dat, V11:V15), na.rm = TRUE),
         V4 = rowMeans(dplyr::select(dat, V16:V20), na.rm = TRUE)) %>%
  dplyr::select(V1:V4)}

datTemp <- dat[,c(5, 7, 9)]
datTemp2 <- cbind(datTemp, dat[,c(1, 3)])


plot(density(datTemp2$V5))

a <- ggplot(datTemp2, aes(x=V5)) + 
  geom_density() + ggtitle("Annual Precipitation") +
  xlab("Precipitation (mm)") + ylab("Density") + 
  theme_bw(base_size = 10) + 
  theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))
b <- ggplot(datTemp2, aes(x=V7)) + 
  geom_density() + ggtitle("Ppt Wet Q") +
  xlab("Precipitation (mm)") + ylab("Density") + 
  theme_bw(base_size = 10) + 
  theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))
c <- ggplot(datTemp2, aes(x=V9)) + 
  geom_density() + ggtitle("Ppt Dry Q") +
  xlab("Precipitation (mm)") + ylab("Density") + 
  theme_bw(base_size = 10) + 
  theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))
d <- ggplot(datTemp2, aes(x=V1)) + 
  geom_density() + ggtitle("Ppt Warm Q") +
  xlab("Precipitation (mm)") + ylab("Density") + 
  theme_bw(base_size = 10) + 
  theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))
e <- ggplot(datTemp2, aes(x=V3)) + 
  geom_density() + ggtitle("Ppt Cold Q") +
  xlab("Precipitation (mm)") + ylab("Density") + 
  theme_bw(base_size = 10) + 
  theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))


(a+b+c)/ (d+e)

# Compare to US 

climUS <- climAll %>%
  filter(years == "19_23") %>%
  filter(scenerio == "lev") %>%
  dplyr::select(X12:X19) 

#climCam<- read.csv("/Users/gabriela.quinlan/Desktop/pptCamDat.csv", row.names = 1)

colnames(datTemp2) <- c("X12", "X16", "X17", "X18", "X19")

climComp <-climCam %>%
  mutate(range = "Cam") %>%
  full_join(climUS) %>%
  mutate(range = ifelse(is.na(range), "US", range)) %>%
  pivot_longer(cols = X12:X19, names_to = "var") %>%
  mutate(var=recode(var, 'X12'='Annual Precipitation','X16'='Ppt Wet Q', 'X17'='Ppt Dry Q',
                    'X18'='Ppt Warm Q', 'X19'='Ppt Cold Q'))


ggplot(climComp, aes(x=range, y=value, fill=(range))) +
  geom_boxplot() + 
  facet_wrap(~var, ncol= 2, scales="free") + 
  xlab("Range") + 
  ylab("Value") + 
  theme_bw(base_size=15) + 
  scale_fill_manual(name = "Range", values = c("lightblue", 'navy')) + 
  theme(axis.text.x = element_text(color="black"), 
        axis.text.y = element_text(color="black"))

summary(lm(value ~ range, climComp[which(climComp$var == "Annual Precipitation"),]))
summary(lm(value ~ range, climComp[which(climComp$var == "Ppt Cold Q"),]))
summary(lm(value ~ range, climComp[which(climComp$var == "Ppt Dry Q"),]))
summary(lm(value ~ range, climComp[which(climComp$var == "Ppt Warm Q"),]))
summary(lm(value ~ range, climComp[which(climComp$var == "Ppt Wet Q"),]))

