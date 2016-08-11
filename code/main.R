# 
#                                                                             
#                 FEDERAL RURAL UNIVERSITY OF RIO DE JANEIRO                  
#                                   and                                       
#                      ISRIC - WORLD SOIL INFORMATION                         
#                                                                             
#         EVALUATION OF BRAZILIAN SOILS WITH HIGH IRON OXIDES CONTENT         
#                                                                             
#                                 WCSS 2014                                   
#                                                                             
#                           ALESSANDRO SAMUEL-ROSA                            
#                                                                             
#               Wageningen, the Netherlands, January of 2014.                 
#                                                                             
# ----------------------------------------------------------------------------
#                                                                             
# Description:                                                                
# Research work aiming at evaluating Brazilian soils with high iron oxides    
# content in order to help defining adequate criteria to classify high iron   
# oxides-containing soils in both WRB and USC.  The results were submited for 
# publication in the WCSS 2014.                                               
#                                                                             
# ----------------------------------------------------------------------------
#                                                                             
# Co-authors:                                                                 
# Dr. Lúcia Helena Cunha dos Anjos (UFRRJ)                                    
#                                                                             
# ----------------------------------------------------------------------------
#                                                                             
#                             Working R script                                
#                                                                             
# Description:                                                                
# Soil data was retrieved from the public database maintained by Esalq in     
# Brazil. Only soil records containing total iron oxides content were used in 
# the analysis. Duplicated soil profiles were authomatically removed from de  
# database. Soil profiles falling outside Brazilian boundaries were also      
# removed. The probability of finding a layer ≥ 30 cm thick and starting      
# ≤ 100 cm of the soil surface, with sulphic acid-extracted iron oxides > 18%,
# was predicted for the whole country using ordinary kriging. The maximum     
# total iron oxides content in Brazilian soils was also predicted using       
# ordinary kriging. A resolution of about 10 km (0.1 degres), or a scale of   
# about 1:40,000,000, was used for the predictions.                           
#                                                                             
# WRB definitions:                                                            
# - Ferritic: having a layer, ≥ 30 cm thick and starting ≤ 100 cm of the soil 
# surface, with Fedith in the fine earth fraction of ≥ 10% and not forming    
# part of a petroplinthic, pisoplinthic or plinthic horizon.                  
# - Hyperferritic: having a layer, ≥ 30 cm thick and starting ≤ 100 cm of the 
# soil surface, with Fedith in the fine earth fraction of ≥ 30% and not       
# forming part of a petroplinthic, pisoplinthic or plinthic horizon.          
#                                                                             
# ----------------------------------------------------------------------------
#                                                                             
# e-mail: alessandrosamuel@yahoo.com.br                                       
# homepage: soil-scientist.net                                                
#                                                                             
# 

# INITIAL CONFIGURATION ########################################################

rm(list = ls())

# Load Packages ================================================================
require(stringr)
require(sp)
require(rgdal)
require(raster)
require(gstat)
require(aqp)
require(GSIF)
require(rgeos)
require(lattice)
require(plotKML)

# Load and Prepare Data ========================================================

load("data/R/wcss2014iron.RData")

brasil <- shapefile("data/gis/brasil.shp")
brasil.grid <- vect2rast(brasil, cell.size = 0.1, method = "SAGA")
proj4string(brasil.grid) <- proj4string(brasil)

meta <- read.table("data/soil/meta.csv", head = TRUE, sep = ";")
data <- read.table("data/soil/data.csv", head = TRUE, sep = "\t", na.strings = c("", "NA"))
str(data)
head(data)

# create unique key
data$key <- paste(data$Source, "_p", data$OrgProfID, sep = "")
head(data$key, 20)

# select records with Fe data
data <- data[is.na(data$TotFe2O3) == FALSE, ]

# select according to wrb: horizon thickness and depth
fedata <- data[data$HzDeFn-data$HzDeIn >= 30 & data$HzDeIn <= 100, ]

# select according total Fe content - SiBCS (Fe >= 180 %)
fedata <- fedata[fedata$TotFe2O3 >= 18, ]

# select records with soil class
fedata <- fedata[is.na(fedata$SoilClass) == FALSE, ]

# remove undesirable soil classes
fedata <- fedata[-grep(c("LATERITA"), fedata$SoilClass, ignore.case = TRUE, value = FALSE), ]
fedata <- fedata[-grep(c("Concre"), fedata$SoilClass, ignore.case = TRUE, value = FALSE), ]
fedata <- fedata[-grep(c("planosol"), fedata$SoilClass, ignore.case = TRUE, value = FALSE), ]
fedata <- fedata[-grep(c("litólico"), fedata$SoilClass, ignore.case = TRUE, value = FALSE), ]
fedata <- fedata[-grep(c("litossolo"), fedata$SoilClass, ignore.case = TRUE, value = FALSE), ]
fedata <- fedata[-grep(c("pedrego"), fedata$SoilClass, ignore.case = TRUE, value = FALSE), ]
length(unique(fedata$key))

# classify according total Fe content - SiBCS (ferric and perferric)
ferric <- fedata[fedata$TotFe2O3 >= 18 & fedata$TotFe2O3 < 36, ]
perferric <- fedata[fedata$TotFe2O3 >= 36, ]

# plot histograms
pdf("res/fig/ferric.pdf")
hist(ferric$TotFe2O3, 20, col = "lightgray", main = "Ferric", xlab = "Total Fe2O3 (sulfuric attack) (%)")
rug(ferric$TotFe2O3)
dev.off()

pdf("res/fig/perferric.pdf")
hist(perferric$TotFe2O3, 20, col = "lightgray", main = "Perferric", xlab = "Total Fe2O3 (sulfuric attack) (%)")
rug(perferric$TotFe2O3)
dev.off()

# PREPARE DATA FOR MODELING ####################################################

# Create spatial objects =======================================================
coordinates(data) <- ~ Long + Lat
proj4string(data) <- proj4string(brasil)
coordinates(fedata) <- ~ Long + Lat
proj4string(fedata) <- proj4string(brasil)
coordinates(ferric) <- ~ Long + Lat
proj4string(ferric) <- proj4string(brasil)
coordinates(perferric) <- ~ Long + Lat
proj4string(perferric) <- proj4string(brasil)

# plot Fe data
pdf("res/fig/ferric-perferric.pdf")
plot(brasil)
plot(ferric, add = TRUE, cex = 0.5, col = "red")
plot(perferric, add = TRUE, cex = 0.8, col = "blue", pch = 20)
title(main = "Total iron content in Brazilian soils", sub = "Red = ferric. Blue = Perferric.")
grid(); axis(1); axis(2)
dev.off()

pdf("res/fig/all-brasil.pdf")
plot(brasil)
plot(data, add = TRUE, cex = 0.3, col = "red")
title(main = "Soil profiles with iron data available")
grid(); axis(1); axis(2)
dev.off()

# MODEL FITTING ################################################################

# Probability of Meeting WRB and SiBCS Criteria ================================

# create indicator variable (proba) --------------------------------------------
proba <- data.frame(coordinates(data), data$key)
colnames(proba) <- c("Long", "Lat", "key")
proba <- proba[!duplicated(proba$key),]
proba$match <- match(proba$key, fedata$key)
proba$proba <- !is.na(proba$match)
head(proba)

# create spatial object and mask to Brazilian boundary -------------------------
coordinates(proba) <- ~ Long + Lat
proj4string(proba) <- proj4string(brasil)
plot(proba)
cut <- over(proba, brasil)[1]
proba <- as.data.frame(proba)
proba <- proba[!is.na(cut),]; rm(cut)
coordinates(proba) <- ~ Long + Lat
proj4string(proba) <- proj4string(brasil)
proba <- remove.duplicates(proba)
spplot(proba, "proba", scales = list(draw = TRUE))

# fit variogram model ----------------------------------------------------------

# first settings
print(plot(variogram(proba ~ 1, loc = proba, cutoff = 350), plot.numbers = TRUE)) # variogram
proba.iv <- variogram(proba ~ 1, loc = proba, cutoff = 350)

# eye fit
print(show.vgms())
proba.ivm <- vgm(psill = 0.030, model = "Gau", range = 300/sqrt(3), nugget = 0.025)
print(plot(proba.iv, main = "Indicator variogram", pch = 20, col = "blue",
           xlim = c(0, max(proba.iv$dist)*1.1),
           ylim = c(0, max(proba.iv$gamma)*1.1),
           plot.numbers = T, model = proba.ivm))

# gstat fit
proba.ivmf <- fit.variogram(object = proba.iv, model = proba.ivm, warn.if.neg = TRUE)
print(plot(proba.iv, main = "Indicator variogram", pch = 20, col = "blue",
           xlim = c(0, max(proba.iv$dist)*1.1),
           ylim = c(0, max(proba.iv$gamma)*1.1),
           plot.numbers = T, model = proba.ivmf))

# make predictions -------------------------------------------------------------
proba.ik <- krige(
  formula = proba ~ 1, locations = proba, newdata = brasil.grid, indicators = TRUE, model = proba.ivmf,
  maxdist = proba.ivmf$range[2]*sqrt(3))

# evaluate predictions ---------------------------------------------------------
summary(proba.ik$var1.pred)

# normalization ----------------------------------------------------------------
proba.ik$var1.pred <- pmin(1, proba.ik$var1.pred)
proba.ik$var1.pred <- pmax(0, proba.ik$var1.pred)
summary(proba.ik$var1.pred)

# plots ------------------------------------------------------------------------

# probability map
pts <- list("sp.points", proba, pch = 20, cex = 0.1,
            col = ifelse(proba$proba, "red", "blue"))
lim <-list("sp.polygons", brasil, col = "black")
dev.off()
pdf(file="res/fig/proba.pdf")
print(spplot(proba.ik, zcol = "var1.pred", at = seq(0, 1, by = 0.1),
             col.regions = rev(heat.colors(10)),
             main = "Probability map", xlab = "Long", ylab = "Lat",
             sub = "Local indicator kriging. Red: TRUE. Blue: FALSE.",
             sp.layout = list(pts, lim), scales = list(draw = TRUE)))
dev.off();rm(pts, lim)

# purity map
purity <- proba.ik
purity$var1.pred[purity$var1.pred < 0.5] <- 1 - purity$var1.pred[purity$var1.pred < 0.5]
summary(purity$var1.pred)
pts <- list("sp.points", proba, pch = 20, cex = 0.1,
            col = ifelse(proba$proba, "red", "blue"))
lim <-list("sp.polygons", brasil, col = "black")
dev.off()
pdf(file = "res/fig/purity.pdf")
print(spplot(purity, zcol = "var1.pred", at = seq(0.5, 1, by = 0.001),
             col.regions = gray(seq(0, 1, by = 0.002)),
             main = "Purity map", xlab = "Long", ylab = "Lat",
             sub = "Local indicator kriging\n Red: TRUE. Blue: FALSE.",
             sp.layout = list(pts, lim), scales = list(draw = TRUE)))
dev.off();rm(pts, lim)

# variance map
pts <- list("sp.points", proba, pch = 20, cex = 0.1,
            col = ifelse(proba$proba, "red", "blue"))
lim <- list("sp.polygons", brasil, col = "black")
dev.off()
pdf(file = "res/fig/variance.pdf")
print(spplot(proba.ik, zcol = "var1.var", col.regions = rev(heat.colors(64)),
             main = "Variance", xlab = "Long", ylab = "Lat",
             sub = "Local indicator kriging. Red: TRUE. Blue: FALSE.",
             sp.layout = list(pts, lim), scales = list(draw = TRUE)))
dev.off();rm(pts, lim)

# Maximum Soil Iron Content ====================================================

# record the maximum Fe content in each profile
high <- data.frame(unique(data$key),
                   as.numeric(by(data$TotFe2O3, data$key, max)),
                   as.numeric(by(data$Lat, data$key, max)),
                   as.numeric(by(data$Long, data$key, max)))
colnames(high) <- c(" key", "TotFe2O3", "Lat", "Long")

# create spatial object and mask to Brazilian boundary -------------------------
coordinates(high) <- ~ Long + Lat
proj4string(high) <- proj4string(brasil)
plot(high)
cut <- over(high, brasil)[1]
high <- as.data.frame(high)
high <- high[!is.na(cut),]; rm(cut)
coordinates(high) <- ~ Long + Lat
proj4string(high) <- proj4string(brasil)
high <- remove.duplicates(high)
spplot(high, "TotFe2O3", scales = list(draw = TRUE))

# fit variogram model ----------------------------------------------------------

# first settings
print(plot(variogram(TotFe2O3 ~ 1, loc = high, cutoff = 400),
           plot.numbers = TRUE)) # variogram
high.v <- variogram(TotFe2O3 ~ 1, loc = high, cutoff = 400)

# eye fit
print(show.vgms())
high.vm <- vgm(psill = 20, model = "Pen", range = 350, nugget = 23)
print(plot(high.v, main = "Variogram", pch = 20, col = "blue",
           xlim = c(0, max(high.v$dist)*1.1),
           ylim = c(0, max(high.v$gamma)*1.1),
           plot.numbers = T, model = high.vm))

# gstat fit
high.vmf <- fit.variogram(object = high.v, model = high.vm,
                            warn.if.neg = TRUE)
print(plot(high.v, main = "Variogram", pch = 20, col = "blue",
           xlim = c(0, max(high.v$dist)*1.1),
           ylim = c(0, max(high.v$gamma)*1.1),
           plot.numbers = T, model = high.vmf))

# make predictions -------------------------------------------------------------
high.k <- krige(formula = TotFe2O3 ~ 1, locations = high, newdata = brasil.grid,
                model = high.vmf, maxdist = high.vmf$range[2])

# evaluate predictions ---------------------------------------------------------
summary(high.k$var1.pred)

# plots ------------------------------------------------------------------------

# total iron oxide content
pts <- list("sp.points", high, pch = 20, cex = 0.1, col = "black")
lim <-list("sp.polygons", brasil, col = "black")
dev.off()
pdf(file="res/fig/high.pdf")
print(spplot(high.k, zcol = "var1.pred", col.regions = rev(heat.colors(64)),
             main = "Predicted maximum Fe2O3 content", xlab = "Long", ylab = "Lat",
             sub = "Local kriging.", sp.layout = list(pts, lim),
             scales = list(draw = TRUE)))
dev.off();rm(pts, lim)

# prediction variance
pts <- list("sp.points", high, pch = 20, cex = 0.1, col = "black")
lim <-list("sp.polygons", brasil, col = "black")
dev.off()
pdf(file="res/fig/high-varia.pdf")
print(spplot(high.k, zcol = "var1.var", col.regions = rev(heat.colors(64)),
             main = "Variance of the predicted maximum Fe2O3 content",
             xlab = "Long", ylab = "Lat",
             sub = "Local kriging.", sp.layout = list(pts, lim),
             scales = list(draw = TRUE)))
dev.off();rm(pts, lim)

# save data
ls()
save(brasil, brasil.grid, data, fedata, ferric, high, high.k, high.v, high.vm,
     high.vmf, meta, perferric, proba, proba.ik,
     proba.iv, proba.ivm, proba.ivmf, purity, file = "wcss2014iron.RData")

# End!
