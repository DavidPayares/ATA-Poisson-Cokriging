#-------------------------------------------------------------------------------
#--------------------- ATA and ATP Poisson Cokriging ---------------------------
#-------------------------------------------------------------------------------


# ---------------------------------- 1 -----------------------------------------
# ------------------------------- Set up ---------------------------------------
# ------------------------------------------------------------------------------

# ----- Libraries
library(sf)
library(gstat)
library(units)
library(lattice)
library(mvtnorm)
library(tidyverse)
library(matrixcalc)

# ----- Set working directory
setwd('C:/Users/PayaresGarciaDE/Documents/Work/PhD/PCK - ATA')

# ----- Functions
source('./atakrig/R/ATA-ATP-Functions.R')

# ---------------------------------- 2 -----------------------------------------
# --------------------------- Data preparation ---------------------------------
# ------------------------------------------------------------------------------

# ----- Population Data

# --- Load geometries (Bogota)
upz = st_read('./data/upz-bogota.shp') %>% st_transform(crs = 3116)
# --- Plot geometry
plot(upz$geometry); axis(1); axis(2)

# --- Load population data
pop = st_read('./data/population-bogota.shp') %>% st_transform(crs = 3116)
# --- Assign block ID
pop <- pop %>% mutate(ID = row_number())

# --- create grid
grid.upz = st_make_grid(upz, cellsize = 500, crs = st_crs(upz))
# --- Assign grid ID
grid.upz = st_sf(grid.upz, 'GRID' = seq(length(grid.upz)), grid.upz)
# --- Intersect grid with Bogota's boundary
grid.upz = st_intersection(grid.upz, st_union(upz))
# --- Intersect grid with Bogota's population data
grid.pop = st_intersection(grid.upz, pop) %>% mutate(area = drop_units(st_area(grid.upz))/1000)
# --- Get area per per cell-block
grid.pop = grid.pop %>% group_by(ID) %>% mutate(percentage = area / sum(area))
# --- Get area per per cell-block
grid.pop$pop = grid.pop$percentage * grid.pop$Pop_Total
# --- Get population per per cell
grid.pop = grid.pop %>% group_by(GRID) %>% summarise(pop.cell = sum(pop)) %>% st_drop_geometry()
# --- Join population data with spatial grid
grid.pop = left_join(grid.upz, grid.pop, by = "GRID") %>% mutate(pop.cell = replace_na(pop.cell,0))
# --- Plot population data
plot(grid.pop["pop.cell"])

# ----- Linear Model of Coregionalization

# --- Size
size = nrow(grid.pop)

# --- Get grid coordinates
coordenadas = st_coordinates(st_centroid(grid.pop))
dist.matrix = as.matrix(dist(coordenadas))

# --- Define covariance structures
Cexp = function(h,a,b){b * exp(-h/a)}
Csph = function(h,a,b){ifelse(h <= a, b * (1-1.5*(h/a)+0.5*(h/a)^3), 0)}

# --- Define coefficient matrices (LMC)
B1 = c(1.9, 0.4, 0.14) %*% t(c(1.9,0.4,0.14))
B2= c(0.09, 2.3, 0.125) %*% t(c(0.09, 2.3, 0.125))
#B2= c(0, 0, 0) %*% t(c(0, 0, 0))
is.positive.semi.definite(B1)
is.positive.semi.definite(B2)

# --- Define processes covariance matrix
range.exp = 2000 
range.sph = 5000
cov11 = B1[1,1] * Cexp(dist.matrix,range.exp,1) + B2[1,1] * Csph(dist.matrix,range.sph,1)
cov22 = B1[2,2] * Cexp(dist.matrix,range.exp,1) + B2[2,2] * Csph(dist.matrix,range.sph,1)
cov33 = B1[3,3] * Cexp(dist.matrix,range.exp,1) + B2[3,3] * Csph(dist.matrix,range.sph,1)
cov12 = cov21 = B1[1,2] * Cexp(dist.matrix,range.exp,1) + B2[1,2] * Csph(dist.matrix,range.sph,1)
cov13 = cov31 = B1[1,3] * Cexp(dist.matrix,range.exp,1) + B2[1,3] * Csph(dist.matrix,range.sph,1)
cov23 = cov32 = B1[2,3] * Cexp(dist.matrix,range.exp,1) + B2[2,3] * Csph(dist.matrix,range.sph,1)
cov.total =rbind(cbind(cov11,cov12, cov13),cbind(cov21,cov22, cov23), cbind(cov31,cov32, cov33))
is.positive.semi.definite(cov.total)

# --- Plot LMC
dum = seq(1, 8000, length.out = 1000)
par(mfrow = c(3,3) ,
    oma = c(5,4,0,0) + 0.1,
    mar = c(2,1,1,1) + 0.1)
plot(dum, (B1[1,1] * Cexp(0,range.exp,1) + B2[1,1] * Csph(0,range.sph,1)) -(B1[1,1] * Cexp(dum,range.exp,1) + B2[1,1] * Csph(dum,range.sph,1)), type = 'l', xlab = 'h' , ylab = 'semivariance', main = expression(~ gamma[1]))
plot.new();plot.new();
plot(dum, (B1[1,2] * Cexp(0,range.exp,1) + B2[1,2] * Csph(0,range.sph,1)) -(B1[1,2] * Cexp(dum,range.exp,1) + B2[1,2] * Csph(dum,range.sph,1)), type = 'l', xlab = 'h' , ylab = 'semivariance', main = expression(~ gamma[12]))
plot(dum, (B1[2,2] * Cexp(0,range.exp,1) + B2[2,2] * Csph(0,range.sph,1)) -(B1[2,2] * Cexp(dum,range.exp,1) + B2[2,2] * Csph(dum,range.sph,1)), type = 'l', xlab = 'h' , ylab = 'semivariance',main = expression(~ gamma[2]))
plot.new();
plot(dum, (B1[1,3] * Cexp(0,range.exp,1) + B2[1,3] * Csph(0,range.sph,1)) -(B1[1,3] * Cexp(dum,range.exp,1) + B2[1,3] * Csph(dum,range.sph,1)), type = 'l', xlab = 'h' , ylab = 'semivariance', main = expression(~ gamma[13]))
plot(dum, (B1[2,3] * Cexp(0,range.exp,1) + B2[2,3] * Csph(0,range.sph,1)) -(B1[2,3] * Cexp(dum,range.exp,1) + B2[2,3] * Csph(dum,range.sph,1)), type = 'l', xlab = 'h' , ylab = 'semivariance', main = expression(~ gamma[23]))
plot(dum, (B1[3,3] * Cexp(0,range.exp,1) + B2[3,3] * Csph(0,range.sph,1)) -(B1[3,3] * Cexp(dum,range.exp,1) + B2[3,3] * Csph(dum,range.sph,1)), type = 'l', xlab = 'h' , ylab = 'semivariance', main = expression(~ gamma[3]))

# --- Generate multivariate gaussian field
seed = 1605
set.seed(seed)
sim = rmvnorm(1,mean=c(rep(11.5, nrow(coordenadas)),rep(49.16, nrow(coordenadas)),rep(3.8, nrow(coordenadas))), sigma=cov.total)
#sim = rmvnorm(1,c(rep(25.2,nrow(coordenadas)), rep(11.9,nrow(coordenadas)), rep(14.5,nrow(coordenadas))), sigma=cov.total)
datos = cbind(coordenadas,Ra = sim[1:size], Rb= sim[(size+1):(2*size)], Rc = sim[(2*size+1):(3*size)])
rho = cor(datos[,3:5])

# --- Convert the matrix into dataframe 
datos = as.data.frame(datos)

# --- Determine shared risk
datos$Rab <- rlk(datos$Ra, datos$Rb, rho[1,2])
datos$Rac <- rlk(datos$Ra, datos$Rc, rho[1,3])
datos$Rbc <- rlk(datos$Rb, datos$Rc, rho[2,3])
datos$GRID <- grid.pop$GRID

# --- Grid with data
grid.data = left_join(grid.pop, datos, by = "GRID")
plot(grid.data["Rc"])

# ----- Upscaling the data

# --- Intersect grid with Bogota's boundary and calculate area per UPZ-grid
grid.data.upz = st_intersection(grid.data, upz)
plot(grid.data.upz["Rc"])
# --- Get area per per cell-upz
grid.data.upz = grid.data.upz %>% mutate(area.cell = drop_units(st_area(grid.data.upz))/1000)
# --- Get percentage area per per cell-upz
grid.data.upz = grid.data.upz %>% group_by(GRID) %>% mutate(percentage.cell = area.cell / sum(area.cell))
# --- Get population per cell-upz
grid.data.upz$pop.cell.upz = grid.data.upz$percentage.cell * grid.data.upz$pop.cell
# --- Get population percentage per cell-upz
upz.data = grid.data.upz %>% group_by(COD_UPZ) %>% mutate(pop.upz = pop.cell.upz / sum(pop.cell.upz))
# --- Get risk per cell-upz
upz.data$Ra.cell.upz = upz.data$Ra * upz.data$pop.upz
upz.data$Rb.cell.upz = upz.data$Rb * upz.data$pop.upz
upz.data$Rc.cell.upz = upz.data$Rc * upz.data$pop.upz
upz.data$Rab.cell.upz = upz.data$Rab * upz.data$pop.upz
upz.data$Rac.cell.upz = upz.data$Rac * upz.data$pop.upz
upz.data$Rbc.cell.upz = upz.data$Rbc * upz.data$pop.upz
upz.data$x.cell.upz = st_coordinates(st_point_on_surface(upz.data))[,1]
upz.data$y.cell.upz = st_coordinates(st_point_on_surface(upz.data))[,2]
upz.data$areaId = as.integer(gsub("\\D", "", upz.data$COD_UPZ))
upz.data =  upz.data %>% arrange(areaId, x.cell.upz, y.cell.upz)

# --- Get population-weighted risk per UPZ
data.upz = upz.data %>% group_by(COD_UPZ) %>% summarise(Ra = sum(Ra.cell.upz),
                                                        Rb = sum(Rb.cell.upz),
                                                        Rc = sum(Rc.cell.upz),
                                                        Rab = sum(Rab.cell.upz),
                                                        Rac = sum(Rac.cell.upz),
                                                        Rbc = sum(Rbc.cell.upz),
                                                        N = int.pop(sum(pop.cell.upz)),
                                                        weightedpop(as.vector(x.cell.upz), 
                                                                    as.vector(y.cell.upz),
                                                                    as.vector(pop.cell.upz))) %>% st_drop_geometry()
# --- Assign spatial units to data
data.upz = left_join(upz, data.upz, by = "COD_UPZ") 
data.upz = na.omit(data.upz) 
plot(data.upz["Rc"])

# --- Generate count data
seed = seed
size = nrow(data.upz)
data.upz$Ya = rpois(size, data.upz$N * data.upz$Ra)
data.upz$Yb = rpois(size, data.upz$N * data.upz$Rb)
data.upz$Yc = rpois(size, data.upz$N * data.upz$Rc)
data.upz$Yab <- rpois(size, data.upz$N * data.upz$Rab)
data.upz$Yac <- rpois(size, data.upz$N * data.upz$Rac)
data.upz$Ybc <- rpois(size, data.upz$N * data.upz$Rbc)
data.upz$areaId <- as.integer(gsub("\\D", "", data.upz$COD_UPZ))
data.upz =  data.upz %>% arrange(areaId)


cor(data.upz$Ra, data.upz$Ya/data.upz$N)
cor(data.upz$Rb, data.upz$Yb/data.upz$N)
cor(data.upz$Rc, data.upz$Yc/data.upz$N)


# ----- Splitting data
set.seed(seed)
sample <- sample.int(n = nrow(data.upz), size = floor(0.7 * nrow(data.upz)), replace = F)
sample.2 <- sample.int(n = nrow(data.upz), size = floor(0.5 * nrow(data.upz)), replace = F)
data.upz.a <- data.upz[sample.2,] %>% arrange(areaId)
data.upz.b <- data.upz[sample,] %>% arrange(areaId)
unknown.a <- data.upz[-sample.2,] %>% arrange(areaId)
data.upz.ab =  data.upz[sample[sample %in% sample.2],] %>% arrange(areaId)

# ------ Plot data
plot(data.upz.a["Ra"])
plot(data.upz.b["Rb"])
plot(data.upz["Rc"])



# ----- Structuring the data

# --- Areal data
risk.a = data.frame(areaId = data.upz.a$areaId, x = data.upz.a$x , y = data.upz.a$y,  ya = data.upz.a$Ya, n = data.upz.a$N) %>% arrange(areaId)
risk.b = data.frame(areaId = data.upz.b$areaId, x = data.upz.b$x , y = data.upz.b$y,  ya = data.upz.b$Yb, n = data.upz.b$N) %>% arrange(areaId)
risk.c = data.frame(areaId = data.upz$areaId, x = data.upz$x , y = data.upz$y,  ya = data.upz$Yc, n = data.upz$N) %>% arrange(areaId)
unknown.risk.a = data.frame(areaId = unknown.a$areaId, x = unknown.a$x , y = unknown.a$y,  ya = unknown.a$Ya, n = unknown.a$N) %>% arrange(areaId)

# --- Population-weighted grid
discretize.points = data.frame(areaId = as.integer(gsub("\\D", "", upz.data$COD_UPZ)), ptx =  upz.data$x.cell.upz,
                               pty = upz.data$y.cell.upz, weight = upz.data$pop.cell.upz) %>% arrange(areaId)

# --- Structure dataset
createDataset <- function(data.test, discretePoints){
  colnames(data.test)[2:5] <- c("centx","centy","counts", "size")
  rslt <- list(areaValues = data.test, discretePoints = discretePoints)
  class(rslt) <- c("list", "discreteArea")
  return(rslt)
}

# --- Final data
data.a <- createDataset(risk.a, discretize.points[discretize.points$areaId %in% risk.a$areaId,])
data.b <- createDataset(risk.b, discretize.points[discretize.points$areaId %in% risk.b$areaId,])
data.c <- createDataset(risk.c, discretize.points)

# --- Validation data
data.unknown.a <- createDataset(unknown.risk.a, discretize.points[discretize.points$areaId %in% unknown.risk.a$areaId,])


# --- Co-infection data
data.co = data.upz %>% na.omit() %>% mutate(areaId = as.integer(gsub("\\D", "", data.upz$COD_UPZ))) %>% dplyr::select(c("areaId","Ya", "Yb","Yc","Yab","Yac", "Ybc", "N")) %>% arrange(areaId) %>% st_drop_geometry()
data.co[!(data.co$areaId %in% data.upz.a$areaId),c("Ya", "Yab", "Yac")] <- NA
data.co[!(data.co$areaId %in% data.upz.b$areaId),c("Yb", "Yab", "Ybc")] <- NA

# --- Means for semivariogram computing
means.co <- c(mean(data.upz$Rab), mean(data.upz$Rac), mean(data.upz$Rbc))

# ---------------------------------- 3 -----------------------------------------
# ------------------------------- ATA Kriging ----------------------------------
# ------------------------------------------------------------------------------

# ----- Poisson semivariogram deconvolution
vg.deconv.a <- deconvPointVgm(data.a, model="Exp", ngroup=15, rd = 0.23, fig = T, fit.method = 2)
vg.deconv.b <- deconvPointVgm(data.b, model="Exp", ngroup=15, rd = 0.23, fig = T, fit.method = 2)
vg.deconv.c <- deconvPointVgm(data.c, model="Exp", ngroup=15, rd = 0.23, fig = T, fit.method = 2)
plotDeconvVgm(vg.deconv.a)
plotDeconvVgm(vg.deconv.b)
plotDeconvVgm(vg.deconv.c)

# --- ATA Poisson kriging 
pred.ata.pk.a <- ataKriging(data.a, unknown = data.unknown.a , vg.deconv.a, nmax = Inf,  showProgress = TRUE)
cor(unknown.a$Ra, pred.ata.pk.a$pred)


# --- ATA Poisson kriging (no observation)
pred.ata.a <- ata.PKriging(data.a, vg.deconv.a, nmax = Inf, showProgress = T)
cor(data.upz.a$Ra, pred.ata.a$pred)

# ---- Predictions
preds.ata <- rbind(pred.ata.pk.a, pred.ata.a[,1:5]) %>% arrange(areaId)
cor(data.upz$Ra, preds.ata$pred)

# --- Map the results
data.upz$Ra.pred = preds.ata$pred
data.upz$Ra.var = preds.ata$var
plot(data.upz["Ra"])
plot(data.upz["Ra.pred"])
plot(data.upz["Ra.var"])

mean((data.upz$Ra - preds.ata$pred)^2)

# ----------------------------------  4 ----------------------------------------
# ------------------------------- ATP Kriging ----------------------------------
# ------------------------------------------------------------------------------
# --- Define discrete points to predict
points = data.frame(x = upz.data$x.cell.upz, y = upz.data$y.cell.upz)

# --- ATP Poisson kriging
pred.atp.ok.a <- atpKriging(data.a, unknown = points , vg.deconv.a, nmax = Inf,  showProgress = TRUE)
cor(upz.data$Ra, pred.atp.ok.a$pred)

# --- ATP Poisson kriging (no observation)
pred.atp.a <- atp.PKriging(data.a, vg.deconv.a, nmax = Inf, showProgress = T)
cor(upz.data[upz.data$areaId %in% data.upz.a$areaId,]$Ra, pred.atp.a$pred)

# --- Map the results
upz.data$Ra.pred = pred.atp.ok.a$pred
upz.data$Ra.var = pred.atp.ok.a$var
plot(upz.data["Ra"])
plot(upz.data["Ra.pred"])
plot(upz.data["Ra.var"])


# ---------------------------------- 5 -----------------------------------------
# ------------------------------- ATA CoKriging --------------------------------
# ------------------------------------------------------------------------------

data.all <- list(data.a = data.a, data.b = data.b, data.c = data.c)


# # --- ATP Poisson kriging
# data.all$data.a$areaValues <- data.all$data.a$areaValues[data.all$data.a$areaValues$areaId %in% data.upz.ab$areaId,]
# data.all$data.b$areaValues <- data.all$data.b$areaValues[data.all$data.b$areaValues$areaId %in% data.upz.ab$areaId,]
# data.all$data.c$areaValues <- data.all$data.c$areaValues[data.all$data.c$areaValues$areaId %in% data.upz.ab$areaId,]
# data.all$data.a$discretePoints <- data.all$data.a$discretePoints[data.all$data.a$discretePoints$areaId %in% data.all$data.a$areaValues$areaId,]
# data.all$data.b$discretePoints <- data.all$data.b$discretePoints[data.all$data.b$discretePoints$areaId %in% data.all$data.a$areaValues$areaId,]
# data.all$data.c$discretePoints <- data.all$data.c$discretePoints[data.all$data.c$discretePoints$areaId %in% data.all$data.a$areaValues$areaId,]


# ----- Poisson cross-semivariogram deconvolution
vg.deconv.cok <- deconvPointVgmForCoKriging(data.all, means.co, model="Exp", fig = T, fit.method = 2, rd = 0.40)

replaceGamma <- function(variogram) {
  
  regularized_data = variogram$regularizedAreaVariogram
  experiential_data = variogram$experientialAreaVariogram
  
  # Ensure both data frames have the same number of rows
  if (nrow(regularized_data) != nrow(experiential_data)) {
    stop("The two data frames must have the same number of rows")
  }
  
  # Replace gamma values in the experiential data with those from the regularized data
  experiential_data$gamma <- regularized_data$gamma
  
  # Return the modified experiential data
  class(experiential_data) <- c("gstatVariogram", "data.frame")
  return(experiential_data)
}

var.a = replaceGamma(vg.deconv.cok$data.a)
var.b = replaceGamma(vg.deconv.cok$data.b)
var.c = replaceGamma(vg.deconv.cok$data.c)
var.ab = replaceGamma(vg.deconv.cok$data.a.data.b)
var.ac = replaceGamma(vg.deconv.cok$data.a.data.c)
var.bc = replaceGamma(vg.deconv.cok$data.b.data.c)


#--- Define initial values
lmc.init <- list( v1 = vgm(B2[1,1], "Sph", range.sph, add.to = vgm(B1[1,1], "Exp", range.exp)),
                  v2 = vgm(B2[2,2], "Sph", range.sph, add.to = vgm(B1[2,2], "Exp", range.exp)),
                  v3 = vgm(B2[3,3], "Sph", range.sph, add.to = vgm(B1[3,3], "Exp", range.exp)),
                  v12 = vgm(B2[1,2], "Sph", range.sph, add.to = vgm(B1[1,2], "Exp", range.exp)),
                  v13 = vgm(B2[1,3], "Sph", range.sph, add.to = vgm(B1[1,3], "Exp", range.exp)),
                  v23 = vgm(B2[2,3], "Sph", range.sph, add.to = vgm(B1[2,3], "Exp", range.exp)))

lmc.poisson.cokrige <- function(var.a, var.b, var.c, crossvar.ab, crossvar.ac, crossvar.bc, data, var.params){
  
  # Change id auxilliary variable
  var.b$id = as.factor("var2")
  var.c$id = as.factor("var3")
  crossvar.ab$id =  as.factor("var1.var2")
  crossvar.ac$id =  as.factor("var1.var3")
  crossvar.bc$id =  as.factor("var2.var3")
  
  names = c(expression(~ gamma[13]), expression(~ gamma[23]), expression(~ gamma[3]), expression(~ gamma[12]), expression(~ gamma[2]), expression(~ gamma[1]))
  my.settings <- list(
    strip.background=list(col="transparent"),
    strip.border=list(col="transparent")
  )
  
  # Integrate variograms
  variograms <- rbind(crossvar.ac, crossvar.bc, var.c, crossvar.ab, var.b, var.a)
  plot(variograms)
  
  # Gstat object
  g <- gstat(id = "var1", formula = c/p ~ 1, loc = ~ x + y, data = data.frame(x = data$data.a$areaValues$centx ,y = data.all$data.a$areaValues$centy, c = data$data.a$areaValues$counts, p = data$data.a$areaValues$size), model = var.params$v1)
  g <- gstat(g,id = "var2", formula = c/p ~ 1, loc = ~ x + y, data = data.frame(x = data$data.b$areaValues$centx ,y = data.all$data.b$areaValues$centy, c = data$data.b$areaValues$counts, p = data$data.b$areaValues$size), model = var.params$v2)
  g <- gstat(g,id = "var3", formula = c/p ~ 1, loc = ~ x + y, data = data.frame(x = data$data.c$areaValues$centx ,y = data.all$data.c$areaValues$centy, c = data$data.c$areaValues$counts, p = data$data.c$areaValues$size), model = var.params$v3)
  g <- gstat(g, id = c("var1","var2"), model = var.params$v12)
  g <- gstat(g, id = c("var1","var3"), model = var.params$v13)
  g <- gstat(g, id = c("var2","var3"), model = var.params$v23)
  fitted.lmc <- fit.lmc(v = variograms, g = g, fit.method = 2, fit.ranges = T)
  print(plot(variograms, fitted.lmc, pch = 20, col = "black", ylab = 'semivariance', par.settings =  my.settings, strip=strip.custom(factor.levels=names)))
  return(fitted.lmc)
}





#--- Linear model of Coregionalization
fitted.lmc <- lmc.poisson.cokrige(var.a, var.b, var.c, var.ab, var.ac, var.bc, data.all, lmc.init)


vg.deconv.cok$data.a$pointVariogram <- fitted.lmc$model$var1
vg.deconv.cok$data.b$pointVariogram <- fitted.lmc$model$var2
vg.deconv.cok$data.c$pointVariogram <- fitted.lmc$model$var3
vg.deconv.cok$data.a.data.b$pointVariogram <- fitted.lmc$model$var1.var2
vg.deconv.cok$data.a.data.c$pointVariogram <- fitted.lmc$model$var1.var3
vg.deconv.cok$data.b.data.c$pointVariogram <- fitted.lmc$model$var2.var3


data.all <- list(data.a = data.a, data.b = data.b, data.c = data.c)

# --- Means for semivariogram computing
means.data <- c(mean(data.upz$Ra), mean(data.upz$Rb), mean(data.upz$Rc), means.co)


#------ ATA Poisson cokriging
pred.atacok.a <- ataCoKriging(data.all, unknownVarId="data.a", unknown = data.unknown.a, comorbidity =  means.data, 
                              ptVgms=vg.deconv.cok, nmax = Inf, oneCondition=FALSE, showProgress = TRUE)
cor(unknown.a$Ra, pred.atacok.a$pred)


pred.atapok.a <-atpCoKriging(data.all, unknownVarId="data.a", unknown = points, comorbidity =  means.data, 
                             ptVgms=vg.deconv.cok, nmax = Inf, oneCondition=FALSE, showProgress = TRUE)
cor(upz.data$Ra, pred.atapok.a$pred)


save.image("./data/simulation-irregular-het.RData")
load("./data/simulation-irregular-het.RData")



# ----- Libraries
library(sf)
library(gstat)
library(units)
library(mvtnorm)
library(tidyverse)
library(matrixcalc)

# ----- Set working directory
setwd('C:/Users/PayaresGarciaDE/Documents/Work/PhD/PCK - ATA')

# ----- Functions
source('./atakrig/R/ATA-ATP-Functions.R')
