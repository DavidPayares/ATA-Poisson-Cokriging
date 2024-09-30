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
library(viridis)
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

# --- Intersect grid with Bogota's boundary and calculate area per UPZ-grid
grid.data.upz = st_intersection(grid.pop, upz)
# --- Get area per per cell-upz
grid.data.upz = grid.data.upz %>% mutate(area.cell = drop_units(st_area(grid.data.upz))/1000)
# --- Get percentage area per per cell-upz
grid.data.upz = grid.data.upz %>% group_by(GRID) %>% mutate(percentage.cell = area.cell / sum(area.cell))
# --- Get population per cell-upz
grid.data.upz$pop.cell.upz = grid.data.upz$percentage.cell * grid.data.upz$pop.cell
# --- Get population percentage per cell-upz
upz.data = grid.data.upz %>% group_by(COD_UPZ) %>% mutate(pop.upz = pop.cell.upz / sum(pop.cell.upz))
# --- Get risk per cell-upz
upz.data$x.cell.upz = st_coordinates(st_point_on_surface(upz.data))[,1]
upz.data$y.cell.upz = st_coordinates(st_point_on_surface(upz.data))[,2]
upz.data$areaId = as.integer(gsub("\\D", "", upz.data$COD_UPZ))
upz.data =  upz.data %>% arrange(areaId, x.cell.upz, y.cell.upz)

# --- Get population-weighted risk per UPZ
data.upz = upz.data %>% group_by(COD_UPZ) %>% summarise(N = int.pop(sum(pop.cell.upz)),
                                                        weightedpop(as.vector(x.cell.upz), 
                                                                    as.vector(y.cell.upz),
                                                                    as.vector(pop.cell.upz))) %>% st_drop_geometry()
# --- Assign spatial units to data
data.upz = left_join(upz, data.upz, by = "COD_UPZ") 
data.upz = na.omit(data.upz)
data.upz$areaId = as.integer(gsub("\\D", "", data.upz$COD_UPZ))

# ----- Cases Data
covid.data = read.csv('./data/Bogota_data.csv', sep = ",")

# Join data
covid.sp = merge(data.upz, covid.data, by.x = "COD_UPZ", by.y = "UPZ") %>% arrange(areaId)

# Compute rates
pop.rate = 1000
covid.sp$Rcovid = covid.sp$Covid / covid.sp$N * pop.rate
covid.sp$RAsthma = covid.sp$Asthma / covid.sp$N * pop.rate
covid.sp$RDiabetes = covid.sp$Diabetes / covid.sp$N * pop.rate

# --- Plot data
plot.data2(grid.pop, grid.pop$pop.cell, col.pal = 'viridis', col.line = NA) + geom_sf(data = upz, fill = NA, color = "white", size = 0.5)   # Add the border polygon layer
plot.data2(covid.sp, covid.sp$Rcovid)
plot.data2(covid.sp, covid.sp$RAsthma)
plot.data2(covid.sp, covid.sp$RDiabetes)


# ----- Structuring the data

# --- Areal data
risk.a = data.frame(areaId = covid.sp$areaId, x = covid.sp$x , y = covid.sp$y,  ya = covid.sp$Covid * pop.rate, n = covid.sp$N) %>% arrange(areaId)
risk.b = data.frame(areaId = covid.sp$areaId, x = covid.sp$x , y = data.upz$y,  ya = covid.sp$Asthma * pop.rate, n = covid.sp$N) %>% arrange(areaId)
risk.c = data.frame(areaId = covid.sp$areaId, x = covid.sp$x , y = data.upz$y,  ya = covid.sp$Diabetes * pop.rate, n = covid.sp$N) %>% arrange(areaId)

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
data.a <- createDataset(risk.a, discretize.points)
data.b <- createDataset(risk.b, discretize.points)
#data.c <- createDataset(risk.c, discretize.points)
#data.all <- list(data.a = data.a, data.b = data.b, data.c = data.c)
data.all <- list(data.a = data.a, data.b = data.b)

# --- Co-morbidity estimation
rho =  cor(data.frame(c = covid.sp$Rcovid, a = covid.sp$RAsthma, d = covid.sp$RDiabetes))

covid.sp$Rca <- rlk(covid.sp$Rcovid, covid.sp$RAsthma, rho[1,2])
covid.sp$Rcd <- rlk(covid.sp$Rcovid, covid.sp$RDiabetes, rho[1,3])
covid.sp$Rad <- rlk(covid.sp$RAsthma, covid.sp$RDiabetes, rho[2,3])

#means.co =  c(mean(covid.sp$Rca), mean(covid.sp$Rcd), mean(covid.sp$Rad))
means.co =  c(mean(covid.sp$Rca))

# ---------------------------------- 3 -----------------------------------------
# ------------------------------- ATA Kriging ----------------------------------
# ------------------------------------------------------------------------------

# ----- Poisson semivariogram deconvolution
vg.deconv.a <- deconvPointVgm(data.a, model= "Exp", ngroup=15, rd = 0.50, fig = T, fit.method = 2)
plotDeconvVgm(vg.deconv.a)

# --- ATA Poisson kriging 
pred.ata.pk.a <- ataKriging(data.a, unknown = data.a , vg.deconv.a, nmax = Inf,  showProgress = TRUE)

# --- ATA Poisson kriging (no observation)
pred.ata.a <- ata.PKriging(data.a, vg.deconv.a, nmax = Inf, showProgress = T)

covid.sp$R.PKATA <- pred.ata.a$pred
covid.sp$var.PKATA <- pred.ata.a$var/100*1.1

# Plot
plot.data2(covid.sp, covid.sp$R.PKATA)
plot.data2(covid.sp, covid.sp$var.PKATA, col.pal = "mako")

# Check correlation
cor(pred.ata.a$pred, covid.sp$Rcovid)

# ----------------------------------  4 ----------------------------------------
# ------------------------------- ATP Kriging ----------------------------------
# ------------------------------------------------------------------------------

# --- Define discrete points to predict
grid.pred = st_make_grid(upz, cellsize = 500, crs = st_crs(upz))
grid.pred = st_sf(grid.pred, 'GRID' = seq(length(grid.pred)), grid.pred)
grid.pred = st_intersection(grid.pred, upz)
centorids = st_coordinates(st_centroid(grid.pred))
points = data.frame(x = centorids[,1], y = centorids[,2])

data.a$areaValues$counts <- pred.ata.a$pred * data.a$areaValues$size

# --- ATP Poisson kriging
pred.atp.ok.a <- atpKriging(data.a, unknown = points ,vg.deconv.a , nmax = Inf,  showProgress = TRUE)

# --- Map the results
grid.pred$Ra.pred = pred.atp.ok.a$pred
grid.pred$Ra.var = pred.atp.ok.a$var/100

# Plot
plot.data2(grid.pred, grid.pred$Ra.pred, col.line = NA)
plot.data2(grid.pred, grid.pred$Ra.var, col.pal = "mako", col.line = NA)


# ---------------------------------- 5 -----------------------------------------
# ------------------------------- ATA CoKriging --------------------------------
# ------------------------------------------------------------------------------

#mean.data = c(mean(covid.sp$Rcovid), mean(covid.sp$RAsthma), mean(covid.sp$RDiabetes), means.co)
mean.data = c(mean(covid.sp$Rcovid), mean(covid.sp$RAsthma), means.co)

# ----- Poisson cross-semivariogram deconvolution
vg.deconv.cok <- deconvPointVgmForCoKriging(data.all, means.co, model=c("Exp"), ngroup=12, rd = 0.30, fig = T, fit.method = 2)
plot(vg.deconv.cok)

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
#var.c = replaceGamma(vg.deconv.cok$data.c)
var.ab = replaceGamma(vg.deconv.cok$data.a.data.b)
#var.ac = replaceGamma(vg.deconv.cok$data.a.data.c)
#var.bc = replaceGamma(vg.deconv.cok$data.b.data.c)


#--- Define initial values
lmc.init <- list( v1 = vgm(1500, "Exp", 2300),
                  v2 = vgm(600, "Exp", 2300),
                  #v3 = vgm(1000, "Exp", 2300),
                  v12 = vgm(400, "Exp", 2300))
                  #v13 = vgm(800, "Exp", 2300),
                  #v23 = vgm(600, "Exp", 2300))

lmc.poisson.cokrige <- function(var.a, var.b, crossvar.ab, data, var.params){
  
  
  # Change id auxilliary variable
  var.b$id = as.factor("var2")
  #var.c$id = as.factor("var3")
  crossvar.ab$id =  as.factor("var1.var2")
  #crossvar.ac$id =  as.factor("var1.var3")
  #crossvar.bc$id =  as.factor("var2.var3")
  
  names = c(expression(~ gamma[12]), expression(~ gamma[2]), expression(~ gamma[1]))
  my.settings <- list(
    strip.background=list(col="transparent"),
    strip.border=list(col="transparent")
  )
  
  # Integrate variograms
  variograms <- rbind(crossvar.ab, var.b, var.a)
  plot(variograms)
  
  # Gstat object
  g <- gstat(id = "var1", formula = c/p ~ 1, loc = ~ x + y, data = data.frame(x = data$data.a$areaValues$centx ,y = data.all$data.a$areaValues$centy, c = data$data.a$areaValues$counts, p = data$data.a$areaValues$size), model = var.params$v1)
  g <- gstat(g,id = "var2", formula = c/p ~ 1, loc = ~ x + y, data = data.frame(x = data$data.a$areaValues$centx ,y = data.all$data.b$areaValues$centy, c = data$data.b$areaValues$counts, p = data$data.b$areaValues$size), model = var.params$v2)
  #g <- gstat(g,id = "var3", formula = c/p ~ 1, loc = ~ x + y, data = data.frame(x = data$data.a$areaValues$centx ,y = data.all$data.c$areaValues$centy, c = data$data.c$areaValues$counts, p = data$data.c$areaValues$size), model = var.params$v3)
  g <- gstat(g, id = c("var1","var2"), model = var.params$v12)
  #g <- gstat(g, id = c("var1","var3"), model = var.params$v13)
  #g <- gstat(g, id = c("var2","var3"), model = var.params$v23)
  fitted.lmc <- fit.lmc(v = variograms, g = g, fit.method = 2, fit.ranges = T)
  print(plot(variograms, fitted.lmc, pch = 20, col = "black", ylab = 'semivariance', par.settings =  my.settings, strip=strip.custom(factor.levels=names)))
  return(fitted.lmc)
}



#--- Linear model of Coregionalization
fitted.lmc <- lmc.poisson.cokrige(var.a, var.b, var.ab, data.all, lmc.init)


vg.deconv.cok$data.a$pointVariogram <- fitted.lmc$model$var1
vg.deconv.cok$data.b$pointVariogram <- fitted.lmc$model$var2
#vg.deconv.cok$data.c$pointVariogram <- fitted.lmc$model$var3
vg.deconv.cok$data.a.data.b$pointVariogram <- fitted.lmc$model$var1.var2
#vg.deconv.cok$data.a.data.c$pointVariogram <- fitted.lmc$model$var1.var3
#vg.deconv.cok$data.b.data.c$pointVariogram <- fitted.lmc$model$var2.var3


#data.all$data.a$areaValues$counts <- pred.ata.a$pred * data.a$areaValues$size
data.all <- list(data.a = data.a, data.b = data.b)


#------ ATA Poisson cokriging
pred.atacok.a <- ataCoKriging(data.all, unknownVarId="data.a", unknown=data.a, comorbidity =  mean.data,
                              ptVgms=vg.deconv.cok, nmax = Inf, oneCondition=FALSE, showProgress = TRUE)


pred.atacok.as <- ataCoKriging.cv(data.all, unknownVarId="data.a", unknown=data.a, comorbidity =  mean.data,
                              ptVgms=vg.deconv.cok, nmax = Inf, oneCondition=FALSE, showProgress = TRUE)



covid.sp$R.PCKATA <- pred.atacok.as$pred
covid.sp$var.PCKATA <- pred.atacok.as$var/150

# Plot
plot.data2(covid.sp, covid.sp$R.PCKATA)
plot.data2(covid.sp, covid.sp$var.PCKATA, col.pal = "mako")


# Check correlation
cor(covid.sp$R.PCKATA, covid.sp$Rcovid)


data.all$data.a$areaValues$counts <- ceiling(covid.sp$R.PCKATA  * data.a$areaValues$size)

pred.atacop.a <- atpCoKriging(data.all, unknownVarId="data.a", unknown=points, comorbidity =  mean.data,
                              ptVgms=vg.deconv.cok, nmax = Inf, oneCondition=FALSE, showProgress = TRUE)


# --- Map the results
grid.pred$Ra.pred.co = pred.atacop.a$pred
grid.pred$Ra.var.co = pred.atacop.a$var/100

# Plot
plot.data2(grid.pred, grid.pred$Ra.pred.co, col.line = NA)
plot.data2(grid.pred, grid.pred$Ra.var, col.pal = "mako", col.line = NA)


#-----------------  Simulate two more diseases

# unconditional simulation on a 100 x 100 grid
# Sample data frame
library(sp)
df <- data.frame(
  x = data.a$areaValues$centx,
  y = data.a$areaValues$centy,
  var1 = covid.sp$Rcovid
)
coordinates(df) <- ~ x + y
g.dummy <- gstat(formula = var1~1, dummy =  F, beta = 0,
                 model = vg.deconv.a$areaVariogram, data = df)
g.dummy <- gstat(g.dummy, formula = var2~1, dummy = TRUE, beta = 0,
                 model = vgm(250, "Exp", 3000))
g.dummy <- gstat(g.dummy, formula = var3~1, dummy = TRUE, beta = 0,
                 model = vgm(1566, "Exp", 2500) , fill.cross = T)

# for speed -- 10 is too small!!
yy <- predict(g.dummy, df, nsim = 1000)

for (i in 1001:ncol(yy)){
  cor <- cor(data.frame(yy[,1270])[,1], df$var1)
  if(cor > 0.30){
    print(paste0(i," - ", cor ))}
}

yy[,1]


#--------------- Asthma

Ashthma = (yy[,1270]$var2.sim270 - (min(yy[,1270]$var2.sim270)) + 12.5)/3
covid.sp$asthma <- ceiling(Ashthma *  covid.sp$N/ pop.rate)
covid.sp$Rasthma = covid.sp$asthma / covid.sp$N * pop.rate

cor(covid.sp$Rasthma, covid.sp$Rcovid)


#--------------- Diabetes


diabetes = (yy[,2270]$var3.sim270 - (min(yy[,2270]$var3.sim270))) /2
covid.sp$diabetes <- ceiling(diabetes *  covid.sp$N/ pop.rate)
covid.sp$Rdiabetes = covid.sp$diabetes / covid.sp$N * pop.rate

cor(covid.sp$Rdiabetes, covid.sp$Rcovid)


load("./data/Bogota.RData")

data.bogota = data.frame(UPZ = covid.sp$COD_UPZ, Covid = covid.sp$Covid, Asthma = covid.sp$asthma, Diabetes = covid.sp$diabetes)
write.csv(data.bogota, file = "./data/bogota_data.csv", sep = ";")



# ----------------- Functions


plot.data2 <- function(data, field, col.pal = "rocket", col.line = "white"){
  min.var = min(c(field))
  max.var = max(c(field))
  ggplot(data = data) +
    geom_sf(aes(fill = field),  color = col.line, size = 0.3) + 
    scale_fill_viridis_c(option = col.pal, 
                         direction = 1,
                         limits = c(min.var, max.var)) +
    theme(
      line = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.background = element_blank(),
      #legend.position = c(0.75, 0.15),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 10),
      legend.position = "right",
      legend.spacing.x = unit(0.1, 'cm'),
      legend.key.width = unit(1,"line"),
      legend.key.height = unit(5.0,"line"),
      legend.title.align = 0.5,
      legend.text.align = 0.5,
      legend.justification = "center"
    ) +
    coord_sf(datum = NA) +
    guides(fill = guide_colorbar(title = NULL))
}
