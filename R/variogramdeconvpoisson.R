#-------------------------------------------------------------------------------
#--------------------- ATA and ATP Poisson Cokriging ---------------------------
#-------------------------------------------------------------------------------


# ---------------------------------- 1 -----------------------------------------
# ------------------------------- Set up ---------------------------------------
# ------------------------------------------------------------------------------

#----- Libraries
library(sf)
library(sp)
library(geoR)
library(Rcpp)
library(raster)
library(gstat)
library(proxy)
library(mvtnorm)
library(viridis)
library(gridExtra)
library(matrixcalc)
library(latticeExtra)
library(RColorBrewer)


#----- Functions
source('C:/Users/PayaresGarciaDE/Documents/Work/PhD/PCK - ATA/atakrig/R/ATA-ATP-Functions.R')

# ---------------------------------- 2 -----------------------------------------
# --------------------------- Data preparation ---------------------------------
# ------------------------------------------------------------------------------

# ----- Linear Model of Coregionalization

# --- Define region
size = 25*25
x = seq(-10,10,length.out = sqrt(size))
y = seq(-10,10,length.out = sqrt(size))
coordenadas = expand.grid(x = x, y = y)
dist.matrix = as.matrix(dist(coordenadas))

# --- Define covariance structures
Cexp = function(h,a,b){b * exp(-h/a)}
Csph = function(h,a,b){ifelse(h <= a, b * (1-1.5*(h/a)+0.5*(h/a)^3), 0)}

# --- Define coefficient matrices (LMC)
B1 = c(1.9, 0.4, 0.14) %*% t(c(1.9,0.4,0.14))
B2= c(0.09, 2.3, 0.125) %*% t(c(0.09, 2.3, 0.125))
is.positive.semi.definite(B1)
is.positive.semi.definite(B2)

# --- Define processes covariance matrix
cov11 = B1[1,1] * Cexp(dist.matrix,1,1) + B2[1,1] * Csph(dist.matrix,5,1)
cov22 = B1[2,2] * Cexp(dist.matrix,1,1) + B2[2,2] * Csph(dist.matrix,5,1)
cov33 = B1[3,3] * Cexp(dist.matrix,1,1) + B2[3,3] * Csph(dist.matrix,5,1)
cov12 = cov21 = B1[1,2] * Cexp(dist.matrix,1,1) + B2[1,2] * Csph(dist.matrix,5,1)
cov13 = cov31 = B1[1,3] * Cexp(dist.matrix,1,1) + B2[1,3] * Csph(dist.matrix,5,1)
cov23 = cov32 = B1[2,3] * Cexp(dist.matrix,1,1) + B2[2,3] * Csph(dist.matrix,5,1)
cov.total =rbind(cbind(cov11,cov12, cov13),cbind(cov21,cov22, cov23), cbind(cov31,cov32, cov33))
is.positive.semi.definite(cov.total)

# --- Plot LMC
dum = seq(1, 10, length.out = 1000)
par(mfrow = c(3,3) ,
    oma = c(5,4,0,0) + 0.1,
    mar = c(2,1,1,1) + 0.1)
plot(dum, (B1[1,1] * Cexp(0,1,1) + B2[1,1] * Csph(0,5,1)) -(B1[1,1] * Cexp(dum,1,1) + B2[1,1] * Csph(dum,5,1)), type = 'l', xlab = 'h' , ylab = 'semivariance', main = expression(~ gamma[1]))
plot.new();plot.new();
plot(dum, (B1[1,2] * Cexp(0,1,1) + B2[1,2] * Csph(0,5,1)) -(B1[1,2] * Cexp(dum,1,1) + B2[1,2] * Csph(dum,5,1)), type = 'l', xlab = 'h' , ylab = 'semivariance', main = expression(~ gamma[12]))
plot(dum, (B1[2,2] * Cexp(0,1,1) + B2[2,2] * Csph(0,5,1)) -(B1[2,2] * Cexp(dum,1,1) + B2[2,2] * Csph(dum,5,1)), type = 'l', xlab = 'h' , ylab = 'semivariance',main = expression(~ gamma[2]))
plot.new();
plot(dum, (B1[1,3] * Cexp(0,1,1) + B2[1,3] * Csph(0,5,1)) -(B1[1,3] * Cexp(dum,1,1) + B2[1,3] * Csph(dum,5,1)), type = 'l', xlab = 'h' , ylab = 'semivariance', main = expression(~ gamma[13]))
plot(dum, (B1[2,3] * Cexp(0,1,1) + B2[2,3] * Csph(0,5,1)) -(B1[2,3] * Cexp(dum,1,1) + B2[2,3] * Csph(dum,5,1)), type = 'l', xlab = 'h' , ylab = 'semivariance', main = expression(~ gamma[23]))
plot(dum, (B1[3,3] * Cexp(0,1,1) + B2[3,3] * Csph(0,5,1)) -(B1[3,3] * Cexp(dum,1,1) + B2[3,3] * Csph(dum,5,1)), type = 'l', xlab = 'h' , ylab = 'semivariance', main = expression(~ gamma[3]))

# --- Generate multivariate gaussian field
seed = 512
set.seed(seed)
sim = rmvnorm(1,c(rep(49.86081,nrow(coordenadas)), rep(11.9,nrow(coordenadas)), rep(5.8,nrow(coordenadas))), sigma=cov.total)
datos = cbind(coordenadas,Ra = sim[1:size], Rb= sim[(size+1):(2*size)], Rc = sim[(2*size+1):(3*size)])
cor(datos[,3:5])

# --- Determine shared risk
datos$Rab <- rlk(datos$Ra, datos$Rb, 0.2572842)
datos$Rac <- rlk(datos$Ra, datos$Rc, 0.7838106)
datos$Rbc <- rlk(datos$Rb, datos$Rc, 0.8017566)

# --- Create dataframe
R = SpatialPixelsDataFrame(points = coordenadas, data = datos[,3:8])

# --- Plot Dataframe
risks <- names(R)
plots <- lapply(risks, function(x){spplot(R, x , col.regions = rev(brewer.pal(n = 11, name = "Spectral")), cuts = 10, par.settings = list(strip.background = list(col =  'transparent')), colorkey = list(space = "right", height = 1), main = x)})
grid.arrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]],plots[[5]],plots[[6]], ncol = 3)

# --- Upscale data
names.risk <- names
R.raster <- raster(R["Rb"])
R.upscale <- aggregate(R.raster, fact = c(2, 2), fun = mean)
R.values <- values(R.upscale)
R.centroid <- coordinates(R.upscale)

ka = SpatialPixelsDataFrame(points = R.centroid, data = data.frame(Rb = R.values))

plot(R.raster)
plot(R.upscale)
plot(ka)

# --- Generate Poisson data
effort <- "poisson"
if(effort == "uniform"){
  n = floor(runif(size, 500,2000))
}else if(effort == "normal"){
  n = floor(rnorm(size, 1000,500))
  n = n - min(n) + 50
}else if(effort == "poisson"){
  n = floor(rpois(size, 2000))
}else if(effort == "beta"){
  n = floor(1000*rbeta(size,0.1,5))
  n = n + 50
}

# --- Generate count data
R$Ya = rpois(size, n * R$Ra)
R$Yb = rpois(size, n * R$Rb)
R$Yc = rpois(size, n * R$Rc)
R$Yab <- rpois(size, n * R$Rab)
R$Yac <- rpois(size, n * R$Rac)
R$Ybc <- rpois(size, n * R$Rbc)

# --- Add data to dataframe
R$N <- n

# --- Structure dataframe
data.test.a <- data.frame(areaId= seq(1,625,1),x = coordenadas$x , y = coordenadas$y,  ya = R$Ya, n = R$N)
data.test.b <- data.frame(areaId= seq(1,625,1),x = coordenadas$x , y = coordenadas$y,  ya = R$Yb, n = R$N)
data.test.c <- data.frame(areaId= seq(1,625,1),x = coordenadas$x , y = coordenadas$y,  ya = R$Yc, n = R$N)

# ----- Discretization points

# --- Parameters
cellsize = abs(x[1]-x[2])/2
blkResoX = abs(x[1]-x[2])
blkResoY = abs(y[1]-y[2])
xResoNum = round(blkResoX / cellsize)
yResoNum = round(blkResoY / cellsize)
w <- 1/(xResoNum*yResoNum)

# --- Add discrete points to dataframe
discretizePoints <- function(data.test){lapply(1:nrow(data.test), function(i) {
  xy <- expand.grid(list(ptx = data.test$x[i] - 0.5 * blkResoX + (1:xResoNum-0.5) * cellsize,
                         pty = data.test$y[i] - 0.5 * blkResoY + (1:yResoNum-0.5) * cellsize))
  data.frame(areaId = data.test$areaId[i], xy, weight=w)
})
}

# --- Create discretized points
discretePoints.a <- do.call(rbind,discretizePoints(data.test.a))
discretePoints.b <- do.call(rbind,discretizePoints(data.test.b))
discretePoints.c <- do.call(rbind,discretizePoints(data.test.c))

# --- Structure dataset

createDataset <- function(data.test, discretePoints){
  colnames(data.test)[2:5] <- c("centx","centy","counts", "size")
  rslt <- list(areaValues = data.test, discretePoints = discretePoints)
  class(rslt) <- c("list", "discreteArea")
  return(rslt)
}

# --- Final data
rslt.a <- createDataset(data.test.a, discretePoints.a)
rslt.b <- createDataset(data.test.b, discretePoints.b)
rslt.c <- createDataset(data.test.b, discretePoints.c)
rslt.combined <- list(data.a = rslt.a, data.b = rslt.b, data.c = rslt.c)

# --- Co-infection data
means.co <- c(sum(R$Yab) /  sum(R$N), sum(R$Yac) / sum(R$N), sum(R$Ybc) / sum((R$N)))

# ---------------------------------- 3 -----------------------------------------
# ------------------------------- ATA Kriging ----------------------------------
# ------------------------------------------------------------------------------

# ----- Poisson semivariogram deconvolution
vg.deconv <- deconvPointVgm(rslt.a, model="Exp", ngroup=15, rd = 0.25, maxSampleNum=400, fig = T)
plotDeconvVgm(vg.deconv)

#------ ATA Poisson kriging
pred.ataok <- ataKriging(rslt.a, unknown = rslt.a , vg.deconv, nmax = 32,  showProgress = TRUE)

# --- Map the results
R$Ra.pred = pred.ataok$pred
R$Ra.var = pred.ataok$var
plot.ra = spplot(R, 'Ra' , col.regions = rev(brewer.pal(n = 11, name = "Spectral")), cuts = 10, par.settings = list(strip.background = list(col =  'transparent')), colorkey = list(space = "right", height = 1))
plot.ra.pred = spplot(R, 'Ra.pred' , col.regions = rev(brewer.pal(n = 11, name = "Spectral")), cuts = 10, par.settings = list(strip.background = list(col =  'transparent')), colorkey = list(space = "right", height = 1))
plot.ra.var = spplot(R, 'Ra.var' , col.regions = rev(brewer.pal(n = 11, name = "PuOr")), cuts = 10, par.settings = list(strip.background = list(col =  'transparent')), colorkey = list(space = "right", height = 1))
grid.arrange(plot.ra, plot.ra.pred, plot.ra.var, ncol = 3)


# ----------------------------------  4 ----------------------------------------
# ------------------------------- ATP Kriging ----------------------------------
# ------------------------------------------------------------------------------

# --- Define discrete points to predict
sizep = 80*80
xp = seq(-10,10,length.out = sqrt(sizep))
yp = seq(-10,10,length.out = sqrt(sizep))
points = expand.grid(x = xp, y = yp)

# --- ATP Poisson kriging
ataStartCluster(2)
pred.atp.ok <- atpKriging(rslt.a, unknown = points , vg.deconv, nmax = 32,  showProgress = TRUE)

# --- Map the results
R.atp = SpatialPixelsDataFrame(points = points, data = pred.atp.ok[4:5])
plot.ra = spplot(R, 'Ra' , col.regions = rev(brewer.pal(n = 11, name = "Spectral")), cuts = 10, par.settings = list(strip.background = list(col =  'transparent')), colorkey = list(space = "right", height = 1))
plot.ra.pred = spplot(R.atp, 'pred' , col.regions = rev(brewer.pal(n = 11, name = "Spectral")), cuts = 10, par.settings = list(strip.background = list(col =  'transparent')), colorkey = list(space = "right", height = 1))
plot.ra.var = spplot(R.atp, 'var' , col.regions = rev(brewer.pal(n = 11, name = "PuOr")), cuts = 10, par.settings = list(strip.background = list(col =  'transparent')), colorkey = list(space = "right", height = 1))
grid.arrange(plot.ra, plot.ra.pred, plot.ra.var, ncol = 3)
 
# ---------------------------------- 5 -----------------------------------------
# ------------------------------- ATA CoKriging --------------------------------
# ------------------------------------------------------------------------------

# ----- Poisson cross-semivariogram deconvolution
vg.deconv.cok <- deconvPointVgmForCoKriging(rslt.combined, means.co, model="Exp", ngroup=12, rd = 0.36, maxSampleNum=500, fig = T)

#------ ATA Poisson cokriging
pred.ataok <- ataKriging(rslt.a, unknown = rslt.a , vg.deconv, nmax = 32,  showProgress = TRUE)

