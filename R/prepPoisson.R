

seed = 123
# ---------------------------------- 1 -----------------------------------------
# ----------- Generating Poisson Multivariate Spatial Data ---------------------
# ------------------------------------------------------------------------------

# ----- Bivariate Linear Model of Coregionalization


# --- Define region
size = 20*20
x = seq(-10,10,length.out = sqrt(size))
y = seq(-10,10,length.out = sqrt(size))
coordenadas = expand.grid(x = x, y = y)
dist.matrix = as.matrix(dist(coordenadas))

# --- Define covariance structures
Cexp = function(h,a,b){b * exp(-h/a)}
Csph = function(h,a,b){ifelse(h <= a, b * (1-1.5*(h/a)+0.5*(h/a)^3), 0)}

# --- Define coefficient matrices (LMC)
B1 = matrix(c(0.81, 0.36, 0.36, 0.16),nrow=2,byrow=T)
B2 = matrix(c(0.0081, 0.054, 0.054 , 0.360),nrow=2,byrow=T)
is.positive.semi.definite(B1)
is.positive.semi.definite(B2)

# --- Define processes covariance matrix
cov11 = B1[1,1] * Cexp(dist.matrix,2,1) + B2[1,1] * Csph(dist.matrix,5,1)
cov22 = B1[2,2] * Cexp(dist.matrix,2,1) + B2[2,2] * Csph(dist.matrix,5,1)
cov12 = B1[1,2] * Cexp(dist.matrix,2,1) + B2[1,2] * Csph(dist.matrix,5,1)
cov21 = B1[2,1] * Cexp(dist.matrix,2,1) + B2[2,1] * Csph(dist.matrix,5,1)
cov.total =rbind(cbind(cov11,cov12),cbind(cov21,cov22))
is.positive.definite(cov.total)

# Plot LMC
dum = seq(1, 10, length.out = 1000)
par(mfrow = c(2,2))
plot(dum, (B1[1,1] * Cexp(0,2,1) + B2[1,1] * Csph(0,5,1)) -(B1[1,1] * Cexp(dum,2,1) + B2[1,1] * Csph(dum,5,1)), type = 'l', xlab = 'h' , ylab = 'semivariance', main = expression(~ Y[alpha]/n[alpha] - Y[alpha]/n[alpha]))
plot(dum, (B1[1,2] * Cexp(0,2,1) + B2[1,2] * Csph(0,5,1)) -(B1[1,2] * Cexp(dum,2,1) + B2[1,2] * Csph(dum,5,1)), type = 'l', xlab = 'h' , ylab = 'semivariance', main = expression(~ Y[alpha]/n[alpha] - Y[beta]/n[beta]))
plot(dum, (B1[2,1] * Cexp(0,2,1) + B2[2,1] * Csph(0,5,1)) -(B1[2,1] * Cexp(dum,2,1) + B2[2,1] * Csph(dum,5,1)), type = 'l', xlab = 'h' , ylab = 'semivariance', main = expression(~ Y[beta]/n[beta] - Y[alpha]/n[alpha]))
plot(dum, (B1[2,2] * Cexp(0,2,1) + B2[2,2] * Csph(0,5,1)) -(B1[2,2] * Cexp(dum,2,1) + B2[2,2] * Csph(dum,5,1)), type = 'l', xlab = 'h' , ylab = 'semivariance',main = expression(~ Y[beta]/n[beta] - Y[beta]/n[beta]))
# --- Generate multivariate gaussian field
set.seed(seed)
sim = rmvnorm(1,c(rep(49.86081,nrow(coordenadas)), rep(11.9,nrow(coordenadas))), sigma=cov.total)
datos = as.data.frame(cbind(coordenadas,Ra = sim[1:size], Rb= sim[(size+1):(2*size)]))
cor(datos$Ra, datos$Rb)

# --- Determine shared risk
rho = 0.6
rab = rho*sqrt(datos$Ra * datos$Rb)
datos$Rab <- rab

# --- Create dataframe
R = SpatialPixelsDataFrame(points = coordenadas, data = datos[,3:5])

# --- Plot Dataframe
plot.ra = spplot(R, 'Ra' , col.regions = rev(brewer.pal(n = 11, name = "Spectral")), cuts = 10, par.settings = list(strip.background = list(col =  'transparent')), colorkey = list(space = "right", height = 1))
plot.rb = spplot(R, 'Rb' ,  col.regions = rev(brewer.pal(n = 11, name = "Spectral"))  , cuts = 10, par.settings = list(strip.background = list(col =  'transparent')), colorkey = list(space = "right", height = 1))
plot.rab = spplot(R, 'Rab' , col.regions = rev(brewer.pal(n = 11, name = "Spectral")) , cuts = 10, par.settings = list(strip.background = list(col =  'transparent')), colorkey = list(space = "right", height = 1))
grid.arrange(plot.ra, plot.rb, plot.rab, ncol = 3)

# --- Generate Poisson data
# Population sizes
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
par(mfrow = c(1,4))
hist(n, breaks = 15, main = '')

# Count data
Ya = rpois(size, n*R$Ra)
Yb = rpois(size, n*R$Rb)
Yab <- rpois(size, n*R$Rab)
# Add data to dataframe
R$Ya = Ya
R$Yb = Yb
R$Yab = Yab
R$Na = n
R$Nb = n

#DATA
data.test <- data.frame(areaId= seq(1,400,1),x = coordenadas$x , y = coordenadas$y,  ya = R$Ya, n = R$Na)


cellsize = 0.5263158
blkResoX = abs(x[1]-x[2])
blkResoY = abs(y[1]-y[2])
xResoNum = round(blkResoX / cellsize)
yResoNum = round(blkResoY / cellsize)
w <- 1/(xResoNum*yResoNum)


discretePoints <- lapply(1:nrow(data.test), function(i) {
  xy <- expand.grid(list(ptx = data.test$x[i] - 0.5 * blkResoX + (1:xResoNum-0.5) * cellsize,
                         pty = data.test$y[i] - 0.5 * blkResoY + (1:yResoNum-0.5) * cellsize))
  data.frame(areaId = data.test$areaId[i], xy, weight=w)
})

discretePoints <- do.call(rbind, discretePoints)


colnames(data.test)[2:5] <- c("centx","centy","counts", "size")
rslt <- list(areaValues = data.test, discretePoints = discretePoints)
class(rslt) <- c("list", "discreteArea")

rslt
