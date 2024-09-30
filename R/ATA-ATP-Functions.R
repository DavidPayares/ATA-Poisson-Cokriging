#-------------------------------------------------------------------------------
#---------------------------- Libraries ----------------------------------------
#-------------------------------------------------------------------------------

library(Rcpp)
sourceCpp('C:/Users/PayaresGarciaDE/Documents/Work/PhD/PCK - ATA/atakrig/src/areaVgm.cpp')
#sourceCpp('C:/Users/PayaresGarciaDE/Documents/Work/PhD/PCK - ATA/atakrig/src/RcppExports.cpp')
source('C:/Users/PayaresGarciaDE/Documents/Work/PhD/PCK - ATA/atakrig/R/parallelsetting.R')

#-------------------------------------------------------------------------------
#---------------------------- Functions ----------------------------------------
#-------------------------------------------------------------------------------
#---------------------- Utils

# --- Determine shared risk (proportion)
rlk = function(Rl, Rk, rho, log = T) {
  if (log) {
    shared = log(rho * sqrt(exp(Rl) * exp(Rk)))
  }
  else{
    shared = rho * sqrt(Rl * Rk)
  }
  return(shared)
}


#--- Calculate weighted centroid based on population
weightedpop <- function(x, y, pop) {
  return(data.frame(x = weighted.mean(x, pop, na.rm = T), y = weighted.mean(y, pop, na.rm = T)))
}

# --- Calculate integer based on population decimals
int.pop <- function(x) {
  integer_part <- floor(x)
  decimal_part <- x - integer_part
  
  if (decimal_part > 0.5) {
    return(ceiling(x))
  } else {
    return(integer_part)
  }
}

# Distance berween observations
spDists <- function(x, y=x, longlat=FALSE) {
  return(spDistsNN(x[,1], x[,2], y[,1], y[,2], longlat))
}

# Covariance structure indices
create_matrix <- function(n) {
  if (n <= 0) {
    stop("Input must be a positive integer")
  }
  
  mat <- matrix(NA, nrow = n, ncol = n)
  
  # Fill diagonal
  for (i in 1:n) {
    mat[i, i] <- i
  }
  
  # Fill upper triangular part
  num <- n + 1
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      mat[i, j] <- num
      mat[j, i] <- num
      num <- num + 1
    }
  }
  
  return(mat)
}

# Bias matrix for cokriging
bias.matrix = function(data, nvars){
  
  # Matrix size
  matrix.size = nrow(data) * nvars
  
  # Empty matrix
  matrix.structured = matrix(0, nrow = matrix.size, ncol = matrix.size)
  
  #Population
  N = data[,ncol(data)]
  
  #--- Fill the values
  
  # Diagonal
  diagonals <- list()
  
  for (var in 2:(ncol(data) - 1)){
    
    # Get variables
    col = data[, var]
    # Select NA
    nas =  ifelse(is.na(col), NA, 1)
    # Select pop
    pop = N * nas
    # Estimate bias
    bias =  sum(col, na.rm = T) / sum(pop, na.rm = T)
    # Calculate the weight
    weight = bias / pop
    # Save value for diagonal
    #diagonals[,var] <- weight
    diag.w =  diag(weight)
    diagonals <- c(diagonals, list(diag.w))
  }
  
  # Get indexes
  result <- create_matrix(nvars)
  
  # Define matrix
  row.list <- list()
  for (row in 1:nvars){
    row.list <- c(row.list, list(do.call(cbind,diagonals[c(result[,row])])))
  }
  
  biases <- do.call(rbind, row.list)
  
  # Append columns to the original matrix
  biases <- cbind(biases, matrix(0, nrow = matrix.size, ncol = nvars))
  
  # Append M rows to the extended matrix
  biases <- rbind(biases, matrix(0, nrow = nvars, ncol = nvars + matrix.size))
  
  # Remove NAs
  biases <- remove.na(biases)
  
  return(biases)
}

remove.na = function(mat){
  
  # Identify rows and columns with NA values in the diagonal
  na.diag.index <- which(is.na(diag(mat)))
  
  if(length(na.diag.index) > 0){
    mat.cleaned <- mat[-na.diag.index, -na.diag.index]
  }else{
    mat.cleaned <- mat
  }
  return(mat.cleaned)
}


sort_sample_ids <- function(sampleIds) {
  
  # Split each element of sampleIds into "data" part and numeric part
  parts <- strsplit(sampleIds, "_")
  
  # Extract "data" part and numeric part
  data_part <- sapply(parts, function(x) x[1])
  numeric_part <- as.numeric(sapply(parts, function(x) x[2]))
  
  # Create a data frame with original sampleIds, data_part, and numeric_part
  df <- data.frame(sampleIds = sampleIds, data_part = data_part, numeric_part = numeric_part)
  
  # Sort data frame first by data_part, then by numeric_part
  df_sorted <- df[order(df$data_part, df$numeric_part), ]
  
  # Return sorted sampleIds vector
  return(df_sorted$sampleIds)
}



#---------------------- Semivariogram deconvolution -------------------------------

#--- Poisson cloud direct semivariogram
calcSvCloud <- function(x, y=x, xy=x, longlat=FALSE) {
  distM <- spDists(as.matrix(x[,1:2]), as.matrix(y[,1:2]), longlat=longlat)
  gammaM <- distM * NA
  wM <- distM * NA
  
  if(!is.null(nrow(xy))){
    bias <- sum(xy[,3]) / sum(xy[,4])
  }else{
    bias = xy
  }

  
  for(i in 1:nrow(x)) {
    wM[i,] <- x[i,4] * y[,4] / (x[i,4] +  y[,4])
    gammaM[i,] <- ( wM[i,] * (x[i,3]/x[i,4]-x[,3]/x[,4]) * (y[i,3]/y[i,4]-y[,3]/y[,4]) - bias) / (2 * wM[i,])
  }
  
  if(identical(x, y)) {
    d <- distM[lower.tri(distM)]
    g <- gammaM[lower.tri(gammaM)]
    dg <- data.frame(d=d, g=g)
  } else {
    dg <- data.frame(d=as.vector(distM), g=as.vector(gammaM))
  }
  return(dg)
}

#--- Poisson empirical direct semivariogram
calcSvDgGroup <- function(dgScatter, ngroup=12, rd=0.66, removeOdd=FALSE) {
  
  if(ngroup > 0) {
    # 		dgGroup <- c()
    dmin <- min(dgScatter$d)
    dmax <- max(dgScatter$d)*rd
    
    dfrom <- dmin
    dinterval <- (dmax-dmin)/ngroup
    dgGroup <- data.frame(np=rep(0,ngroup), dist=rep(0,ngroup),
                          gamma=rep(NA,ngroup), dir.hor=rep(0,ngroup),
                          dir.ver=rep(0,ngroup), id=rep('var1',ngroup))
    ncnt <- 0
    for(i in 1:ngroup) {
      dto <- dfrom+dinterval
      indx <- (dgScatter$d >= dfrom & dgScatter$d < dto)
      
      m <- sum(indx)
      if(m >= 1) {
        ncnt <- ncnt+1
        dgGroup[ncnt,] <- data.frame(np=m, dist=mean(dgScatter$d[indx]),
                                     gamma=mean(dgScatter$g[indx]), dir.hor=0, dir.ver=0, id='var1')
      }
      dfrom <- dto
    }
    if(ncnt < ngroup) dgGroup <- dgGroup[1:ncnt,]
  } else {
    m <- length(dgScatter$d)
    dgGroup <- data.frame(np=rep(1,m), dist=dgScatter$d, gamma=dgScatter$g,
                          dir.hor=rep(0,m), dir.ver=rep(0,m), id=rep('var1',m))
  }
  
  dgGroup$np <- as.numeric(dgGroup$np)
  
  row.names(dgGroup) <- NULL
  class(dgGroup) <- c("gstatVariogram", "data.frame")
  return(dgGroup)
}

#--- Fitting theoterical direct semivariogram to empirical one
fitPointVgm <- function(vgmdef, model=c("Sph","Exp"), fit.nugget=FALSE, fixed.range=NA, fig=FALSE, ...) {
  init_nugget <- ifelse(!fit.nugget, 0, min(vgmdef$gamma))
  init_range <- 0.1*(max(vgmdef$dist)-min(vgmdef$dist))   # 0.10 times the length of the central axis through the area
  init_sill <- mean(c(max(vgmdef$gamma), median(vgmdef$gamma)))
  if(identical(model,''))
    model <- c("Sph","Exp") # as.character(vgm()$short)
  
  fitModels <- function(psill, range, nugget, ...) {
    serror <-Inf
    fitmodel <- NULL
    mfit <- NULL
    for(m in model) {
      if(!is.na(fixed.range)) range <- fixed.range
      initSv <- if(!fit.nugget) vgm(psill,m,range,...) else vgm(psill,m,range,nugget,...)
      
      fit.method <- list(...)$fit.method
      if(is.null(fit.method)) fit.method <- 7
      tryCatch(suppressWarnings(mfit <- fit.variogram(vgmdef, initSv, fit.ranges = is.na(fixed.range), fit.method = fit.method)),
               error=function(e) mfit <<- NULL)
      if(is.null(mfit) || (!is.null(mfit) && any(c(mfit$psill, mfit$range) < 0))) {
        # remove Odds
        gmn <- mean(vgmdef$gamma)
        gst <- sd(vgmdef$gamma)
        indx <- (vgmdef$gamma > gmn+3*gst) | (vgmdef$gamma < gmn-3*gst)
        if(sum(indx) > 0) {
          tryCatch(suppressWarnings(mfit <- fit.variogram(vgmdef[!indx,], initSv, fit.ranges = is.na(fixed.range), fit.method = fit.method)),
                   error=function(e) mfit <<- NULL)
        }
      }
      if(is.null(mfit)) next
      
      sserr <- attr(mfit, "SSErr")
      if(is.na(sserr)) next
      
      if(sserr < serror) {
        serror <- attr(mfit, "SSErr")
        fitmodel <- mfit
      }
      # 			if(!attr(mfit, "singular")) break
    }
    # 		return(list(bins=vgmdef, model=fitmodel, sserr=serror))
    return(list(model=fitmodel, sserr=serror))
  }
  
  # try with initial parameters
  psill <- init_sill-init_nugget
  range <- init_range
  nugget <- init_nugget
  fModel <- fitModels(psill, range, nugget, ...)
  if(!is.null(fModel) && !is.null(fModel$model)) {
    fModel$bins <- vgmdef
    if(fig) print(plot(fModel$bins, fModel$model))
    return(fModel)
  }
  
  # try with different parameter combinations
  serror <-Inf
  # comb <- c(1/3, 3, 1/9, 9, 1/27, 27, 1/100, 100, 1/500, 500)
  comb <- c(1/3, 3, 1/9, 9, 1/27, 27)
  for(n in nugget*comb) {
    for(r in range*comb) {
      for(p in psill*comb) {
        mfit <- fitModels(p, r, n, ...)
        if(!is.null(attr(mfit, "SSErr")) && attr(mfit, "SSErr") < serror) {
          serror <- attr(mfit, "SSErr")
          fModel <- mfit
          if(!attr(fModel, "singular")) {
            fModel$bins <- vgmdef
            if(fig) plot(fModel$bins, fModel$model)
            return(fModel)
          }
        }
      }
    }
  }
  
  fModel$bins <- vgmdef
  if(fig) plot(fModel$bins, fModel$model)
  return(fModel)
}

#--- Computing and fitting Poisson direct semivariogram
autofitVgm <- function(x, y=x, xy=x, ngroup=c(12,15), rd=seq(0.3,0.9,by=0.1), model=c("Sph","Exp"),
                       fit.nugget=FALSE, fixed.range=NA, longlat=FALSE, fig=FALSE, ...) {
  
  if(!is.null(nrow(xy))){
    xy <- xy[,2:5]
  }
  
  areaSvCloud <- calcSvCloud(x[,2:5], y[,2:5], xy, longlat=longlat)
  mSel <- NULL
  for (i in 1:length(ngroup)) {
    for(j in 1:length(rd)) {
      areaVgm <- calcSvDgGroup(areaSvCloud, ngroup=ngroup[i], rd=rd[j])
      #print(areaVgm)
      if(nrow(areaVgm) < 10) {
        #message("too few points to fit variogram model!")
        next
      }
      
      m <- fitPointVgm(vgmdef=areaVgm, model=model, fit.nugget=fit.nugget, fixed.range=fixed.range, ...)
      if(is.null(m$model)) next
      
      m$ngroup <- ngroup[i]
      m$rd <- rd[j]
      if(is.null(mSel)) {
        mSel <- m
      } else {
        if(m$sserr < mSel$sserr) mSel <- m
      }
    }
  }
  if(fig && !is.null(mSel)) print(plot(mSel$bins, mSel$model))
  return(mSel)
}

#--- Sampling observations for semivariogram deconvolution
sysSampling <- function(pop, n) {
  if(length(pop) == 1 && is.numeric(pop)) {
    pop <- 1:pop
  }
  
  N <- length(pop)
  ss <- pop[seq(1, N, ceiling(N / n))]
  if(length(ss) < n) {
    ss <- sort(c(ss, sysSampling(setdiff(pop, ss), n-length(ss))))
  }
  return(ss)
}


#--- Plot deconvoluted semivariogram
plotDeconvVgm <- function(v, main=NULL, posx=NULL, posy=NULL, lwd=2, showRegVgm=FALSE) {
  plotDVgm <- function(v, main) {
    xlim <- c(0, 1.1 * max(v$experientialAreaVariogram$dist))
    xx <- seq(0, xlim[2], length=100)
    xx[1] <- 1.0e-3
    yy <- variogramLine(v$areaVariogram, covariance=FALSE, dist_vector=xx)$gamma
    yy2 <- variogramLine(v$pointVariogram, covariance=FALSE, dist_vector=xx)$gamma
    
    ylim <- c(min(0, min(v$experientialAreaVariogram$gamma)),
              1.1 * max(c(v$experientialAreaVariogram$gamma, yy, yy2)))
    plot(v$experientialAreaVariogram$dist, v$experientialAreaVariogram$gamma,
         xaxs="i", yaxs="i", col="black", xlim=xlim, ylim=ylim, ann = FALSE, pch = 3,
         xlab="distance", ylab="semivariance", main=main)
    lines(xx, yy, col="black", lty=3)
    lines(xx, yy2, col="black", lwd=1)
    if(showRegVgm)
      lines(v$regularizedAreaVariogram$dist, v$regularizedAreaVariogram$gamma, col="black", type="l", lty=2, lwd=lwd)
    mtext(main, side = 3, line = 0.2, cex = 0.8)
  }
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  par(mar=c(2.5, 3.5, 1, .5)) # reduce the margins around the figure
  par(mgp=c(1.5, .5, 0)) #  reduce the spacing between the figure plotting region and the axis labels
  par(oma=c(1, 1, 2.5, 0)) #  add an outer margin to the top of the graph
  
  if (hasName(v, "regularizedAreaVariogram")) {
    plotDVgm(v, main = "")
    if(is.null(posx)) posx <- "bottomright"
  } else {
    m0 <- (sqrt(1+8*length(v))-1)/2
    n <- sort(names(v)[1:m0])
    m1 <- m0*(m0-1)/2 + m0
    m <- matrix(0, nrow=m0, ncol=m0)
    m[lower.tri(m, diag = TRUE)] <- seq(m1)
    m[1,2:m0] <- m1 + 1
    layout(mat=m)
    
    for (i in 1:m0) {
      plotDVgm(v[[n[i]]], main = n[i])
      for (j in 1:m0) {
        if(j > i) {
          n2 <- crossName(n[i],n[j])
          plotDVgm(v[[n2]], main = n2)
        }
      }
    }
    plot.new()
    if(is.null(posx)) posx <- "left"
  }
  
  if(showRegVgm) {
    legend(posx, posy, bty="n",
           legend=c("Empirical area variogram", "Fitted area variogram",
                    "Deconvoluted point variogram", "Regularized area variogram"),
           pch=c(3,NA,NA,NA), lty=c(NA,1,1,2), lwd=c(NA,lwd,lwd,lwd),
           col=c("black","black","black","black"))
  } else {
    legend(posx, posy, bty="n",
           legend=c("Empirical area variogram", "Fitted area variogram",
                    "Deconvoluted point variogram"),
           pch=c(3,NA,NA), lty=c(NA,3,1), lwd=c(NA,1,1),
           col=c("black","black","black","black"))
  }
  
  if(is.null(main)) main <- "Deconvoluted variogram"
  #mtext(main, outer = TRUE,  side = 3, cex = 1.2, line = 0)
  mtext("distance", outer = TRUE,  side = 1, cex = 0.8, line = 0)
  mtext("semivariance", outer = TRUE,  side = 2, cex = 0.8, line = -1)
}

#--- Deconvolution of the direct semivariogram
deconvPointVgm <- function(x, model="Exp", maxIter=100, fixed.range=NA, longlat=FALSE,
                           maxSampleNum=100, fig=TRUE, rd= rd, ...) {
  
  # Pierre Goovaerts suggested a different weight (similar to fit.method = 2) in fitting the experimental variogram.
  fit.nugget <- FALSE
  
  if(!hasName(x, "discretePoints"))
    x$discretePoints <- cbind(x$areaValues[,1:3], data.frame(weight=rep(1,nrow(x$areaValues))))
  
  if(nrow(x$areaValues) > maxSampleNum) {
    # warnings("Too many points for generating point-pairs. Only a sample will be used.")
    indx <- sort(sysSampling(nrow(x$areaValues), maxSampleNum))
    x <- subsetDiscreteArea(x, x$areaValues[indx,1])
  }
  
  ## 1. area scale semivariogram
  # areaSvmodel <- autofitVgm(x$areaValues, fit.nugget=fit.nugget, fixed.range=fixed.range, longlat=longlat, model=model, ...)
  areaSvmodel <- autofitVgm(x$areaValues, fit.nugget=fit.nugget, longlat=longlat, rd =  rd) #model=model, ngroup = ngroup, rd = rd, fig = TRUE)
  if(is.null(areaSvmodel)) {
    message("Fitting area-scale variogram model failed!\n")
    return(NULL)
  }
  if(nrow(x$areaValues) == nrow(x$discretePoints)) {
    vgms <- list(pointVariogram = areaSvmodel$model,
                 areaVariogram = areaSvmodel$model,
                 experientialAreaVariogram = areaSvmodel$bins,
                 regularizedAreaVariogram = NULL,
                 status = 1,
                 sserr = areaSvmodel$sserr)
    class(vgms) <- c("list", "ataKrigVgm")
    return(vgms)
  }
  ngroup <- areaSvmodel$ngroup
  rd <- areaSvmodel$rd
  
  ## 0. prepare
  x$areaValues <- x$areaValues[order(x$areaValues$areaId),]
  
  # distance between area centroids
  areaDistByCentroid <- spDists(as.matrix(x$areaValues[,2:3]), longlat=longlat)
  
  ## 2. initialize
  pointSvmodel <- areaSvmodel$model
  dgGroup <- areaSvmodel$bins
  gamaExp <- variogramLine(areaSvmodel$model, covariance=FALSE, dist_vector=dgGroup$dist)$gamma
  s2Exp <- areaSvmodel$model$psill
  
  svAreaCloudByPointVgmInit(x$discretePoints, areaDistByCentroid, longlat)
  
  ## 3. regularize
  gamaReg <- calcSvDgGroup(svAreaCloudByPointVgm(pointSvmodel), ngroup = ngroup, rd = rd)$gamma
  # gamaReg <- calcSvDgGroup(calcSvAreaCloudByPointVgm(x$discretePoints, pointSvmodel, areaDistByCentroid, areaDistByPts),
  #                           ngroup = ngroup, rd = rd)$gamma
  
  ## 4. difference of areal regularized and experimental semivariogram
  D0 <- mean(abs(gamaReg-gamaExp)/gamaExp)
  
  ## 5. backup for comparison
  pointSvmodelOpt <- pointSvmodel
  gamaOpt <- gamaReg
  DOpt <- D0
  
  nvib <- 0
  bNewW <- FALSE
  status <- -1
  for(iter in 1:maxIter) {
    message(sprintf("\riterating: %d", iter))
    
    ## 6. update regularized semivariance
    if(!bNewW)
      wl <- 1 + (gamaExp-gamaOpt)/(s2Exp*sqrt(iter))
    # gPoint <- variogramLine(pointSvmodelOpt, covariance=FALSE, dist_vector=dgGroup$dist)$gamma * wl
    gPoint <- variogramLineSimple(pointSvmodelOpt, dgGroup$dist, FALSE)$gamma * wl
    
    ## 7. fit new point scale semivariogram
    dgGroup$gamma <- gPoint
    pointSvmodel <- fitPointVgm(dgGroup, pointSvmodel[nrow(pointSvmodel),1],
                                fit.nugget=fit.nugget, fixed.range=fixed.range)$model
    if(is.null(pointSvmodel)) {
      pointSvmodel <- pointSvmodelOpt
      status <- 0
      # warning("fitting point-scale variogram failed!")
      break
    }
    
    ## 8. regularize
    # dgScatter <- calcSvAreaCloudByPointVgm(x$discretePoints, pointSvmodel, areaDistByCentroid, areaDistByPts)
    dgScatter <- svAreaCloudByPointVgm(pointSvmodel)
    gamaReg <- calcSvDgGroup(dgScatter, ngroup = ngroup, rd = rd)$gamma
    
    ## 9. difference of areal regularized and experimental semivariogram
    D1 <- mean(abs(gamaReg-gamaExp)/gamaExp)
    
    # 10. stop criterion
    if(max(abs(wl-1)) < 0.001) {
      status <- 1
      break
    }
    if(abs(D1-DOpt)/DOpt <= 0.01) {
      nvib <- nvib+1
      if(nvib >= 5) {
        status <- 2
        break
      }
    } else {
      nvib <- 0
    }
    if(D1/D0 <= 0.01) {
      status <- 3
      break
    }
    
    if(D1 < DOpt) {
      pointSvmodelOpt <- pointSvmodel
      gamaOpt <- gamaReg
      DOpt <- D1
      bNewW <- FALSE
    } else {
      pointSvmodel <- pointSvmodelOpt
      wl <- 1 + (wl-1)/2
      bNewW <- TRUE
    }
  } # end of iteration
  svAreaCloudByPointVgmEnd()
  
  if(iter == maxIter) {
    status <- 4
  }
  # message("\r",rep(" ",15),"\r")
  
  vgms <- list(pointVariogram = pointSvmodel,
               areaVariogram = areaSvmodel$model,
               experientialAreaVariogram = areaSvmodel$bins,
               regularizedAreaVariogram = data.frame(dist=areaSvmodel$bins$dist, gamma=gamaOpt),
               status = status,
               sserr = areaSvmodel$sserr)
  class(vgms) <- c("list", "ataKrigVgm")
  
  if(fig) try(plotDeconvVgm(vgms), silent = TRUE)
  
  return(vgms)
}



#---------------------- ATA Poisson Kriging ------------------------------------

#--- Extract area centrods
calcAreaCentroid <- function(discretePoints) {
  uId <- sort(unique(discretePoints[,1]))
  xys <- matrix(nrow=length(uId), ncol=2)
  for(i in 1:length(uId)) {
    indx <- (discretePoints[,1] == uId[i])
    ws <- sum(discretePoints[indx,4])
    x <- sum(discretePoints[indx,2] * discretePoints[indx,4]/ws)
    y <- sum(discretePoints[indx,3] * discretePoints[indx,4]/ws)
    xys[i,] <- c(x,y)
  }
  
  colnames(xys) <- c("centx","centy")
  xys <- as.data.frame(xys)
  xys$areaId <- uId
  xys <- xys[,c("areaId","centx","centy")]
  return(xys)
}

#--- ATA Poisson kriging (local neighborhood)
ataKriging.local <- function(x, unknown, ptVgm, nmax=100, longlat=FALSE, showProgress=FALSE, nopar=FALSE, clarkAntiLog=FALSE) {
  
  
  if(is(unknown, "discreteArea")) unknown <- unknown$discretePoints
  
  unknown <- unknown[sort.int(unknown[,1], index.return = TRUE)$ix,]
  unknownCenter <- calcAreaCentroid(unknown)
  
  nb <- FNN::get.knnx(as.matrix(x$areaValues[,2:3,drop=FALSE]), as.matrix(unknownCenter[,2:3,drop=FALSE]), nmax)
  nb$nn.index <- matrix(x$areaValues[nb$nn.index,1], ncol = ncol(nb$nn.index))
  
  unknownAreaIds <- unknownCenter[,1]
  
  krigOnce <- function(k) {
    curUnknown <- unknown[unknown[,1] == unknownAreaIds[k],]
    curAreaPts <- x$discretePoints[x$discretePoints[,1] %in% nb$nn.index[k,],]
    curAreaVals <- x$areaValues[x$areaValues[,1] %in% nb$nn.index[k,],]
    
    estResult <- ataKriging(x=list(discretePoints=curAreaPts, areaValues = curAreaVals),
                            curUnknown, ptVgm, nmax=Inf, longlat=longlat, showProgress=FALSE, nopar=TRUE, clarkAntiLog)
    return(estResult)
  }
  
  hasCluster <- ataIsClusterEnabled()
  if(showProgress) pb <- txtProgressBar(min=0, max=length(unknownAreaIds), width = 50, style = 3)
  
  if(!hasCluster || nopar || length(unknownAreaIds) == 1) {
    estResults <- c()
    for (k in 1:length(unknownAreaIds)) {
      estResults <- rbind(estResults, krigOnce(k))
      if(showProgress) setTxtProgressBar(pb, k)
    }
  } else {
    progress <- function(k) if(showProgress) setTxtProgressBar(pb, k)
    estResults <-
      foreach(k = 1:length(unknownAreaIds), .combine = rbind, .options.snow=list(progress=progress),
              .export = c("ataKriging","ataCov","calcAreaCentroid","extractPointVgm", "spDistsNN"),
              .packages = c("sp","gstat")) %dopar% {
                krigOnce(k)
              }
    ataClusterClearObj()
  }
  if(showProgress) close(pb)
  
  return(estResults)
}

#--- ATA Covariance matrix
ataCov <- function(areaPts1, areaPts2, ptVgm, longlat=FALSE) {
  # disM <- spDists(as.matrix(areaPts1[,1:2,drop=FALSE]), as.matrix(areaPts2[,1:2,drop=FALSE]), longlat=longlat)
  # mCov <- variogramLine(ptVgm, covariance=TRUE, dist_vector=disM)
  # return(sum(outer(areaPts1[,3], areaPts2[,3]) * mCov))
  disM <- spDistsNN(areaPts1[,1], areaPts1[,2], areaPts2[,1], areaPts2[,2], longlat=longlat)
  mCov <- variogramLineSimple(ptVgm, disM, bCov = TRUE)
  # ------------------------- IMPORTANT
  return(1/sum(outerProd(areaPts1[,3], areaPts2[,3])) * (sum(outerProd(areaPts1[,3], areaPts2[,3]) * mCov)))
}

#--- ATA Poisson kriging
ataKriging <- function(x, unknown, ptVgm, nmax=10, longlat=FALSE, showProgress=FALSE, nopar=FALSE, clarkAntiLog = FALSE) {
  
  stopifnot(nmax > 0)
  if(nmax < Inf) { # local neigbourhood Kriging.
    return(ataKriging.local(x, unknown, ptVgm, nmax, longlat, showProgress, nopar, clarkAntiLog))
  }
  
  if(is(unknown, "discreteArea")) unknown <- unknown$discretePoints
  if(is(ptVgm, "ataKrigVgm")) ptVgm <- extractPointVgm(ptVgm)
  
  sampleIds <- x$areaValues[,1]
  nSamples <- length(sampleIds)		# number of samples
  
  ## Kriging system: C * wmu = D
  # C matrix
  C <- matrix(1, nrow=nSamples+1, ncol=nSamples+1)
  sampleIndex <- list()
  for(i in 1:nSamples) {
    if(length(sampleIndex) < i) sampleIndex[[i]] <- x$discretePoints[,1] == sampleIds[i]
    sampleI <- x$discretePoints[sampleIndex[[i]],]
    for(j in i:nSamples) {
      if(length(sampleIndex) < j) sampleIndex[[j]] <- x$discretePoints[,1] == sampleIds[j]
      sampleJ <- x$discretePoints[sampleIndex[[j]],]
      C[i,j] <- ataCov(sampleI[,2:4], sampleJ[,2:4], ptVgm, longlat = longlat)
      C[j,i] <- C[i,j]
    }
  }
  
  C[nSamples+1,nSamples+1] <- 0
  
  # Population weighted bias
  W <- diag((sum(x$areaValues[,4])/sum(x$areaValues[,5]))/x$areaValues[,5])
  W <- rbind(cbind(W, rep(0, nSamples)), rep(0, nSamples + 1))
  
  # Add bias
  C = C + W
  #print(W)
  
  unknownAreaIds <- sort(unique(unknown[,1]))
  
  
  krigOnce <- function(k) {
    cur <- unknown[unknown[,1] == unknownAreaIds[k], 2:4, drop=FALSE]
    
    # D matrix
    D <- matrix(1, nrow=nSamples+1, ncol=1)
    for(i in 1:nSamples) {
      sampleI <- x$discretePoints[sampleIndex[[i]],]
      D[i] <- ataCov(sampleI[,2:4,drop=FALSE], cur, ptVgm, longlat = longlat)
    }
    
    # solving
    solvedByGInv <- FALSE
    wmu <- try(solve(C, D), TRUE)
    if(is(wmu, "try-error")) {
      wmu <- MASS::ginv(C) %*% D
      solvedByGInv <- TRUE
    }
    w <- wmu[1:nSamples]
    mu <- wmu[(nSamples+1):nrow(wmu)]
    
    # estimation
    yest <- sum(w * (x$areaValues[,4]/x$areaValues[,5]))
    yvar <- ataCov(cur, cur, ptVgm, longlat = longlat) - sum(wmu * D)
    
    if(!clarkAntiLog)
      return(data.frame(areaId=unknownAreaIds[k], pred=yest, var=yvar))
    else
      return(data.frame(areaId=unknownAreaIds[k], pred=yest, var=yvar,
                        pred.Clark=yest + sum(wmu * D)))
  }
  
  hasCluster <- ataIsClusterEnabled()
  if(showProgress) pb <- txtProgressBar(min=0, max=length(unknownAreaIds), width = 50, style = 3)
  
  if(!hasCluster || nopar || length(unknownAreaIds) == 1) {
    estResults <- c()
    for (k in 1:length(unknownAreaIds)) {
      estResults <- rbind(estResults, krigOnce(k))
      if(showProgress) setTxtProgressBar(pb, k)
    }
  } else {
    progress <- function(k) if(showProgress) setTxtProgressBar(pb, k)
    estResults <-
      foreach(k = 1:length(unknownAreaIds), .combine = rbind, .options.snow=list(progress=progress),
              .export = c("ataCov","calcAreaCentroid"),
              .packages = c("sp","gstat")) %dopar% {
                krigOnce(k)
              }
    ataClusterClearObj()
  }
  
  if(showProgress) close(pb)
  
  unknownCenter <- calcAreaCentroid(unknown)
  estResults <- merge(unknownCenter, estResults)
  
  return(estResults)
}


#--- select discretized data by area id
subsetDiscreteArea <- function(x, selAreaId, revSel=FALSE) {
  if(!revSel) {
    rslt <- list(areaValues = x$areaValues[x$areaValues[,1] %in% selAreaId,,drop=FALSE],
                 discretePoints = x$discretePoints[x$discretePoints[,1] %in% selAreaId,,drop=FALSE])
  } else {
    rslt <- list(areaValues = x$areaValues[!x$areaValues[,1] %in% selAreaId,,drop=FALSE],
                 discretePoints = x$discretePoints[!x$discretePoints[,1] %in% selAreaId,,drop=FALSE])
  }
  
  class(rslt) <- c("list", "discreteArea")
  return(rslt)
}

#---------------------- ATP Poisson Kriging ------------------------------------

#--- ATP Poisson kriging
atpKriging <- function(x, unknown0, ptVgm, nmax=10, longlat=FALSE, showProgress=FALSE, nopar=FALSE) {
  unknown <- cbind(areaId=1:nrow(unknown0), unknown0, weight=1)
  return(ataKriging(x, unknown, ptVgm, nmax, longlat, showProgress, nopar))
}



#---------------------- ATP Poisson Kriging ------------------------------------
atpCoKriging <- function(x, unknownVarId, unknown0, comorbidity, ptVgms, nmax=10, longlat=FALSE, oneCondition=FALSE,
                         meanVal=NULL, auxRatioAdj=TRUE, showProgress=FALSE, nopar=FALSE) {
  unknown <- cbind(areaId=1:nrow(unknown0), unknown0, weight=1)
  return(ataCoKriging(x, unknownVarId, unknown, comorbidity, ptVgms, nmax, longlat, oneCondition, meanVal, auxRatioAdj, showProgress, nopar))
}
#---------------------- Cross-semivariogram ----------------------------------------


deconvPointCrossVgm <- function(x, y, m, xPointVgm, yPointVgm, model="Exp", maxIter=100,
                                fixed.range=NA, longlat=FALSE, maxSampleNum=100, fig=TRUE, rd = rd, ...) {
  fit.nugget <- FALSE
  if(!hasName(x, "discretePoints"))
    x$discretePoints <- cbind(x$areaValues[,1:3], data.frame(weight=rep(1,nrow(x$areaValues))))
  if(!hasName(y, "discretePoints"))
    y$discretePoints <- cbind(y$areaValues[,1:3], data.frame(weight=rep(1,nrow(y$areaValues))))
  
  if(nrow(x$areaValues) > maxSampleNum) {
    # warnings("Too many points for generating point-pairs. Only a sample will be used.")
    indx <- sort(sysSampling(nrow(x$areaValues), maxSampleNum))
    x <- subsetDiscreteArea(x, x$areaValues[indx,1])
  }
  if(nrow(y$areaValues) > maxSampleNum) {
    # warnings("Too many points for generating point-pairs. Only a sample will be used.")
    indx <- sort(sysSampling(nrow(y$areaValues), maxSampleNum))
    y <- subsetDiscreteArea(y, y$areaValues[indx,1])
  }
  
  ## 1. area scale cross-semivariogram
  # areaCrossVgm <- autofitVgm(x$areaValues, y$areaValues, fit.nugget=fit.nugget, fixed.range=fixed.range, longlat=longlat, model=model, ...)
  areaCrossVgm <- autofitVgm(x$areaValues, y$areaValues, m, fit.nugget=fit.nugget, longlat=longlat, model=model, rd = rd, ...)
  if(is.null(areaCrossVgm$model)) {
    warning("Fitting area-scale cross-variogram model failed!\n")
    return(NULL)
  }
  if(nrow(x$areaValues) == nrow(x$discretePoints) && nrow(y$areaValues) == nrow(y$discretePoints)) {
    vgms <- list(pointVariogram = areaCrossVgm$model,
                 areaVariogram = areaCrossVgm$model,
                 experientialAreaVariogram = areaCrossVgm$bins,
                 regularizedAreaVariogram = NULL,
                 status = 1,
                 sserr = areaCrossVgm$sserr)
    class(vgms) <- c("list", "ataKrigVgm")
    return(vgms)
  }
  ngroup <- areaCrossVgm$ngroup
  rd <- areaCrossVgm$rd
  
  ## 0. prepare
  x$areaValues <- x$areaValues[order(x$areaValues$areaId),]
  y$areaValues <- y$areaValues[order(y$areaValues$areaId),]
  
  # distance between area centroids
  areaDistByCentroidXY <- spDists(as.matrix(x$areaValues[,2:3]), as.matrix(y$areaValues[,2:3]), longlat=longlat)
  
  # entry index for each area
  uId1 <- sort(unique(x$discretePoints[,1]))
  uId2 <- sort(unique(y$discretePoints[,1]))
  
  # distances between discretized area points
  hasCluster <- FALSE #ataIsClusterEnabled() # parallel version has some problem currently!
  if(!hasCluster) {
    areaDistByPtsX <- list()
    areaWeightByPtsX <- list()
    for(i in 1:length(uId1)) {
      areaIndexI <- which(x$discretePoints[,1] == uId1[i])
      areaDistByPtsX[[i]] <- spDists(as.matrix(x$discretePoints[areaIndexI,2:3]),
                                     as.matrix(x$discretePoints[areaIndexI,2:3]), longlat=longlat)
      # ------------------------- IMPORTANT
      areaWeightByPtsX[[i]] <- 1/sum(outer(x$discretePoints[areaIndexI,4], x$discretePoints[areaIndexI,4])) * outer(x$discretePoints[areaIndexI,4], x$discretePoints[areaIndexI,4])
    }
    
    areaDistByPtsY <- list()
    areaWeightByPtsY <- list()
    for(i in 1:length(uId2)) {
      areaIndexI <- which(y$discretePoints[,1] == uId2[i])
      areaDistByPtsY[[i]] <- spDists(as.matrix(y$discretePoints[areaIndexI,2:3]),
                                     as.matrix(y$discretePoints[areaIndexI,2:3]), longlat=longlat)
      # ------------------------- IMPORTANT
      areaWeightByPtsY[[i]] <- 1/sum(outer(y$discretePoints[areaIndexI,4], y$discretePoints[areaIndexI,4])) * outer(y$discretePoints[areaIndexI,4], y$discretePoints[areaIndexI,4])
    }
    
    areaDistByPtsXY <- list()
    areaWeightByPtsXY <- list()
    for(i in 1:length(uId1)) {
      areaDistByPtsXY[[i]] <- list()
      areaWeightByPtsXY[[i]] <- list()
      areaIndexI <- which(x$discretePoints[,1] == uId1[i])
      for(j in 1:length(uId2)) {
        areaIndexJ <- which(y$discretePoints[,1] == uId2[j])
        areaDistByPtsXY[[i]][[j]] <- spDists(as.matrix(x$discretePoints[areaIndexI,2:3]),
                                             as.matrix(y$discretePoints[areaIndexJ,2:3]), longlat=longlat)
        # ------------------------- IMPORTANT
        areaWeightByPtsXY[[i]][[j]] <- 1/sum(outer(x$discretePoints[areaIndexI,4], y$discretePoints[areaIndexJ,4])) * outer(x$discretePoints[areaIndexI,4], y$discretePoints[areaIndexJ,4])
      }
    }
  } else {
    ll <- foreach(i = 1:length(uId1), .packages = "sp") %dopar% {
      areaIndexI <- which(x$discretePoints[,1] == uId1[i])
      areaDistByPtsX <- spDists(as.matrix(x$discretePoints[areaIndexI,2:3]),
                                as.matrix(x$discretePoints[areaIndexI,2:3]), longlat=longlat)
      # ------------------------- IMPORTANT
      areaWeightByPtsX <- 1/sum(outer(x$discretePoints[areaIndexI,4], x$discretePoints[areaIndexI,4])) * outer(x$discretePoints[areaIndexI,4], x$discretePoints[areaIndexI,4])
      return(list(areaDistByPtsX, areaWeightByPtsX))
    }
    areaDistByPtsX <- ll[[1]]
    areaWeightByPtsX <- ll[[2]]
    
    ll <- foreach(i = 1:length(uId2), .packages = "sp") %dopar% {
      areaIndexI <- which(y$discretePoints[,1] == uId2[i])
      areaDistByPtsY <- spDists(as.matrix(y$discretePoints[areaIndexI,2:3]),
                                as.matrix(y$discretePoints[areaIndexI,2:3]), longlat=longlat)
      # ------------------------- IMPORTANT
      areaWeightByPtsY <- 1/sum(outer(y$discretePoints[areaIndexI,4], y$discretePoints[areaIndexI,4])) * outer(y$discretePoints[areaIndexI,4], y$discretePoints[areaIndexI,4])
      
      return(list(areaDistByPtsY, areaWeightByPtsY))
    }
    areaDistByPtsY <- ll[[1]]
    areaWeightByPtsY <- ll[[2]]
    
    ll <- foreach(i = 1:length(uId1), .packages = "sp") %dopar% {
      areaDistByPtsXY <- list()
      areaWeightByPtsXY <- list()
      areaIndexI <- which(x$discretePoints[,1] == uId1[i])
      for(j in 1:length(uId2)) {
        areaIndexJ <- which(y$discretePoints[,1] == uId2[j])
        areaDistByPtsXY[[j]] <- spDists(as.matrix(x$discretePoints[areaIndexI,2:3]),
                                        as.matrix(y$discretePoints[areaIndexJ,2:3]), longlat=longlat)
        # ------------------------- IMPORTANT
        areaWeightByPtsXY[[j]] <- 1/sum(outer(x$discretePoints[areaIndexI,4], y$discretePoints[areaIndexJ,4])) * outer(x$discretePoints[areaIndexI,4], y$discretePoints[areaIndexJ,4])
      }
      return(list(areaDistByPtsXY, areaWeightByPtsXY))
    }
    areaDistByPtsXY <- list()
    areaWeightByPtsXY <- list()
    for (i in 1:length(ll)) {
      areaDistByPtsXY[[i]] <- ll[[i]][[1]]
      areaWeightByPtsXY[[i]] <- ll[[i]][[2]]
    }
    rm(ll)
    ataClusterClearObj()
  }
  
  ## 2. initialize
  pointCrossVgm <- areaCrossVgm$model
  dgGroup <- areaCrossVgm$bins
  gamaExp <- variogramLine(areaCrossVgm$model, covariance=FALSE, dist_vector=dgGroup$dist)$gamma
  s2Exp <- areaCrossVgm$model$psill
  
  crossSvAreaCloudByPointVgmInit(x$discretePoints, y$discretePoints, xPointVgm, yPointVgm,
                                 areaDistByCentroidXY, areaDistByPtsX, areaDistByPtsY, areaDistByPtsXY,
                                 areaWeightByPtsX, areaWeightByPtsY, areaWeightByPtsXY)
  
  ## 3. regularize
  # crossDgScatter <- calcCrosssvAreaCloudByPointVgm(x$discretePoints, y$discretePoints, xPointVgm, yPointVgm,
  #                                                  xyPointCrossVgm = pointCrossVgm, areaDistByCentroidXY,
  #                                                  areaDistByPtsX, areaDistByPtsY, areaDistByPtsXY)
  crossDgScatter <- crossSvAreaCloudByPointVgm(pointCrossVgm)
  
  gamaReg <- calcSvDgGroup(crossDgScatter, ngroup = ngroup, rd = rd)$gamma
  
  ## 4. difference of areal regularized and experimental semivariogram
  D0 <- mean(abs(gamaReg-gamaExp)/gamaExp)
  
  ## 5. backup for comparison
  pointCrossVgmOpt <- pointCrossVgm
  gamaOpt <- gamaReg
  DOpt <- D0
  
  nvib <- 0
  bNewW <- FALSE
  status <- -1
  for(iter in 1:maxIter) {
    message(sprintf("\riterating: %d", iter))
    
    ## 6. update regularized semivariance
    if(!bNewW)
      wl <- 1 + (gamaExp-gamaOpt)/(s2Exp*sqrt(iter))
    gPoint <- variogramLine(pointCrossVgmOpt, covariance=FALSE, dist_vector=dgGroup$dist)$gamma * wl
    
    ## 7. fit new point scale semivariogram
    dgGroup$gamma <- gPoint
    pointCrossVgm <- fitPointVgm(dgGroup, pointCrossVgm[nrow(pointCrossVgm),1], fit.nugget=fit.nugget, fixed.range=fixed.range, ...)$model
    if(is.null(pointCrossVgm)) {
      pointCrossVgm <- pointCrossVgmOpt
      status <- 0
      # warning("fitting point-scale cross-variogram failed!")
      break
    }
    
    ## 8. regularize
    # crossDgScatter <- calcCrosssvAreaCloudByPointVgm(x$discretePoints, y$discretePoints, xPointVgm, yPointVgm,
    #                                                  xyPointCrossVgm = pointCrossVgm, areaDistByCentroidXY,
    #                                                  areaDistByPtsX, areaDistByPtsY, areaDistByPtsXY)
    crossDgScatter <- crossSvAreaCloudByPointVgm(pointCrossVgm)
    gamaReg <- calcSvDgGroup(crossDgScatter, ngroup = ngroup, rd = rd)$gamma
    
    ## 9. difference of areal regularized and experimental semivariogram
    D1 <- mean(abs(gamaReg-gamaExp)/gamaExp)
    
    # 10. stop criterion
    if(max(abs(wl-1)) < 0.001) {
      status <- 1
      break
    }
    if(abs(D1-DOpt)/DOpt <= 0.01) {
      nvib <- nvib+1
      if(nvib >= 5) {
        status <- 2
        break
      }
    } else {
      nvib <- 0
    }
    if(D1/D0 <= 0.01) {
      status <- 3
      break
    }
    
    if(D1 < DOpt) {
      pointCrossVgmOpt <- pointCrossVgm
      gamaOpt <- gamaReg
      DOpt <- D1
      bNewW <- FALSE
    } else {
      pointCrossVgm <- pointCrossVgmOpt
      wl <- 1 + (wl-1)/2
      bNewW <- TRUE
    }
  }
  crossSvAreaCloudByPointVgmEnd()
  
  if(iter == maxIter) {
    status <- 4
  }
  # message("\r",rep(" ",15),"\r")
  
  vgms <- list(pointVariogram = pointCrossVgm,
               areaVariogram = areaCrossVgm$model,
               experientialAreaVariogram = areaCrossVgm$bins,
               regularizedAreaVariogram = data.frame(dist=areaCrossVgm$bins$dist, gamma=gamaOpt),
               status = status,
               sserr = areaCrossVgm$sserr)
  class(vgms) <- c("list", "ataKrigVgm")
  
  if(fig) try(plotDeconvVgm(vgms), TRUE)
  return(vgms)
}

crossName <- function(id1, id2) {
  if(id1 == id2) {
    return(id1)
  } else {
    id <- sort(c(id1,id2))
    return(paste(id[1], id[2], sep = "."))
  }
}

comb = function(n, x) {
  factorial(n) / factorial(n-x) / factorial(x)
}

deconvPointVgmForCoKriging <- function(x, means, model="Exp", maxIter=100, fixed.range=NA, maxSampleNum=100, fig=TRUE, ...) {
  
  varnames <- sort(names(x))
  vgms <- list()
  
  if(!(comb(length(varnames),2) == length(means))){
    stop("The amount of co-morbidity means does not match the number of variables")
  }
  
  # vargioram
  for (i in 1:length(varnames)) {
    id <- varnames[i]
    if(!hasName(x[[id]], "discretePoints")) {
      x[[id]]$discretePoints <- cbind(x[[id]]$areaValues[,1:3], data.frame(weight=rep(1,nrow(x[[id]]$areaValues))))
      names(x[[id]]$discretePoints)[2:3] <- c("ptx","pty")
    }
    
    if(nrow(x[[id]]$areaValues) > maxSampleNum) {
      # warnings("Too many points for generating point-pairs. Only a sample will be used.")
      indx <- sort(sysSampling(nrow(x[[id]]$areaValues), maxSampleNum))
      x[[id]] <- subsetDiscreteArea(x[[id]], x[[id]]$areaValues[indx,1])
    }
    
    message(sprintf("Deconvoluting variogram of %s ...\n", id))
    suppressMessages(
      vgms[[id]] <- deconvPointVgm(x[[id]], model = model, maxIter = maxIter, fixed.range=fixed.range, fig=FALSE, ...)
    )
    # stopifnot(!is.null(vgms[[id]]))
    if(is.null(vgms[[id]]) || vgms[[id]]$status == 0) {
      warning(sprintf("deconvolution failed for %s!", id))
      return(NULL)
    }
  }
  
  # cross-variogram
  meanpos = 1
  for (i in 1:(length(varnames)-1)) {
    id1 <- varnames[i]
    for (j in (i+1):length(varnames)) {
      id2 <- varnames[j]
      message(sprintf("Deconvoluting cross-variogram between %s and %s ...\n", id1, id2))
      id <- crossName(id1,id2)
      vgms[[id]] <- deconvPointCrossVgm(x[[id1]], x[[id2]], means[[meanpos]],
                                        vgms[[id1]]$pointVariogram,
                                        vgms[[id2]]$pointVariogram,
                                        model = model, maxIter = maxIter, fixed.range=fixed.range,
                                        fig=FALSE, ...)
      meanpos =  meanpos + 1
      stopifnot(!is.null(vgms[[id]]))
      if(vgms[[id]]$status == 0) {
        warning(sprintf("deconvolution failed for %s!", id))
      }
    }
  }
  
  
  if (!is.na(fixed.range)) {
    posdef = function(X) {
      q = eigen(X)
      d = q$values
      d[d < 0] = 0
      q$vectors %*% diag(d, nrow = length(d)) %*% t(q$vectors)
    }
    
    psill = matrix(NA, nrow = length(varnames), ncol = length(varnames))
    for (i in 1:length(varnames)) {
      for (j in i:length(varnames)) {
        id = ifelse(i == j, varnames[i], crossName(varnames[i], varnames[j]))
        psill[i, j] = psill[j, i] = vgms[[id]]$pointVariogram[1,"psill"]
      }
    }
    psill = posdef(psill)
    for (i in 1:length(varnames)) {
      for (j in i:length(varnames)) {
        id = ifelse(i == j, varnames[i], crossName(varnames[i], varnames[j]))
        vgms[[id]]$pointVariogram[1,"psill"] = psill[i, j]
      }
    }
  }
  
  
  if(fig) plotDeconvVgm(vgms, main = "Deconvoluted variograms/cross-variograms")
  
  class(vgms) <- c("list", "ataKrigVgm")
  return(vgms)
}


#---------------------- ATA Poisson Cokriging ----------------------------------------

ataCoKriging.local <- function(x, unknownVarId, unknown, comorbidity, ptVgms, nmax=10, longlat=FALSE,
                               oneCondition=FALSE, meanVal=NULL, auxRatioAdj=TRUE,
                               showProgress=FALSE, nopar=FALSE, clarkAntiLog=FALSE) {
  
  if(is(unknown, "discreteArea")) unknown <- unknown$discretePoints
  
  # sort areaId in ascending order.
  unknown <- unknown[sort.int(unknown[,1], index.return = TRUE)$ix,]
  unknownCenter <- calcAreaCentroid(unknown)
  
  # neighbor indexes for each unknown point.
  varIds <- sort(names(x))
  nb <- list()
  for (id in varIds) {
    nb[[id]] <- FNN::get.knnx(as.matrix(x[[id]]$areaValues[,2:3,drop=FALSE]), as.matrix(unknownCenter[,2:3,drop=FALSE]), nmax)
    nb[[id]]$nn.index <- matrix(x[[id]]$areaValues[,1][nb[[id]]$nn.index], ncol = nmax)
  }
  # only consider covariables within the radius of unknownVarId
  for (id in varIds[varIds != unknownVarId]) {
    indx <- nb[[id]]$nn.dist > matrix(rep(nb[[unknownVarId]]$nn.dist[,nmax] * 1.5, nmax), ncol = nmax)
    nb[[id]]$nn.dist[indx] <- NA
    nb[[id]]$nn.index[indx] <- NA
  }
  
  unknownAreaIds <- sort(unique(unknown[,1]))
  
  krigOnce <- function(k) {

    curUnknown <- unknown[unknown[,1] == unknownAreaIds[k], ]
    
    curx <- list()
    for (id in varIds) {
      if(!hasName(x[[id]], "discretePoints")) {
        x[[id]]$discretePoints <- cbind(x[[id]]$areaValues[,1:3], data.frame(weight=rep(1,nrow(x[[id]]$areaValues))))
        names(x[[id]]$discretePoints)[2:3] <- c("ptx","pty")
      }
      
      curVals <- x[[id]]$areaValues[x[[id]]$areaValues[,1] %in% nb[[id]]$nn.index[k,],]
      curPts <- x[[id]]$discretePoints[x[[id]]$discretePoints[,1] %in% nb[[id]]$nn.index[k,],]
      if(nrow(curVals) > 0) {
        curx[[id]] <- list(areaValues=curVals, discretePoints=curPts)
        #como <- comorbidity[comorbidity$areaId %in% nb[[id]]$nn.index[k,],]
      }
    }
    
    estResult <- ataCoKriging(curx, unknownVarId, curUnknown, comorbidity, ptVgms, nmax=Inf, longlat, oneCondition,
                              meanVal, auxRatioAdj, showProgress=FALSE, nopar=TRUE, clarkAntiLog)
    return(estResult)
  }
  
  hasCluster <- ataIsClusterEnabled()
  if(showProgress) pb <- txtProgressBar(min=0, max=length(unknownAreaIds), width = 50, style = 3)
  
  if(!hasCluster || nopar || length(unknownAreaIds) == 1) {
    estResults <- c()
    for (k in 1:length(unknownAreaIds)) {
      estResults <- rbind(estResults, krigOnce(k))
      if(showProgress) setTxtProgressBar(pb, k)
    }
  } else {
    progress <- function(k) if(showProgress) setTxtProgressBar(pb, k)
    estResults <-
      foreach(k = 1:length(unknownAreaIds), .combine = rbind, .options.snow=list(progress=progress),
              .export = c("x","ataCoKriging","crossName","ataCov","calcAreaCentroid"),
              .packages = c("sp","gstat")) %dopar% {
                krigOnce(k)
              }
    ataClusterClearObj()
  }
  if(showProgress) close(pb)
  
  return(estResults)
}


ataCoKriging <- function(x, unknownVarId, unknown, comorbidity, ptVgms, nmax=10, longlat=FALSE, oneCondition=FALSE,
                         meanVal=NULL, auxRatioAdj=TRUE, showProgress=FALSE, nopar=FALSE, clarkAntiLog=FALSE) {
  
  # x = data.all
  # unknownVarId="data.a"
  # unknown=data.a
  # comorbidity =  means.data
  # ptVgms=vg.deconv.cok
  # nmax = Inf
  # oneCondition=FALSE
  # showProgress = TRUE
  # nopar=FALSE
  # clarkAntiLog=FALSE
  # longlat=FALSE

  stopifnot(nmax > 0)
  if(nmax < Inf) {
    return(ataCoKriging.local(x, unknownVarId, unknown, comorbidity, ptVgms, nmax, longlat, oneCondition,
                              meanVal, auxRatioAdj, showProgress, nopar, clarkAntiLog))
  }
  
  if(is(unknown, "discreteArea")) unknown <- unknown$discretePoints
  if(is(ptVgms, "ataKrigVgm")) ptVgms <- extractPointVgm(ptVgms)
  
  # sort areaId in ascending order.
  for (i in 1:length(x)) {
    x[[i]]$areaValues <- x[[i]]$areaValues[sort.int(x[[i]]$areaValues[,1], index.return = TRUE)$ix,]
  }
  
  # combine all data together.
  varIds <- sort(names(x))
  xAll <- list(areaValues=NULL, discretePoints=NULL)
  for (id in varIds) {
    if(!hasName(x[[id]], "discretePoints")) {
      x[[id]]$discretePoints <- cbind(x[[id]]$areaValues[,1:3], data.frame(weight=rep(1,nrow(x[[id]]$areaValues))))
      names(x[[id]]$discretePoints)[2:3] <- c("ptx","pty")
    }
    
    x[[id]]$areaValues$varId <- id
    x[[id]]$areaValues$var_areaId <- paste(id, x[[id]]$areaValues[,1], sep = "_")
    x[[id]]$discretePoints$varId <- id
    x[[id]]$discretePoints$var_areaId <- paste(id, x[[id]]$discretePoints[,1], sep = "_")
    
    xAll$areaValues <- rbind(xAll$areaValues, x[[id]]$areaValues)
    xAll$discretePoints <- rbind(xAll$discretePoints, x[[id]]$discretePoints)
  }
  
  sampleIds <- sort_sample_ids(sort(unique(xAll$discretePoints$var_areaId)))
  nSamples <- length(sampleIds)		# number of all samples
  nVars <- length(x) # number of variables
  
  ## Kriging system: C * wmu = D
  if(oneCondition) {
    C <- matrix(0, nrow=nSamples+1, ncol=nSamples+1)
    D <- matrix(0, nrow=nSamples+1, ncol=1)
  } else {
    C <- matrix(0, nrow = nSamples + nVars,  ncol= nSamples + nVars)
    W <- matrix(0, nrow = nSamples + nVars, ncol= nSamples + nVars)
    D <- matrix(0, nrow=nSamples+nVars, ncol=1)
  }
  
  
  # C matrix
  sampleIndex <- list()
  for(i in 1:nSamples) {
    if(length(sampleIndex) < i) sampleIndex[[i]] <- xAll$discretePoints$var_areaId == sampleIds[i]
    sampleI <- xAll$discretePoints[sampleIndex[[i]],]
    for(j in i:nSamples) {
      if(length(sampleIndex) < j) sampleIndex[[j]] <- xAll$discretePoints$var_areaId == sampleIds[j]
      sampleJ <- xAll$discretePoints[sampleIndex[[j]],]
      ptVgm <- ptVgms[[crossName(sampleI$varId[1], sampleJ$varId[1])]]
      C[i,j] <- ataCov(sampleI[,2:4], sampleJ[,2:4], ptVgm, longlat = longlat)
      C[j,i] <- C[i,j]
    }
  }

  
  #is.symmetric.matrix(C[1:339,1:339])
  #print(C)
  
  if(oneCondition) {
    C[nSamples+1, ] <- 1
    C[, nSamples+1] <- 1
    C[nSamples+1, nSamples+1] <- 0
    D[nSamples+1] <- 1
  } else {
    for (i in 1:nVars) {
      indx <- xAll$areaValues$varId == varIds[i]
      C[nSamples+i, (1:nSamples)[indx]] <- 1
      C[(1:nSamples)[indx], nSamples+i] <- 1
    }
    D[nSamples + which(unknownVarId == varIds)] <- 1
  }
  
  # varsIds
  crossVarIds <- c()
  for (i in 1:(length(varIds)-1)) {
    name1 <- varIds[i]
    for (j in (i+1):length(varIds)) {
      name2 <- varIds[j]
      names.id <- crossName(name1,name2)
      crossVarIds <- c(crossVarIds, names.id)
    }
  }
  
  # Assign names to mean values
  names(comorbidity) <- c(varIds,crossVarIds)
  
  # Population weighted bias
  for(i in 1:nSamples) {
    for(j in i:nSamples) {
      if (xAll$areaValues$areaId[i] == xAll$areaValues$areaId[j]){
        W[i,j] <- comorbidity[[crossName(xAll$areaValues$varId[i], xAll$areaValues$varId[j])]] / xAll$areaValues$size[j]
        W[j,i] <- W[i,j]
      }
    }
  }
  
  # Add bias term
  C = C + W
  
  unknownAreaIds <- sort(unique(unknown[,1]))
  
  krigOnce <- function(k) {
    
    #k = 1
    curUnknown <- unknown[unknown[,1] == unknownAreaIds[k], 2:4]
    
    # D matrix
    for(i in 1:nSamples) {
      sampleI <- xAll$discretePoints[xAll$discretePoints$var_areaId == sampleIds[i],]
      ptVgm <- ptVgms[[crossName(sampleI$varId[1], unknownVarId)]]
      D[i] <- ataCov(curUnknown, sampleI[,2:4], ptVgm, longlat = longlat)
    }
    
    # solving
    solvedByGInv <- FALSE
    wmu <- try(solve(C, D), TRUE)
    if(is(wmu, "try-error")) {
      wmu <- MASS::ginv(C) %*% D
      solvedByGInv <- TRUE
    }
    
    # estimation
    if(oneCondition) {
      if(is.null(meanVal)) {
        for (id in varIds) {
          meanVal <- rbind(meanVal, data.frame(varId=id, value=mean(x[[id]]$areaValues[,4]/x[[id]]$areaValues[,5])))
        }
      }
      rownames(meanVal) <- meanVal$varId
      
      w <- wmu[1:nSamples]
      w1 <- w[unknownVarId == xAll$areaValues$varId]
      yest <- sum(w1 * x[[unknownVarId]]$areaValues[,4]/x[[unknownVarId]]$areaValues[,5])
      for (id in varIds[varIds != unknownVarId]) {
        w2 <- w[id == xAll$areaValues$varId]
        if (auxRatioAdj && abs(meanVal[id,2]) > 1e-6) {
          yest <- yest + sum(w2 * ((x[[id]]$areaValues[,4]/x[[id]]$areaValues[,5] - meanVal[id, 2])*(meanVal[unknownVarId, 2]/meanVal[id,2]) + meanVal[unknownVarId, 2]))
        } else {
          yest <- yest + sum(w2 * (x[[id]]$areaValues[,4]/x[[id]]$areaValues[,5] - meanVal[id, 2] + meanVal[unknownVarId, 2]))
        }
      }
    } else {
        #w <- wmu[1:nSamples][unknownVarId == xAll$areaValues$varId]
        #yest <- sum(w * (x[[unknownVarId]]$areaValues[,4]/x[[unknownVarId]]$areaValues[,5]))
        w <- wmu[1:nSamples]
        yest <- sum(w * (xAll$areaValues[,4]/xAll$areaValues[,5]) )
    }
    yvar <- ataCov(curUnknown, curUnknown, ptVgms[[unknownVarId]], longlat = longlat) - sum(wmu * D)

    if(!clarkAntiLog)
      return(data.frame(areaId=unknownAreaIds[k], pred=yest, var=yvar))
    else
      return(data.frame(areaId=unknownAreaIds[k], pred=yest, var=yvar,
                        pred.Clark=yest + sum(wmu * D)))
  }
  
  hasCluster <- ataIsClusterEnabled()
  if(showProgress) pb <- txtProgressBar(min=0, max=length(unknownAreaIds), width = 50, style = 3)
  
  if(!hasCluster || nopar || length(unknownAreaIds) == 1) {
    estResults <- c()
    for (k in 1:length(unknownAreaIds)) {
      estResults <- rbind(estResults, krigOnce(k))
      if(showProgress) setTxtProgressBar(pb, k)
    }
  } else {
    progress <- function(k) if(showProgress) setTxtProgressBar(pb, k)
    estResults <-
      foreach(k = 1:length(unknownAreaIds), .combine = rbind, .options.snow=list(progress=progress),
              .export = c("D","meanVal","crossName","ataCov","calcAreaCentroid"),
              .packages = c("sp","gstat")) %dopar% {
                krigOnce(k)
              }
    ataClusterClearObj()
  }
  
  if(showProgress) close(pb)
  
  unknownCenter <- calcAreaCentroid(unknown)
  estResults <- merge(unknownCenter, estResults)
  
  return(estResults)
}


#--- Getting empirical direct semivariogram
extractPointVgm <- function(g) {
  if(!is(g, "ataKrigVgm")) return(NULL)
  
  if(hasName(g, "pointVariogram")) {
    return(g$pointVariogram)
  } else {
    for(id in names(g)) {
      if(hasName(g[[id]], "pointVariogram"))
        g[[id]] <- g[[id]]$pointVariogram
    }
    return(g)
  }
}


## ataKriging.cv: ataKriging cross-validation ----
#   nfold: integer; n-fold cross validation.
ata.PKriging <- function(x, ptVgm, nmax=32,longlat=FALSE, showProgress=FALSE, nopar=FALSE, clarkAntiLog=FALSE) {
  
  N <- nrow(x$areaValues)
  nfold <- N
  indexM <- matrix(sort(x$areaValues[,1]), ncol = 1)
  hasCluster <- ataIsClusterEnabled()
  if(showProgress) pb <- txtProgressBar(min=0, max=nrow(indexM), width = 50, style = 3)
  
  
  if(!hasCluster || nopar || nrow(indexM) == 1) {
    estResults <- c()
    for (k in 1:nrow(indexM)) {
      xknown <- subsetDiscreteArea(x, indexM[k,], revSel = TRUE)
      unknown <- subsetDiscreteArea(x, indexM[k,])$discretePoints
      estResults <- rbind(estResults, ataKriging(xknown, unknown, ptVgm, nmax, longlat, showProgress = FALSE, nopar = TRUE, clarkAntiLog))
      if(showProgress) setTxtProgressBar(pb, k)
    }
  } else {
    bInnerParallel <- ncol(indexM) > 2*nrow(indexM)
    if(bInnerParallel) {
      estResults <- c()
      for (k in 1:nrow(indexM)) {
        xknown <- subsetDiscreteArea(x, indexM[k,], revSel = TRUE)
        unknown <- subsetDiscreteArea(x, indexM[k,])$discretePoints
        estResults <- rbind(estResults, ataKriging(xknown, unknown, ptVgm, nmax, longlat, showProgress = FALSE, nopar = FALSE, clarkAntiLog))
        if(showProgress) setTxtProgressBar(pb, k)
      }
    } else {
      progress <- function(k) if(showProgress) setTxtProgressBar(pb, k)
      estResults <-
        foreach(k = 1:nrow(indexM), .combine = rbind, .options.snow=list(progress=progress),
                .export = c("subsetDiscreteArea","ataCov","calcAreaCentroid","ataKriging","ataKriging.local"),
                .packages = c("sp","gstat","FNN")) %dopar% {
                  xknown <- subsetDiscreteArea(x, indexM[k,], revSel = TRUE)
                  unknown <- subsetDiscreteArea(x, indexM[k,])$discretePoints
                  ataKriging(xknown, unknown, ptVgm, nmax, longlat, showProgress = FALSE, nopar = TRUE, clarkAntiLog)
                }
      ataClusterClearObj()
    }
  }
  if(showProgress) close(pb)
  
  estResults <- estResults[order(estResults$areaId),]
  indx <- match(estResults[,1], x$areaValues[,1])
  # estResults$diff <- x$areaValues[indx,4] - estResults[,4]
  estResults$value <- x$areaValues[indx,4]/x$areaValues[indx,5]
  
  return(estResults)
}


## ataKriging.cv: ataKriging cross-validation ----
#   nfold: integer; n-fold cross validation.
atp.PKriging<- function(x, ptVgm, nmax=10, ATP = FALSE,longlat=FALSE, showProgress=FALSE, nopar=FALSE) {
  
  N <- nrow(x$areaValues)
  nfold <- N
  indexM <- matrix(sort(x$areaValues[,1]), ncol = 1)
  
  hasCluster <- ataIsClusterEnabled()
  if(showProgress) pb <- txtProgressBar(min=0, max=nrow(indexM), width = 50, style = 3)
  
  
  if(!hasCluster || nopar || nrow(indexM) == 1) {
    estResults <- c()
    for (k in 1:nrow(indexM)) {
      xknown <- subsetDiscreteArea(x, indexM[k,], revSel = TRUE)
      unknown <- subsetDiscreteArea(x, indexM[k,])$discretePoints
      unknown <- data.frame(x = unknown$ptx, y = unknown$pty)
      estResults <- rbind(estResults, atpKriging(xknown, unknown, ptVgm, nmax, longlat, showProgress = FALSE, nopar = TRUE))
      if(showProgress) setTxtProgressBar(pb, k)
    }
  } else {
    bInnerParallel <- ncol(indexM) > 2*nrow(indexM)
    if(bInnerParallel) {
      estResults <- c()
      for (k in 1:nrow(indexM)) {
        xknown <- subsetDiscreteArea(x, indexM[k,], revSel = TRUE)
        unknown <- subsetDiscreteArea(x, indexM[k,])$discretePoints
        estResults <- rbind(estResults, ataKriging(xknown, unknown, ptVgm, nmax, longlat, showProgress = FALSE, nopar = FALSE, clarkAntiLog))
        if(showProgress) setTxtProgressBar(pb, k)
      }
    } else {
      progress <- function(k) if(showProgress) setTxtProgressBar(pb, k)
      estResults <-
        foreach(k = 1:nrow(indexM), .combine = rbind, .options.snow=list(progress=progress),
                .export = c("subsetDiscreteArea","ataCov","calcAreaCentroid","ataKriging","ataKriging.local"),
                .packages = c("sp","gstat","FNN")) %dopar% {
                  xknown <- subsetDiscreteArea(x, indexM[k,], revSel = TRUE)
                  ataKriging(xknown, unknown, ptVgm, nmax, longlat, showProgress = FALSE, nopar = TRUE, clarkAntiLog)
                }
      ataClusterClearObj()
    }
  }
  if(showProgress) close(pb)
  
   estResults$areaId <- x$discretePoints$areaId
  # indx <- match(estResults[,1], x$areaValues[,1])
  # estResults$value <- x$areaValues[indx,4]/x$areaValues[indx,5]
  
  return(estResults)
}


## ataCoKriging.cv: ataCoKriging cross validation. ----
#   nfold: integer; n-fold cross validation.
ataCoKriging.cv <- function(x, unknownVarId, unknownVarCountId, comorbidity, ptVgms, nmax=10, longlat=FALSE, oneCondition=FALSE,
                            meanVal=NULL, auxRatioAdj=TRUE, showProgress=FALSE, nopar=FALSE, clarkAntiLog=FALSE) {
  
  # x = data.all
  # unknownVarId="data.a"
  # unknownVarCountId = "Ya"
  # comorbidity = data.co
  # ptVgms = vg.deconv.cok
  # nmax=Inf
  # longlat=FALSE
  # oneCondition=FALSE
  # meanVal=NULL
  # auxRatioAdj=TRUE
  # showProgress=FALSE
  # nopar=FALSE
  # clarkAntiLog=FALSE
  
  
  N <- nrow(x[[unknownVarId]]$areaValues)
  nfold <- N
  indexM <- matrix(sort(x[[unknownVarId]]$areaValues[,1]), ncol = 1)
  

  hasCluster <- ataIsClusterEnabled()
  if(showProgress) pb <- txtProgressBar(min=0, max=nrow(indexM), width = 50, style = 3)

  xknown <- x

  if(!hasCluster || nopar || nrow(indexM) == 1) {
    estResults <- c()
    for (k in 1:nrow(indexM)) {
      # k=  1 
      xknown[[unknownVarId]] <- subsetDiscreteArea(x[[unknownVarId]], indexM[k,], revSel = TRUE)
      unknown <- subsetDiscreteArea(x[[unknownVarId]], indexM[k,])$discretePoints
      comor <- comorbidity
      estResults <- rbind(estResults,
                          ataCoKriging(xknown, unknownVarId, unknown, comor, ptVgms, nmax, longlat, oneCondition,
                                       meanVal, auxRatioAdj, showProgress = FALSE, nopar = TRUE, clarkAntiLog))
      if(showProgress) setTxtProgressBar(pb, k)
    }
  } else {
    bInnerParallel <- ncol(indexM) > 2*nrow(indexM)
    if(bInnerParallel) {
      estResults <- c()
      for (k in 1:nrow(indexM)) {
        xknown[[unknownVarId]] <- subsetDiscreteArea(x[[unknownVarId]], indexM[k,], revSel = TRUE)
        unknown <- subsetDiscreteArea(x[[unknownVarId]], indexM[k,])$discretePoints
        comor <- comorbidity
        #comor[k, unknownVarCountId] <- NA
        estResults <- rbind(estResults,
                            ataCoKriging(xknown, unknownVarId, unknown, comor, ptVgms, nmax, longlat, oneCondition,
                                         meanVal, auxRatioAdj, showProgress = FALSE, nopar = FALSE, clarkAntiLog))
        if(showProgress) setTxtProgressBar(pb, k)
      }
    } else {
      progress <- function(k) if(showProgress) setTxtProgressBar(pb, k)
      estResults <-
        foreach(k = 1:nrow(indexM), .combine = rbind, .options.snow=list(progress=progress),
                .export = c("crossName","ataCov","calcAreaCentroid","subsetDiscreteArea","ataCoKriging","ataCoKriging.local"),
                .packages = c("sp","gstat","FNN")) %dopar% {
                  xknown[[unknownVarId]] <- subsetDiscreteArea(x[[unknownVarId]], indexM[k,], revSel = TRUE)
                  unknown <- subsetDiscreteArea(x[[unknownVarId]], indexM[k,])$discretePoints
                  ataCoKriging(xknown, unknownVarId, unknown, ptVgms, nmax, longlat, oneCondition, meanVal,
                               auxRatioAdj, showProgress = FALSE, nopar = TRUE, clarkAntiLog)
                }
      ataClusterClearObj()
    }
  }
  if(showProgress) close(pb)

  estResults <- estResults[order(estResults$areaId),]
  indx <- match(estResults[,1], x[[unknownVarId]]$areaValues[,1])
  # estResults$diff <- x[[unknownVarId]]$areaValues[indx,4] - estResults[,4]
  estResults$value <- x[[unknownVarId]]$areaValues[indx,4]

  return(estResults)
}



#------ Plot data
plot.data <- function(data, field, breaks, title, subtitle){
  
  cols = RColorBrewer::brewer.pal(9, "Spectral")
  ggplot(data = data, aes(fill = cut(field, breaks = breaks, include.lowest = TRUE, dig.lab = 5))) + geom_sf(col =  'white', size = 0.1) +
    # scale_fill_brewer(type = "qual",
    #                   palette = "Spectral",
    #                   direction = -1,
    #                   name = "rate/10,000 persons-year") +
    scale_fill_viridis(
      discrete = TRUE,
      direction = 1 ,
      option = 'magma',
      name = '',
      labels = function(breaks) {breaks[is.na(breaks)] <- "No Data"; breaks},
      na.value = "lightgray",
      guide = guide_legend(
        direction = "horizontal",
        reverse = FALSE,
        label.position = "bottom",
        title.position = 'top',
        title.hjust = 0.9,
        label.hjust = 0.9,
        nrow = 1
      )
    ) +
    labs(
      x = NULL,
      y = NULL,
      title = title,
      subtitle = subtitle
    ) +
    theme(
      line = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.background = element_blank(),
      #legend.position = c(0.75, 0.15),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      legend.position = "bottom",
      legend.spacing.x = unit(0.0, 'cm'),
      legend.key.width = unit(2.5,"line"),
      legend.key.height = unit(1,"line"),
      legend.title.align = 0.5,
      legend.text.align = 0.5,
      legend.justification = "center"
    ) +
    coord_sf(datum = NA)
}
