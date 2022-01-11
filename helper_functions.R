################################################################################################  Libraries  #####################################################################################################
library(ggplot2)
library(mgcv)
library(dplyr)
library(reshape2)
library(devtools)
library(matrixStats)
library(ebbr)
library(BiocParallel)
library(furniture)
################################################################################################    #####################################################################################################
FindIndices <- function(data, trainvar, trainlevel) {
  #
  # Constructs vector of rownames or indices based on subsetting parameters
  # Helpful for GroupModel if interested in training models based
  #  on subset of data
  #
  # Args:
  #  data: full dataframe
  #  trainvar: string of variable name (must be factor variable) 
  #  trainlevel: string of factor level of interest
  #
  # Returns:
  #  Vector of rownames or indices
  #
  index      <- which(data[[trainvar]] == trainlevel, 
                      arr.ind = TRUE)
  index      <- as.vector(index)
  return(index)
}

################################################################################################    #####################################################################################################
BuildDict <- function(covar.data, ref.batch=NULL) {
  #
  # Extract and store batch characteristics for later use
  # Returns dictionary of batch characteristics
  #
  # Args:
  #  covar.data: dataframe of covariates/batch assosciation 
  #
  # Returns:
  #  batches: row indicies split by batch assosciation 
  #  n.batch: number of batches
  #  n.array: number of observations
  #  n.batches: number of observations per batch
  #
  dict    <- list()
  
  batches <- lapply(levels(covar.data[["STUDY"]]), 
                    function(x)which(covar.data[["STUDY"]] == x))

  dict[["batches"]]     <- batches 
  dict[["n.batch"]]     <- nlevels(covar.data[["STUDY"]])
  dict[["n.array"]]     <- nrow(covar.data)
  dict[["n.batches"]]   <- sapply(batches, length)
  return(dict)
}

 ################################################################################################    #####################################################################################################

BuildFormula <- function(covar.data, smooth.terms = NULL, k.val = NULL) {
  #
  # Constructs formula for GAM model
  # Allows specification of smooth predictors along with smoothing parameters
  #
  # Args:
  #  covar.data: Dataframe of covariate data
  #  smooth.terms: Vector of column names for smoothing
  #  k.val:  vector of knot values for smooth column
  # (note) length of k.val must equal length of smooth.terms
  #
  # Returns:
  #  A string to be coerced into formula class for GAM regression
  #
  formstring <- " ~  -1 + STUDY"
  cols <- colnames(covar.data)
  cols <- cols[! cols == "STUDY" ]
  if(!is.null(smooth.terms)) {
    for(i in 1:length(smooth.terms)) {
      smth         <- paste("s", "(",  smooth.terms[i], ",", "k=", 
                            toString(k.val[i]),  ")", sep = "" )
      formstring   <- paste(formstring, smth, sep  = " + ")
    }
    
    cols <- cols[!cols %in% smooth.terms]
  }
  for(i in 1:length(cols)) {
    formstring <- paste(formstring, cols[i], sep= " + " )
  }
  return(formstring)
}

################################################################################################    #####################################################################################################


FitModel <- function(feature.data, covar.data, model.formula, verbose = FALSE) {
  #
  # Fits a GAM model for each feature
  #
  # Args:
  #  feature.data: Imaging data
  #  covar.data: Corresponding covariate data
  #  training.indices: (OPTIONAL) Vector of indices to train models
  #  verbose: (OPTIONAL) Extra output
  #
  # Returns:
  #  A list of GAM models 
  #
  imcols  <- colnames(feature.data)
  covcols <- colnames(covar.data)
  data       <- as.data.frame(cbind(feature.data, covar.data))
  colnames(data) <- c(imcols, covcols)
  feature.data <- as.data.frame(feature.data)
  covar.data <-as.data.frame(covar.data)
  checkzero  <- which(feature.data <= 0, arr.ind = TRUE)[,1]
  
  if(length(checkzero) > 0  ) {
    if(verbose) {
      cat( "[ComGam] ", "Removed", length(checkzero), 
           "rows with non-positive response values from training" , "\n")
    }
    data     <- data[ -checkzero, ]
    
  }
  modlist    <- list()
  for (j in 1:length(imcols)) {
    if(verbose) {
      cat("[ComGam] ", "Fit model", j, "out of", length(imcols), "\n")
    }
    gammod              <- gam(as.formula(paste(imcols[j], model.formula, sep = "")) ,   
                               method = "REML", data = data)
    modlist[[j]]        <- gammod
  }
  names(modlist)        <- imcols

  return(modlist)
}

################################################################################################    #####################################################################################################

StanAcrossFeatures <- function(feature.data, covar.data, models.list, data.dict) {
  #
  # Takes data to be harmonized and standardizes the values
  # Removes and preserves biological variance to isolate just batch variance
  # 
  # Args:
  #  dat: Just features to be harmonized (e.x Imaging data/Biomarker data)
  #  covar.data: Dataframe of covariate data
  #  models.list: List of fitted GAM models
  #  data.dict: Dictionary with batch characteristics
  #
  # Returns:
  #  std.data: standardized features data
  #  design:   design matrix
  #  mod.names: names of harmonization features
  #  stand.mean: mean of each feature weighted by batches
  #  var.pooled: variance of each feature
  #
  n.batch   <- data.dict[["n.batch"]]
  n.batches <- data.dict[["n.batches"]]
  n.array   <- data.dict[["n.array"]]
  feature.data       <- t(feature.data)
  design    <- predict.gam(models.list[[1]], 
                           type = "lpmatrix", 
                           newdata = covar.data)
  
  ## dat feature names
  mod.names  <- names(models.list)
  
  ## covariate names
  coef.names <- names(models.list[[1]][["coefficients"]])
  
  B.hat <- data.frame() 
  for(i in 1:length(models.list)) {
    #extract coeffs
    B.hat.row <- models.list[[i]][["coefficients"]]
    B.hat     <- rbind(B.hat, B.hat.row)
  }
  ## build coeff matrix
  colnames(B.hat) <- coef.names
  rownames(B.hat) <- mod.names
  B.hat           <- as.matrix(B.hat)
  B.hat           <- t(B.hat)
  predicted       <- design %*% B.hat
  predicted       <- t(predicted)
  ## find mean value for features using weighted average
  
  grand.mean <- crossprod(n.batches / n.array, B.hat[1 : n.batch, ])
  stand.mean <- crossprod(grand.mean, t(rep(1, n.array)))
  stand.mean1 <- stand.mean
  
  # add bio variance to stand.mean
  design.tmp <- design
  design.tmp[ ,c(1 : n.batch)] <- 0
  bio.var    <- design.tmp %*% B.hat
  bio.var    <- t(bio.var)
  stand.mean <- stand.mean + bio.var
  
  ## variance
  var.pooled <- (feature.data - predicted) ^ 2
  var.adj    <- as.matrix(rep(1, n.array)) / n.array
  var.pooled <- var.pooled %*% var.adj
  var.pooled <- sqrt(var.pooled)
  var.pooled1 <- var.pooled
  var.pooled <- var.pooled %*% as.matrix(t(rep(1, n.array)))
  std.data   <- (feature.data - stand.mean) / var.pooled
  return(list("std.data"   = std.data, 
              "design"     = design, 
              "mod.names"  = mod.names,
              "stand.mean" = stand.mean,
              "var.pooled" = var.pooled,
              "grand.mean" = grand.mean,
              "var.pooled1" = var.pooled1))
}



StanAcrossFeaturesREFBATCH <- function(feature.data, covar.data, models.list, data.dict) {
  #
  # Takes data to be harmonized and standardizes the values
  # Removes and preserves biological variance to isolate just batch variance
  # 
  # Args:
  #  dat: Just features to be harmonized (e.x Imaging data/Biomarker data)
  #  covar.data: Dataframe of covariate data
  #  models.list: List of fitted GAM models
  #  data.dict: Dictionary with batch characteristics
  #
  # Returns:
  #  std.data: standardized features data
  #  design:   design matrix
  #  mod.names: names of harmonization features
  #  stand.mean: mean of each feature weighted by batches
  #  var.pooled: variance of each feature
  #
  n.batch   <- data.dict[["n.batch"]]
  n.batches <- data.dict[["n.batches"]]
  n.array   <- data.dict[["n.array"]]
  feature.data       <- t(feature.data)
  design    <- predict.gam(models.list[[1]], 
                           type = "lpmatrix", 
                           newdata = covar.data)
  
  ## dat feature names
  mod.names  <- names(models.list)
  
  ## covariate names
  coef.names <- names(models.list[[1]][["coefficients"]])
  
  B.hat <- data.frame() 
  for(i in 1:length(models.list)) {
    #extract coeffs
    B.hat.row <- models.list[[i]][["coefficients"]]
    B.hat     <- rbind(B.hat, B.hat.row)
  }
  ## build coeff matrix
  colnames(B.hat) <- coef.names
  rownames(B.hat) <- mod.names
  B.hat           <- as.matrix(B.hat)
  B.hat           <- t(B.hat)
  predicted       <- design %*% B.hat
  predicted       <- t(predicted)
  ## find mean value for features using weighted average
  
  grand.mean <- crossprod(n.batches / n.array, B.hat[1 : n.batch, ])
  stand.mean <- crossprod(grand.mean, t(rep(1, n.array)))
  
  # add bio variance to stand.mean
  design.tmp <- design
  design.tmp[ ,c(1 : n.batch)] <- 0
  bio.var    <- design.tmp %*% B.hat
  bio.var    <- t(bio.var)
  stand.mean <- stand.mean + bio.var
  
  ## variance
  var.pooled <- (feature.data - predicted) ^ 2
  var.adj    <- as.matrix(rep(1, n.array)) / n.array
  var.pooled <- var.pooled %*% var.adj
  var.pooled <- sqrt(var.pooled)
  var.pooled <- var.pooled %*% as.matrix(t(rep(1, n.array)))
  std.data   <- (feature.data - stand.mean) / var.pooled
  return(list("std.data"   = std.data, 
              "design"     = design, 
              "mod.names"  = mod.names,
              "stand.mean" = stand.mean,
              "var.pooled" = var.pooled))
}







################################################################################################    #####################################################################################################

CalcGammaDelta <- function(stan.dict, data.dict) {
  #
  # Calculates the Mean/Variance (Shift/Scale) of batch effect
  # 
  # Args:
  # stan.dict: Dictionary with standardized data and other values
  # Output of StanAcrossFeatures
  #
  # Returns:
  # gamma.hat: Shift value per batch per feature
  # delta.hat: Scale value per batch per feature
  #
  batches   <- data.dict[["batches"]]
  n.batch   <- data.dict[["n.batch"]]
  design    <- stan.dict[["design"]]
  std.data  <- stan.dict[["std.data"]]
  mod.names <- stan.dict[["mod.names"]]
  batch.mod <- design[, 1:n.batch]

    ## shift (mean)
  gamma.hat <- tcrossprod(solve(crossprod(batch.mod, batch.mod)), batch.mod)
  gamma.hat <- tcrossprod(gamma.hat, std.data)
  
  ## scale (variance)
  delta.hat <- data.frame()
  for(i in batches) {
    delta.hat <- rbind(delta.hat, 
                       rowVars(std.data, cols = i, na.rm = TRUE))
  }
  colnames(delta.hat) <- mod.names
  rownames(delta.hat) <- colnames(batch.mod)
  return(list("gamma.hat" = gamma.hat,
              "delta.hat" = delta.hat,
              "stan.data" = std.data))
}


################################################################################################    #####################################################################################################

ModelDiagnostics <- function(mod.list) {
  #
  # Pulls information on model fitting for each features smooth terms 
  # Useful for checking model fit accuracy in harmonizationn
  #
  # Args:
  #  mod.list: Output from FitModel
  #
  # Returns:
  #  matrix with diagnostic information 
  #
  diag.list  <-list()
  mod.names  <- names(mod.list)
  diagnos.df <- data.frame(matrix(ncol = 7))
  for (i in 1:length(mod.list)) {
    gam.o <- capture.output(gam.check(mod.list[[i]]))
    dev.off()
    for(j in 13:length(gam.o)) {
      o.vals.split <- c(strsplit(gam.o[j], split = " "))
      o.vals <- o.vals.split[[1]]
      o.vals <- o.vals[o.vals != ""]
      if("---" %in% o.vals | "Signif." %in% o.vals) {
        #skip
      } else {
        if(length(o.vals) == 6) {
          row.val <- c(mod.names[i], o.vals)
        } else {
          row.val <- c(mod.names[i], o.vals, "")
        }
      }
      
      diagnos.df <- rbind(diagnos.df, as.character(row.val))
    }
  }
  diagnos.df <- as.data.frame(diagnos.df)
  colnames(diagnos.df) <- c("feature",
                            "smooth",
                            "k",
                            "edf",
                            "k.index",
                            "p.value",
                            "signif code (alpha = 0.05")
  diagnos.df <- unique(diagnos.df[2:nrow(diagnos.df), ])
  return(diagnos.df)
}
################################################################################################    #####################################################################################################

ApplyGammaDelta <- function(stan.dict, site.params, data.dict) {
  #
  # Takes gamma and delta values and applies them to the data
  # Reintroduces biological variance after gamma/delta correction
  #
  # Args:
  #  stan.dict: StanAcrossFeatures output
  #  site.params: CalcGammaDelta output
  #  data.dict: BuildDict Output
  #
  # Returns:
  #  Harmonized data
  #
  ## extract all variables
  std.data   <- stan.dict[["std.data"]]
  design     <- stan.dict[["design"]]
  mod.names  <- stan.dict[["mod.names"]]
  stand.mean <- stan.dict[["stand.mean"]]
  var.pooled <- stan.dict[["var.pooled"]]
  batches    <- data.dict[["batches"]]
  n.batches  <- data.dict[["n.batches"]]
  n.batch    <- data.dict[["n.batch"]]
  gamma.hat  <- site.params[["gamma.hat"]]
  delta.hat  <- as.matrix(site.params[["delta.hat"]])
  batch.mod  <- design[, 1:n.batch]
  ## apply gamma/delta
  std.adj <- std.data
  j <- 1 
  for(i in batches){ 
    shft <- std.adj[ ,i] - t(batch.mod[i, ] %*% gamma.hat)
    scl  <- tcrossprod(sqrt(delta.hat[j, ]), rep(1, n.batches[j]))
    std.adj[,i] <- shft / scl
    j <- j + 1
  }
  
  #reintroduce biological mean and variance
  std.adj <- (std.adj * var.pooled) + stand.mean
  return(std.adj)
}

################################################################################################    #####################################################################################################

NLAdjustment  <- function(covar.data, stan.dict, ref.cohort, k.val.nlt) {
  t.stan            <- as.data.frame(t(stan.dict[["std.data"]]))
  t.std.mean        <- as.data.frame(t(stan.dict[["stand.mean"]]))
  t.std.var         <- as.data.frame(t(stan.dict[["var.pooled"]]))
  feats             <- colnames(t.stan)
  batch.var         <- factor(covar.data[["STUDY"]], levels = unique(covar.data[["STUDY"]]))
  t.std.mean$batch.var   <- batch.var
  t.std.var$batch.var    <- batch.var
  t.stan[["batch.var"]]  <- batch.var
  split.batch       <- split(t.stan, t.stan[["batch.var"]])
  batch.names       <- names(split.batch)
  batch.names       <- batch.names[!batch.names %in% ref.cohort]
  pair.list         <- list()
  ref.df            <- split.batch[[ref.cohort]]
  len.ref           <- nrow(ref.df)
  for(i in 1:length(batch.names)) {
    element.name <- batch.names[i]
    batch.df.red <- split.batch[[batch.names[i]]]
    ref.df.red   <- ref.df
    if(nrow(batch.df.red) >= len.ref) {
      batch.df.red <- batch.df.red[1:len.ref, ]
    } else {
      ref.df.red   <- ref.df.red[1:nrow(batch.df.red),]  
    }
     batch.df.red[["batch.var"]] <- NULL
     ref.df.red[["batch.var"]]   <- NULL
     colnames(batch.df.red) <- paste(colnames(batch.df.red), batch.names[i], sep="_")
     colnames(ref.df.red)   <- paste(colnames(ref.df.red), ref.cohort, sep="_")
     full.map.df            <- cbind(ref.df.red, batch.df.red)
     pair.list[[element.name]] <- full.map.df
     
  }
  all.models <- list()
  all.data   <- list()
  for(i in 1:length(pair.list)) {
    dfname    <- names(pair.list[i])
    clnames   <- c()
    mod.list   <- list()
    dat <- as.data.frame(pair.list[[i]])
    newdata <- as.data.frame(split.batch[[batch.names[i]]])
    pred.frame <- data.frame(matrix(nrow = nrow(newdata)))
    for(j in 1:length(feats)) {
      feat <- feats[j]
      clnames <- append(clnames, feat)
      outcome.name      <- paste(feat, ref.cohort, sep = "_")
      pred.name         <- paste(feat, batch.names[i], sep = "_")
      gamdf <- data.frame("outcome" = dat[[outcome.name]],
                           "pred" = dat[[pred.name]])
      gam.mod      <- gam(outcome ~ s(pred, k = k.val.nlt), data = gamdf)
      pred.data    <- data.frame(newdata[feat])
      colnames(pred.data) <- c("pred")
      pred.feat    <- predict.gam(gam.mod, newdata = pred.data, type = "response")
      pred.frame   <- cbind(pred.frame, pred.feat)
      mod.list[[feat]] <- gam.mod
    }
    pred.frame[,1] <- NULL
    colnames(pred.frame) <- clnames
    list.name <- paste(ref.cohort, batch.names[i], sep = " ~ ")
    all.data[[dfname]] <- pred.frame
    all.models[[list.name]] <- mod.list
  }
  refco <- split.batch[[ref.cohort]]
  refco$batch.var <-NULL
  t.std.var$batch.var <-NULL
  t.std.mean$batch.var <-NULL
  all.data[[ref.cohort]] <- refco
  batch.order <- as.character(unique(covar.data[["STUDY"]]))
  all.data <- all.data[match(batch.order, names(all.data))]
  fulldata <- do.call(rbind, all.data)
  fulldata <- fulldata * t.std.var
  fulldata <- fulldata + t.std.mean
  return(fulldata)
}


getEbEstimators <- function(naiveEstimators,
                            s.data, 
                            dataDict,
                            parametric=TRUE, 
                            mean.only=FALSE,
                            BPPARAM=bpparam("SerialParam")
){
  gamma.hat <- naiveEstimators[["gamma.hat"]]
  delta.hat <- naiveEstimators[["delta.hat"]]
  batches   <- dataDict$batches
  n.batch   <- dataDict$n.batch
  ref.batch <- dataDict$ref.batch
  ref <- dataDict$ref
  
  .getParametricEstimators <- function(){
    gamma.star <- delta.star <- NULL
    for (i in 1:n.batch){
      if (mean.only){
        gamma.star <- rbind(gamma.star, postmean(gamma.hat[i,], gamma.bar[i], 1, 1, t2[i]))
        delta.star <- rbind(delta.star, rep(1, nrow(s.data)))
      } else {
        temp <- it.sol(s.data[,batches[[i]]],
                       gamma.hat[i,],
                       delta.hat[i,],
                       gamma.bar[i],
                       t2[i],
                       a.prior[i],
                       b.prior[i])
        gamma.star <- rbind(gamma.star,temp[1,])
        delta.star <- rbind(delta.star,temp[2,])
      }
    }
    rownames(gamma.star) <- rownames(delta.star) <- names(batches)
    out <- list(gamma.star=gamma.star, delta.star=delta.star)
    return(out)
  }
  
  .getNonParametricEstimators <- function(BPPARAM=bpparam("SerialParam")){
    gamma.star <- delta.star <- NULL
    
    results <- bplapply(1:n.batch, function(i){
      if (mean.only){
        delta.hat[i, ] = 1
      }
      temp <- int.eprior(as.matrix(s.data[, batches[[i]]]),
                         gamma.hat[i,],
                         delta.hat[i,])
      return(temp)
    }, BPPARAM = BPPARAM)
    gamma.star <- lapply(results, function(x) x[1,])
    delta.star <- lapply(results, function(x) x[2,])
    gamma.star <- do.call("rbind",gamma.star)
    delta.star <- do.call("rbind",delta.star)
    #for (i in 1:n.batch){
    #    gamma.star <- rbind(gamma.star,temp[1,])
    #    delta.star <- rbind(delta.star,temp[2,])
    #}
    rownames(gamma.star) <- rownames(delta.star) <- names(batches)
    out <- list(gamma.star=gamma.star, delta.star=delta.star)
    return(out)
  }
  
  gamma.bar <- rowMeans(gamma.hat, na.rm=TRUE)
  t2 <- rowVars(gamma.hat, na.rm=TRUE)
  names(t2) <- rownames(gamma.hat)
  a.prior <- apriorMat(delta.hat)
  b.prior <- bpriorMat(delta.hat)
  if (parametric){
    temp <- .getParametricEstimators()
  } else {
    temp <- .getNonParametricEstimators(BPPARAM=BPPARAM)
  }
  if(!is.null(ref.batch)){
    temp[["gamma.star"]][ref,] <- 0  ## set reference batch mean equal to 0
    temp[["delta.star"]][ref,] <- 1  ## set reference batch variance equal to 1
  }
  out <- list()
  out[["gamma.star"]] <- temp[["gamma.star"]]
  out[["delta.star"]] <- temp[["delta.star"]]
  out[["gamma.bar"]] <- gamma.bar
  out[["t2"]] <- t2
  out[["a.prior"]] <- a.prior
  out[["b.prior"]] <- b.prior
  return(out)
}




 ################################################################################################    #####################################################################################################

################################################################################################    #####################################################################################################

LongSelect <- function(cs, fulldata, id = "RID") {
  cs.id <- as.list(cs[[id]])
  fd.id <- as.list(fulldata[[id]])
  newlist <- as.numeric(intersect(unique(cs.id), fd.id))
  rowvec <- c()
  for(i in 1:nrow(fulldata)) {
    if(fulldata["RID"][i,] %in% newlist)
      rowvec <- append(rowvec, i)
  }
  fullsubset <- fulldata[rowvec,]
  fullsubset <- OrderData(fullsubset)
  return(fullsubset)
}

aprior <- function(delta.hat){
m=mean(delta.hat)
s2=var(delta.hat)
return((2*s2+m^2)/s2)
}

bprior <- function(delta.hat){
  m=mean(delta.hat)
  s2=var(delta.hat)
  return((m*s2+m^3)/s2)
}

apriorMat <- function(delta.hat) {
  m  <- rowMeans2(delta.hat)
  s2 <- rowVars(delta.hat)
  out <- (2*s2+m^2)/s2
  names(out) <- rownames(delta.hat)
  return(out)
}

bpriorMat <- function(delta.hat) {
  m <- rowMeans2(delta.hat)
  s2 <- rowVars(delta.hat)
  out <- (m*s2+m^3)/s2
  names(out) <- rownames(delta.hat)
  return(out)
}

postmean <- function(g.hat, g.bar, n, d.star, t2){
  (t2*n*g.hat+d.star*g.bar)/(t2*n+d.star)
}

postvar <- function(sum2, n, a, b){
  (.5*sum2+b)/(n/2+a-1)
}

# Helper function for parametric adjustements:
it.sol  <- function(sdat, g.hat, d.hat, g.bar, t2, a, b, conv=.0001){
  #n <- apply(!is.na(sdat),1,sum)
  n <- rowSums(!is.na(sdat))
  g.old  <- g.hat
  d.old  <- d.hat
  change <- 1
  count  <- 0
  ones <- rep(1,ncol(sdat))
  
  while(change>conv){
    g.new  <- postmean(g.hat,g.bar,n,d.old,t2)
    sum2   <- rowSums2((sdat-tcrossprod(g.new, ones))^2, na.rm=TRUE)
    d.new  <- postvar(sum2,n,a,b)
    change <- max(abs(g.new-g.old)/g.old,abs(d.new-d.old)/d.old)
    g.old <- g.new
    d.old <- d.new
    count <- count+1
  }
  adjust <- rbind(g.new, d.new)
  rownames(adjust) <- c("g.star","d.star")
  return(adjust)
}



# Helper function for non-parametric adjustements:
int.eprior <- function(sdat, g.hat, d.hat){
  g.star <- d.star <- NULL
  r <- nrow(sdat)
  for(i in 1:r){
    g <- g.hat[-i]
    d <- d.hat[-i]		
    x <- sdat[i,!is.na(sdat[i,])]
    n <- length(x)
    j <- numeric(n)+1
    dat <- matrix(as.numeric(x), length(g), n, byrow=TRUE)
    resid2 <- (dat-g)^2
    sum2 <- resid2 %*% j
    LH <- 1/(2*pi*d)^(n/2)*exp(-sum2/(2*d))
    LH[LH=="NaN"]=0
    g.star <- c(g.star, sum(g*LH)/sum(LH))
    d.star <- c(d.star, sum(d*LH)/sum(LH))
  }
  adjust <- rbind(g.star,d.star)
  rownames(adjust) <- c("g.star","d.star")
  return(adjust)
} 


getNaiveEstimators <- function(s.data,
                               dataDict,
                               hasNAs,
                               mean.only
){
  batch.design <- dataDict$batch.design
  batches <- dataDict$batches
  if (!hasNAs){
    gamma.hat <- tcrossprod(solve(crossprod(batch.design, batch.design)), batch.design)
    gamma.hat <- tcrossprod(gamma.hat, s.data)
  } else{
    gamma.hat <- apply(s.data, 1, .betaNA, batch.design) 
  }
  delta.hat <- NULL
  for (i in dataDict$batches){
    if (mean.only){
      delta.hat <- rbind(delta.hat,rep(1,nrow(s.data))) 
    } else {
      delta.hat <- rbind(delta.hat,rowVars(s.data, cols=i, na.rm=TRUE))
    }    
  }
  colnames(gamma.hat)  <- colnames(delta.hat) <- rownames(s.data)
  rownames(gamma.hat)  <- rownames(delta.hat) <- names(batches)
  delta.hat[delta.hat==0] <- 1
  out <- list(gamma.hat=gamma.hat, delta.hat=delta.hat)
  return(out)
}




getNonEbEstimators <- function(naiveEstimators, dataDict){
  out <- list()
  out[["gamma.star"]] <- naiveEstimators[["gamma.hat"]]
  out[["delta.star"]] <- naiveEstimators[["delta.hat"]]
  out[["gamma.bar"]]  <- NULL
  out[["t2"]] <- NULL
  out[["a.prior"]] <- NULL
  out[["b.prior"]] <- NULL
  ref.batch=dataDict$ref.batch
  ref=dataDict$ref
  if(!is.null(ref.batch)){
    out[["gamma.star"]][ref,] <- 0  ## set reference batch mean equal to 0
    out[["delta.star"]][ref,] <- 1  ## set reference batch variance equal to 1
  }
  return(out)
}


getEbEstimators <- function(naiveEstimators,
                            s.data, 
                            dataDict,
                            parametric=TRUE, 
                            mean.only=FALSE,
                            BPPARAM=bpparam("SerialParam")
){
  gamma.hat <- naiveEstimators[["gamma.hat"]]
  delta.hat <- naiveEstimators[["delta.hat"]]
  batches   <- dataDict$batches
  n.batch   <- dataDict$n.batch
  ref.batch <- dataDict$ref.batch
  ref <- dataDict$ref
  
  .getParametricEstimators <- function(){
    gamma.star <- delta.star <- NULL
    for (i in 1:n.batch){
      if (mean.only){
        gamma.star <- rbind(gamma.star, postmean(gamma.hat[i,], gamma.bar[i], 1, 1, t2[i]))
        delta.star <- rbind(delta.star, rep(1, nrow(s.data)))
      } else {
        temp <- it.sol(s.data[,batches[[i]]],
                       gamma.hat[i,],
                       delta.hat[i,],
                       gamma.bar[i],
                       t2[i],
                       a.prior[i],
                       b.prior[i])
        gamma.star <- rbind(gamma.star,temp[1,])
        delta.star <- rbind(delta.star,temp[2,])
      }
    }
    rownames(gamma.star) <- rownames(delta.star) <- names(batches)
    out <- list(gamma.star=gamma.star, delta.star=delta.star)
    return(out)
  }
  
  .getNonParametricEstimators <- function(BPPARAM=bpparam("SerialParam")){
    gamma.star <- delta.star <- NULL
    
    results <- bplapply(1:n.batch, function(i){
      if (mean.only){
        delta.hat[i, ] = 1
      }
      temp <- int.eprior(as.matrix(s.data[, batches[[i]]]),
                         gamma.hat[i,],
                         delta.hat[i,])
      return(temp)
    }, BPPARAM = BPPARAM)
    gamma.star <- lapply(results, function(x) x[1,])
    delta.star <- lapply(results, function(x) x[2,])
    gamma.star <- do.call("rbind",gamma.star)
    delta.star <- do.call("rbind",delta.star)
    #for (i in 1:n.batch){
    #    gamma.star <- rbind(gamma.star,temp[1,])
    #    delta.star <- rbind(delta.star,temp[2,])
    #}
    rownames(gamma.star) <- rownames(delta.star) <- names(batches)
    out <- list(gamma.star=gamma.star, delta.star=delta.star)
    return(out)
  }
  
  gamma.bar <- rowMeans(gamma.hat, na.rm=TRUE)
  t2 <- rowVars(gamma.hat, na.rm=TRUE)
  names(t2) <- rownames(gamma.hat)
  a.prior <- apriorMat(delta.hat)
  b.prior <- bpriorMat(delta.hat)
  if (parametric){
    temp <- .getParametricEstimators()
  } else {
    temp <- .getNonParametricEstimators(BPPARAM=BPPARAM)
  }
  if(!is.null(ref.batch)){
    temp[["gamma.star"]][ref,] <- 0  ## set reference batch mean equal to 0
    temp[["delta.star"]][ref,] <- 1  ## set reference batch variance equal to 1
  }
  out <- list()
  out[["gamma.star"]] <- temp[["gamma.star"]]
  out[["delta.star"]] <- temp[["delta.star"]]
  out[["gamma.bar"]] <- gamma.bar
  out[["t2"]] <- t2
  out[["a.prior"]] <- a.prior
  out[["b.prior"]] <- b.prior
  return(out)
}


ApoeIndADNI <- function(df) {
  indvar <- c()
  drop.rows <- c()
  na.drop <- which(is.na(df$APGEN1), arr.ind = TRUE)
  if(length(na.drop) != 0) {
    df <- df[-na.drop,]
  }
  for( i in 1:nrow(df)) {
    a <- df["APGEN1"][i,]
    b <- df["APGEN2"][i,]
    geneval <- c(a, b)
    #if(a == 2 | b == 2) {
    # drop.rows <- append(drop.rows, i)
    #}
    if(sum(geneval == 4) == 0) { 
      indvar <- append(indvar, 0)
    } else if(sum(geneval == 4) == 1){
      indvar <- append(indvar, 1)
    } else {
      indvar <- append(indvar, 2)
    }
  }
  df["ApoE4Ind"] <- indvar
  return(df)
}

ApoeInd2ADNI <- function(df) {
  indvar <- c()
  drop.rows <- c()
  na.drop <- which(is.na(df$APGEN1), arr.ind = TRUE)
  if(length(na.drop) != 0) {
    df <- df[-na.drop,]
  }
  for( i in 1:nrow(df)) {
    a <- df["APGEN1"][i,]
    b <- df["APGEN2"][i,]
    geneval <- c(a, b)
    if(sum(geneval == 2) == 0) { 
      indvar <- append(indvar, 0)
    } else if(sum(geneval == 2) == 1){
      indvar <- append(indvar, 1)
    } else {
      indvar <- append(indvar, 2)
    }
  }
  df["ApoE2Ind"] <- indvar
  return(df)
}
