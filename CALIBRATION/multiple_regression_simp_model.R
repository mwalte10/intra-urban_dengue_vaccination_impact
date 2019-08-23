library(fields)
sp9.h <- c()
sp9 <- c()
sp9.l <- c()

setwd('~/Desktop/MANUSCRIPT/CALIBRATION/')
##select which calibration is being run
load('sp9_no_travel.RData')
# load('sp9_travel.RData')
# load('sp9_skew_h.RData')
# load('sp9_skew_l.RData')
sp9 <- sp9.list[[1]]
sp9.h <- sp9.list[[2]]
sp9.l <- sp9.list[[3]]


#set up the betas which are being run
beta_h.vec <- seq(0,0.2, length.out = 10)
beta_l.vec <- seq(0.2,0.6, length.out = 10)
beta.mat <- as.vector(sapply(beta_h.vec, function(i)rep(i,10)))
beta.mat <- cbind(beta.mat, rep(beta_l.vec, 10))
parms <- cbind(beta.mat, rep(0.05,100), sp9, sp9.h, sp9.l)
parm_mat <- as.data.frame(parms)
travel_levels <- unique(parm_mat[,3])


outputs.list <- list()
for(i in 1:length(travel_levels)){
  parms <- parm_mat.list[[i]]
  sp9 <- parms[,4]
  sp9.h <- parms[,5]
  sp9.l <- parms[,6]
  beta_h <- parms[,1]
  beta_l <- parms[,2]
  
  ##USE TPS TO FIND SP9 OF POP FROM BOTH BETAS
  travel_1 <- data.frame(cbind(beta_h, beta_l, sp9))
  names(travel_1) <- c('beta_h', 'beta_l', 'sp9')
  spline_sp9 = Tps(x = as.matrix(travel_1[c("beta_h","beta_l")]),
                   Y = matrix(travel_1["sp9"][,1]),
                   scale.type = "unscaled",lambda = 0)
  x <- seq(0, 0.3, length.out = 50)
  list.x <- list()
  for(j in 1:length(x)){
    list.x[[j]] <- rep(x[j], 50)
  }
  list.x <- unlist(list.x)
  y <- seq(0.4, 0.9, length.out = 50)
  new.data <- data.frame(beta_h = list.x, beta_l = rep(y, 50))
  sp9.predict <- predict(spline_sp9,  new.data)

  ##RESTRICT TO ONLY BETAS WHICH GIVE POP SP9 OF ~0.5
  restrict <- which(sp9.predict > 0.4 & sp9.predict < 0.6)
  sp9.predict <- sp9.predict[which(sp9.predict > 0.4 & sp9.predict < 0.6)]
  new.data <- new.data[restrict,]
  
  ##USE BETAS TO PREDICT SP9.H
  travel_1 <- data.frame(cbind(beta_h, beta_l, sp9.h))
  names(travel_1) <- c('beta_h', 'beta_l', 'sp9.h')
  spline_sp9.h = Tps(x = as.matrix(travel_1[c("beta_h","beta_l")]),
                     Y = matrix(travel_1["sp9.h"][,1]),
                     scale.type = "unscaled",lambda = 0)
  sp9.h.predict <- predict(spline_sp9.h,  new.data)

  ##USE BETAS TO PREDICT SP9.L
  travel_1 <- data.frame(cbind(beta_h, beta_l, sp9.l))
  names(travel_1) <- c('beta_h', 'beta_l', 'sp9.l')
  spline_sp9.l = Tps(x = as.matrix(travel_1[c("beta_h","beta_l")]),
                     Y = matrix(travel_1["sp9.l"][,1]),
                     scale.type = "unscaled",lambda = 0)
  sp9.l.predict <- predict(spline_sp9.l,  new.data)
 
  ##FIND BETAS WHICH MINIMIZE THE DIFFERENCES BETWEEN ALL WANTED SP9'S
  sums.mat <- cbind((sp9.predict - 0.5), (sp9.h.predict - 0.2), (sp9.l.predict - 0.8))
  sums.mat <- data.frame(abs(sums.mat))
  sums <- sums.mat[,1] + sums.mat[,2] + sums.mat[,3]
  output <- c(new.data[which.min(sums),], sp9.predict[which.min(sums)], 
              sp9.l.predict[which.min(sums),], sp9.h.predict[which.min(sums),])
  output <- unlist(output)
  names(output) <- c('beta_h', 'beta_l', 'sp9', 'sp9.l', 'sp9.h')
  output <- as.vector(output)
  outputs.list[[i]] <- output
}




