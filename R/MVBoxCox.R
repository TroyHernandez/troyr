#MVBoxCox.R

############################################################
##Optimal lambda values as $lambda
MVBoxCox <- function(xmat, lamseq = seq(-2, 2, .05), NumInstances = 4,
                     Instance = 1, output = T) {
  lambda <- c(rep(0, ncol(xmat)))
  #vector telling if small positive value (1.01) was added
  smpv <- rep(0, length(lamseq))
  bcscore <- rep(0, length(lamseq))
  n <- nrow(xmat)
  xt <- xmat / Inf
  
  for(j in 1:ncol(xmat)){
    xtemp <- xmat[, j]
    if(sd(xtemp) != 0 && sd(xtemp + (-min(xtemp) + 1.01)) != 0){
      if(min(xtemp) <= 1){
        smpv[j] = (-min(xtemp) + 1.01)
        xtemp <- xtemp + smpv[j]
      }
      for(i in 1:length(lamseq)){
        if(abs(lamseq[i]) < 1e-10){
          bct <- log(xtemp)
        }else{
          bct <- (xtemp ^ (lamseq[i]) - 1) / lamseq[i]
        }
        bcscore[i] <- -n / 2 * log(1 / n * sum((bct - mean(bct)) ^ 2 )) + 
                        (lamseq[i] - 1) * sum(log(xtemp))
      }
      badscores <- which(bcscore == Inf | bcscore == -Inf)
      if(length(badscores) > 0){
        lamseq <- lamseq[-badscores]
        bcscore <- bcscore[-badscores]
      }
      #       cat("length bc seq:",length(bcscore),"\n")
      #plot(bcscore)
      if((bcscore[which.max(bcscore)] == bcscore[1] |
            bcscore[which.max(bcscore)] == bcscore[length(bcscore)]) &
           Instance < NumInstances & length(bcscore) > 20){
        TempInstance <- Instance + 1
        tmp <- MVBoxCox(xmat = cbind(xmat[, j]),
                        lamseq = seq(min(lamseq * 2), max(lamseq * 2),
                                     2 * (lamseq[2] - lamseq[1])),
                        NumInstances, Instance = TempInstance)
        lambda[j] <- tmp$lambda
        xt[, j] <- tmp$xt
        #Go backwards if it results in single value for all entries
        if(sd(xt[, j]) == 0){
          if(Instance == 1){
            lambda[j] <- 2
            xt[, j] <- (xtemp ^ (lambda[j]) - 1) / lambda[j]
          }else{
            TempInstance <- NumInstances
            tmp <- MVBoxCox(cbind(xmat[, j]),
                            lamseq = seq(min(lamseq / 2), max(lamseq / 2),
                                         (lamseq[2] - lamseq[1]) / 2),
                            NumInstances, Instance = TempInstance)
            lambda[j] <- tmp$lambda
            xt[, j] <- tmp$xt
          }
        }
        #Only output when it's the final iteration
        if(Instance == 1 & output == T & round(j / 1000) == j / 1000){
          cat("Segment ", j, ": ", lambda[j], "\n")
        }
      }else{

        lambda[j] <- lamseq[which.max(bcscore)]
        if(abs(lambda[j]) < 1e-10){
          xt[, j] <- log(xtemp)
        }else{
          xt[, j] <- (xtemp ^ (lambda[j]) - 1) / lambda[j]
        }
        if(sd(xt[, j]) == 0){
          if(Instance == 1){
            lambda[j] <- 1
            xt[, j] <- (xtemp ^ (lambda[j]) - 1) / lambda[j]
          }else{
            TempInstance <- NumInstances
            tmp <- MVBoxCox(cbind(x[, j]),
                            lamseq = seq(min(lamseq / 2), max(lamseq / 2),
                                         (lamseq[2] - lamseq[1]) / 2),
                            NumInstances, Instance = TempInstance)
            lambda[j] <- tmp$lambda
            xt[, j] <- tmp$xt
          }
        }
        #Only output when it's the final iteration
        if(Instance == 1 & output == T & round(j / 1000) == j / 1000){
          cat("Segment ", j, ": ", lambda[j], "\n")
        }
      }
    }else{
      if(output == T & round(j / 1000) == j / 1000){
        cat("Segment ", j, "is uniform.\n")
      }
    }
  }
  list(xt = xt, lambda = lambda, smpv = smpv)
}
