require("glmnet")

covTestBinom <- function (x, y, weights, maxp = "num.nonzero"){
  # First build the model
  set.seed(1)
  glm.mod <- cv.glmnet(x = x, y = y, family = "binomial", weights = weights)
  
  glmobj = glm.mod$glmnet.fit
  n = nrow(x)
  p = ncol(x)
  my = mean(y)
  
  # Determine non-zero betas and at which betas they enter
  beta.mat <- glmobj$beta != 0
  beta.mat.enter <- rep(0, nrow(beta.mat))
  for(i in 1:nrow(beta.mat)){
    beta.mat.enter[i] <- which(beta.mat[i, ] == TRUE)[1]
  }
  lambda.ind <- as.numeric(names(table(beta.mat.enter)))
  lamlist = c(glmobj$lambda[lambda.ind], 0)
  
  # jlist is the active set of variables
  jlist = vector("list", length = length(lambda.ind))#unlist(fitobj$act) # Predictors added at each step.
  for(i in 1:length(lambda.ind)){
    jlist[[i]] <- which(beta.mat.enter == lambda.ind[i])
  }
  
  # maxp is limiting the number of things to check
  if(is.character(maxp)){
    which.lam.min <- which(glm.mod$lambda == glm.mod$lambda.min)
    num.nonzero <- sum(glmobj$beta[, which.lam.min] != 0)
    maxp = num.nonzero
  } else {
    maxp = min(c(length(jlist), which(lamlist == 0), nrow(x), ncol(x)), maxp)
  }
  
  jlist = jlist[1:maxp]
  cov0 = cov = sig = rep(NA, maxp)
  yy = y - my

  # Iterate through each variable entrance
  for (j in 1:maxp) {
    cat("Calculating covTest ", j, " of ", maxp, "\n")
    lambda = lamlist[j + 1]
    yhat = as.vector(predict(glmobj, x, type = "link", s = lambda/n))
    cov[j] = sum(yy * yhat)
    
    if (j == 1) {
      cov0[j] = 0
    }
    if (j > 1) {
      tt0 = unlist(jlist[1:j])
      #       This is the shortcut that does not work.  It needs to rerun path to 
      #       yhat0 = 1 / (1 + exp(x[, tt0, drop = F] %*%
      #                              glmobj$beta[tt0,lambda.ind[j]] +
      #                              glmobj$a0[lambda.ind[j]]))
      #       yhat0 = as.vector(predict(glmobj, x, type = "link", s = lambda/n))
      glmobj0 = glmnet(x[, tt0, drop = F], y, family = "binomial", 
                       lambda = glmobj$lambda, weights = weights)
      yhat0 = as.vector(predict(glmobj0, x[, tt0, drop = F],
                                type = "link", s = lambda/n))
      cov0[j] = sum(yy * yhat0)
    }
  }
  
  # compile results
  tt = cov - cov0
  
  results = cbind(tt, pexp(tt, 1, lower.tail = FALSE))
  colnames(results) = c("Drop_in_Cov", "P-value")

  which.lam.min <- which(glm.mod$lambda == glm.mod$lambda.min)
  which.lam.1se <- which(glm.mod$lambda == glm.mod$lambda.1se)

  var.pvalue <- data.frame(p.value = rep(1, dim(x)[2]),
                           beta.min = glmobj$beta[, which.lam.min],
                           beta.1se = glmobj$beta[, which.lam.1se])
  for(i in 1:length(jlist)){
    for(j in 1:length(jlist[i])){
      var.pvalue[jlist[i][[1]][j], "p.value"] <- results[i, 2]
    }
  }
  
  list(results = results, var.list = jlist,
       var.pvalue = var.pvalue, glm.mod = glm.mod)
}
