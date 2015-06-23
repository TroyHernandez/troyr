require("glmnet")

#' @export
covTestBinom <- function (x, y, weights, maxp = "num.nonzero",
                          prop.sample = F, nfolds = 10){
  # First build the model
  
  if(prop.sample == TRUE){
    foldid <- rep(0, length(y))
    foldid[which(y == 0)] <- rep(1:nfolds,
                                 ceiling(length(y) / nfolds))[1:sum(y == 0)]
    foldid[which(y == 1)] <- rep(1:nfolds,
                                 ceiling(length(y) / nfolds))[1:sum(y == 1)]
  }

  set.seed(1)
  glm.mod <- cv.glmnet(x = x, y = y, family = "binomial", weights = weights,
                       foldid = foldid, nfolds = nfolds)
  
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
  
  # var.list is the active set of variables
  var.list = vector("list", length = length(lambda.ind))#unlist(fitobj$act) # Predictors added at each step.
  for(i in 1:length(lambda.ind)){
    var.list[[i]] <- which(beta.mat.enter == lambda.ind[i])
  }
  
  # maxp is limiting the number of things to check
  if(is.character(maxp)){
    which.lam.min <- which(glm.mod$lambda == glm.mod$lambda.min)
    # num.nonzero <- sum(glmobj$beta[, 1:which.lam.min] != 0)
    maxp = which.lam.min #num.nonzero
  } else {
    maxp = min(c(length(var.list), which(lamlist == 0), nrow(x), ncol(x)), maxp)
  }
  
  var.list = var.list[1:maxp]
  cov0 = cov = sig = rep(NA, maxp)
  yy = y - my

  # Iterate through each lambda step with a variable entrance
  for (j in 1:maxp) {
    cat("Calculating covTest ", j, " of ", maxp, "\n")
    lambda = lamlist[j + 1]
    yhat = as.vector(predict(glmobj, x, type = "link", s = lambda/n))
    cov[j] = sum(yy * yhat)
    
    if (j == 1) {
      cov0[j] = 0
    }
    if (j > 1) {
      tt0 = unlist(var.list[1:j])
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
  
  cov.stat = cov - cov0

  which.lam.min <- which(glm.mod$lambda == glm.mod$lambda.min)
  which.lam.1se <- which(glm.mod$lambda == glm.mod$lambda.1se)

  # compile results
  var.list.pvalue = cbind(cov.stat = cov.stat,
                          p.value = pexp(cov.stat, 1, lower.tail = FALSE),
                          beta.min = rep(NA, length(cov.stat)),
                          beta.1se = rep(NA, length(cov.stat)),
                          log.min = rep(NA, length(cov.stat)),
                          log.1se = rep(NA, length(cov.stat)))
  # Set rownames so we can rewrite them later.
  rownames(var.list.pvalue) <- as.character(c(1:length(var.list)))
  
  for(i in 1:length(var.list)){
    ind <- var.list[[i]]
    # if(length(ind) > 1){}
    rownames(var.list.pvalue)[i] <- paste(row.names(glmobj$beta)[ind],
                                          collapse = "|")
    var.list.pvalue[i, "beta.min"] <- mean(glmobj$beta[ind, which.lam.min])
    var.list.pvalue[i, "beta.1se"] <- mean(glmobj$beta[ind, which.lam.1se])
    var.list.pvalue[i, "log.min"] <- mean(exp(glmobj$beta[ind, which.lam.min]))
    var.list.pvalue[i, "log.1se"] <- mean(exp(glmobj$beta[ind, which.lam.1se]))
  }

  var.pvalue <- data.frame(cov.stat = rep(NA, dim(x)[2]),
                           p.value = rep(NA, dim(x)[2]),
                           beta.min = glmobj$beta[, which.lam.min],
                           beta.1se = glmobj$beta[, which.lam.1se],
                           log.min = exp(glmobj$beta[, which.lam.min]),
                           log.1se = exp(glmobj$beta[, which.lam.1se]))

  for(i in 1:length(var.list)){
    indlen <- length(var.list[[i]])
    for(j in 1:length(var.list[[i]])){
      temp.cov.stat <- var.list.pvalue[i, 1] / indlen
      var.pvalue[var.list[[i]][j], "cov.stat"] <- temp.cov.stat
      var.pvalue[var.list[[i]][j], "p.value"] <- pexp(temp.cov.stat, 1,
                                                      lower.tail = FALSE)
    }
  }

  list(var.list.pvalue = var.list.pvalue, var.list = var.list,
       var.pvalue = var.pvalue, glm.mod = glm.mod)
}
