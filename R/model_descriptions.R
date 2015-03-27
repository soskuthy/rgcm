bayesian.model.text.training.test <- "
model {
  # outcome probabilities

  for (i in 1:ntest) {
    outcome_test[i] ~ dcat(probcat[i,1:ncat])
  }

  for (i in 1:ntest) {
    for (k in 1:ntrain) {
      target_matrix[i,k,1:npredictors] <- test_predictors[i,]
    }
  }

  for (i in 1:ntrain) {
    for (l in 1:ncat) {
      outcome_matrix[i,l] <- (l==outcome_train[i])
    }
    for (o in 1:npredictors) {
      matrix_mask[i,o] <- vector_mask[o]
    }
  }

  for (i in 1:ntest) {  
    for (j in 1:ncat) {

# defining category probabilites for specific stimuli
# g: determinism

      probcat[i,j] <- b[j] * pow(s_overall[i,j],g) / inprod(b, pow(s_overall[i,], g))
      s_overall[i,j] <- inprod(pow(freq, m)*(outcome_matrix[,j]), pow(2.718282,-c*pow(pow(pow(matrix_mask * (target_matrix[i,1:ntrain,1:npredictors]!=train_predictors) + abs((1-matrix_mask) * (target_matrix[i,1:ntrain,1:npredictors]-train_predictors)),r) %*% w[1:npredictors], 1/r), alpha)))
    }
  }

# priors
  c <- sqrt(csquared)
  csquared ~ dgamma(0.001,0.001)
  for (i in 1:npredictors) {
    dirichlet_alphas_w[i] <- 1
  }
  for (i in 1:ncat) {
    dirichlet_alphas_b[i] <- 1
  }
  w[1:npredictors] ~ ddirch(dirichlet_alphas_w[1:npredictors])
  b[1:ncat] ~ ddirch(dirichlet_alphas_b[1:ncat])
  alpha <- 1
  g <- 1
  r <- 1
  m <- 0
}
"



bayesian.model.text.leave.one.out <- "
model {
  # outcome probabilities

  for (i in 1:n) {
    outcome[i] ~ dcat(probcat[i,1:ncat])
  }

  for (i in 1:n) {
    for (k in 1:n) {
      target_matrix[i,k,1:npredictors] <- predictors[i,]
      identity_matrix[i,k] <- (sum(predictors[i,]==predictors[k,]) != npredictors)
    }
    for (l in 1:ncat) {
      outcome_matrix[i,l] <- (l==outcome_temp[i])
    }
    for (o in 1:npredictors) {
      matrix_mask[i,o] <- vector_mask[o]
    }
  }

  for (i in 1:n) {  
    for (j in 1:ncat) {

# defining category probabilites for specific stimuli
# g: determinism

      probcat[i,j] <- b[j] * pow(s_overall[i,j],g) / inprod(b, pow(s_overall[i,], g))
      s_overall[i,j] <- inprod(pow(freq, m)*outcome_matrix[,j]*identity_matrix[i,], pow(2.718282,-c*pow(pow(pow(matrix_mask * (target_matrix[i,1:n,1:npredictors]!=predictors) + abs((1-matrix_mask) * (target_matrix[i,1:n,1:npredictors]-predictors)),r) %*% w[], 1/r), alpha)))
    }
  }

# priors
  c <- sqrt(csquared)
  csquared ~ dgamma(0.001,0.001)
  for (i in 1:npredictors) {
    dirichlet_alphas_w[i] <- 1
  }
  for (i in 1:ncat) {
    dirichlet_alphas_b[i] <- 1
  }
  w[1:npredictors] ~ ddirch(dirichlet_alphas_w[1:npredictors])
  b[1:ncat] ~ ddirch(dirichlet_alphas_b[1:ncat])

  alpha <- 1
  g <- 1
  r <- 1
  m <- 0
}
"




bayesian.model.text.2 <- "
model {
  # outcome probabilities

  for (i in 1:n) {
    outcome[i] ~ dcat(probcat[i,1:ncat])
  }

  for (i in 1:n) {
    for (k in 1:n) {
      identity_matrix[i,k] <- ((sum(predictors[i,]==predictors[k,]) / npredictors) == 1)*i
    }
    for (l in 1:ncat) {
      outcome_matrix[i,l] <- (l==outcome_temp[i])*outcome_temp[i]
    }
  }

  for (i in 1:n) {  
    for (j in 1:ncat) {

# defining category probabilites for specific stimuli
# g: determinism

      probcat[i,j] <- b[j] * pow(s_overall[i,j],g) / inprod(b, pow(s_overall[i,], g))
      s_overall[i,j] <- inprod(pow(freq, m)*(outcome_temp==outcome_matrix[,j])*(1 - (indices==identity_matrix[i,])),pow(2.718282,-c*pow(pow(pow(diff_matrix[i,1:n,1:npredictors],r) %*% w[], 1/r), alpha)))
    }
  }

# priors
  c <- sqrt(csquared)
  csquared ~ dgamma(0.001,0.001)
  for (i in 1:npredictors) {
    phi[i] ~ dnorm(0,1)
    w[i] <- pow(2.718282, phi[i]) / sum(pow(2.718282, phi[1:npredictors]))
  }
  for (i in 1:ncat) {
    dirichlet_alphas_b[i] <- 1
  }
  #w[1:npredictors] ~ ddirch(dirichlet_alphas_w[1:npredictors])
  b[1:ncat] ~ ddirch(dirichlet_alphas_b[1:ncat])
  alpha <- 1
  g <- 1
  r <- 1
  m <- 1
}
"

bayesian.model.text.training.test.single <- "
model {
  # outcome probabilities

  for (i in 1:ntest) {
    outcome_test[i] ~ dcat(probcat[i,1:ncat])
  }

  for (i in 1:ntest) {
    for (k in 1:ntrain) {
      target_matrix[i,k,1:npredictors] <- test_predictors[i,]
    }
  }

  for (i in 1:ntrain) {
    for (l in 1:ncat) {
      outcome_matrix[i,l] <- (l==outcome_train[i])
    }
    for (o in 1:npredictors) {
      matrix_mask[i,o] <- vector_mask[o]
    }
  }

  for (i in 1:ntest) {  
    for (j in 1:ncat) {

# defining category probabilites for specific stimuli
# g: determinism

      probcat[i,j] <- b[j] * pow(s_overall[i,j],g) / inprod(b, pow(s_overall[i,], g))
      s_overall[i,j] <- inprod(pow(freq, m)*(outcome_matrix[,j]), pow(2.718282,-c*pow(pow(pow(matrix_mask * (target_matrix[i,1:ntrain,1:npredictors]!=train_predictors) + abs((1-matrix_mask) * (target_matrix[i,1:ntrain,1:npredictors]-train_predictors)),r) * w, 1/r), alpha)))
    }
  }

# priors
  c <- sqrt(csquared)
  csquared ~ dgamma(0.001,0.001)
  
  for (i in 1:ncat) {
    dirichlet_alphas_b[i] <- 1
  }
  w <- 1
  b[1:ncat] ~ ddirch(dirichlet_alphas_b[1:ncat])
  alpha <- 1
  g <- 1
  r <- 1
  m <- 0
}
"


bayesian.model.text.leave.one.out.single <- "
model {
  # outcome probabilities

  for (i in 1:n) {
    outcome[i] ~ dcat(probcat[i,1:ncat])
  }

  for (i in 1:n) {
    for (k in 1:n) {
      target_matrix[i,k,1:npredictors] <- predictors[i,]
      identity_matrix[i,k] <- (sum(predictors[i,]==predictors[k,]) != npredictors)
    }
    for (l in 1:ncat) {
      outcome_matrix[i,l] <- (l==outcome_temp[i])
    }
    for (o in 1:npredictors) {
      matrix_mask[i,o] <- vector_mask[o]
    }
  }

  for (i in 1:n) {  
    for (j in 1:ncat) {

# defining category probabilites for specific stimuli
# g: determinism

      probcat[i,j] <- b[j] * pow(s_overall[i,j],g) / inprod(b, pow(s_overall[i,], g))
      s_overall[i,j] <- inprod(pow(freq, m)*outcome_matrix[,j]*identity_matrix[i,], pow(2.718282,-c*pow(pow(pow(matrix_mask * (target_matrix[i,1:n,1:npredictors]!=predictors) + abs((1-matrix_mask) * (target_matrix[i,1:n,1:npredictors]-predictors)),r) * w, 1/r), alpha)))
    }
  }

# priors
  c <- sqrt(csquared)
  csquared ~ dgamma(0.001,0.001)
  for (i in 1:ncat) {
    dirichlet_alphas_b[i] <- 1
  }
  w <- 1
  b[1:ncat] ~ ddirch(dirichlet_alphas_b[1:ncat])

  alpha <- 1
  g <- 1
  r <- 1
  m <- 0
}
"


bayesian.model.text.single.pred <- "
model {
  # outcome probabilities

  for (i in 1:n) {
    outcome[i] ~ dcat(probcat[i,1:ncat])
  }

  for (i in 1:n) {
    for (k in 1:n) {
      target_matrix[i,k,1:npredictors] <- predictors[i,]
      identity_matrix[i,k] <- ((sum(predictors[i,]==predictors[k,]) / npredictors) == 1)*i
    }
    for (l in 1:ncat) {
      outcome_matrix[i,l] <- (l==outcome_temp[i])*outcome_temp[i]
    }
  }

  for (i in 1:n) {  
    for (j in 1:ncat) {

# defining category probabilites for specific stimuli
# g: determinism

      probcat[i,j] <- b[j] * pow(s_overall[i,j],g) / inprod(b, pow(s_overall[i,], g))
      s_overall[i,j] <- inprod(pow(freq, m)*(outcome_temp==outcome_matrix[,j])*(1 - (indices==identity_matrix[i,])),pow(2.718282,-c*pow(pow(pow(diff_matrix[i,1:n,1:npredictors],r) * w, 1/r), alpha)))
    }
  }

# priors
  c <- sqrt(csquared)
  csquared ~ dgamma(0.001,0.001)
  for (i in 1:ncat) {
    dirichlet_alphas_b[i] <- 1
  }
  w <- 1
  b[1:ncat] ~ ddirch(dirichlet_alphas_b[1:ncat])
  alpha <- 1
  g <- 1
  r <- 1
  m <- 1
}

"