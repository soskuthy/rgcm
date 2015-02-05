library( rjags )
library( doParallel )
library( foreach )
library( random )


bayesian.model.text <- "
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
  m <- 1
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


rjags.distributions <- c("dbeta", "dchisqr", "ddexp", "dexp", "df", "dgamma", "dgen.gamma", "dlogis", "dlnorm", "dnchisqr", "dnorm", "dpar", "dt", "dunif", "dweib", "dbetabin", "dbern", "dbin", "dcat", "dhyper", "dnegbin", "dpois", "ddirch", "dmnorm", "dwish", "dmt", "dmulti")
one.dim.priors <- c("c", "alpha", "g", "r", "m")
multi.dim.priors <- c("w", "b")

default.priors <- list(c=quote(sqrt(dgamma(0.001, 0.001))),
                       w=quote(ddirch(1)),
                       b=quote(ddirch(1)),
                       alpha=quote(dgamma(3,3)),
                       g=quote(1),
                       r=quote(1),
                       m=quote(dgamma(0.5,1)))

recursive.replace.dist <- function (obj, replacement="") {
  obj.temp <- as.list(obj)
  output <- list()
  output.2 <- NULL
  for (i in 1:length(obj.temp)) {
    if (toString(obj.temp[[i]]) %in% rjags.distributions) {
      output <- replacement
      output.2 <- obj
      return(list(as(output, "name"), output.2))
    } else if (is.call(obj.temp[[i]])) {
      temp.output <- recursive.replace.dist(obj.temp[[i]], replacement=replacement)
      output[[i]] <- temp.output[[1]]
      if (!is.null(temp.output[[2]])) {
        output.2 <- temp.output[[2]]
      }
    } else {
      output[[i]] <- obj.temp[[i]]
    }
  }
  if (length(obj.temp)==1 & class(obj)!="call") {
    return(list(obj, output.2))
  } else {
    return(list(as(output, "call"), output.2))
  }
}

prior.list.to.string <- function (prior.list, single=F) {
  output <- list()
  for (p in 1:length(prior.list)) {
    if (!(names(prior.list)[p] %in% one.dim.priors) & !(names(prior.list)[p] %in% multi.dim.priors)) {
      stop(paste(names(prior.list)[p], "is not a GCM parameter."))
    }
    prior.def <- recursive.replace.dist(prior.list[[p]], replacement=paste(names(prior.list)[p], ".temp", sep=""))
    if ((names(prior.list)[p] %in% one.dim.priors) | (names(prior.list)[p]=="w" & single==T)) {
      if (is.null(prior.def[[2]])) {
        output[[p]] <- paste(names(prior.list)[p], " <- ", deparse(prior.def[[1]]), sep="")
      } else {
        if (deparse(prior.def[[1]])==paste(names(prior.list)[p], ".temp", sep="")) {
          output[[p]] <- paste(names(prior.list)[p], " ~ ", deparse(prior.def[[2]]), sep="")
        } else {
          output[[p]] <- paste(names(prior.list)[p], " <- ", deparse(prior.def[[1]]), "\n", paste(names(prior.list)[p], ".temp", sep=""), " ~ ", deparse(prior.def[[2]]), sep="")
        }
      }
    } else {
      dim.num.name <- ifelse(names(prior.list)[p]=='w', "npredictors", "ncat")
      pars.name <- ifelse(names(prior.list)[p]=='w', "dirichlet_alphas_w", "dirichlet_alphas_b")
      if (is.null(prior.def[[2]])) {
        output[[p]] <- paste("for (i in 1:", dim.num.name, ") {", names(prior.list)[p], "[i] <- ", deparse(prior.def[[1]]), "}", sep="")
      } else {
        output[[p]] <- paste("for (i in 1:", dim.num.name, ") {", pars.name, "[i] <- ", deparse(prior.def[[2]][[2]]), "}\n", 
                             names(prior.list)[p], "[1:", dim.num.name, "] ~ ", deparse(prior.def[[2]][[1]]), "(", 
                             pars.name, "[1:", dim.num.name, "])",sep="")
      }
    }
  }
  return(paste(output, collapse="\n"))
}

jinits <- function( ) {
  return(list(.RNG.name = "lecuyer::RngStream", 
              .RNG.seed = runif(1,1,1e+06) ) )
}

mcmc.combine <- function( ... ) {
  return( as.mcmc.list( sapply( list( ... ), mcmc ) ) )
}


run.model.parallel <- function (model.string, formula, data, weights, adaptSteps=500, burnInSteps=500, nChains=4, iterations=10000, nCores=1) {
  outcome <- data[,1]
  predictors <- data.frame(data[,-1])
  datalist <- create.datalist(predictors, outcome, weights)
  parameters = c( "b" , "w" , "c", "m", "alpha", "g", "r")
  nPerChain = ceiling( iterations / nChains )
  cat("Sorry, you're entering parallel-land, so progress monitoring is not available.\n")
  cat("Generating", as.character(nPerChain), "posterior samples per chain in", as.character(nChains), "chains.\n")
  
  if (nCores==1) {
    coresToUse <- min(nChains, detectCores())
  } else {
    coresToUse <- min(nChains, nCores)
  }
  cat("Using", as.character(coresToUse), "processor cores.\n")
  cat("This might take a while...\n")
  
  cl <- makePSOCKcluster(coresToUse)
  clusterSetRNGStream(cl)
  registerDoParallel(cl)
  
  time.1 <- system.time(jags.parsamples <- foreach(i=1:nChains, .inorder = FALSE, 
                             .packages = c( 'rjags', 'random' ), .combine = "mcmc.combine", 
                             .multicombine = TRUE, .export=c("jinits") ) 
  
  %dopar% {
    
    load.module( "lecuyer" )
    jagsModel <- jags.model( textConnection(model.string) , data=datalist , 
                             inits = jinits,
                             n.adapt=adaptSteps ) # this really ought to work, but it doesn't...
                             #n.adapt=adaptSteps )
    update( jagsModel, n.iter=burnInSteps )
    result <- coda.samples( jagsModel, variable.names = parameters,
                            n.iter = nPerChain)
    
    return( result )
  })
  stopCluster(cl)
  output <- list(posterior.samples=jags.parsamples,
                 model=model.string, 
                 formula=formula,
                 data=data,
                 weights=datalist$freq, 
                 adaptSteps=adaptSteps,
                 burnInSteps=burnInSteps,
                 nChains=nChains,
                 iterations=iterations,
                 time=time.1)
  attr(output, "class") <- "gcm"
  return(output)
}

run.model.parallel.custom <- function (model.string, datalist, parameters, adaptSteps=500, burnInSteps=500, nChains=4, iterations=10000) {
  nPerChain = ceiling( iterations / nChains )
  cat("Sorry, you're entering parallel-land, so progress monitoring is not available.\n")
  cat("Generating", as.character(nPerChain), "posterior samples per chain in", as.character(nChains), "chains.\n")
  cat("This might take a while...\n")
  
  registerDoParallel(min(nChains, detectCores()))
  time.1 <- system.time(jags.parsamples <- foreach(i=1:nChains, .inorder = FALSE, 
                                                   .packages = c( 'rjags', 'random' ), .combine = "mcmc.combine", 
                                                   .multicombine = TRUE, .export=c("jinits") ) 
                        
                        %dopar% {
                          
                          load.module( "lecuyer" )
                          
                          jagsModel <- jags.model( textConnection(model.string) , data=datalist , 
                                                   inits = jinits, n.adapt=adaptSteps )
                          update( jagsModel, n.iter=burnInSteps )
                          result <- coda.samples( jagsModel, variable.names = parameters,
                                                  n.iter = nPerChain)
                          
                          return( result )
                        })
  output <- list(posterior.samples=jags.parsamples,
                 model=model.string,
                 datalist=datalist,
                 adaptSteps=adaptSteps,
                 burnInSteps=burnInSteps,
                 nChains=nChains,
                 iterations=iterations,
                 time=time.1)
  attr(output, "class") <- "gcm.custom"
  return(output)
}





run.model.single <- function (model.string, formula, data, weights, adaptSteps=500, burnInSteps=500, nChains=4, iterations=10000) {
  outcome <- data[,1]
  predictors <- data.frame(data[,-1])
  datalist <- create.datalist(predictors, outcome, weights)
  parameters <- c( "b" , "w" , "c", "m", "alpha", "g", "r")
  nPerChain <- ceiling( iterations / nChains )
  cat("Running rjags on a single processor core.\n")
  cat("Initialising model parameters...\n")
  #print(model.string) # 
  jagsModel <- jags.model(textConnection(model.string), data=datalist , 
                          n.adapt=adaptSteps, n.chains=nChains)
  cat("Burn-in...\n")
  update(jagsModel, n.iter=burnInSteps)
  cat("Generating", as.character(iterations), "posterior samples...\n")
  jags.parsamples <- coda.samples( jagsModel, variable.names = parameters,
                                                       n.iter = nPerChain)
  output <- list(posterior.samples=jags.parsamples, 
                 model=jagsModel, 
                 formula=formula,
                 data=data,
                 weights=datalist$freq, 
                 adaptSteps=adaptSteps,
                 burnInSteps=burnInSteps,
                 nChains=nChains,
                 iterations=iterations)
  attr(output, "class") <- "gcm"
  return(output)
}

calculate.diffs <- function (x,n) {
  if (class(x)=="factor" | class(x)=="character") {
    return(as.numeric(x!=n))
  } 
  if (class(x)=="numeric") {
    return(abs(x-n))
  }
}

create.diff.matrix <- function (preds) {
  matrices <- list()
  for (i in 1:nrow(preds)) {
    row <- preds[i,]
    m <- list()
    for (j in 1:ncol(preds)) {
      m[[j]] <- calculate.diffs(preds[,j], row[[j]])
    }
    matrices[[i]] <- do.call(cbind, m)
  }
  return(aperm(simplify2array(matrices),c(3,1,2)))
}

create.datalist <- function (predictors, outcome, weights) {
  if (is.null(weights)) {
    freq <- rep(1, length(outcome))
  } else {
    freq <- weights
  }
  datalist = list(
    predictors = as.matrix(apply(predictors, 2, FUN=function (x) {as.numeric(factor(x))})),
    diff_matrix = create.diff.matrix(predictors),
    outcome = as.numeric(factor(outcome)),
    outcome_temp = as.numeric(factor(outcome)),
    npredictors = ncol(predictors),
    n = nrow(predictors),
    ncat = length(unique(outcome)),
    indices = as.matrix(1:nrow(predictors)),
    freq = freq
  )
  return(datalist)
}

gcm <- function (formula, data, iterations=10000, priors=NULL, weights=NULL, parallel=FALSE,
                          adaptSteps=500, burnInSteps=500, nChains=4) {
  if (length(all.vars(formula))==2) {
    model.string.lines <- strsplit(bayesian.model.text.single.pred, "\n")[[1]]
    single <- T
  } else if (length(all.vars(formula)) > 2) {
    model.string.lines <- strsplit(bayesian.model.text, "\n")[[1]] #readLines("bayesian_GCM_model.txt")
    single <- F
  } else {
    stop("Formula should be of the form X ~ A + B + ..., where at least X and A have to be specified.")
  }
  if (!is.null(priors)) {
    if (is.character(priors)) {
      model.string <- paste(paste(model.string.lines[1:which(model.string.lines=="# priors")], collapse="\n"), "\n", priors, "\n}", sep="")
    } else if (is.list(priors)) {
      model.string <- paste(paste(model.string.lines[1:which(model.string.lines=="# priors")], collapse="\n"), "\n", paste(prior.list.to.string(priors, single=single), collapse="\n"), "\n}", sep="")
    } else {
      stop(paste("class '", class(priors), "' is not allowed as a prior.", sep=""))
    }
  } else {
    model.string <- paste(model.string.lines, collapse="\n")
  }
  if (parallel) {
    return(run.model.parallel(model.string, formula, data[,all.vars(formula)], weights, iterations=iterations, adaptSteps=adaptSteps, burnInSteps=burnInSteps,
                              nChains=nChains, nCores=parallel)
    )  
  } else {
    return(run.model.single(model.string, formula, data[,all.vars(formula)], weights, iterations=iterations, adaptSteps=adaptSteps, burnInSteps=burnInSteps,
                              nChains=nChains)
    )
  }
}

#gcm.simple <- function (formula, parameters) {
#  return(formula)
#}


summary.gcm <- function (object) {
  cat("GCM formula:\n")
  print(object$formula)
  summary(object$posterior.samples)
}

plot.gcm <- function (object, ...) {
  plot(object$posterior.samples, ...)
}



predict.gcm <- function (object, user.samples=NULL, newdata=NULL, weights=NULL, type=c("response", "probabilities"), sample.type=c("random", "user", "mean"), iterations=1000) {
  if (is.null(user.samples)) {
    posterior.samples <- do.call("rbind", object$posterior.samples)
  } else {
    if (class(user.samples) %in% c("mcmc", "matrix", "data.frame")) {
      if (prod(colnames(object$posterior.samples[[1]]) %in% colnames(user.samples))) {
        posterior.samples <- user.samples[,colnames(object$posterior.samples[[1]])]
      } else {
        stop("Missing parameter values or parameter names not specified.")
      }
    } else if (class(user.samples) == c("numeric")) {
        if (prod(colnames(object$posterior.samples[[1]]) %in% names(user.samples))) {
          posterior.samples <- user.samples[colnames(object$posterior.samples[[1]])]
          cnames <- names(posterior.samples)
          dim(posterior.samples) <- c(1, length(posterior.samples))
          colnames(posterior.samples) <- cnames
        } else {
          stop("Missing parameter values or parameter names not specified.")
        }
    } else {
      stop(paste("'", class(user.samples), "' is not a valid class for user samples", sep=""))
    }
  }
  if (is.null(newdata)) {
    predictors <- object$data[,-1]
    outcome <- object$data[,1]
    if (is.null(weights)) {
      freq <- object$weights
    } else {
      freq <- weights
    }
  } else {
    predictors <- newdata[,all.vars(object$formula)][,-1]
    outcome <- newdata[,all.vars(object$formula)][,1]
    if (is.null(weights)) {
      freq <- rep(1, nrow(predictors))
    } else {
      freq <- weights
    }
  }
  datalist <- create.datalist(predictors, outcome, freq)
  n <- datalist$n
  ncat <- datalist$ncat
  if (sample.type[1] == "mean") {
    posterior.samples <- apply(posterior.samples, 2, mean)
    cnames <- names(posterior.samples)
    dim(posterior.samples) <- c(1, length(posterior.samples))
    colnames(posterior.samples) <- cnames
    iterations <- 1
  }
  output.list <- vector(mode="list", length=iterations)
  preds.string <- apply(datalist$predictors,1, FUN = function(x) {return(paste(x, collapse=''))})
  if (sample.type[1] == "user" & !is.null(user.samples)) {
    iterations <- nrow(user.samples)
  }
  pb <- txtProgressBar(min=0, max=iterations, style=3)
  for (s in 1:iterations) {
    if (sample.type[1] == "user" & !is.null(user.samples)) {
      ss <- s
    } else {
      ss <- sample(1:nrow(posterior.samples), 1)
    }
    step <- posterior.samples[ss,]
    if (type[1]=="response") {
      predicted.outcomes <- rep("", n)
    } else if (type[1]=="probabilities") {
      predicted.outcomes <- matrix(rep(0, n*ncat), nrow=n)
      colnames(predicted.outcomes) <- levels(factor(outcome))
    }
    diff.matrix <- create.diff.matrix(datalist$predictors)
    for (i in 1:n) {
      s_overall <- rep(0, ncat)
      for (j in 1:ncat) {
        observation.weights <-(freq**step["m"])*(datalist$outcome==j)*(preds.string!=preds.string[i])
        observation.similarities <- exp(-step["c"]*(( (diff.matrix[i,,]**step["r"]) %*% step[paste("w[", 1:ncol(datalist$predictors), "]", sep="")])**(1/step["r"]))**step["alpha"])
        s_overall[j] <- sum(observation.weights*observation.similarities)
      }
      probs <- step[paste("b[", 1:ncat, "]", sep="")] * (s_overall**step["g"]) / (sum(step[paste("b[", 1:ncat, "]", sep="")] * (s_overall**step["g"])))
      if (type[1]=="response") {
        predicted.outcomes[i] <- levels(factor(outcome))[which.max(probs)]
      } else if (type[1]=="probabilities") {
        predicted.outcomes[i,] <- probs
      }
    }
    output.list[[s]] <- list(index=ss, parameters=step, outcome=predicted.outcomes)
    if (s %% 10 == 0) {setTxtProgressBar(pb, s)}
  }
  close(pb)
  output <- list(data.type=type[1], predictions=output.list, correct=as.character(object$data[,1]), formula=object$formula)
  attr(output, "class") <- "gcm.predictions"
  return(output)
}

summary.gcm.predictions <- function (object) {
  if (object$data.type=="response") {
    m <- do.call("rbind", lapply(object[[2]], function (x) return(x$outcome)))
    ans <- list()
    ans$formula <- object$formula
    ans$accuracies <- apply(m, 1, function (x) {return(mean(x==object$correct))})
    confusion.data.frame <- data.frame(correct=c(t(m)), predicted=rep(object$correct, nrow(m)))
    ans$confusion.matrix <- xtabs(~predicted+correct, data=confusion.data.frame)
    attr(ans, "class") <- "summary.gcm.predictions"
  } else if (object$data.type=="probabilities") {
    ans <- list()
    ans$formula <- object$formula
    ans$likelihoods <- unlist(lapply(object[[2]], function (x) return(sum(log(t(x$outcome)[seq(0, (nrow(x$outcome)-1)*length(levels(factor(object$correct))), length(levels(factor(object$correct)))) + as.numeric(factor(object$correct))])))))
    attr(ans, "class") <- "summary.gcm.predictions"
  }
  ans
}

print.summary.gcm.predictions <- function (object) {
  cat("GCM formula:\n")
  print(object$formula)
  if ("accuracies" %in% attributes(object)$names) {
    cat("\n", "Prediction accuracy\n", sep="")
    cat("   Min. =", round(min(object$accuracies), 3), "\n")
    cat("   1st Qu. =", round(quantile(object$accuracies, 0.25), 3), "\n")
    cat("   Median =", round(mean(object$accuracies), 3), "\n")
    cat("   3rd Qu. =", round(quantile(object$accuracies, 0.75), 3), "\n")
    cat("   Max. =", round(max(object$accuracies), 3), "\n\n")
    cat("Confusion matrix:\n")
    print(object$confusion.matrix)
  } else {
    cat("\n", "Log Likelihoods\n", sep="")
    cat("   Min. =", round(min(object$likelihoods), 3), "\n")
    cat("   1st Qu. =", round(quantile(object$likelihoods, 0.25), 3), "\n")
    cat("   Median =", round(mean(object$likelihoods), 3), "\n")
    cat("   3rd Qu. =", round(quantile(object$likelihoods, 0.75), 3), "\n")
    cat("   Max. =", round(max(object$likelihoods), 3), "\n\n")
  }
}

plot.gcm.predictions <- function (object) {
  if (object$data.type=="response") {
    accuracies <- summary(object)$accuracies
    par(mfrow=c(1,1))
    plot(density(accuracies), type="n", xlab="Accuracy", ylab="Density", main="Distribution of accuracy values")
    lines(density(accuracies), lw=2)
  } else if (object$data.type=="probabilities"){
    likelihoods <- summary(object)$likelihoods
    par(mfrow=c(1,1))
    plot(density(likelihoods), type="n", xlab="Log Likelihood", ylab="Density", main="Distribution of Log Likelihoods")
    lines(density(likelihoods), lw=2)
  }
}

