\name{predict.gcm}
\alias{predict.gcm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for generating GCM model predictions
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
predict.gcm(object, user.samples = NULL, newdata = NULL, weights = NULL, type = c("response", "probabilities"), sample.type = c("random", "user", "mean"), iterations = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
  \item{user.samples}{
%%     ~~Describe \code{user.samples} here~~
}
  \item{newdata}{
%%     ~~Describe \code{newdata} here~~
}
  \item{weights}{
%%     ~~Describe \code{weights} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{sample.type}{
%%     ~~Describe \code{sample.type} here~~
}
  \item{iterations}{
%%     ~~Describe \code{iterations} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (object, user.samples = NULL, newdata = NULL, weights = NULL, 
    type = c("response", "probabilities"), sample.type = c("random", 
        "user", "mean"), iterations = 1000) 
{
    if (is.null(user.samples)) {
        posterior.samples <- do.call("rbind", object$posterior.samples)
    }
    else {
        if (class(user.samples) \%in\% c("mcmc", "matrix", "data.frame")) {
            if (prod(colnames(object$posterior.samples[[1]]) \%in\% 
                colnames(user.samples))) {
                posterior.samples <- user.samples[, colnames(object$posterior.samples[[1]])]
            }
            else {
                stop("Missing parameter values or parameter names not specified.")
            }
        }
        else if (class(user.samples) == c("numeric")) {
            if (prod(colnames(object$posterior.samples[[1]]) \%in\% 
                names(user.samples))) {
                posterior.samples <- user.samples[colnames(object$posterior.samples[[1]])]
                cnames <- names(posterior.samples)
                dim(posterior.samples) <- c(1, length(posterior.samples))
                colnames(posterior.samples) <- cnames
            }
            else {
                stop("Missing parameter values or parameter names not specified.")
            }
        }
        else {
            stop(paste("'", class(user.samples), "' is not a valid class for user samples", 
                sep = ""))
        }
    }
    if (is.null(newdata)) {
        predictors <- object$data[, -1]
        outcome <- object$data[, 1]
        if (is.null(weights)) {
            freq <- object$weights
        }
        else {
            freq <- weights
        }
    }
    else {
        predictors <- newdata[, all.vars(object$formula)][, -1]
        outcome <- newdata[, all.vars(object$formula)][, 1]
        if (is.null(weights)) {
            freq <- rep(1, nrow(predictors))
        }
        else {
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
    output.list <- vector(mode = "list", length = iterations)
    preds.string <- apply(datalist$predictors, 1, FUN = function(x) {
        return(paste(x, collapse = ""))
    })
    if (sample.type[1] == "user" & !is.null(user.samples)) {
        iterations <- nrow(user.samples)
    }
    pb <- txtProgressBar(min = 0, max = iterations, style = 3)
    for (s in 1:iterations) {
        if (sample.type[1] == "user" & !is.null(user.samples)) {
            ss <- s
        }
        else {
            ss <- sample(1:nrow(posterior.samples), 1)
        }
        step <- posterior.samples[ss, ]
        if (type[1] == "response") {
            predicted.outcomes <- rep("", n)
        }
        else if (type[1] == "probabilities") {
            predicted.outcomes <- matrix(rep(0, n * ncat), nrow = n)
            colnames(predicted.outcomes) <- levels(factor(outcome))
        }
        for (i in 1:n) {
            s_overall <- rep(0, ncat)
            for (j in 1:ncat) {
                observation.weights <- (freq^step["m"]) * (datalist$outcome == 
                  j) * (preds.string != preds.string[i])
                observation.similarities <- exp(-step["c"] * 
                  (((((rep(1, n) \%*\% t(datalist$predictors[i, 
                    ])) != datalist$predictors)^step["r"]) \%*\% 
                    step[paste("w[", 1:ncol(datalist$predictors), 
                      "]", sep = "")])^(1/step["r"]))^step["alpha"])
                s_overall[j] <- sum(observation.weights * observation.similarities)
            }
            probs <- step[paste("b[", 1:ncat, "]", sep = "")] * 
                (s_overall^step["g"])/(sum(step[paste("b[", 1:ncat, 
                "]", sep = "")] * (s_overall^step["g"])))
            if (type[1] == "response") {
                predicted.outcomes[i] <- levels(factor(outcome))[which.max(probs)]
            }
            else if (type[1] == "probabilities") {
                predicted.outcomes[i, ] <- probs
            }
        }
        output.list[[s]] <- list(index = ss, parameters = step, 
            outcome = predicted.outcomes)
        if (s\%\%10 == 0) {
            setTxtProgressBar(pb, s)
        }
    }
    close(pb)
    output <- list(data.type = type[1], predictions = output.list, 
        correct = as.character(object$data[, 1]), formula = object$formula)
    attr(output, "class") <- "gcm.predictions"
    return(output)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
