\name{gcm}
\alias{gcm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bayesian parameter estimation for the Generalised Context Model
}
\description{
Performs Bayesian parameter estimation for the Generalised Context Model.
}
\usage{
gcm(formula, data, iterations = 10000, priors = NULL, 
    weights = NULL, parallel = FALSE, adaptSteps = 500, 
    burnInSteps = 500, nChains = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
an object of class "formula", with the outcome variable on the left hand slide and a list of the predictors on the right hand side, separated by +'s. Interactions are not supported.
}
  \item{data}{
a data frame containing the variables in the model.
}
  \item{iterations}{
an optional integer specifying the overall number of MCMC steps. Each chain runs for \code{iterations / nChains} steps.
}
  \item{priors}{
an optional list specifying the prior distributions for the model parameters. See \link{default.priors} for more detail.
}
  \item{weights}{
an optional vector containing the frequencies of all items in the data set. If it is set to \code{NULL}, a vector of \code{1}'s is used.
}
  \item{parallel}{
an optional logical or integer, which specifies whether the model parameters should be estimated in parallel on several processor cores. If set to \code{FALSE} or \code{0}, a single processor core is used. If set to \code{TRUE} or \code{1}, \code{gcm} attempts to automatically determine the number of available processor cores and runs the estimation in parallel. If set to an integer \code{n > 1}, \code{gcm} runs the estimation on \code{n} cores. The number of cores used is never greater than \code{nChains}.
}
  \item{adaptSteps}{
an optional integer specifying the number of adaptation steps to be run at the beginning of the MCMC chain.
}
  \item{burnInSteps}{
an optional integer specifying the number of burn-in steps to be run before the estimation of the posterior distribution begins.
}
  \item{nChains}{
an optional integer specifying the number of chains to be run in parallel.
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
# load example data set
data("rgcm.example")

# run gcm in parallel mode on 2 cores

model <- gcm(category~size + colour + shape + symmetry, 
             data=rgcm.example, parallel=2, adaptSteps=1000)
summary(model)
plot(model)

# run gcm with frequency data included on single processor core

model.2 <- gcm(category~size + colour + shape + symmetry, 
               data=rgcm.example, weights=rgcm.example$freq, 
               adaptSteps=1000)
summary(model.2)
plot(model.2)

# run gcm with frequency data included and custom priors

new.priors <- default.priors
new.priors$m <- quote(dgamma(0.5, 1))
new.priors$b <- quote(0.5)
model.3 <- gcm(category~size + colour + shape + symmetry, 
               data=rgcm.example, weights=rgcm.example$freq, 
               priors=new.priors, adaptSteps=1000)
summary(model.3)
plot(model.3)

}