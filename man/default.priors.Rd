\name{default.priors}
\alias{default.priors}
\docType{data}
\title{
Default priors used by the GCM
}
\description{
A list of the \link{rgcm} default priors used during the Bayesian estimation of model parameters for \link{gcm} models. This help page offers some suggestions relating to rgcm priors.
}
\usage{data("default.priors")}
\format{
The GCM models implemented in rgcm use the following priors:
\describe{
    \item{\code{c}}{scale: the precision of the similarity metric}
    \item{\code{w}}{feature weights: the importance of the kth feature in calculating
the distance metric}
    \item{\code{b}}{baseline: the baseline probabilities of the outcome categories}
    \item{\code{alpha}}{abruptness: the extent to which the similarity metric behaves
as a window functiona}
    \item{\code{g}}{determinism: the extent to which the model makes deterministic
predictions}
    \item{\code{r}}{metric type: a parameter that interpolates between different
distance metrics (city-block vs. Euclidean)}
    \item{\code{m}}{token weight: the importance of token frequency in calculating
the similarity metric}
  }
}  
   
\details{
The names of the prior distributions are taken from \link[pkg]{rjags}. Priors specified as constants are not estimated from the data. Priors can be transformed in various ways. Note that both constants and priors have to be quoted using \code{quote()} when they are specified. Here is a list of suggested priors:
\itemize{
  \item \code{c}: \code{quote(sqrt(dgamma(0.001, 0.001)))} (square root of gamma distribution, following Vanpaemel 2009)
  \item \code{w}: \itemize{ \item \code{quote(ddirch(1))} (dirichlet distribution)
                            \item \code{quote(1/npredictors)} (use this exact format, including \code{npredictors}; this produces equal weighting for all features)}
  \item \code{b}: \itemize{ \item \code{quote(ddirch(1))} (dirichlet distribution)
                            \item \code{quote(1/ncat)} (use this exact format, including \code{ncat}; this produces equal baseline probabilities for all outcome categories)}
  \item \code{alpha}: \itemize{ \item \code{quote(1)}
                                \item \code{quote(dgamma(1, 1))} (gamma distribution)}
  \item \code{g}: \itemize{ \item \code{quote(1)} 
                            \item \code{quote(dgamma(1, 1))} (gamma distribution)}
  \item \code{r}: \itemize{ \item \code{quote(1)} (city-block)
                            \item \code{quote(2)} (Euclidean)
                            \item \code{quote(dgamma(1, 1))} (gamma distribution)}
  \item \code{m}: \itemize{ \item \code{quote(0)} (no frequency weighting)
                            \item \code{quote(1)} (linear frequency weighting)
                            \item \code{quote(dgamma(0.5, 1))} (gamma distribution)}
}
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
Vanpaemel, W. (2008) BayesGCM: Software for Bayesian inference with the generalized context model. Behavior Research Methods 41 (4), 1111-1120.
}
\examples{

# estimating frequency weighting from data, fixed baseline probabilities

new.priors <- default.priors
new.priors$m <- quote(dgamma(0.5, 1))
new.priors$b <- quote(1/ncat)

data("rgcm.example")
model <- gcm(category~size + colour + shape + symmetry, 
             data=rgcm.example, weights=rgcm.example$freq, 
             priors=new.priors, adaptSteps=1000)
             
# estimating r from data, c fixed 

new.priors <- default.priors
new.priors$r <- quote(dgamma(1, 1))
new.priors$c <- quote(20)

data("rgcm.example")
model <- gcm(category~size + colour + shape + symmetry, 
             data=rgcm.example, weights=rgcm.example$freq, 
             priors=new.priors, adaptSteps=1000)

}
\keyword{datasets}
