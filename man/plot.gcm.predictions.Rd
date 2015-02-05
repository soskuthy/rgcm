\name{plot.gcm.predictions}
\alias{plot.gcm.predictions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visual summary of GCM predictions 
}
\description{
Creates a density plot of model accuracy / log-likelihood based on the output of \link{predict.gcm}.
}
\usage{
plot.gcm.predictions(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an object of class \code{gcm.predictions}, created by \link{predict.gcm}.
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
# generating posterior samples

model <- gcm(category~size + colour + shape + symmetry, 
             data=rgcm.example, adaptSteps=1000)

# model predictions based on 1000 samples from the posterior 
# distribution

preds <- predict(model)
plot(preds)

# categorisation probabilities & model likelihoods based 
# on 1000 samples from the posterior distribution

preds <- predict(model, type="probabilities")
plot(preds)
}