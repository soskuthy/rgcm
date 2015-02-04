\name{plot.gcm.predictions}
\alias{plot.gcm.predictions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visual summary of GCM predictions 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.gcm.predictions(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
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
function (object) 
{
    if (object$data.type == "response") {
        accuracies <- summary(object)$accuracies
        par(mfrow = c(1, 1))
        plot(density(accuracies), type = "n", xlab = "Accuracy", 
            ylab = "Density", main = "Distribution of accuracy values")
        lines(density(accuracies), lw = 2)
    }
    else if (object$data.type == "probabilities") {
        likelihoods <- summary(object)$likelihoods
        par(mfrow = c(1, 1))
        plot(density(likelihoods), type = "n", xlab = "Log Likelihood", 
            ylab = "Density", main = "Distribution of Log Likelihoods")
        lines(density(likelihoods), lw = 2)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
