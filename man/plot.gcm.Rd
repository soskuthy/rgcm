\name{plot.gcm}
\alias{plot.gcm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Trace and density plots for the GCM posterior samples
}
\description{
This function is a wrapper around the \link[coda]{plot.mcmc} command from the \link{coda} package.
}
\usage{
plot.gcm(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{gcm}.}
  \item{\dots}{further parameters that can be passed to \link[coda]{plot.mcmc}. See the help page for \link[coda]{plot.mcmc}.}
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
model <- gcm(category~size + colour + shape + symmetry, 
             data=rgcm.example, parallel=2, adaptSteps=1000)
summary(model)
plot(model)
}