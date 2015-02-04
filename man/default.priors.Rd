\name{default.priors}
\alias{default.priors}
\docType{data}
\title{
Default priors used by the GCM
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("default.priors")}
\format{
  The format is:
List of 7
 $ c    : language sqrt(dgamma(0.001, 0.001))
 $ w    : language ddirch(1)
 $ b    : language ddirch(1)
 $ alpha: language dgamma(3, 3)
 $ g    : num 1
 $ r    : num 1
 $ m    : language dgamma(0.5, 1)
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(default.priors)
## maybe str(default.priors) ; plot(default.priors) ...
}
\keyword{datasets}
