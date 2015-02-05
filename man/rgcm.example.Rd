\name{rgcm.example}
\alias{rgcm.example}
\docType{data}
\title{
Example data set to illustrate the use of the rgcm package
}
\description{
An example data set with two numeric and two categorical predictor variables and a single categorical outcome variable to illustrate the use of the rgcm package. The variables \code{colour} and \code{shape} are strongly associated with the outcome variable. The variables \code{size} and \code{symmetry} were generated randomly.
}
\usage{data("rgcm.example")}
\format{
  A data frame with 40 observations on the following 5 variables.
  \describe{
    \item{\code{colour}}{a numeric predictor variable}
    \item{\code{size}}{a numeric predictor variable}
    \item{\code{shape}}{a categorical predictor variable with levels \code{a} \code{b} \code{c}}
    \item{\code{symmetry}}{a categorical predictor variable with levels \code{no} \code{yes}}
    \item{\code{category}}{the outcome variable, with levels \code{X} \code{Y}}
  }
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
data(rgcm.example)

# plotting the relationship between the predictors and the outcome variable

# colour
hist(rgcm.example[rgcm.example$category=="X",]$colour, 10, col="red")
hist(rgcm.example[rgcm.example$category=="Y",]$colour, 5, col="blue", add=T)

# size
boxplot(size~category, data=rgcm.example)

# shape
mosaicplot(table(rgcm.example$shape,rgcm.example$category))

# symmetry
mosaicplot(table(rgcm.example$symmetry,rgcm.example$category))
}
\keyword{datasets}
