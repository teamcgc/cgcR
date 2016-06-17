\name{isbApp}
\alias{isbApp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
shiny App for oncoPrint from TCGA in ISB bigquery store
}
\description{
shiny App for oncoPrint from TCGA in ISB bigquery store
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
isbApp(bq)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bq}{
a src_bigquery instance
}
}
\details{
user supplies study code and gene set to generate an oncoPrint display
}
\value{
interactive
}
\references{
depends completely on the \code{\link[ComplexHeatmap]{oncoPrint}} function and documentation of the ComplexHeatmap package
}
\author{
VJ Carey <stvjc@channing.harvard.edu>
}
\note{
demonstration prototype, has only a few ad hoc selections of gene sets
}

\examples{
\dontrun{
 # get a src_bigquery in bq
 isbApp(bq)
 }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graphics }% use one of  RShowDoc("KEYWORDS")
