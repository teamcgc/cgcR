\name{rainfallBQ}
\alias{rainfallBQ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
kataegis plot for a single TCGA sample 
}
\description{
kataegis plot for a single TCGA sample: mutation types are scatterplotted in genomic coordinates with
altitude of points measuring proximity in genomic coordinates to neighbor
}
\usage{
rainfallBQ(bq, studytag = "LUAD", VariantType = "SNP", id = "TCGA-05-4398", 
   colmap = kataColors(), 
   ymax = 9, inylab = "-log10 dbp", inxlab = "chromosome", 
   intitle, dotrim = TRUE, addbounds = TRUE, cnames = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bq}{
a src\_bigquery instance
}
  \item{studytag}{
alphabetic short code for TCGA tumor type, e.g, LUAD for lung adenocarcinoma, OV for ovarian cystadenocarcinoma, etc.
}
  \item{VariantType}{
code for type of variant, currently SNP is all that is expected
}
  \item{id}{
a ParticipantBarcode value for selected study
}
  \item{colmap}{
vector map from mutation code to colors
}
  \item{ymax}{
maximum altitude of points
}
  \item{inylab}{
label for y axis -- will be -log10dbp to indicate -log10 distance  in base pairs to nearest mutation
}
  \item{inxlab}{
label for x axis
}
  \item{intitle}{
graph title
}
  \item{addbounds}{
logical: if TRUE add some space at ends of x axis
}
  \item{cnames}{
logical: if TRUE add chromosome names
}
}
\details{
uses ggplot2 to plan a visualization of mutation density, genome-wide

note that if addbounds and cnames are set to false, the resulting
gg instance can in principle be made interactive with ggplotly
}
\value{
ggplot instance
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Alexandrov, Ludmil B
Nik-Zainal, Serena
Wedge, David C et al., Signatures of mutational processes in human cancer, Nature 2013.
}
\author{
VJ Carey <stvjc@channing.harvard.edu>
}
\examples{
 \dontrun{
  # get a bigquery instance in bq
  rainfallBQ(bq)
 }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graphics }% use one of  RShowDoc("KEYWORDS")
