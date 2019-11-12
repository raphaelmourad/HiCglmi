\name{createHiCDataset}
\alias{createHiCDataset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a dataset for further processing with the R package.
}
\description{
This is an internal function. 
}
\usage{
createHiCDataset(HTCL, distInter = NULL, verbose = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{HTCL}{
An HTClist object from HiTC R package.
}
  \item{distInter}{
A numeric vector of two values containing the distance interval used for Hi-C data analysis.
}
  \item{verbose}{
If verbose is True, then the different processing steps will be displayed.
}
}
\value{
Return a list of objects.
}
\author{
Raphael Mourad
}
