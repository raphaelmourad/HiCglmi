\name{HiCglmiProcData}
\alias{HiCglmiProcData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Process data for further generalized linear regression with interactions (GLMI).
}
\description{
This function is used to annotate Hi-C interaction bin pairs with genomic feature information such as ChIP-seq peak intervals.
}
\usage{
HiCglmiProcData(genomicFeatureList.GR, annotNames, HTCList, distInter, 
		filters.GR = NULL, filterMode = NULL, verbose = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genomicFeatureList.GR}{
A list of GRanges objects. Each GRanges object has been built from either coordinate data using readGFBed function (for instance ChIP-seq peak coordinates or functional element coordinates). 
}
  \item{annotNames}{
A character vector defining the name of each genomic feature. 
Names should not comprise any special character such as ":/+-*^,;!?" because an internal R formula object is created inside the function borderAnalysisFun.

}
  \item{HTCList}{
A HTClist object (HiTC R package). This object contains Hi-C intrachromosomal matrices.
}
  \item{distInter}{
Optional. A vector of two numerical values to set a distance range for Hi-C analysis.
}
  \item{filters.GR}{
Optional. A GRanges object. This is a set of genomic intervals to filter Hi-C data.
}
  \item{filterMode}{
Optional. If filterMode=="in", only Hi-C pairs that are within intervals from filters.GR are kept. If filterMode=="out", only Hi-C pairs that are outside intervals from filters.GR are kept. 
}
  \item{verbose}{
If verbose is True, then the different processing steps will be displayed.
}
}
\value{
An object further used for HiCglmiModel function.
}
\author{
Raphael
}
