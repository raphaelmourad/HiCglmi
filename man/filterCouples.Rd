\name{filterCouples}
\alias{filterCouples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Filter Hi-C data based on knowledge.
}
\description{
This is an internal function. Function to filter Hi-C data for processing. Hi-C data are filtered according to a set of regions to focus on.
}
\usage{
filterCouples(Left_IC_GR, Right_IC_GR, HiC_data, Filters_GR, FilterMode)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Left_IC_GR}{
GRanges object.
}
  \item{Right_IC_GR}{
GRanges object.
}
  \item{HiC_data}{
Matrix.
}
  \item{Filters_GR}{
GRanges object.
}
  \item{FilterMode}{
Character object.
}
}
\value{
A list of objects.
}
\author{
Raphael Mourad
}
