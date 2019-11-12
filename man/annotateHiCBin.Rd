\name{annotateHiCBin}
\alias{annotateHiCBin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Annotate Hi-C bins with genomic feature data.
}
\description{
This is an internal function. Function to compute the average value of a genomic feature over a set of non-overlapping bins (Hi-C bins). Genomic features can be ChIP-seq peak intervals.
}
\usage{
annotateHiCBin(HiC_bin.GR, GenomicFeature.GR)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{HiC_bin.GR}{
A GRanges object containing the Hi-C bins.
}
  \item{GenomicFeature.GR}{
A GRanges object containing the genomic feature intervals.
}
}
\value{
A vector object is returned.
}
\author{
Raphael Mourad
}
