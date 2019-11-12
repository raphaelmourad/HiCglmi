\name{dataExample}
\alias{dataExample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hi-C matrix and ChIP-seq peaks.
}
\description{
This dataExample object comprises four objects: GenomicFeatureList.GR, AnnotNames, SeqInfoChr and HTCL. All data are from Drosophila melanogaster Kc167 cells. Only data for chromosome 2L are provided.

GenomicFeatureList.GR is a list of GRanges objects, one for each genomic feature to analyze. ChIP-seq data (object GenomicFeatureList.GR) were downloaded from Gene Expression Omnibus (GEO) accessions GSE30740, GSE42085 and GSE54529. AnnotNames is a vector of the genomic feature names. SeqInfoChr is a seqinfo object that contains chromosome 2L information. 

HTCL is a HTClist object (HiTC R package) that contains a 20 kb contact matrix Hi-C experiments processed from Gene Expression Omnibus (GEO) accession GSE62904. This is a low resolution contact matrix provided as an example.
}
\usage{
data(dataExample)
}
\references{
Kevin Van Bortle, Michael H. Nichols, Li Li, Chin-Tong Ong, Naomi Takenaka, Zhaohui S. Qin, and Victor G. Corces. Insulator function and topological domain border strength scale with architectural protein occupancy. Genome Biology, 15(5):R82+, June 2014.

Li Li, Xiaowen Lyu, Chunhui Hou, Naomi Takenaka, Huy Q. Nguyen, Chin-Tong Ong, Caelin Cubenas-Potts, Ming Hu, Elissa P. Lei, Giovanni Bosco, Zhaohui S. Qin and Victor G. Corces. Widespread rearrangement of 3D chromatin organization underlies Polycomb-mediated stress-induced silencing. Molecular Cell, 15:S1097-2765, March 2015.
}
\author{
Raphael Mourad
}
