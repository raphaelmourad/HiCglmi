# Raphael MOURAD
# LBME lab
# 07/10/2016




# Function to process Hi-C and genomic feature data (such as ChIP-seq peaks). CHECKED 11/06/2016.
HiCglmiProcData<-function(genomicFeatureList.GR,annotNames,HTCList,
		      distInter,filters.GR=NULL,filterMode=NULL,verbose=F){


# CHECK INPUT DATA ------------------------------------------------------------ 

if(verbose){print("Data checking")}

if(class(genomicFeatureList.GR)!="list"){print("genomicFeatureList.GR is not a list object!"); return(0)}
for(i in 1:length(genomicFeatureList.GR)){
 if(class(genomicFeatureList.GR[[i]])!="GRanges"){print("i-th object of genomicFeatureList.GR is not a GenomicRanges object!"); return(0)}
}
if(class(annotNames)!="character"){print("annotNames is not a character object!"); return(0)}
if(class(HTCList)!="HTClist"){print("HTCList is not a HiTClist object!"); return(0)}
if(class(distInter)!="numeric" & length(distInter)==2){print("distInter is not a vector of 2 numerical values!"); return(0)}
if(!is.null(filters.GR) & class(filters.GR)!="GRanges"){print("filters.GR is not a GRanges object!"); return(0)}



#  DATA TREATMENT---------------------------------------------------------
if(verbose){print("Data parsing")}

# Create HiC_dataset # CHECKED!
Chr=seqlevels(HTCList)
HiC_dataset=createHiCDataset(HTCList,distInter,verbose)
HiC_binList=lapply(HTCList,x_intervals)
names(HiC_binList)=Chr
binSize=width(HiC_binList[[1]][1])
rm(HTCList)

# Filter HiC bin couples based on specific filters.GR (e.g. TAD intervals) # CHECKED!
if(!is.null(filters.GR)){
 HiC_filterRes=filterCouples(HiC_dataset$left.GR,HiC_dataset$right.GR,HiC_dataset$data,filters.GR,filterMode)
 HiC_left.GR=HiC_filterRes$left_filtered
 HiC_right.GR=HiC_filterRes$right_filtered
 HiC_data=HiC_filterRes$data_filtered
 rm(HiC_filterRes)
}else{
 HiC_left.GR=HiC_dataset$left.GR
 HiC_right.GR=HiC_dataset$right.GR
 HiC_data=HiC_dataset$data
}
Idx_IC=1:length(HiC_left.GR)
HiC_left.GR$idxPair=Idx_IC
HiC_right.GR$idxPair=Idx_IC
rm(HiC_dataset,Idx_IC) 

# Biases # CHECKED!
HiC_left.Bias=mcols(HiC_left.GR)[,1:3]
HiC_right.Bias=mcols(HiC_right.GR)[,1:3]

# Count overlapping with annotations # CHECKED!
HiC_left.Fac=NULL
HiC_right.Fac=NULL
for(i in 1:length(genomicFeatureList.GR)){
 annot_righti=NULL
 annot_lefti=NULL
 for(j in 1:length(Chr)){
  annot_binij=annotateHiCBin(HiC_binList[[Chr[j]]],genomicFeatureList.GR[[i]])
  #annot_binij=countOverlaps(HiC_binList[[Chr[j]]],genomicFeatureList.GR[[i]])
  idx_rightij=HiC_right.GR[seqnames(HiC_right.GR)==Chr[j]]$idxBin
  idx_leftij=HiC_left.GR[seqnames(HiC_left.GR)==Chr[j]]$idxBin
  annot_rightij=unlist(annot_binij[idx_rightij])
  annot_leftij=unlist(annot_binij[idx_leftij])
  annot_righti=c(annot_righti,annot_rightij)
  annot_lefti=c(annot_lefti,annot_leftij)
  rm(annot_leftij,annot_rightij)
 }
 if(is.null(HiC_left.Fac)){
  HiC_left.Fac=as(annot_lefti,"Matrix")
  HiC_right.Fac=as(annot_righti,"Matrix")
 }else{
  HiC_left.Fac=cbind(HiC_left.Fac,annot_lefti)
  HiC_right.Fac=cbind(HiC_right.Fac,annot_righti)
 }
 rm(annot_lefti,annot_righti)

 if(verbose){print(paste0(annotNames[i]," annotated"))}
}
rm(idx_leftij,idx_rightij,HiC_binList)
colnames(HiC_left.Fac)=paste0(annotNames,"L")
colnames(HiC_right.Fac)=paste0(annotNames,"R")


ProcData=list(left.Bias=HiC_left.Bias,right.Bias=HiC_right.Bias,left.Fac=HiC_left.Fac,right.Fac=HiC_right.Fac,data=HiC_data,annotNames=annotNames,binSize=binSize)
return(ProcData)

}# End of function HiCglmiProcData


