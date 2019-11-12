# Raphael MOURAD
# LBME lab
# 10/10/2016





# Function to compute generalized linear model with interactions.
HiCglmiModel<-function(hrpd,model,facModel,regressionMode="NB",includeBias=T,orient=F,sampleSize=NULL,distInter=NULL){


# MODEL VARIABLE COMPUTATION ---------------------------------------------------------------
print("Model variable computation")

# Left and right factors and cofactors # CHECKED
annotNames=hrpd$annotNames
facModelS=annotNames[annotNames%in%facModel]
if(!is.null(distInter)){ # filter by distance
testDist=hrpd$data[,1]>=distInter[1] & hrpd$data[,1]<=distInter[2]
}else{
testDist=1:nrow(hrpd$data)
}
HiC_left.Fac=as(hrpd$left.Fac[testDist,annotNames%in%facModelS],"Matrix")
HiC_right.Fac=as(hrpd$right.Fac[testDist,annotNames%in%facModelS],"Matrix")
HiC_left.Bias=as(hrpd$left.Bias[testDist,],"Matrix")
HiC_right.Bias=as(hrpd$right.Bias[testDist,],"Matrix")
HiC_data=hrpd$data[testDist,]

if(!is.null(sampleSize)){
 if(sampleSize<nrow(HiC_data)){ # sampling
  totSize=nrow(HiC_data)
  idxS=sample(1:totSize,sampleSize)
 }else{
  totSize=nrow(HiC_data)
  idxS=1:totSize
 }
}else{
 totSize=nrow(HiC_data)
 idxS=1:totSize
}
if(ncol(HiC_left.Fac)==1){
HiC_left.Fac=as(HiC_left.Fac[idxS],"Matrix")
HiC_right.Fac=as(HiC_right.Fac[idxS],"Matrix")
}else{
HiC_left.Fac=HiC_left.Fac[idxS,]
HiC_right.Fac=HiC_right.Fac[idxS,]
}
HiC_left.Bias=HiC_left.Bias[idxS,]
HiC_right.Bias=HiC_right.Bias[idxS,]
HiC_data=HiC_data[idxS,]
binSize=hrpd$binSize
rm(hrpd)
print(paste0("Data size: ",nrow(HiC_data)," rows"))


# Marginal effects of factors and cofactors # CHECKED 
if(!orient){
 HiC_mat.FacMarg=as((HiC_left.Fac+HiC_right.Fac)/2,"Matrix")
 colnames(HiC_mat.FacMarg)=facModelS
}else if(orient){
 HiC_mat.FacMarg=cbind(HiC_left.Fac,HiC_right.Fac)
}
HiC_mat.FacMarg=as(HiC_mat.FacMarg,"dgCMatrix")


# Marginal couple effects of factors and cofactors 
HiC_mat.FacMarg2=NULL
colnames_HiC_mat.FacMarg2=NULL
if(length(strsplit(as.character(model[3]),"_m_")[[1]])>1){
for(i in 1:length(facModelS)){
 for(j in 1:length(facModelS)){
  HiC_mat.FacMarg2ij=((HiC_left.Fac[,i]*HiC_left.Fac[,j])+(HiC_right.Fac[,i]*HiC_right.Fac[,j]))/2
  colnames_HiC_mat.FacMarg2=c(colnames_HiC_mat.FacMarg2,paste0(facModelS[i],"_m_",facModelS[j]))
  if(is.null(HiC_mat.FacMarg2)){
   HiC_mat.FacMarg2=as(HiC_mat.FacMarg2ij,"Matrix")
  }else{
   HiC_mat.FacMarg2=cbind(HiC_mat.FacMarg2,as(HiC_mat.FacMarg2ij,"Matrix"))
  }
 }
}
colnames(HiC_mat.FacMarg2)=colnames_HiC_mat.FacMarg2
}


# Interaction effects of factors # CHECKED!
HiC_mat.FacInter=NULL
colnames_HiC_mat.FacInter=NULL
for(i in 1:length(facModelS)){
 for(j in 1:length(facModelS)){
  if(!orient){
   HiC_mat.FacInterIJ=as(((HiC_left.Fac[,i]*HiC_right.Fac[,j])+(HiC_left.Fac[,j]*HiC_right.Fac[,i]))/2,"Matrix")
   colnames_HiC_mat.FacInter=c(colnames_HiC_mat.FacInter,paste(facModelS[i],facModelS[j],sep="_"))
  }else{
   HiC_mat.FacInterIJ=as(HiC_left.Fac[,i]*HiC_right.Fac[,j],"Matrix")
   colnames_HiC_mat.FacInter=c(colnames_HiC_mat.FacInter,paste0(facModelS[i],"L_",facModelS[j],"R"))
  }
  HiC_mat.FacInterIJ=as(HiC_mat.FacInterIJ,"dgCMatrix")
  if(is.null(HiC_mat.FacInter)){
   HiC_mat.FacInter=HiC_mat.FacInterIJ
  }else{
   HiC_mat.FacInter=cbind(HiC_mat.FacInter,HiC_mat.FacInterIJ)
  }
  
 }
 #print(i)
}
colnames(HiC_mat.FacInter)=colnames_HiC_mat.FacInter
HiC_mat.FacInter=as(HiC_mat.FacInter,"dgCMatrix")


# Remove Left and Right # CHECKED!
HiC_mat.bias=as(log(HiC_left.Bias*HiC_right.Bias),"dgCMatrix")
rm(HiC_left.Fac,HiC_right.Fac,HiC_left.Bias,HiC_right.Bias)


# Distance # CHECKED!
dist=HiC_data[,1]
dist[HiC_data[,1]==0]=binSize/2
HiC_logDist=log(dist)

# HiC data # CHECKED!
HiC_count=HiC_data[,2]

# All data: Matrix format # CHECKED!
HiC_vecbind=as(cbind(HiC_count,HiC_logDist),"dgCMatrix")
if(includeBias){
 if(is.null(HiC_mat.FacMarg2)){
  HiC_mat.All=cbind(HiC_vecbind,HiC_mat.bias,HiC_mat.FacMarg,HiC_mat.FacInter)
  rm(HiC_vecbind,HiC_mat.bias,HiC_mat.FacMarg,HiC_mat.FacInter)
 }else{
  HiC_mat.All=cbind(HiC_vecbind,HiC_mat.bias,HiC_mat.FacMarg,HiC_mat.FacMarg2,HiC_mat.FacInter)
  rm(HiC_vecbind,HiC_mat.bias,HiC_mat.FacMarg,HiC_mat.FacMarg2,HiC_mat.FacInter)
 }
}else{
 if(is.null(HiC_mat.FacMarg2)){
  HiC_mat.All=cbind(HiC_vecbind,HiC_mat.FacMarg,HiC_mat.FacInter)
  rm(HiC_vecbind,HiC_mat.FacMarg,HiC_mat.FacInter)
 }else{
  HiC_mat.All=cbind(HiC_vecbind,HiC_mat.FacMarg,HiC_mat.FacMarg2,HiC_mat.FacInter)
  rm(HiC_vecbind,HiC_mat.FacMarg,HiC_mat.FacMarg2,HiC_mat.FacInter)
 }
}
colnames(HiC_mat.All)[1:2]=c("Count","logDist")




# REGRESSION ANALYSIS ---------------------------------------------------------------
print("Regression analysis")
print(paste0(regressionMode," regression"))

# GLM # CHECKED!
if(regressionMode!="PoissonLasso" & regressionMode!="RF"){
dataGLM=as.data.frame(as.matrix(HiC_mat.All))
dataGLM=dataGLM[!is.na(dataGLM[,5]),]
dataGLM=dataGLM[dataGLM[,3]!=-Inf,]
if(regressionMode=="Poisson"){
 GLM=suppressWarnings(glm(model,data=dataGLM, family=poisson()))
 #print(dispersiontest(GLM,trafo=1))
}else if(regressionMode=="QP"){
 GLM=suppressWarnings(glm(model,data=dataGLM, family=quasipoisson()))
}else if(regressionMode=="NB"){
 GLM=suppressWarnings(glm.nb(model,data=dataGLM))
}
toReturn=summary(GLM)
}

# Poisson Lasso # CHECKED!
if(regressionMode=="PoissonLasso"){
#varsLasso=c(c("logDist","len","GC","map"),facModelS,paste0(facModelS,"_",facModelS),apply(t(combn(facModelS,2)),1,paste,collapse="_"))
varsLasso=strsplit(gsub("\n    ","",as.character(model[3]))," [+] ")[[1]]
HiC_mat.All=HiC_mat.All[!is.na(HiC_mat.All[,5]),]
HiC_mat.All=HiC_mat.All[HiC_mat.All[,3]!=-Inf,]

CVLasso=cv.glmnet(HiC_mat.All[,varsLasso],HiC_mat.All[,1],family="poisson",parallel=F)
lambda=CVLasso$lambda.min # CVLasso$lambda.min or CVLasso$lambda.1se
CVLassoError=CVLasso$cvm[which(CVLasso$lambda==lambda)]
devLasso=deviance.glmnet(CVLasso$glmnet.fit)[which(CVLasso$lambda==lambda)]
coefLasso=CVLasso$glmnet.fit$beta[,which(CVLasso$lambda==lambda)]
coefLassoMat=data.frame(Variable=names(coefLasso),Coefficient=round(coefLasso,5))
toReturn=coefLassoMat
}


return(toReturn)
}# End of function HiCglmiModel


