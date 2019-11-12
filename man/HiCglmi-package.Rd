\name{HiCglmi-package}
\alias{HiCglmi-package}
\alias{HiCglmi}
\docType{package}
\title{
\packageTitle{HiCglmi}
}
\description{
\packageDescription{HiCglmi}
}
\author{
Raphael Mourad

Maintainer: Raphael Mourad, raphael.mourad@ibcg.biotoul.fr
}
\references{
Raphael Mourad, Lang Li and Olivier Cuvier. Uncovering direct and indirect molecular determinants of chromatin loops using a computational integrative approach, PLOS Computational Biology, 2017, 13(5):e1005538. 
}

\examples{
# Load data
# The Hi-C matrix is at 20kb resolution (low resolution only for example)
data(dataExample)
genomicFeatureList.GR=dataExample$GenomicFeatureList.GR
annotNames=dataExample$AnnotNames
HTCList=dataExample$HTC
distInter=c(100e3,140e3)

# Annotate Hi-C data with genomic features
HRPD=HiCglmiProcData(genomicFeatureList.GR, annotNames, HTCList, distInter)

# Compute generalized linear model with interactions (GLMI)
model1=as.formula(paste0("Count~logDist+len+GC+map+BEAF32+BEAF32_BEAF32"))
facModel1="BEAF32"
testGLMI1=HiCglmiModel(HRPD, model1, facModel1)
print(testGLMI1$coefficients[7,]) # Homologous interaction variable

if(FALSE){ # Replace FALSE by TRUE to run these examples.
# Annotate Hi-C data with genomic features
distInter=c(40e3,1e6) # more data
HRPD2=HiCglmiProcData(genomicFeatureList.GR, annotNames, HTCList, distInter)

model2=as.formula(paste0("Count~logDist+len+GC+map+BEAF32+dCTCF+BEAF32_dCTCF"))
facModel2=c("BEAF32","dCTCF")
testGLMI2=HiCglmiModel(HRPD2, model2, facModel2)
print(testGLMI2$coefficients[8,]) # Heterologous interaction variable

model3=as.formula(paste0("Count~logDist + len + GC + map + BEAF32 + CP190 + BEAF32_m_CP190 + 
	BEAF32_BEAF32 + BEAF32_CP190 + CP190_CP190 + 
	BEAF32_BEAF32:CP190 + CP190_CP190:BEAF32 + 
	BEAF32_BEAF32:CP190_CP190"))
facModel3=c("BEAF32","CP190")
testGLMI3=HiCglmiModel(HRPD2, model3, facModel3)
print(testGLMI3$coefficients[14,]) # Homologous cofactor interaction variable

model4=as.formula(paste0("Count~logDist + len + GC + map + BEAF32 + GAF + CP190 + BEAF32_m_CP190 + 
    GAF_m_CP190 + BEAF32_GAF + BEAF32_CP190 + GAF_CP190 + CP190_CP190 + 
    BEAF32_CP190 + BEAF32_GAF:CP190 + CP190_CP190:BEAF32 + CP190_CP190:GAF + 
    BEAF32_GAF:CP190_CP190"))
facModel4=c("BEAF32","GAF","CP190")
testGLMI4=HiCglmiModel(HRPD2, model4, facModel4)
print(testGLMI4$coefficients[18,]) # Heterologous cofactor interaction variable

model5=as.formula(paste0("Count~logDist+len+GC+map+",paste(annotNames,collapse="+"),"+",
	paste(paste0(annotNames,"_",annotNames),collapse="+"),"+",
	paste(apply(t(combn(annotNames,2)),1,paste,collapse="_"),collapse='+')))
facModel5=annotNames
testGLMI5=HiCglmiModel(HRPD2, model5, facModel5,regressionMode="PoissonLasso")
print(testGLMI5) # Full model with all interaction variables


}
}
