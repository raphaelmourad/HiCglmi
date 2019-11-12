\name{testDistancePowerLawDistrib}
\alias{testDistancePowerLawDistrib}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Check log-log linear relation between distance and Hi-C count.
}
\description{
Function to check the log-log linear relation between distance and Hi-C count. At a particular distance range, the log-log function should be linear. The generalized linear regression with interactions (GLMI) can be used only if this log-log linear relation holds. One could consider the log-log linear relation holds if R2>0.95.
}
\usage{
testDistancePowerLawDistrib(HTCL, distInter)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{HTCL}{
A HTClist object (HiTC R package).
}
  \item{distInter}{
Optional. A vector of two numerical values to set a distance range for Hi-C analysis.
}
}
\value{
A list containing the log-log regression information. 
}
\author{
Raphael Mourad
}


