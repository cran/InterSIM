\name{InterSIM-package}
\alias{InterSIM-package}
\docType{package}
\title{
Simulation of inter-related genomic datasets 
}
\description{
InterSIM is an R package that generates three inter-related genomic datasets : methylation, gene expression and protein expression.
}
\details{
\tabular{ll}{
Package: \tab InterSIM\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-07-23\cr
License: \tab GPL>=2\cr
}
}
\author{
Prabhakar Chalise, Rama Raghavan, Brooke Fridley
Maintainer: Rama Raghavan <rraghavan@kumc.edu>
}
\examples{
prop <- c(0.30,0.40,0.30)
effect <- 2.0
InterSIM(n.sample=500,cluster.sample.prop=prop,delta.methyl=effect,delta.expr=effect,
delta.protein=effect,p.DMP=0.2,p.DEG=NULL,p.DEP=NULL,do.plot=FALSE, sample.cluster=TRUE, 
feature.cluster=TRUE)
}
