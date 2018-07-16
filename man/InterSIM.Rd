\name{InterSIM}
\alias{InterSIM}
\title{
InterSIM
}
\description{
This function simulates three inter-related genomic datasets : DNA methylation, gene expression and protein expression.
}
\usage{
InterSIM(n.sample=500, cluster.sample.prop=c(0.30,0.30,0.40), delta.methyl=2.0,
delta.expr=2.0, delta.protein=2.0, p.DMP=0.2,
p.DEG=NULL, p.DEP=NULL, sigma.methyl=NULL, sigma.expr=NULL, sigma.protein=NULL,
cor.methyl.expr=NULL, cor.expr.protein=NULL, do.plot=FALSE, sample.cluster=TRUE,
feature.cluster=TRUE)
}
\arguments{
  \item{n.sample}{
Number of subjects to simulate
}
  \item{cluster.sample.prop}{
Proportion of samples in the clusters. The number of proportions entered is used to determine the number of clusters in the
simulated data. e.g. if (0.3,0.4,0.3) is entered then the number of clusters will be 3.
}
  \item{delta.methyl}{
Cluster mean shift for methylation data
}
  \item{delta.expr}{
Cluster mean shift for expression data
}
  \item{delta.protein}{
Cluster mean shift for protein data
}
  \item{p.DMP}{
proportion of DE CpGs (DE = Differentially Expressed)
}
  \item{p.DEG}{
proportion of DE mRNA, if NULL (default) mRNAs mapped by DE CpGs will be selected
}
  \item{p.DEP}{
proportion of DE protein, if NULL (default) proteins mapped by DE mRNAs will be selected
}
  \item{sigma.methyl}{
Covariance structure methylation data, if NULL (default) precomputed values will be used.
"indep" gives covariance structure with diagonal elements only (Independent features)
}
  \item{sigma.expr}{
Covariance structure mRNA data, if NULL (default) precomputed values will be used.
"indep" gives covariance structure with diagonal elements only (Independent features)
}
  \item{sigma.protein}{
Covariance structure Protein data, if NULL (default) precomputed values will be used.
"indep" gives covariance structure with diagonal elements only (Independent features)
}
  \item{do.plot}{
TRUE to generate heatmap, default is FALSE
}
  \item{sample.cluster}{
TRUE (default), if clustering should be done on samples for heatmap. This option will be applicable only if do.plot=TRUE.
}
  \item{feature.cluster}{
TRUE (default), if clustering should be done on genomic features for heatmap.This option will be applicable only if do.plot=TRUE.
}
  \item{cor.methyl.expr}{
Correlation between methylation and mRNA, if NULL (default) precomputed values will be used
}
  \item{cor.expr.protein}{
Correlation between mRNA and protein, if NULL (default) precomputed values will be used
}
}
   \alias{CpG.gene.map.for.DEG}
   \alias{cov.M}
   \alias{cov.expr}
   \alias{cov.protein}
   \alias{mean.M}
   \alias{mean.expr}
   \alias{mean.expr.with.mapped.protein}
   \alias{mean.protein}
   \alias{methyl.gene.level.mean}
   \alias{protein.gene.map.for.DEP}
   \alias{rho.expr.protein}
   \alias{rho.methyl.expr}
   \alias{rev.logit}
   \alias{logit}
\value{
This function returns three datasets as matrices - DNA methylation, gene expression and protein expression. It also returns a vector that has true cluster assignment for each subject in the generated data.
}
\author{
Prabhakar Chalise <pchalise@kumc.edu>, Rama Raghavan <rraghavan@kumc.edu>, Brooke Fridley <bfridley@kumc.edu>
}
\examples{
#
prop <- c(0.20,0.30,0.27,0.23)
effect <- 5
sim.data <- InterSIM(n.sample=500, cluster.sample.prop = prop,
delta.methyl=effect, delta.expr=effect, delta.protein=effect,
p.DMP=0.2, p.DEG=NULL, p.DEP=NULL,
sigma.methyl=NULL, sigma.expr=NULL, sigma.protein=NULL,
cor.methyl.expr=NULL, cor.expr.protein=NULL,
do.plot=FALSE, sample.cluster=TRUE, feature.cluster=TRUE)
sim.methyl <- sim.data$dat.methyl
sim.expr <- sim.data$dat.expr
sim.protein <- sim.data$dat.protein

}
