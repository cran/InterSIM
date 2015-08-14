\name{InterSIM}
\alias{InterSIM}
\title{
InterSIM
}
\description{
This function simulates three inter-related genomic datasets : DNA methylation, gene expression and protein expression.
}
\usage{
InterSIM(n.sample = 500, cluster.sample.prop = c(0.3, 0.3, 0.4), delta.methyl = 2, 
delta.expr = 2, delta.protein = 2, p.DMP = 0.2, p.DEG = NULL, p.DEP = NULL, 
sigma.methyl = NULL, sigma.expr = NULL, sigma.protein = NULL, do.plot = F, 
sample.cluster = F, feature.cluster = F)
}
\arguments{
  \item{n.sample}{
Number of subjects to generate
}
  \item{cluster.sample.prop}{
Proportion of number of samples in each cluster
}
  \item{delta.methyl}{
Effect size for DNA methylation
}
  \item{delta.expr}{
Effect size for gene expression
}
  \item{delta.protein}{
Effect size for protein expression
}
  \item{p.DMP}{
Proportion of differentially methylated CpG probes
}
  \item{p.DEG}{
Proportion of differentially expressed genes
}
  \item{p.DEP}{
Proportion of differentially expressed proteins
}
  \item{sigma.methyl}{
Correlation structure of true methylation data
}
  \item{sigma.expr}{
Correlation structure of true gene expression data
}
  \item{sigma.protein}{
Correlation structure of true protein expression data
}
  \item{do.plot}{
Binary (Y/N) variable for plotting - for data visualization purpose only
}
  \item{sample.cluster}{
Binary (Y/N) variable for clustering samples - for data visualization purpose only
}
  \item{feature.cluster}{
Binary (Y/N) variable for clustering genomic features - for data visualization purpose only
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
Function returns three datasets as matrices - DNA methylation, gene expression and protein expression. It also returns a vector that has true cluster assignment
}
\author{
Prabhakar Chalise <pchalise@kumc.edu>
}
\examples{

InterSIM(n.sample = 500, cluster.sample.prop = c(0.3, 0.3, 0.4), delta.methyl = 2, 
delta.expr = 2, delta.protein = 2, p.DMP = 0.2, p.DEG = NULL, p.DEP = NULL, 
sigma.methyl = NULL, sigma.expr = NULL, sigma.protein = NULL, do.plot = FALSE, 
sample.cluster = FALSE, feature.cluster = FALSE)

## The function is currently defined as
function (n.sample = 500, cluster.sample.prop = c(0.3, 0.3, 0.4), 
    delta.methyl = 5, delta.expr = 5, delta.protein = 5, p.DMP = 0.2, 
    p.DEG = NULL, p.DEP = NULL, sigma.methyl = NULL, sigma.expr = NULL, 
    sigma.protein = NULL, do.plot = F, sample.cluster = F, feature.cluster = F) 
{
    library(MASS)
    library(NMF)
    if (sum(cluster.sample.prop) != 1) 
        stop("The proportions must sum up to 1")
    if (p.DMP < 0 | p.DMP > 1) 
        stop("p.DMP must be between 0 to 1")
    if (!is.null(p.DEG) && (p.DEG < 0 | p.DEG > 1)) 
        stop("p.DEG must be between 0 and 1")
    if (!is.null(p.DEP) && (p.DEP < 0 | p.DEP > 1)) 
        stop("p.DEP must be between 0 and 1")
    n.cluster <- length(cluster.sample.prop)
    cluster.id <- do.call(c, sapply(1:n.cluster, function(x) rep(x, 
        cluster.sample.prop[x] * n.sample)))
    n.CpG <- ncol(cov.M)
    if (!is.null(sigma.methyl)) {
        if (sigma.methyl == "indep") 
            cov.str <- diag(diag(cov.M))
    }
    else cov.str <- cov.M
    DMP <- sapply(1:n.cluster, function(x) rbinom(n.CpG, 1, prob = p.DMP))
    rownames(DMP) <- names(mean.M)
    d <- lapply(1:n.cluster, function(i) {
        effect <- mean.M + DMP[, i] * delta.methyl
        mvrnorm(n = cluster.sample.prop[i] * n.sample, mu = effect, 
            Sigma = cov.str)
    })
    sim.methyl <- do.call(rbind, d)
    sim.methyl <- rev.logit(sim.methyl)
    n.gene <- ncol(cov.expr)
    if (!is.null(sigma.expr)) {
        if (sigma.expr == "indep") 
            cov.str <- diag(diag(cov.expr))
    }
    else cov.str <- cov.expr
    if (!is.null(p.DEG)) {
        DEG <- sapply(1:n.cluster, function(x) rbinom(n.gene, 
            1, prob = p.DEG))
        rownames(DEG) <- names(mean.expr)
    }
    else {
        DEG <- sapply(1:n.cluster, function(x) {
            cg.name <- rownames(subset(DMP, DMP[, x] == 1))
            gene.name <- as.character(CpG.gene.map.for.DEG[cg.name, 
                ]$tmp.gene)
            as.numeric(names(mean.expr) \%in\% gene.name)
        })
        rownames(DEG) <- names(mean.expr)
    }
    d <- lapply(1:n.cluster, function(i) {
        effect <- (rho.methyl.expr * methyl.gene.level.mean + 
            sqrt(1 - rho.methyl.expr^2) * mean.expr) + DEG[, 
            i] * delta.expr
        mvrnorm(n = n.sample * cluster.sample.prop[i], mu = effect, 
            Sigma = cov.str)
    })
    sim.expr <- do.call(rbind, d)
    n.protein <- ncol(cov.protein)
    if (!is.null(sigma.protein)) {
        if (sigma.protein == "indep") 
            cov.str <- diag(diag(cov.protein))
    }
    else cov.str <- cov.protein
    if (!is.null(p.DEP)) {
        DEP <- sapply(1:n.cluster, function(x) rbinom(n.protein, 
            1, prob = p.DEP))
        rownames(DEP) <- names(mean.protein)
    }
    else {
        DEP <- sapply(1:n.cluster, function(x) {
            gene.name <- rownames(subset(DEG, DEG[, x] == 1))
            protein.name <- rownames(protein.gene.map.for.DEP[protein.gene.map.for.DEP$gene \%in\% 
                gene.name, ])
            as.numeric(names(mean.protein) \%in\% protein.name)
        })
        rownames(DEP) <- names(mean.protein)
    }
    d <- lapply(1:n.cluster, function(i) {
        effect <- (rho.expr.protein * mean.expr.with.mapped.protein + 
            sqrt(1 - rho.expr.protein^2) * mean.protein) + DEP[, 
            i] * delta.protein
        mvrnorm(n = n.sample * cluster.sample.prop[i], mu = effect, 
            Sigma = cov.str)
    })
    sim.protein <- do.call(rbind, d)
    indices <- sample(1:n.sample)
    cluster.id <- cluster.id[indices]
    sim.methyl <- sim.methyl[indices, ]
    sim.expr <- sim.expr[indices, ]
    sim.protein <- sim.protein[indices, ]
    rownames(sim.methyl) <- paste("subject", 1:nrow(sim.methyl), 
        sep = "")
    rownames(sim.expr) <- paste("subject", 1:nrow(sim.expr), 
        sep = "")
    rownames(sim.protein) <- paste("subject", 1:nrow(sim.protein), 
        sep = "")
    d.cluster <- data.frame(rownames(sim.methyl), cluster.id)
    colnames(d.cluster)[1] <- "subjects"
    if (do.plot) {
        hmcol <- colorRampPalette(c("blue", "deepskyblue", "white", 
            "orangered", "red3"))(100)
        if (dev.interactive()) 
            dev.off()
        if (sample.cluster && feature.cluster) {
            windows(width = 15, height = 5)
            par(mfrow = c(1, 3))
            aheatmap(t(sim.methyl), col = hmcol, Rowv = F, Colv = F, 
                labRow = NA, labCol = NA, annLegend = T, main = "Methylation", 
                fontsize = 8, breaks = 0.5)
            aheatmap(t(sim.expr), col = hmcol, Rowv = F, Colv = F, 
                labRow = NA, labCol = NA, annLegend = T, main = "Gene expression", 
                fontsize = 8, breaks = 0.5)
            aheatmap(t(sim.protein), col = hmcol, Rowv = F, Colv = F, 
                labRow = NA, labCol = NA, annLegend = T, main = "Protein", 
                fontsize = 8, breaks = 0.5)
        }
        else if (sample.cluster) {
            windows(width = 15, height = 5)
            par(mfrow = c(1, 3))
            aheatmap(t(sim.methyl), col = hmcol, Rowv = NA, Colv = F, 
                labRow = NA, labCol = NA, annLegend = T, main = "Methylation", 
                fontsize = 8, breaks = 0.5)
            aheatmap(t(sim.expr), col = hmcol, Rowv = NA, Colv = F, 
                labRow = NA, labCol = NA, annLegend = T, main = "Gene expression", 
                fontsize = 8, breaks = 0.5)
            aheatmap(t(sim.protein), col = hmcol, Rowv = NA, 
                Colv = F, labRow = NA, labCol = NA, annLegend = T, 
                main = "Protein", fontsize = 8, breaks = 0.5)
        }
        else if (feature.cluster) {
            windows(width = 15, height = 5)
            par(mfrow = c(1, 3))
            aheatmap(t(sim.methyl), col = hmcol, Rowv = F, Colv = NA, 
                labRow = NA, labCol = NA, annLegend = T, main = "Methylation", 
                fontsize = 8, breaks = 0.5)
            aheatmap(t(sim.expr), col = hmcol, Rowv = F, Colv = NA, 
                labRow = NA, labCol = NA, annLegend = T, main = "Gene expression", 
                fontsize = 8, breaks = 0.5)
            aheatmap(t(sim.protein), col = hmcol, Rowv = F, Colv = NA, 
                labRow = NA, labCol = NA, annLegend = T, main = "Protein", 
                fontsize = 8, breaks = 0.5)
        }
        else {
            windows(width = 15, height = 5)
            par(mfrow = c(1, 3))
            aheatmap(t(sim.methyl), col = hmcol, Rowv = NA, Colv = NA, 
                labRow = NA, labCol = NA, annLegend = T, main = "Methylation", 
                fontsize = 8, breaks = 0.5)
            aheatmap(t(sim.expr), col = hmcol, Rowv = NA, Colv = NA, 
                labRow = NA, labCol = NA, annLegend = T, main = "Gene expression", 
                fontsize = 8, breaks = 0.5)
            aheatmap(t(sim.protein), col = hmcol, Rowv = NA, 
                Colv = NA, labRow = NA, labCol = NA, annLegend = T, 
                main = "Protein", fontsize = 8, breaks = 0.5)
        }
    }
    return(list(dat.methyl = sim.methyl, dat.expr = sim.expr, 
        dat.protein = sim.protein, clustering.assignment = d.cluster))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
