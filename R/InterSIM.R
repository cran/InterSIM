InterSIM <-
function(n.sample=500,cluster.sample.prop=c(0.30,0.30,0.40),delta.methyl=2.0,delta.expr=2.0,delta.protein=2.0,
						p.DMP=0.2,p.DEG=NULL,p.DEP=NULL,sigma.methyl=NULL,sigma.expr=NULL,sigma.protein=NULL,
						do.plot=F, sample.cluster=F, feature.cluster=F)
	{
	#---------------------------------------------------------------------------------------------------------------
	# n.sample = Number of samples to simulate
	# cluster.sample.prop = Proportion of samples in the clusters. Also, the number of such proportions entered is
	#                       used to determine the number of clusters in the simulated data. e.g. if (0.6,0.4) is
	#						entered then the number of clusters will be 2. 
	# delta.methyl = Cluster mean shift for methylation data
	# delta.expr = Cluster mean shift for expression data
	# delta.protein = Cluster mean shift for protein data
	# p.DMP = proportion of differentially expressed CpGs
	# p.DEG = proportion of differentially expressed expression 
	# p.DEP = proportion of differentially expressed protein 
	# sigma = Covariance structure, covariance structure from the real data is the default
	#         sigma="indep" gives covariance structure with diagonal elements only (Independent features)  
	#---------------------------------------------------------------------------------------------------------------
	if (sum(cluster.sample.prop)!=1) stop("The proportions must sum up to 1")
	if (p.DMP<0 | p.DMP>1) stop("p.DMP must be between 0 to 1")
	if (!is.null(p.DEG) && (p.DEG<0 | p.DEG>1)) stop("p.DEG must be between 0 and 1")
	if (!is.null(p.DEP) && (p.DEP<0 | p.DEP>1)) stop("p.DEP must be between 0 and 1")
	
	n.cluster <- length(cluster.sample.prop) 			# Number of clusters
	cluster.id <- do.call(c,sapply(1:n.cluster, function(x) rep(x,cluster.sample.prop[x]*n.sample))) 
	
	#-----------------
	# Methylation
	#-----------------
	n.CpG <- ncol(cov.M) 								# Number of CpG probes  in the data	
	# Covariance structure
	if (!is.null(sigma.methyl)){                           
		if (sigma.methyl=="indep") cov.str <- diag(diag(cov.M))  # Independent features
		}
	else cov.str <- cov.M   							    # Dependent features	
	# Differenatially methylated CpGs (DMP)
	DMP <- sapply(1:n.cluster,function(x) rbinom(n.CpG, 1, prob = p.DMP))
	rownames(DMP) <- names(mean.M)
	
	d <- lapply(1:n.cluster,function(i) {
			effect <- mean.M + DMP[,i]*delta.methyl
			mvrnorm(n=cluster.sample.prop[i]*n.sample, mu=effect, Sigma=cov.str)})
	sim.methyl <- do.call(rbind,d)		
	sim.methyl <- rev.logit(sim.methyl) 						# Transform back to beta values between (0,1)
	
	#-----------------
	# Gene expression
	#-----------------
	n.gene <- ncol(cov.expr) 										# Number of genes in the data
	# Covariance structure
	if (!is.null(sigma.expr)){                           
		if (sigma.expr=="indep") cov.str <- diag(diag(cov.expr)) 	# Independent features
		}
	else cov.str <- cov.expr   							   	# Dependent features	
	# Differenatially expressed genes (DEG)
	if (!is.null(p.DEG)){
		DEG <- sapply(1:n.cluster,function(x) rbinom(n.gene, 1, prob = p.DEG))
		rownames(DEG) <- names(mean.expr)}
	else {DEG <- sapply(1:n.cluster,function(x){
			cg.name <- rownames(subset(DMP,DMP[,x]==1))
			gene.name <- as.character(CpG.gene.map.for.DEG[cg.name,]$tmp.gene)
			as.numeric(names(mean.expr) %in% gene.name)})
		rownames(DEG) <- names(mean.expr)	
		}
		
	d <- lapply(1:n.cluster,function(i) {
			effect <- (rho.methyl.expr*methyl.gene.level.mean+sqrt(1-rho.methyl.expr^2)*mean.expr) + DEG[,i]*delta.expr
			mvrnorm(n=n.sample*cluster.sample.prop[i], mu=effect, Sigma=cov.str)})
	sim.expr <- do.call(rbind,d)
	
	#-----------------
	# Protein
	#-----------------
	n.protein <- ncol(cov.protein) 								   # Number of genes in the data
	# Covariance structure
	if (!is.null(sigma.protein)){                           
		if (sigma.protein=="indep") cov.str <- diag(diag(cov.protein))  # Independent features
		}
	else cov.str <- cov.protein   							   # Dependent features	
	# Differenatially expressed proteins (DEP)
	if (!is.null(p.DEP)){	
		DEP <- sapply(1:n.cluster,function(x) rbinom(n.protein, 1, prob = p.DEP))
		rownames(DEP) <- names(mean.protein)}
	else {DEP <- sapply(1:n.cluster,function(x){
			gene.name <- rownames(subset(DEG,DEG[,x]==1))
			protein.name <- rownames(protein.gene.map.for.DEP[protein.gene.map.for.DEP$gene %in% gene.name,])
			as.numeric(names(mean.protein) %in% protein.name)})
		rownames(DEP) <- names(mean.protein)
	}
		
	d <- lapply(1:n.cluster,function(i) {
			effect <- (rho.expr.protein*mean.expr.with.mapped.protein+sqrt(1-rho.expr.protein^2)*mean.protein) + DEP[,i]*delta.protein
			mvrnorm(n=n.sample*cluster.sample.prop[i], mu=effect, Sigma=cov.str)})
	sim.protein <- do.call(rbind,d)	

	# Randomly order the samples in the data 	
	indices <- sample(1:n.sample)
	cluster.id <- cluster.id[indices]
	sim.methyl <- sim.methyl[indices,]   	
	sim.expr <- sim.expr[indices,]   	
	sim.protein <- sim.protein[indices,] 
	rownames(sim.methyl) <- paste("subject",1:nrow(sim.methyl),sep="")
	rownames(sim.expr) <- paste("subject",1:nrow(sim.expr),sep="")
	rownames(sim.protein) <- paste("subject",1:nrow(sim.protein),sep="")
	d.cluster <- data.frame(rownames(sim.methyl),cluster.id)
	colnames(d.cluster)[1] <- "subjects"
	
	if(do.plot){
		hmcol <- colorRampPalette(c("blue","deepskyblue","white","orangered","red3"))(100)
		if (dev.interactive()) dev.off() 
		if(sample.cluster && feature.cluster) {
			#windows(width=15, height=5)
			par(mfrow=c(1,3))
			aheatmap(t(sim.methyl),color=hmcol,Rowv=F, Colv=F, labRow=NA, labCol=NA,annLegend=T,main="Methylation",fontsize=8,breaks=0.5)
			aheatmap(t(sim.expr),color=hmcol,Rowv=F, Colv=F, labRow=NA, labCol=NA,annLegend=T,main="Gene expression",fontsize=8,breaks=0.5)
			aheatmap(t(sim.protein),color=hmcol,Rowv=F, Colv=F, labRow=NA, labCol=NA,annLegend=T,main="Protein",fontsize=8,breaks=0.5)}
		else if(sample.cluster) {
			#windows(width=15, height=5)
			par(mfrow=c(1,3))
			aheatmap(t(sim.methyl),color=hmcol,Rowv=NA, Colv=F, labRow=NA, labCol=NA,annLegend=T,main="Methylation",fontsize=8,breaks=0.5)
			aheatmap(t(sim.expr),color=hmcol,Rowv=NA, Colv=F, labRow=NA, labCol=NA,annLegend=T,main="Gene expression",fontsize=8,breaks=0.5)
			aheatmap(t(sim.protein),color=hmcol,Rowv=NA, Colv=F, labRow=NA, labCol=NA,annLegend=T,main="Protein",fontsize=8,breaks=0.5)}
		else if(feature.cluster){
			#windows(width=15, height=5)
			par(mfrow=c(1,3))
			aheatmap(t(sim.methyl),color=hmcol,Rowv=F, Colv=NA, labRow=NA, labCol=NA,annLegend=T,main="Methylation",fontsize=8,breaks=0.5)
			aheatmap(t(sim.expr),color=hmcol,Rowv=F, Colv=NA, labRow=NA, labCol=NA,annLegend=T,main="Gene expression",fontsize=8,breaks=0.5)
			aheatmap(t(sim.protein),color=hmcol,Rowv=F, Colv=NA, labRow=NA, labCol=NA,annLegend=T,main="Protein",fontsize=8,breaks=0.5)}
		else {
			#windows(width=15, height=5)
			par(mfrow=c(1,3))
			aheatmap(t(sim.methyl),color=hmcol,Rowv=NA, Colv=NA, labRow=NA, labCol=NA,annLegend=T,main="Methylation",fontsize=8,breaks=0.5)
			aheatmap(t(sim.expr),color=hmcol,Rowv=NA, Colv=NA, labRow=NA, labCol=NA,annLegend=T,main="Gene expression",fontsize=8,breaks=0.5)
			aheatmap(t(sim.protein),color=hmcol,Rowv=NA, Colv=NA, labRow=NA, labCol=NA,annLegend=T,main="Protein",fontsize=8,breaks=0.5)}
		}
	return(list(dat.methyl=sim.methyl,dat.expr=sim.expr,dat.protein=sim.protein,clustering.assignment=d.cluster))
	}
