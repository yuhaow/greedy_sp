library(pcalg)
library(graph)

set.seed(1)

#input to Greedy SP: mat -- precision matrix, order -- initial permutation
sp.greedy.alg <- function(mat, lambda){
	#get the initial parameter
	p <- nrow(mat)
	cormat <- cov2cor(solve(mat))
	order <- sample(1:p, p, replace=F)

	#get the initial dag
	init.dag <- function(order){
		return(sapply(1:p, function(j) sapply(1:p, function(i) if(i < j) abs(pcorOrder(i, j, c(1:(j-1))[-i], 
			cormat[order, order])) > lambda else 0)))
	}

	#get new dag based on edge flip
	get.newdag <- function(dag, order, edge){
		#intialize parameters
		neworder <- c(0:(edge[1]-1), edge[2], edge[1]:(edge[2]-1), (edge[2]+1):(p+1))[2:(p+1)]
		newdag <- dag[neworder, neworder]
		newdag[edge[1], edge[1] + 1] <- 1
		newdag[edge[1] + 1, edge[1]] <- 0
		neworder <- order[neworder]
		newcormat <- cormat[neworder, neworder]
		newdag[, c(edge[1], edge[1]+1)] <- sapply(c(edge[1], edge[1]+1), function(x) sapply(1:p, function(i) if(newdag[i,x] != 0) 
			abs(pcorOrder(i, x, c(1:(x-1))[-i], cormat[neworder, neworder])) > lambda else 0))
		return(list(dag=newdag, order=neworder))
	}

	#the stack for visited orders
	vorders <- list()
	vtrace <- list()
	vdags <- list()
	dag <- init.dag(order)
	while(TRUE){
		cov.edge <- which(dag != 0, arr.ind = TRUE)
		cov.edge <- data.frame(subset(cov.edge, apply(cov.edge, 1, function(x) all.equal(c(dag[-x[1], x[1]]), c(dag[-x[1], x[2]])) == TRUE)))
		rdags <- if(nrow(cov.edge) > 0) apply(cov.edge, 1, function(edge) get.newdag(dag, order, edge)) else list()
		if(length(rdags) > 0){
			unvisited <- sapply(rdags, function(rdag) !(list(rdag$order) %in% vorders))
			rdags <- subset(rdags, unvisited)
		}
		if(length(rdags) > 0){
			select <- which(sapply(rdags, function(rdag) sum(rdag$dag != 0) < sum(dag != 0)))
			if(length(select) != 0){
				vorders <- list()
				vtrace <- list()
				vdags <- list()
				order <- rdags[[select[1]]]$order
				dag <- rdags[[select[1]]]$dag
				#print("reduction")
			}else{
				vorders <- append(vorders, list(order))
				vtrace <- append(vtrace, list(order))
				vdags <- append(vdags, list(dag))
				order <- rdags[[1]]$order
				dag <- rdags[[1]]$dag
				#print("forward")
			}
		}else{
			if(length(vtrace) == 0)
				break
			vorders <- append(vorders, list(order))
			order <- tail(vtrace, 1)[[1]]
			vtrace <- head(vtrace, -1)
			dag <- tail(vdags, 1)[[1]]
			vdags <- head(vdags, -1)
			#print("backward")
		}
	}
	dag[order, order] <- dag
	return(dag)
}

pc.alg <- function(mat, lambda){
	partest <- function(x, y, S, suffStat){
		test <- abs(pcorOrder(x, y, S, suffStat$C))
		return(if(test <= lambda) 1 else 0)
	}
	p <- nrow(mat)
	cormat <- cov2cor(solve(mat))
	pc.fit <- pc(suffStat=list(C=cormat, n=10^9), indepTest=partest, alpha=0.01, labels=sapply(1:p, toString), verbose=TRUE)
	return(pc.fit@graph)
}

#set up base parameters
simulation <- function(dag.list, lambda){

	#simulation of single random DAG
	dag.alg <- function(dag){
		#PC algorithm
		pcdag <- pc.alg(dag$dagmat, lambda)

		#greedy SP algorithm
		gspdag <- sp.greedy.alg(dag$dagmat, lambda)
		gspdag <- dag2cpdag(as(gspdag, "graphNEL"))
		#SHD value
		a <- c(shd(dag$cpdag, gspdag), shd(dag$cpdag, pcdag))
		#return the proportion of consistent estimation
		cpdag <- as(dag$cpdag, "matrix")
		b <- sapply(list(gspdag, pcdag), function(dag) all.equal(cpdag, as(dag, "matrix")) == TRUE)
		return(c(b,a))
	}

	#get table of faithful DAGs
	dagnum <- length(dag.list)
	result <- sapply(dag.list, dag.alg)
	faithprop <- rowSums(result[1:2,]) / dagnum
	aveshd <- rowSums(result[3:4,]) / dagnum
	return(list(prop=faithprop, aveshd=aveshd, result=result))
}

#draw pictures
args <- commandArgs()
load(args[4])
lambda <- as.numeric(args[5])

curve <- lapply(dag.list, function(x) simulation(x, lambda))
suffix <- strsplit(basename(args[4]), ".rda")[[1]][1]
save(curve, file=paste("result/lambda_", toString(lambda), "_", suffix, "_methods.rda", sep=""))
