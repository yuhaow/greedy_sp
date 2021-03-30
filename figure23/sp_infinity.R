library(pcalg)
library(graph)

set.seed(1)

#input to Greedy SP with random restarts: mat -- precision matrix, order -- initial permutation
sp.restart.alg <- function(mat, lambda, ordernum, depth){
	#get the initial parameter
	p <- nrow(mat)
	cormat <- cov2cor(solve(mat))

	#get the initial dag
	init.dag <- function(order){
		cormat <- cormat[order, order]
		return(sapply(1:p, function(j) sapply(1:p, function(i) if(i < j) abs(pcorOrder(i, j, c(1:(j-1))[-i], cormat)) > lambda else 0)))
	}

	#get new dag based on edge flip
	get.newdag <- function(dag, order, edge, vorders){
		#intialize parameters
		neworder <- c(0:(edge[1]-1), edge[2], edge[1]:(edge[2]-1), (edge[2]+1):(p+1))[2:(p+1)]
		final.order <- order[neworder]
		if(list(final.order) %in% vorders) return(NULL)
		newdag <- dag[neworder, neworder]
		newdag[edge[1], edge[1] + 1] <- 1
		newdag[edge[1] + 1, edge[1]] <- 0
		newcormat <- cormat[final.order, final.order]
		newdag[, c(edge[1], edge[1]+1)] <- sapply(c(edge[1], edge[1]+1), function(x) sapply(1:p, function(i) if(newdag[i,x] != 0) 
			abs(pcorOrder(i, x, c(1:(x-1))[-i], newcormat)) > lambda else 0))
		return(list(dag=newdag, order=final.order))
	}

	#the stack for visited orders
	sing.restart <- function(order){
		vorders <- list()
		vtrace <- list()
		vdags <- list()
		dag <- init.dag(order)
		while(TRUE){
			cov.edge <- which(dag != 0, arr.ind = TRUE)
			cov.edge <- data.frame(subset(cov.edge, apply(cov.edge, 1, function(x) sum(dag[-x[1], x[1]] != dag[-x[1], x[2]]) == 0)))
			rdags <- if(nrow(cov.edge) > 0) apply(cov.edge, 1, function(edge) get.newdag(dag, order, edge, vorders)) else list()
			if(length(rdags) > 0) rdags <- subset(rdags, sapply(rdags, function(t) !is.null(t)))
			select <- which(sapply(rdags, function(rdag) sum(rdag$dag != 0) < sum(dag != 0)) == TRUE)
			if(length(rdags) > 0){
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

	#main part of the algorithm
	start.order <- lapply(1:ordernum, function(x) sample(1:p, p, replace=F))
	dag.list <- lapply(start.order, function(order) sing.restart(order))
	edgenum.list <- sapply(dag.list, function(dag) sum(dag != 0))
	return(dag.list[[which.min(edgenum.list)]])
}

#set up base parameters
simulation <- function(dag.list, lambda){

	#simulation of single random DAG
	dag.alg <- function(dag){

		#greedy SP with random restarts
		grspdag.list <- lapply(c(1, 5, 10), function(ordernum) lapply(c(-1), function(depth) sp.restart.alg(dag$dagmat, lambda, ordernum, depth)))
		grspdag.list <- lapply(grspdag.list, function(t) lapply(t, function(grspdag) dag2cpdag(as(grspdag, "graphNEL"))))

		#SHD value
		a <- sapply(grspdag.list, function(t) sapply(t, function(grspdag) shd(dag$cpdag, grspdag)))
		#return the proportion of consistent estimation
		cpdag <- as(dag$cpdag, "matrix")
		b <- sapply(grspdag.list, function(t) sapply(t, function(grspdag) all.equal(cpdag, as(grspdag, "matrix")) == TRUE))
		return(c(b, a))
	}

	#get table of faithful DAGs
	dagnum <- length(dag.list)
	result <- sapply(dag.list, dag.alg)
	faithprop <- rowSums(result[1:3,]) / dagnum
	aveshd <- rowSums(result[4:6,]) / dagnum
	return(list(prop=faithprop, aveshd=aveshd, result=result))
}

#draw pictures
args <- commandArgs()
load(args[4])
lambda <- as.numeric(args[5])
curve <- lapply(dag.list, function(x) simulation(x, lambda))
suffix <- strsplit(basename(args[4]), ".rda")[[1]][1]
save(curve, file=paste("result/lambda_", toString(lambda), "_", suffix, "_sp_infinity.rda", sep=""))
