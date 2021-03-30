library(pcalg)
library(graph)
library(glasso)

set.seed(1)

#input to Greedy SP with random restarts: mat -- precision matrix, order -- initial permutation
sp.restart.alg <- function(suffstat, alpha){
	#set up the initial parameters for all functions
	p <- nrow(suffstat$C)

	#get new dag based on edge flip
	get.newdag <- function(dag, order, edge){
		#set up new orders as well as dag corresponding to the new order
		neworder <- c(0:(edge[1]-1), edge[2], edge[1]:(edge[2]-1), (edge[2]+1):(p+1))[2:(p+1)]
		newdag <- dag[neworder, neworder]
		newdag[edge[1], edge[1] + 1] <- 1
		newdag[edge[1] + 1, edge[1]] <- 0
		neworder <- order[neworder]
		#parent set of the flipped components
		par <- subset(1:p, newdag[,edge[1]] == 1)
		if(length(par) != 0){
			suffstat$C <- suffstat$C[neworder, neworder]
			newdag[par, edge[1]] <- sapply(1:length(par), function(i) gaussCItest(par[i], edge[1], par[-i], suffstat) < alpha)
			newdag[par, edge[1] + 1] <- sapply(1:length(par), function(i) gaussCItest(par[i], edge[1] + 1, c(par[-i], edge[1]), suffstat) < alpha)
		}
		return(list(dag=newdag, order=neworder))
	}

	#get the initial dag
	init.dag <- function(order){
		suffstat$C <- suffstat$C[order, order]
		return(sapply(1:p, function(j) sapply(1:p, function(i) if(i < j) gaussCItest(i, j, c(1:(j-1))[-i], suffstat) < alpha else 0)))
	}

	#the stack for visited orders
	sing.restart <- function(order){
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
			select <- which(sapply(rdags, function(rdag) sum(rdag$dag != 0) < sum(dag != 0)) == TRUE)
			if((length(rdags) > 0 && length(vtrace) != 4) || length(select) != 0){
				if(length(select) != 0){
					vorders <- list()
					vtrace <- list()
					vdags <- list()
					order <- rdags[[select[1]]]$order
					dag <- rdags[[select[1]]]$dag
				}else{
					vorders <- append(vorders, list(order))
					vtrace <- append(vtrace, list(order))
					vdags <- append(vdags, list(dag))
					order <- rdags[[1]]$order
					dag <- rdags[[1]]$dag
				}
			}else{
				if(length(vtrace) == 0)
					break
				vorders <- append(vorders, list(order))
				order <- tail(vtrace, 1)[[1]]
				vtrace <- head(vtrace, -1)
				dag <- tail(vdags, 1)[[1]]
				vdags <- head(vdags, -1)
			}
		}
		dag[order, order] <- dag
		return(dag)
	}

	#main part of the algorithm
	start.order <- lapply(1:10, function(x) sample(1:p, p, replace=F))
	dag.list <- lapply(start.order, function(order) sing.restart(order))
	edgenum.list <- sapply(dag.list, function(dag) sum(dag != 0))
	return(dag.list[[which.min(edgenum.list)]])
}

#get data for simulation
args <- commandArgs()
load(args[4])
print(title)
alpha <- as.numeric(args[5])

#get result for each dag
grspdag.list <- lapply(1:dagnum, function(i) {
	suffstat <- list(C=cor(data.list[[i]]), n=nrow(data.list[[i]]))
	grspdag <- sp.restart.alg(suffstat, alpha)
	dag2cpdag(as(grspdag, "graphNEL"))
})
save(grspdag.list, dag.list, file=paste("result/", basename(args[4]), "_alpha_", toString(alpha), ".rda", sep=""))
