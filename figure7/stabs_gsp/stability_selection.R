library(pcalg)
library(graph)
library(stabs)

#input to Greedy SP with random restarts: mat -- precision matrix, order -- initial permutation
#input to Greedy SP with random restarts: mat -- precision matrix, order -- initial permutation
sp.restart.alg <- function(suffstat, alpha){
	#set up the initial parameters for all functions
	p <- nrow(suffstat$C)

	#minimum degree algorithm to propose initial dag
	mindeg <- function(g.undir){
		#get initial undirected graph
		nodeset <- 1:p
		order <- c()
		dag <- matrix(0, p, p)

		#start the iteration for construction
		while(length(nodeset) > 1){
			#select the node to be removed
			t <- sapply(1:length(nodeset), function(i) sum(g.undir[, i]))
			minnode <- which(t == min(t))
			if(length(minnode) > 1) minnode <- sample(minnode, 1)
			#add to the new permutation
			order <- c(nodeset[minnode], order)
			#update undirected graph
			neighset <- which(g.undir[,minnode] != 0)
			if(length(neighset > 1)){
				g.undir[neighset, neighset] <- sapply(neighset, function(i) sapply(neighset, function(j) if(i != j) 
					{if(g.undir[i, j] == 0) 1 else gaussCItest(i, j, (1:length(nodeset))[-c(i, j, minnode)], suffstat) < alpha} else 0))
			}
			#set up incoming edges for the removed node
			dag[nodeset, nodeset[minnode]] <- g.undir[, minnode]
			g.undir <- g.undir[-minnode, -minnode]
			suffstat$C <- suffstat$C[-minnode, -minnode]
			nodeset <- nodeset[-minnode]
		}
		order <- c(nodeset, order)
		dag <- dag[order, order]
		print(dag)
		return(list(dag=dag, order=order))
	}

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
	sing.restart <- function(dag, order){
		vorders <- list()
		vtrace <- list()
		vdags <- list()
		#dag <- init.dag(order)
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
	g.undir <- sapply(1:p, function(j) sapply(1:p, function(i) if(i != j) gaussCItest(i, j, (1:p)[-c(i, j)], suffstat) < alpha else 0))
	start.order <- lapply(1:20, function(x) mindeg(g.undir))
	dag.list <- lapply(start.order, function(x) sing.restart(x$dag, x$order))
	edgenum.list <- sapply(dag.list, function(dag) sum(dag != 0))
	return(dag.list[[which.min(edgenum.list)]])
}

#get the test significance for estimated edges
conf.mat <- function(dag, order, suffstat){
	p <- ncol(suffstat$C)
	revorder <- sapply(1:p, function(t) which(order==t))
	t <- sapply(1:p, function(j) sapply(1:p, function(i) if(dag[i,j] != 0) 1 - gaussCItest(i, j, order[c(1:(revorder[j]-1))[-revorder[i]]], suffstat) else 0))
	t <- sapply(1:p, function(j) sapply(1:p, function(i) if(dag[i,j] != 0) t[i, j] else t[j, i]))
}

stabs.gsp <- function(x, y, q, ...){
	alphas <- c(0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05)
	dt <- x[, 1:p]
	suffstat <- list(C=cor(dt), n=nrow(dt))
	model.alpha <- function(alpha){
		dag <- sp.restart.alg(suffstat, alpha)
		dag <- as(dag2essgraph(as(dag, "graphNEL")), "matrix")
		as.vector(dag != 0)
	}
	path <- sapply(alphas, model.alpha)
	selected <- rowSums(path) != 0
	return(list(selected=selected, path=path))

}

#get data for simulation
args <- commandArgs()
load(args[4])
print(title)

set.seed(1)

#construct x
x <- data.list[[1]]
p <- ncol(x)
x <- cbind(x, matrix(0, nrow=nrow(x), ncol=p * (p-1)))

y <- rnorm(nrow(x))

stab.result <- stabsel(x = x, y = y, fitfun = stabs.gsp, cutoff = 0.55, PFER = 1)
save(stab.result, file=paste("result/", basename(args[4]), ".rda", sep=""))
