library(pcalg)
library(graph)
library(parallel)

set.seed(1)

#input to Greedy SP with random restarts: mat -- precision matrix, order -- initial permutation
sp.restart.alg <- function(suffstat){
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
		return(list(dag=dag, order=order))
	}

	#get new dag based on edge flip
	get.newdag <- function(dag, order, edge, vorders, cov.edge){

		#get the new orders
		a <- which(order == edge[1])
		b <- which(order == edge[2])
		order <- order[c(0:(a-1), b, a:(b-1), (b+1):(p+1))[2:(p+1)]]
		#check if the new order has been visited
		if(list(order) %in% vorders) return(NULL)
		#if it has not been visited, check if this edge is an I-covered edge
		par <- subset(1:p, dag[,edge[1]] == 1)
		dag[edge[1], edge[2]] <- 0
		dag[edge[2], edge[1]] <- 1
		#parent set of the flipped components
		if(length(par) != 0){
			dag[par, edge[1]] <- sapply(1:length(par), function(i) gaussCItest(par[i], edge[1], c(par[-i], edge[2]), suffstat) < alpha)
			dag[par, edge[2]] <- sapply(1:length(par), function(i) gaussCItest(par[i], edge[2], par[-i], suffstat) < alpha)
		}
		#get updates of the number of contradicting edges
		edgenum <- sum(dag[par, c(edge[1], edge[2])]) - 2 * length(par)
		#cov.edge <- subset(cov.edge, row != edge[1] & row != edge[2] & col != edge[1] & col != edge[2])
		#cov.edge <- data.frame()
		return(list(dag=dag, order=order, edgenum=edgenum))
	}

	#get the initial dag
	init.dag <- function(order){
		revorder <- sapply(1:p, function(t) which(order==t))
		return(sapply(1:p, function(j) sapply(1:p, function(i) if(revorder[i] < revorder[j]) gaussCItest(i, j, order[c(1:(revorder[j]-1))[-revorder[i]]], suffstat) < alpha else 0)))
	}

	#the stack for visited orders
	sing.restart <- function(dag, order){
		vorders <- list()
		vtrace <- list()
		vdags <- list()
		dag <- init.dag(order)
		while(TRUE){
			#get the list of covered edges
			cov.edge <- which(dag != 0, arr.ind = TRUE)
			cov.edge <- data.frame(subset(cov.edge, apply(cov.edge, 1, function(x) sum(dag[-x[1], x[1]] != dag[-x[1], x[2]]) == 0)))
			#get the list of DAGs after I-covered edge reversals
			rdags <- if(nrow(cov.edge) > 0) apply(cov.edge, 1, function(edge) get.newdag(dag, order, edge, vorders)) else list()
			if(length(rdags) > 0) rdags <- subset(rdags, sapply(rdags, function(t) !is.null(t)))
			select <- which(sapply(rdags, function(rdag) rdag$edgenum < 0) == TRUE)
			#start the searching
			if((length(rdags) > 0 && length(vtrace) != 1) || length(select) != 0){
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
		return(dag)
	}

	#main part of the algorithm
	g.undir <- sapply(1:p, function(j) sapply(1:p, function(i) if(i != j) gaussCItest(i, j, (1:p)[-c(i, j)], suffstat) < alpha else 0))
	start.order <- lapply(1:50, function(x) mindeg(g.undir))
	est.dag.list <- lapply(start.order, function(x) sing.restart(x$dag, x$order))
	edgenum.list <- sapply(est.dag.list, function(dag) sum(dag != 0))
	return(est.dag.list[[which.min(edgenum.list)]])
}

#get data for simulation
args <- commandArgs()
load(args[4])
print(title)
alpha <- as.numeric(args[5])

#set up parallels
cor.num <- min(detectCores(), length(dag.list))
cl <- makeCluster(getOption("cl.cores", cor.num))
clusterExport(cl, c("alpha"))
clusterEvalQ(cl, c(library(pcalg), library(graph)))

#start parallel computation
suffstat.list <- lapply(1:dagnum, function(i) list(C=cor(data.list[[i]]), n=nrow(data.list[[i]])))
grspdag.list <- parLapply(cl, suffstat.list, sp.restart.alg)
grspdag.list <- lapply(grspdag.list, function(dag) dag2cpdag(dag))
stopCluster(cl)
save(grspdag.list, dag.list, file=paste("result/", basename(args[4]), "_alpha_", toString(alpha), "_restart.rda", sep=""))
