library(pcalg)
library(graph)

#get results for SP algorithm
tprate <- function(prefix){
	load(paste(prefix, sep=""))
	#mean(sapply(1:length(dag.list), function(i) shd(grspdag.list[[i]], dag2essgraph(dag.list[[i]], t.list[[i]]))))
	sum(sapply(1:length(dag.list), function(i) {
		true_A = t(as(dag.list[[i]], "matrix"))
		true_DAG = true_A
		true_DAG[true_DAG != 0] = 1
		true_skel = true_DAG + t(true_DAG)
		est_skel =  as(grspdag.list[[i]], "matrix")
		est_skel = est_skel | t(est_skel)
		sum(est_skel != true_skel) == 0})) / length(dag.list)
	#mean(sapply(1:length(dag.list), function(i) shd(grspdag.list[[i]], dag.list[[i]])))
}

ps <- c(8)
ns <- c(10000, 1000)
neighs <- c(0.2, 1, 2, 3, 4, 5, 6, 7)
dagnum <- 100
alphas <- c(0.001)

for(p in ps) for(n in ns) for(alpha in alphas){
	#get the prefix of data
	prefix <- paste("p_", toString(p), "_neigh_", sep="")
	suffix <- paste("_n_", toString(n), ".rda", "_alpha_", toString(alpha), ".rda", sep="")

	print(paste(prefix, suffix, sep=""))
	a <- sapply(neighs, function(neigh) tprate(paste("sp/result/", prefix, toString(neigh), suffix, sep="")))
	print(a)
}
