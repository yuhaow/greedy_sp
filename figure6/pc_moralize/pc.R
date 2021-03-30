library(pcalg)
library(graph)

set.seed(1)

#PC algorithm 
pc.alg <- function(suffstat, alpha, data.dag){
	p <- nrow(suffstat$C)
	g.undir <- as(data.dag, "matrix")
	g.undir <- ((diag(rep(1, p)) - g.undir) %*% t(diag(rep(1, p)) - g.undir)) == 0
	g.undir[cbind(1:p, 1:p)] <- TRUE
	pc.fit <- pc(suffStat=suffstat, indepTest=gaussCItest, alpha=alpha, labels=sapply(1:p, toString), verbose=TRUE, fixedGaps=g.undir)
	return(pc.fit@graph)
}

#get data for simulation
args <- commandArgs()
load(args[4])
print(title)
alpha <- as.numeric(args[5])

#get result for each dag
pcdag.list <- lapply(1:dagnum, function(i) {
	suffstat <- list(C=cor(data.list[[i]]), n=nrow(data.list[[i]]))
	pcdag <- pc.alg(suffstat, alpha, dag.list[[i]])
})
save(pcdag.list , dag.list, file=paste("result/", basename(args[4]), "_alpha_", toString(alpha), ".rda", sep=""))
