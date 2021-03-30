library(pcalg)
library(graph)

set.seed(1)

#PC algorithm 
pc.alg <- function(suffstat, alpha){
	p <- nrow(suffstat$C)
	pc.fit <- pc(suffStat=suffstat, indepTest=gaussCItest, alpha=alpha, labels=sapply(1:p, toString), verbose=TRUE)
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
	pcdag <- pc.alg(suffstat, alpha)
})
save(pcdag.list , dag.list, file=paste("result/", strsplit(basename(args[4]), ".", fixed=TRUE)[[1]][1], "_alpha_", toString(alpha), ".rda", sep=""))
