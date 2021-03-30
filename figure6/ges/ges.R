library(pcalg)
library(graph)

set.seed(1)

#get data for simulation
args <- commandArgs()
load(args[4])
print(title)
const <- as.numeric(args[5])

ges.alg <- function(data){
	#set up input
	score <- new("GaussL0penObsScore", data=data, intercept = FALSE, lambda = const * log(ncol(data)), use.cpp = TRUE)

	#estimate essential graph
	ges.fit <- ges(score)
	return(as(ges.fit$essgraph, "graphNEL"))
}

gesdag.list <- lapply(1:dagnum, function(i) ges.alg(data.list[[i]]))
save(gesdag.list, dag.list, file=paste("result/", basename(args[4]), "_const_", toString(const), ".rda", sep=""))
