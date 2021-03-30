load("result/Cebpb__excluded_data.rda.rda")

p <- 24
cutoff <- 0.6
dag <- matrix(as.numeric(as.vector(stab.result$max > cutoff)), nrow=p, ncol=p)
gene_names = read.csv('genes', stringsAsFactors = F, header = F)
gene_names = gene_names$V1

colnames(dag) <- gene_names
rownames(dag) <- gene_names

print(dag)
