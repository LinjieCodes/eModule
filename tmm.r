library(edgeR)

tmm_nor <- function(count, output){
	y <- DGEList(counts=count)
	y <- calcNormFactors(y)
	write.table(y$counts, sep =",", file=output)
}