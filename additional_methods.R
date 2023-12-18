## This Script contains some additional methods that could be alternately, according to the requirements of the data.

## 1.2 Other Filtering Methods ----
# 1.2 Method 2: Filtering low quality cells

set.seed(100)
limit = 100
out.1 <- emptyDrops(raw_mtx, lower = limit, test.ambient = TRUE)
summary(out.1$FDR <= 0.001)    # Tested by FDR
table(Sig = out.1$FDR <= 0.001, Limited = out.1$Limited)    #Concordance by testing with FDR and limited
hist(out$PValue[out.1$Total <= 100 & out.1$Total > 0], xlab = "P-value", main = "", col = "grey80")    #Plot the distribution of significance of non-empty reads.
#filt_mtx1 = filt_mtx[, 1:1000]


# 1.3 Method 3: Filtering low quality cells
counts_per_cell = colSums(filt_mtx)
counts_per_gene = rowSums(filt_mtx)
genes_per_cell = colSums(filt_mtx > 0)    # count gene only if it has non-zero reads mapped.
cells_per_gene = rowSums(filt_mtx > 0)    # only count cells where the gene is expressed.

hist(log10(counts_per_cell + 1), main = "counts per cell", col = "lightblue")    # +1 is used to avoid log10(0) which is not defined
hist(log10(genes_per_cell + 1) ,main = "genes per cell", col = "lightblue")
plot(counts_per_cell, genes_per_cell, log = 'xy',main = " counts per cell v/s genes per cell", col = "lightblue")
plot(sort(genes_per_cell), xlab = "cell", log='y', main = "genes per cell (ordered)")
#here we rank each cells according to the number of genes detected per cell. The lower end outliers correspond to the failed libraries while the higher end outliers could be cell doublets.
