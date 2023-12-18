# Quality Assessment (QA) and Analysis of single cell RNA-Seq data

## Import Data into R and filter out empty droplets.
getwd()
setwd("/home/data/scRNAseq-learning/")

# Install Libraries
library(tidyverse)
library(DropletUtils)
library(Seurat)    #install.packages("Seurat")
library(Matrix)
library(scales)
library(rjson)
library(DT)

raw_mtx = readMM('counts_unfiltered/cellranger/matrix.mtx') #load the generated unfiltered count matrix. Library - Matrix
raw_cells = read.csv('counts_unfiltered/cellranger/barcodes.tsv', header = FALSE, sep = "\t")[,1]
genes = read.csv('counts_unfiltered/cellranger/genes.tsv', sep = '\t', header = F)    #load genes for the same

rownames(raw_mtx) = genes[,1]    #add ensemble gene_ids to the data matrix as rownames
colnames(raw_mtx) = raw_cells[,1]    # add cell barcodes as column names
raw_mtx[1:25,1:3]
dim(raw_mtx)    #how big is the matrix. Report the number of genes (rows) and number of cells (column)

## 1.1 Method 1: Filtering low quality cells ----
# use DropletUtils package to get probability that each barcode is a cell
out = emptyDrops(raw_mtx)  #emptyDrops function will estimate  the probability that a given barcode comes from an empty droplet or the one containing the cell.
table(out$FDR <= 0.001)
table(out$Limited)

keep = out$FDR <= 0.001 #set threshold probability for calling a cell. This can also be set to 0.05
keep[is.na(keep)] = FALSE    #use threshold to remove empty drops
filt_mtx = raw_mtx[,keep]
dim(filt_mtx) 
object.size(filt_mtx) ##300MB. Almost 10 times as much memory.   #size in bytes. How much memory does a sparse matrix take up relative to a dense matrix

write10xCounts('counts_filtered', gene.symbol = genes[,2], filt_mtx, overwrite=T)    # write out filtered results.

## 2. Creating Barcode Rank Plot  ---- 

source('./createBarcodeplot.R')     # Should be present in the same working directory
bcrank = barcodeRanks(raw_mtx)
filt_cells = read.csv('counts_filtered/barcodes.tsv', sep = "\t",header = FALSE)[,1]

# Generate the plot. After sourcing the file, function will be added to the Global Environment. 
createBarcodeRankPlot(bcrank = bcrank, all_cells = raw_cells, trim_cells = filt_cells, save = 'counts_filtered/barcode_rank.png') 

## 3. Knowing some statistics of data ----

# load run info from JSON files produced by Kb
kb_stats = c(fromJSON(file = 'inspect.json'), 
             fromJSON(file = 'run_info.json'))

# determine chemistry version
tech = grep('10X(.*)', strsplit(kb_stats$call, '\\s')[[1]], value=T) 

# make a nice/simple table that summarizes that stats
seq_stats = data.frame(stat = c('Sequencing technology', 'Number of reads processed', '% reads pseudoaligned', # get sequencing/alignment stats 
                                '% reads on whitelist'), 
                       value = prettyNum(c(tech, kb_stats$n_processed, kb_stats$p_pseudoaligned, 
                                           round(kb_stats$percentageReadsOnWhitelist,2)), big.mark = ','))
print(seq_stats)

# calculate cell stats and save to df
p_cnts_in_cells = round((sum(filt_mtx)/sum(raw_mtx))*100, 2) 
med_cnts_cell = median(colSums(filt_mtx))    # Here the function calculates the sum of all values present in the column.
med_genes_cell = median(apply(filt_mtx, 2, function(x) sum(x >= 1)))    ## Here the function counts how many values in each column are greater than or equal to one. i.e. Any zero value will correspond to some gene.
tot_genes_detected = sum(rowSums(filt_mtx)>=1)
cell_stats = data.frame(stat = c('Estimated number of cells', '% counts in cells', 
                                 'Median counts per cell', 'Median genes per cell', 'Total genes detected'), 
                        value = prettyNum(c(ncol(filt_mtx), p_cnts_in_cells, med_cnts_cell,
                                            med_genes_cell, tot_genes_detected), big.mark = ','))
print(cell_stats)

## 4. Create seurat Object ----








## ----



