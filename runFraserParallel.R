#!/usr/bin/env Rscript

library(tidyverse)
library(dplyr)
library(FRASER)
library(batchtools)
library(BiocParallel)

# Set the working directory to the directory housing the RNASeq bams
setwd("/uufs/chpc.utah.edu/common/home/pezzolesi-group2/RNASeq_data/Pezzolesi_RNASeq_data/fraser/RNASeqBams")

# Number of workers to divide tasks between and the slurm template
tmpl <- system.file(package="batchtools", "templates", "slurm-simple.tmpl")
#param <- BatchtoolsParam(workers=5, cluster="slurm", template=tmpl, resources=list(ncpus=1, walltime=04:00:00, memory=7680))
# walltime is accepted as seconds so 1sec * 60 (for minute) * 60 (for hour) * 24 (for day) * 3 (for num days)
param <- BatchtoolsParam(workers=5, cluster="slurm", template=tmpl, resources=list(ncpus=1, walltime=3 * 24 * 60 * 60, memory=7680, partition="notchpeak", account="pezzolesi"))
register(param)
# Usage of this will look like: BPPARAM = param

# Read in your file with cols of sampleIDs, bamFiles, etc per FRASER doc...
sampleTable <- fread("sampleTable.tsv")

bamFiles <- sampleTable[,bamFile]
sampleTable[,bamFile:=bamFiles]

# Make FRASER settings object
settings <- FraserDataSet(colData=sampleTable,
                          workingDir="/uufs/chpc.utah.edu/common/home/pezzolesi-group2/RNASeq_data/Pezzolesi_RNASeq_data/fraser/working")

# Set whether or not the RNASeq data was sequenced for strand specificity
strandSpecific(settings) <- TRUE

# Do counts on splice sites and splice junctions
fds <- countRNAData(settings, BPPARAM = param)

# Filtering
fds <- calculatePSIValues(fds, BPPARAM = param)

# filtering junction with low expression
fds <- filterExpressionAndVariability(fds, minExpressionInOneSample=20,
                                      minDeltaPsi=0.0, filter=TRUE)

fds <- filterExpressionAndVariability(fds, minDeltaPsi=0.0, filter=FALSE)
plotFilterExpression(fds, bins=100)

fds_filtered <- fds[mcols(fds, type="j")[,"passed"],]

png("my_plot.png", height = 800, width = 600)
plotCountCorHeatmap(fds, type="psi5", logit=TRUE, normalized=FALSE)
dev.off()

png("my_plot2.png", height = 800, width = 600)
plotCountCorHeatmap(fds, type="psi5", logit=TRUE, normalized=FALSE,
                    plotType="junctionSample", topJ=100, minDeltaPsi = 0.01)

dev.off()


# START ANALYSIS ==========


# This is computationally heavy on real size datasets and can take awhile
fds <- FRASER(fds, q=c(psi5=3, psi3=5, theta=2), BPPARAM = param)

png("my_plot3.png", height = 800, width = 600)
plotCountCorHeatmap(fds, type="psi5", normalized=TRUE, logit=TRUE)
dev.off()

# annotate introns with the HGNC symbols of the corresponding gene
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
orgDb <- org.Hs.eg.db
#fds <- annotateRangesWithTxDb(fds, txdb=txdb, orgDb=orgDb)
fds <- annotateRanges(fds) # alternative way using biomaRt
# retrieve results with default and recommended cutoffs (padj <= 0.05 and
# |deltaPsi| >= 0.3)
res <- results(fds)
res
data2 <- as.data.table(res)
write.table(data, "unfilteredResultsFrom_runFraser3")

res <- results(fds, zScoreCutoff=2, padjCutoff=NA, deltaPsiCutoff=0.1)
res
data <- as.data.table(res)
write.table(data, "resultsFrom_runFraser3")

# Edit these plots below for your sample names
png("my_plot_19250X1.png", height = 800, width = 600)
plotVolcano(fds, type="psi5", "19250X1")
dev.off()

png("my_plot_19250X2.png", height = 800, width = 600)
plotVolcano(fds, type="psi5", "19250X2")
dev.off()

png("my_plot_19250X3.png", height = 800, width = 600)
plotVolcano(fds, type="psi5", "19250X3")
dev.off()

png("my_plot_19250X4.png", height = 800, width = 600)
plotVolcano(fds, type="psi5", "19250X4")
dev.off()

png("my_plot_19250X5.png", height = 800, width = 600)
plotVolcano(fds, type="psi5", "19250X5")
dev.off()

png("my_plot_19250X6.png", height = 800, width = 600)
plotVolcano(fds, type="psi5", "19250X6")
dev.off()

png("my_plot_19250X7.png", height = 800, width = 600)
plotVolcano(fds, type="psi5", "19250X7")
dev.off()

png("my_plot_19250X7.png", height = 800, width = 600)
plotVolcano(fds, type="psi5", "19250X7")
dev.off()

png("my_plot_19250X8.png", height = 800, width = 600)
plotVolcano(fds, type="psi5", "19250X8")
dev.off()

png("my_plot_19250X9.png", height = 800, width = 600)
plotVolcano(fds, type="psi5", "19250X9")
dev.off()

png("my_plot_19250X10.png", height = 800, width = 600)
plotVolcano(fds, type="psi5", "19250X10")
dev.off()

# saving a fds
#workingDir(fds) <- "/uufs/chpc.utah.edu/common/home/pezzolesi-group2/RNASeq_data/Pezzolesi_RNASeq_data/fraser/working"
#name(fds) <- "C5orf42_ASAH2_Analysis"
#saveFraserDataSet(fds, dir=workingDir(fds), name=name(fds))
