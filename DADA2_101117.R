#!/usr/bin/Rscript

multiplot <- function(..., plotlist=NULL, file, cols=1, rows=1) {
  require(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
	
	i = 1
	while (i < numPlots) {
		numToPlot <- min(numPlots-i+1, cols*rows)
		# Make the panel
		# ncol: Number of columns of plots
		# nrow: Number of rows needed, calculated from # of cols
		layout <- matrix(seq(i, i+cols*rows-1), ncol = cols, nrow = rows, byrow=T)
		if (numToPlot==1) {
		  print(plots[[i]])
		} else {
		  # Set up the page
		  grid.newpage()
		  pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
		  # Make each plot, in the correct location
		  for (j in i:(i+numToPlot-1)) {
		    # Get the i,j matrix positions of the regions that contain this subplot
		    matchidx <- as.data.frame(which(layout == j, arr.ind = TRUE))
		    print(plots[[j]], vp = viewport(layout.pos.row = matchidx$row,
		                                    layout.pos.col = matchidx$col))
		  }
		}
		i <- i+numToPlot
  }
}

library(dada2)
library(ShortRead)
library(ggplot2)

fnFs <- list.files("/Lab_Share/Assorted_101117/fastq/R1/split/fastq/")
fnRs <- list.files("/Lab_Share/Assorted_101117/fastq/R2/split/fastq/")
sample.names <- sapply(strsplit(fnFs, ".fastq"), `[`, 1)
fnFs <- file.path("/Lab_Share/Assorted_101117/fastq/R1/split/fastq/", fnFs)
fnRs <- file.path("/Lab_Share/Assorted_101117/fastq/R2/split/fastq/", fnRs)


# plot quality profiles
pdf("/Lab_Share/Assorted_101117/DADA2/qualityProfiles.pdf")
for (sn in sample.names) {
	plist <- list()
 	fn <- sprintf("/Lab_Share/Assorted_101117/fastq/R1/split/fastq/%s.fastq", sn)
 	plist[[length(plist)+1]] <- plotQualityProfile(fn)
 	fn <- sprintf("/Lab_Share/Assorted_101117/fastq/R2/split/fastq/%s.fastq", sn)
 	plist[[length(plist)+1]] <- plotQualityProfile(fn)
 	print(multiplot(plotlist=plist, cols=1, rows=2))
}
dev.off()

filt_path <- file.path("/Lab_Share/Assorted_101117/fastq/DADA2_filtered")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

for (i in seq_along(fnFs)){
	fastqPairedFilter(c(fnFs[i], fnRs[i]), c(filtFs[i], filtRs[i]),
			trimLeft=c(10, 10), truncLen=c(150, 140),
			maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE,
			compress=TRUE, verbose=TRUE, matchIDs=TRUE)
}

pdf("/Lab_Share/Assorted_101117/DADA2/qualityProfiles.filt.pdf")
for (sn in sample.names) {
	plist <- list()
 	fn <- sprintf("/Lab_Share/Assorted_101117/fastq/DADA2_filtered/%s_F_filt.fastq.gz", sn)
 	plist[[length(plist)+1]] <- plotQualityProfile(fn)
 	fn <- sprintf("/Lab_Share/Assorted_101117/fastq/DADA2_filtered/%s_R_filt.fastq.gz", sn)
 	plist[[length(plist)+1]] <- plotQualityProfile(fn)
 	print(multiplot(plotlist=plist, cols=1, rows=2))
}
dev.off()


derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs.lrn <- dada(derepFs, err=NULL, selfConsist = TRUE, multithread=TRUE)
errF <- dadaFs.lrn[[1]]$err_out
dadaRs.lrn <- dada(derepRs, err=NULL, selfConsist = TRUE, multithread=TRUE)
errR <- dadaRs.lrn[[1]]$err_out

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

seqtab <- makeSequenceTable(mergers)

seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(230, 236)]

seqtab.nochim <- removeBimeraDenovo(seqtab2, verbose=TRUE, multithread=TRUE)

saveRDS(seqtab.nochim, "/Lab_Share/Assorted_101117/DADA2/merged_seqtab.rds")

save.image(file=sprintf("/Lab_Share/Assorted_101117/fastq/DADA2_workspace.RData"))


taxa <- assignTaxonomy(seqtab.nochim, "/Lab_Share/DADA2/rdp_train_set_14.fa.gz", multithread=TRUE)
colnames(taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")


out <- data.frame(id=sprintf(">seq%06d", 1:nrow(taxa)), seq=rownames(taxa))
write.table(out, file="/Lab_Share/Assorted_101117/fastq/DADA2_sequences.fasta", quote=F, sep="\n", row.names=F, col.names=F)
