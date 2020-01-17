
library(dada2)
library(phyloseq)
library(useful)
library(stringi)
library(reshape2)

source("/Lab_Share/fanli/code/Neotrans/utils.R")

ps.project <- readRDS("/Lab_Share/Neotrans/phyloseq/Neotrans_phyloseq_object.rds")
mapping.project <- read.table("/Lab_Share/Neotrans/Neotrans_Mapping.011620.txt", header=T, as.is=T, sep="\t", row.names=1)
project <- "Neotrans"

# contaminant filtering
ps.sel <- subset_samples(ps.project, SampleType %in% c("BufferControl", "PCRWater"))
mapping.sel <- as(sample_data(ps.sel), "data.frame")
otu_table <- as.matrix(as.data.frame(otu_table(ps.sel)))
## distribution of read counts by blank vs. sample for each SV
negids <- rownames(subset(mapping.sel, SampleType %in% c("BufferControl", "PCRWater")))
rs <- rowSums(otu_table[,negids])
rs.true <- rowSums(otu_table(ps.project)[, setdiff(sample_names(ps.project), negids)])
pct_blank <- 100* (rs / (rs + rs.true))
hist(pct_blank, breaks=40)
otus_to_exclude <- names(which(pct_blank > 10))

ps.project <- prune_taxa(setdiff(taxa_names(ps.project), otus_to_exclude), ps.project)

otu.filt <- normalizeByCols(as.data.frame(otu_table(ps.project)))
otu.filt$Family <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(ps.project), level="Family")
agg <- aggregate(. ~ Family, otu.filt, sum)
families <- agg$Family
agg <- agg[,-1, drop=F]
agg <- normalizeByCols(agg)
inds_to_grey <- which(rowMeans(agg)<0.005 & !(apply(agg, 1, function(x) any(x>0.05))))
families[inds_to_grey] <- "Other"
agg$Family <- families
df <- melt(agg, variable.name="SampleID")
df2 <- aggregate(as.formula("value ~ Family + SampleID"), df, sum)
df2$SampleType <- mapping.project[as.character(df2$SampleID), "SampleType"]
df2$Run <- mapping.project[as.character(df2$SampleID), "Run"]
p <- ggplot(df2, aes_string(x="SampleID", y="value", fill="Family", order="Family"))  + geom_bar(stat="identity", position="stack") + ylim(c(-0.1, 1.0))  + theme_classic() + theme(legend.position="right", axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8)) + scale_fill_manual(values=coloring) + ggtitle(sprintf("%s", project)) + guides(col = guide_legend(ncol = 3))
print(p)

for (sampletype in unique(mapping.project$SampleType)){
  ps.sampletype <- subset_samples(ps.project, SampleType == sampletype)
  mapping.sampletype <- as(sample_data(ps.sampletype), "data.frame")
  otu.filt <- normalizeByCols(as.data.frame(otu_table(ps.sampletype)))
  otu.filt$Family <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(ps.sampletype), level="Family")
  agg <- aggregate(. ~ Family, otu.filt, sum)
  families <- agg$Family
  agg <- agg[,-1, drop=F]
  agg <- normalizeByCols(agg)
  inds_to_grey <- which(rowMeans(agg)<0.005 & !(apply(agg, 1, function(x) any(x>0.05))))
  families[inds_to_grey] <- "Other"
  agg$Family <- families
  df <- melt(agg, variable.name="SampleID")
  df2 <- aggregate(as.formula("value ~ Family + SampleID"), df, sum)
  df2$SampleType <- mapping.sampletype[as.character(df2$SampleID), "SampleType"]
  df2$Run <- mapping.sampletype[as.character(df2$SampleID), "Run"]
  df2$SampleID <- as.character(df2$SampleID)
  out <- dcast(df2, Family ~ SampleID, value.var="value")
  write.table(out, file=sprintf("/Lab_Share/Neotrans/phyloseq/counts.Family.%s.txt", sampletype), quote=F, sep="\t", row.names=F, col.names=T)
  p <- ggplot(df2, aes_string(x="SampleID", y="value", fill="Family", order="Family"))  + geom_bar(stat="identity", position="stack") + ylim(c(-0.1, 1.0))  + theme_classic() + theme(legend.position="right", axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8)) + scale_fill_manual(values=coloring) + ggtitle(sprintf("%s", sampletype)) + guides(col = guide_legend(ncol = 3))
  print(p)
}



# overall taxa composition at various levels
out_prefix <- "/Lab_Share/Neotrans/phyloseq"
for (lvl in c("Genus", "Species")) {
	otu.filt <- as.data.frame(otu_table(ps.project))
	otu.filt[["taxa"]] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(ps.project), level=lvl)
	agg <- aggregate(. ~ taxa, otu.filt, sum)
	write.table(agg, file=sprintf("%s/raw_abundance.%s.txt", out_prefix, lvl), quote=F, sep="\t", row.names=F, col.names=T)
}



dev.off()

