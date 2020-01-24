#!/bin/bash

# "This document is a batch file that contains all Mothur commands (v 1.39.5) that were used for analysis of Canadian cis women cohort for Birse et al. Script is meant to be run using Unix command line. To begin, need .fasta file containing sequences for all samples merged (merged.fasta), and a .file matching sequence IDs to sample IDs (mergegroups)."

#Find total number of seqs in merged.fasta
echo "summary.seqs(fasta=merged.fasta, processors=1)"|srun -p NMLResearch --mem 5G -c 1 mothur

#Keep only reads that are <=427 bp
echo "screen.seqs(fasta=merged.fasta, group=mergegroups, maxambig=0, maxlength=427, processors=8)"|srun -p NMLResearch --mem 20G -c 8 mothur

#Keep unique seqs to reduce comp time
echo "unique.seqs(fasta=merged.good.fasta)"|srun -p NMLResearch --mem 20G -c 8 mothur

#Check number of unique seqs vs total seqs
echo "summary.seqs(fasta=merged.good.unique.fasta, processors=1)"|srun -p NMLResearch --mem 5G -c 1 mothur

#Get count_table of seqs in groups
echo "count.seqs(name=merged.good.names, group=mergegroupsgood)"|srun -p NMLResearch --mem 10G -c 2 mothur

#Get summary of reads (optional)
echo "summary.seqs(fasta=merged.good.unique.fasta, count=merged.good.count_table, processors=1)"|srun -p NMLResearch --mem 5G -c 1 mothur

#Align reads to SILVA reference fasta, which has been premade to be used to align the V3-V4 region
#Keep flip=T, allows for better alignment 
echo "align.seqs(fasta=merged.good.unique.fasta, reference=silva.v3v4.fasta, flip=T, processors=32)"|srun -p NMLResearch --mem 300G -c 32 mothur

#Get summary of aligned reads (optional)
echo "summary.seqs(fasta=merged.good.unique.align, count=merged.good.count_table, processors=2)"|srun -p NMLResearch --mem 5G -c 2 mothur

#Keep only reads that aligned to correct region of SILVA reference (V3-V4)
echo "screen.seqs(fasta=merged.good.unique.align, count=merged.good.count_table, summary=merged.good.unique.summary, start=6428, end=23440, maxhomop=8, processors=32)"|srun -p NMLResearch --mem 50G -c 32 mothur

#Get summary of aligned reads that were retained
echo "summary.seqs(fasta=merged.good.unique.good.align, count=merged.good.good.count_table, processors=2)"|srun -p NMLResearch --mem 5G -c 2 mothur

#Filter out gaps, etc. to get a smaller .align file
echo "filter.seqs(fasta=merged.good.unique.good.align, trump=., vertical=T, processors=8)"|srun -p NMLResearch --mem 40G -c 8 mothur

#Keep unique seqs to reduce comp time
echo "unique.seqs(fasta=merged.good.unique.good.filter.fasta, count=merged.good.good.count_table)"|srun -p NMLResearch --mem 10G -c 2 mothur

#Precluster: sequences within a sample (aka group) that differ by 2 nt or less are clustered together for later steps
echo "pre.cluster(fasta=merged.good.unique.good.filter.unique.fasta, count=merged.good.unique.good.filter.count_table, diffs=2, processors=32)"|srun -p NMLResearch --mem 100G -c 32 mothur

#Identify chimeras. MAKE SURE that dereplicate=T !!! or else 1 of the output files will not be made  
echo "chimera.uchime(fasta=merged.good.unique.good.filter.unique.precluster.fasta, count=merged.good.unique.good.filter.unique.precluster.count_table, dereplicate=T, processors=32)"| srun -p NMLResearch --mem 300G -c 32 mothur

#Remove chimeras
echo "remove.seqs(fasta=merged.good.unique.good.filter.unique.precluster.fasta, accnos=merged.good.unique.good.filter.unique.precluster.denovo.uchime.accnos)"|srun -p NMLResearch --mem 5G -c 2 mothur

#Get summary of reads (optional)
echo "summary.seqs(fasta=merged.good.unique.good.filter.unique.precluster.pick.fasta, count=merged.good.unique.good.filter.unique.precluster.denovo.uchime.pick.count_table, processors=8)"| srun -p NMLResearch --mem 10G -c 8 mothur

#Classify seqs by comparing to RDP taxa database for 16S
echo "classify.seqs(fasta=merged.good.unique.good.filter.unique.precluster.pick.fasta, count=merged.good.unique.good.filter.unique.precluster.denovo.uchime.pick.count_table, reference=trainset16_022016.rdp.fasta, taxonomy=trainset16_022016.rdp.tax, cutoff=60, processors=32)"|srun -p NMLResearch --mem 100G -c 32 mothur

#Get summary file of taxonomy (.tax.summary file) #output from this step is the final mothur output you need for R
echo "summary.tax(taxonomy=merged.good.unique.good.filter.unique.precluster.pick.rdp.wang.taxonomy,count=merged.good.unique.good.filter.unique.precluster.denovo.uchime.pick.count_table)"|srun -p NMLResearch --mem 10G -c 2 mothur

#Get summary of 'phylotypes' (sim to OTUs - see MiSeq SOP wiki) (extra info about OTUs in each sample - not needed for R)
echo "phylotype(taxonomy=merged.good.unique.good.filter.unique.precluster.pick.rdp.wang.taxonomy)"|srun -p NMLResearch --mem 5G -c 2 mothur
echo "make.shared(list=merged.good.unique.good.filter.unique.precluster.pick.rdp.wang.tx.list, count=merged.good.unique.good.filter.unique.precluster.denovo.uchime.pick.count_table, label=1)"|srun -p NMLResearch --mem 5G -c 2 mothur
echo "classify.otu(list=merged.good.unique.good.filter.unique.precluster.pick.rdp.wang.tx.list, count=merged.good.unique.good.filter.unique.precluster.denovo.uchime.pick.count_table, taxonomy=merged.good.unique.good.filter.unique.precluster.pick.rdp.wang.taxonomy, label=1)"|srun -p NMLResearch --mem 5G -c 2 mothur

echo "move the final *tax.summary file (NOT cons.tax.summary) and *cons.taxomony file to a final output folder"
echo "done"
