#Functions for Neovagina Paper_16S analysis_R Markdown_Jan.14.2020.RMD
#Date: Jan 14 2020


# Remove '_unclassified' taxa from your phylo level of interest
remove.unclassified.taxa <- function (tax.summary.raw, genus.level.subset) {
    all.taxa <- genus.level.subset$taxon                         
    taxa.unclassified <- grep("unclassified", 
                              all.taxa, 
                              value = TRUE)  #selects only taxons with 'unclassified' in name. value=T returns the values in vector
    
    #Finds values from taxa.unclassified vector in taxon column of genus.level.subset, and returns everything that does NOT match the taxa.unclassified vector
    result <- filter(genus.level.subset, 
                     !grepl(paste(taxa.unclassified, 
                                  collapse="|"), 
                            taxon))  
    
    num.taxa.start <- length(genus.level.subset$taxon)     #number of taxa in genus.level.subset before removal
    num.taxa.removed <- length(taxa.unclassified)          #number of taxa removed 
    num.taxa.result <- length(result$taxon)                #number of taxa left after removal
    
    
    if((num.taxa.start-num.taxa.removed) == num.taxa.result) {
     cat("\n", num.taxa.start,"taxa to start",
         "\n", num.taxa.removed, "taxa unclassified at genus-level",
         "\n", num.taxa.result, "taxa retained")
      }  
    
    return(result)
  }

# Calculate fold difference of taxa in Samples vs. Negative Controls
taxa.rm.calc.table <- function (df_samples,      
                                df_neg,          
                                genus.ordered) {         

    # df_samples is a numeric df with sample names as colnames, taxa as rownames
    # df_neg is a numeric df with negative controls as colnames, taxa as rownames
    # genus.ordered is a numeric df with all samples and all controls
  
    #///REMOVAL OF TAXA FROM NEGATIVE CONTROLS TABLE/// 
    calc.table <- as.data.frame(matrix(ncol=10, 
                                       nrow=nrow(genus.ordered)))
    
    rownames(calc.table) <- rownames(genus.ordered)
    colnames(calc.table) <- c("Keep.or.Lose.Taxa",
                              "Mean.Reads.across.Samples",
                              "Total.Reads.in.Samples",
                              "Total.Reads.in.Samples.as.proportion",
                              "Total.Reads.in.Samples.as.percentage",
                              "Mean.Reads.across.Negative.Controls",
                              "Total.Reads.in.Negative.Controls",
                              "Total.Reads.in.Negative.Controls.as.proportion",
                              "Total.Reads.in.Negative.Controls.as.percentage",
                              "Fold.Change.(Sample - Negative)/Negative")
    # /// SAMPLES /// 
      # Avg reads per taxa across samples
      calc.table$Mean.Reads.across.Samples <- apply(df_samples, 1, mean) 
      
      # Total reads per taxa across samples
      calc.table$Total.Reads.in.Samples <- apply(df_samples, 1, sum) 
      
      # Relative abundance of each taxa across samples as a whole
        # i.e. 80% of reads from Samples in this experiment belong to Lactobacillus
      for(k in 1:nrow(calc.table)){
        sum.of.all.sample.reads <- sum(calc.table$Total.Reads.in.Samples)
        calc.table$Total.Reads.in.Samples.as.proportion[k] <- (calc.table$Total.Reads.in.Samples[k] / sum.of.all.sample.reads)
        calc.table$Total.Reads.in.Samples.as.percentage[k] <- (calc.table$Total.Reads.in.Samples.as.proportion[k]*100)
      }
    
    # /// NEGATIVE CONTROLS /// 
      # Avg reads per taxa across negative controls
      calc.table$Mean.Reads.across.Negative.Controls <- apply(df_neg, 1, mean)   
      
      # Total reads per taxa across negative controls
      calc.table$Total.Reads.in.Negative.Controls <- apply(df_neg, 1, sum) 
      
      # Relative abundance of each taxa across negative controls as a whole
        # i.e. 42% of reads from Negative Controls in this experiment belong to Delftia
      for(k in 1:nrow(calc.table)){
        sum.of.all.neg.reads <- sum(calc.table$Total.Reads.in.Negative.Controls)
        calc.table$Total.Reads.in.Negative.Controls.as.proportion[k] <- (calc.table$Total.Reads.in.Negative.Controls[k] / sum.of.all.neg.reads)
        calc.table$Total.Reads.in.Negative.Controls.as.percentage[k] <- (calc.table$Total.Reads.in.Negative.Controls.as.proportion[k]*100)
      }
  
    # /// METHOD FOR REMOVING NEG READS /// 
      
      # Calculate *fold change* of mean reads per taxa in samples
      # relative to means reads per taxa in negative controls
      # Note: if mean reads in negative controls == 0, result will be 'Inf', 
      # because cannot divide by 0.
      
      for (k in 1:nrow(calc.table)){
        s <- calc.table$Mean.Reads.across.Samples[k]
        n <- calc.table$Mean.Reads.across.Negative.Controls[k]
        fold.change <- (s-n)/n
        calc.table$`Fold.Change.(Sample - Negative)/Negative`[k] <- fold.change
      }
      
      # If fold change is <0, this means there were more reads of this taxa in 
      # negative controls than in real samples
      # Therefore, remove these taxa from analysis!
      
      for (k in 1:nrow(calc.table)) {
        calc.table$Keep.or.Lose.Taxa[k] <- ifelse(calc.table$`Fold.Change.(Sample - Negative)/Negative`[k] >=0, "Keep","Lose")
      }
  
      return(calc.table)
  }

# Pool the 3 triplicates for each biological sample
pool.samples <- function(table.to.pool) {
    # table.to.pool is a numeric df, where taxa are rownames and ptid are colnames
    
    samples1 <- colnames(table.to.pool)
    samples2 <- gsub("\\_.*","",samples1)              #trim '_SXXX'
    samples3 <- substr(samples2, 0, nchar(samples2)-1) #remove last char
    colnames(table.to.pool) <- samples3
    pooled.samples <- unique(samples3)
    pool.table <- as.data.frame(matrix(0, 
                                       ncol=length(pooled.samples), 
                                       nrow=nrow(table.to.pool),
                                       dimnames = list(rownames(table.to.pool),
                                                       pooled.samples)))
    
    for (k in 1:length(pooled.samples)) {
      current.sample <- pooled.samples[k]
      cols.to.sum <- grep(current.sample, colnames(table.to.pool))
      pool.table[,current.sample] <- rowSums(table.to.pool[,cols.to.sum])
    }
    
    return(pool.table)
  } 

# Bin low abundant taxa to 'Other'
binning <- function(table.to.bin) {
  #table.to.bin is a numeric df with taxa as rownames, samples as colnames
  
    binning.table <- as.data.frame(table.to.bin)
    binning.table$abundance <- rowSums(binning.table)
    binning.table <- binning.table[order(binning.table$abundance,decreasing = T),]
    
    # Taxa with abundance <0.0515% of total abundance are binned to 'Other'
    sum.of.all.sample.reads <- sum(binning.table$abundance)
    taxa.in.other <- rownames(binning.table)[which(binning.table$abundance < (sum.of.all.sample.reads*0.000515))]
    bin.indices <- which(binning.table$abundance < (sum.of.all.sample.reads*0.000515))
    
    # Identify individual low abundant taxa, remove 
    in.other <- binning.table[bin.indices,] 
    other <- colSums(in.other)
    binned.table <- binning.table[-bin.indices,]
    
    # Add as single 'Other' group
    binned.table <- rbind(binned.table,other)
    row.names(binned.table)[nrow(binned.table)] <- "Other"
    binned.table <- binned.table[order(binned.table$abundance, decreasing = T),]
    
    return(binned.table[, 1:ncol(binned.table)-1]) 
  } 


##Burgener lab Painter's plot functions:
# Convert raw counts to proportions
getProp <- function(data){
  #for each column, 
  # each element is divided by the sum of that column
  for(i in 1:ncol(data)){
    data[ , i] <- (data[ , i] / sum(data[ , i]))
  }
  return(data)
}

# Convert data to long format and determine hclust branches 
meltData <- function(data, n_branches = 1, subject_numbers = TRUE, 
                     dist = "euclidean"){
  
  #reverse each column so that bar colors are stacked in same order as legend
  #data <- apply(data, 2, rev)  #Jan 31 edit - try not reversing for 16S plot, since already in order
  
  #add subject number IDs
  if(subject_numbers == TRUE){
    colnames(data) <- 1:ncol(data)
  }
  
  #transpose data
  data <- t(data)
  
  
  ##main code to melt data into long format
  #subject IDs are replicated
  s <- rep(cbind(rownames(data)), ncol(data))
  
  #Protein or bug names are replicated
  p <- rep(cbind(colnames(data)), each = nrow(data))
  
  #convert data matrix into vector
  prop <- as.numeric(c(data))
  
  #combine data into dataframe and add column names
  dl <- data.frame(s, p, prop)
  colnames(dl) <- c("Subject", "Protein", "Proportion")
  
  #select clustering method
  
  if(dist == "euclidean"){
    hc <- hclust(dist(data))
  }
  if(dist == "pearson"){
    hc <- hclust(dist(cor(t(data), method = "pearson")))
  }
  if(dist == "spearman"){
    hc <- hclust(dist(cor(t(data), method = "spearman")))
  }
  
  #add branch information to long format data
  b <- paste("Branch", as.numeric(cutree(hc, n_branches)))
  b <- rep(b, ncol(data))
  dl <- cbind(dl, b)
  colnames(dl)[4] <- "Branch"
  
  #row order
  ro <- rownames(data)[hc$order]
  
  #use row order to change factor ordering
  # needed for corerct ordering of subjects in plots
  dl$Subject <- factor(dl$Subject, levels = ro)
  
  #add clustering order to factors of long data
  # dl$Subject <- factor(dl$Subject, levels = hc$order)
  
  
  return(dl)
}

# Generate Painter's Plot (clustered bar plot)
PPlot <- function(data, n_branches = 1, subject_numbers = FALSE, 
                  prop = FALSE, omit = "none", 
                  colors = "default", dist = "euclidean"){
  
  #default high contrast colors
  # clrs <- c("navy", "greenyellow", "salmon", "seagreen1", "#A30059","hotpink", "purple2", "orange1",
  #         "orange3", "orange4",
  #           "#5A0007", "red3", "#0000A6", "#63FFAC", "lightgoldenrod3", "royalblue4", "forestgreen", "yellow1",
  #           "#1CE6FF", "#809693", "#FEFFE6", "black", "steelblue4", "#3B5DFF", "#4A3B53", "#FF2F80",
  #           "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
  #           "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
  #           colors())
  
  #if user defines any other colors it will copy over default colors
  # while retaining the rest of the list
  if(any(colors != "default")){
    clrs[1:length(colors)] <- colors
  }
  
  #remove omitted bugs
  if(any(omit != "none")){
    rn <- rownames(data)         #data rownames
    x <- which(rn %in% omit)     #the row numbers of the omitted bugs
    data <- data[-(x), ]         #remove the rows of the omitted bugs
    clrs <- clrs[-(x)]           #remove the corresponding colors
  }
  
  #get proportion data
  if(prop == FALSE){
    data <- getProp(data)
  }
  
  
  #use meltData function
  dl <- meltData(data, n_branches, subject_numbers, dist)
  
  #Keep bacteria in order loaded in
  mo = 1:nrow(data)
  co = rownames(data)[mo]
  dl$Protein <- factor(dl$Protein, levels = co)
  
  #use ggplot2 to make basic bar plot
  pp <- ggplot(dl, aes(x = as.factor(Subject), y = Proportion, fill = Protein)) +
    geom_bar(stat = "identity") + xlab("Subject") + ylab("Microbial proportions") + 
    scale_fill_manual(values = clrs, guide = guide_legend(title = " "))
  
  #if branches is not equal to one, this adds the different groupings
  # to the already existing plot
  if(n_branches != 1) {
    pp <- pp + facet_grid(~ Branch, scale="free_x", space = "free_x")
  }
  
  #if default subject IDs are used, this makes the subject IDs in the vertical 
  # so they don't get too crowded in the graph
  if(subject_numbers == FALSE){
    pp <- pp + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  }
  
  #return the painter's plot
  return(pp)
}

# View hierarchical clustering dendrograms
plotClust <- function(data, subject_numbers = TRUE, 
                      prop = FALSE, dist = "euclidean"){
  
  #convert raw data into proportions
  if(prop == FALSE){
    data <- getProp(data)
  }
  
  #add subject number IDs
  if(subject_numbers == TRUE){
    colnames(data) <- 1:ncol(data)
  }
  
  #select clustering distance metric
  if(dist == "euclidean"){
    hc <- hclust(dist(t(data)))
  }
  if(dist == "pearson"){
    hc <- hclust(dist(cor(data, method = "pearson")))
  }
  if(dist == "spearman"){
    hc <- hclust(dist(cor(data, method = "spearman")))
  }
  
  #widen margins
  par(mar = c(5,5,5,10))
  
  #return plot
  return(plot(as.dendrogram(hc),horiz=FALSE))
  
  
  #turn graphics device off (resets margins)
  dev.off()
  
}

# Get hierarchical clutering branches
branches <- function(data, n_branches = 1, prop = FALSE, 
                     subject_numbers = TRUE, omit = "none", 
                     dist = "euclidean"){
  
  #add subject number IDs
  if(subject_numbers == TRUE){
    colnames(data) <- 1:ncol(data)
  }
  
  #remove omitted bugs
  if(any(omit != "none")){
    rn <- rownames(data)         #data rownames
    x <- which(rn %in% omit)     #the row numbers of the omitted bugs
    data <- data[-(x), ]         #remove the rows of the omitted bugs
  }
  
  if(prop == FALSE){
    data <- getProp(data)
  }
  
  #transpose data so that hclust clusters by subjects, not microbes
  data <- t(data)
  
  if(dist == "euclidean"){
    hc <- hclust(dist(data))
  }
  if(dist == "pearson"){
    hc <- hclust(dist(cor(t(data), method = "pearson")))
  }
  if(dist == "spearman"){
    hc <- hclust(dist(cor(t(data), method = "spearman")))
  }
  
  #add branch information to long format data
  b <- as.character(cutree(hc, n_branches))
  
  #column bind branch information and transposed data
  r <- cbind(b, data)
  
  #re-transpose data so it in same structure as original
  r <- t(r)
  r <- as.data.frame(r)
  
  rownames(r)[1] <- "Branch"
  
  #reorder data so that it is the same as the painter's plot
  r <- r[ , hc$order]
  
  #reorder data so that branches are order 1..n
  r <- r[ , order(r[1,])]
  
  return(r)
}


#### Colours vector for pplot ####

clrs <- c("Atopobium" = rgb(8/255, 20/255, 109/255),
          "Atopobium vaginae" = rgb(20/255, 16/255, 114/255),
          "Atopobium rimae" = rgb(12/255, 25/255, 113/255),
          "Atopobium parvulum" = rgb(8/255, 20/255, 109/255),
          "Achromobacter" = "bisque4",
          "Acidaminococcus" = rgb(38/255, 100/255, 90/255),
          "Acidaminococcus fermentans" = rgb(38/255, 100/255, 90/255),
          "Acidovorax" = rgb(40/255, 102/255, 92/255),
          "Acinetobacter gerneri" = rgb(49/255, 1/255, 120/255), 
          "Acinetobacter" = "#efae9a", 
          "Acinetobacter lwoffii" = "#c44236", 
          "Acinetobacter pittii" = "#da7258", 
          "Actinobacillus" = rgb(206/255, 88/255, 67/255),
          "Actinobacteria" = rgb(206/255, 82/255, 69/255),
          "Actinotignum" = "chocolate3",
          "Actinomyces" = rgb(200/255, 85/255, 60/255),
          "Actinomyces europaeus" = rgb(200/255, 85/255, 60/255),
          "Aerococcus" = rgb (205/255,173/255,0/255),
          "Aerococcus christensenii" = "gold3",
          "Agathobacter" = rgb(200/255, 2/255, 30/255),
          "Agathobacter rectalis" = rgb(200/255, 2/255, 30/255),
          "Alcaligenaceae_unclassified" = rgb(233/255, 235/255, 33/255),
          "Alloscardovia" = rgb(0/255,0/255,205/255),
          "Alloscardovia omnicolens" = "mediumblue",
          "Alicyclobacillus" = rgb(80/255, 100/255, 120/255),
          "Alicyclobacillus acidocaldarius" = rgb(80/255, 100/255, 120/255),
          "Anaerococcus" = "#072301",
          "Anaerococcus lactolyticus" = "green",
          "Anaerococcus hydrogenalis" = "green2",
          "Anaeroglobus" = rgb(190/255, 130/255, 20/255),
          "Anaeroglobus geminatus" = rgb(190/255, 130/255, 20/255),
          "Anaerosphaera" = rgb(200/255, 200/255, 170/255),
          "Aneurinibacillus" = rgb(20/255, 50/255, 30/255),
          "Azospirillum" = "#1f75fe", 
          "Bacillus" = rgb(59/255, 201/255, 70/255),
          #"Bacteroides" = rgb(225/255, 180/255, 116/255),
          "Bacteroides" = rgb(225/255, 154/255, 116/255),   
          "Bacteroides caccae" = rgb(225/255, 180/255, 116/255),
          "Bacteroidetes" = rgb(215/255, 181/255, 114/255), 
          "Bifidobacterium" = rgb(133/255, 237/255, 73/255),
          "Bifidobacterium adolescentis" = rgb(125/255, 237/255, 73/255), 
          "Bifidobacterium animalis" = rgb(120/255, 236/255, 72/255),
          "Bifidobacterium angulatum" = rgb(132/255, 236/255, 72/255),
          "Bifidobacterium asteroides" = rgb(131/255, 235/255, 71/255),
          "Bifidobacterium biavatii" = rgb(130/255, 234/255, 70/255),
          "Bifidobacterium bifidum" = rgb(129/255, 233/255, 69/255),
          "Bifidobacterium bohemicum" = rgb(128/255, 234/255, 70/255),
          "Bifidobacterium boum" = rgb(122/255, 234/255, 70/255),
          "Bifidobacterium boum" = rgb(122/255, 234/255, 70/255),
          "Bifidobacterium commune" = rgb(117/255, 230/255, 62/255),
          "Bifidobacterium cuniculi" = rgb(130/255, 230/255, 62/255),
          "Bifidobacterium gallicum" = rgb(117/255, 215/255, 61/255),
          "Bifidobacterium kashiwanohense" = rgb(120/255, 212/255, 65/255),
          "Bifidobacterium longum" = rgb(116/255, 213/255, 60/255),
          "Bifidobacterium magnum" = rgb(115/255, 212/255, 59/255),
          "Bifidobacterium moukalabense" = rgb(111/255, 230/255, 62/255),
          "Bifidobacterium mongoliense" = rgb(50/255, 212/255, 59/255),
          "Bifidobacterium pseudocatenulatum" = rgb(119/255, 209/255, 54/255),
          "Bifidobacterium psychraerophilum" = rgb(111/255, 210/255, 55/255),
          "Bifidobacterium ruminantium" = rgb(112/255, 210/255, 55/255),
          "Bifidobacterium stellenboschense" = rgb(109/255, 208/255, 53/255),
          "Bifidobacterium subtile" = rgb(103/255, 209/255, 43/255),
          "Bifidobacterium thermophilum" = rgb(200/255, 230/255, 62/255),
          "Bifidobacterium tsurumiense" = rgb(102/255, 208/255, 42/255),
          "Blautia" = rgb(99/255, 228/255, 113/255),
          "Bordetella" = rgb(138/255, 20/255, 14/255),
          "Bordetella ansorpii" = rgb(138/255, 20/255, 14/255),
          "Bradyrhizobium" = rgb(153/255, 71/255, 123/255),
          "Brevundimonas" =  rgb(30/255,30/255, 90/255),             
          "Brevibacterium" =  rgb(30/255,30/255, 90/255),              
          "Brucellaceae_unclassified" = rgb(158/255, 81/255, 123/255), 
          "Bulleidia" = rgb(78/255, 71/255, 114/255),
          "Burkholderia" = rgb(190/255, 87/255, 200/255),
          "Butyrivibrio" = rgb(192/255, 173/255, 219/255),
          "Butyrivibrio crossotus" = rgb(193/255, 172/255, 220/255),
          "Butyrivibrio fibrisolvens" = rgb(192/255, 171/255, 217/255),
          "Butyrivibrio proteoclasticus" = rgb(191/255, 170/255, 216/255),
          "Campylobacter" = "chocolate4",
          "Catenibacterium" = rgb(20/255, 60/255, 40/255),
          "Catenibacterium mitsuokai" = rgb(20/255, 60/255, 40/255),
          "Collinsella" = rgb(48/255, 127/255, 140/255),
          "Collinsella aerofaciens" = rgb(52/255, 140/255, 155/255),
          "Clostridium" = rgb(252/255, 196/255, 217/255),
          "Clostridium " = rgb(252/255, 196/255, 217/255),
          "Clostridium_XlVa" = rgb(252/255, 196/255, 217/255),
          "Clostridium nexile" = rgb(251/255, 194/255, 219/255),
          "Clostridium magnum" = rgb(241/255, 198/255, 220/255),
          "Christensenella" = rgb(23/255, 190/255, 22/255),
          "Christensenella minuta" = rgb(23/255, 190/255, 22/255),
          "Chlamydia" = "grey",
          "Chlamydia trachomatis" = "grey",
          "Comamonadaceae_unclassified" = rgb(67/255, 99/255, 89/255),
          "Coprococcus" = rgb(66/255, 101/255, 24/255),
          "Coriobacteriaceae_unclassified" = rgb(99/255, 99/255, 19/255),
          "Corynebacterium" = rgb(92/255, 92/255, 92/255),
          "Corynebacteriaceae_unclassified" = rgb(94/255, 92/255, 92/255),
          "Cricetibacter" = rgb(88/255, 90/255, 149/255),
          "Criibacterium bergeronii" = rgb(90/255, 90/255, 150/255),
          "Criibacterium" = rgb(90/255, 90/255, 150/255),
          "Dialister" = "seagreen3",
          "Dialister invisus" = "seagreen2", 
          "Dialister succinatiphilus" = "seagreen4",
          "Dialister microaerophilus" = "seagreen1",
          "Dialister succinatiphilus" = "seagreen2",
          "Delftia" = "steelblue1",
          "Enterobacteriaceae_unclassified" = "#41f4c7",
          "Enterococcaceae_unclassified" = "yellow4",
          "Enterococcus" = "yellow3",
          "Escherichia" = rgb(163/255, 0/255, 89/255),
          "Escherichia/Shigella" = rgb(163/255, 0/255, 89/255),
          "Escherichia coli" = rgb(178/255, 1/255, 98/255),
          "Escherichia fergusonii" = rgb(177/255, 2/255, 97/255),
          "Escherichia hermannii" = rgb(160/255, 1/255, 88/255),
          "Escherichia vulneris" = rgb(186/255, 1/255, 102/255),
          "Eubacterium" = rgb(76/255, 111/255, 249/255),
          "[Eubacterium]" = rgb(76/255, 111/255, 249/255),
          "[Eubacterium] eligens" = rgb(79/255, 113/255, 250/255),
          "Eubacterium rectale" = rgb(76/255, 111/255, 249/255),
          "Ezakiella" = rgb(176/255, 104/255, 119/255),
          "Faecalibacterium" = rgb(202/255, 202/255, 3/255),
          "Faecalibacterium prausnitzii" = rgb(202/255, 202/255, 3/255),
          "Fastidiosipila" = rgb(128/255, 150/255, 160/255),
          "Finegoldia" = rgb(121/255, 47/255, 99/255),
          "Finegoldia magna" = rgb(121/255, 47/255, 99/255),
          "Firmicutes" = rgb(5/255, 90/255, 220/255),
          "Firmicutes bacterium" = rgb(5/255, 90/255, 220/255),
          "Flavobacterium" = "royalblue3",
          "Fusobacterium" = "khaki",
          "Fusobacterium gonidiaformans" = "khaki1",
          "Fusobacterium hwasookii" = "khaki2",
          "Fusobacterium nucleatum" = "khaki3",
          "Fusobacterium periodonticum" = "khaki1",
          "Gardnerella" = "purple2",
          "Gardnerella vaginalis" = "purple2",
          "Gemella" = rgb(205/255, 199/255, 172/255),
          "Haemophilus" = rgb(76/255, 112/255, 10/255),
          "Helcococcus" = rgb(246/255, 168/255, 179/255),
          "Howardella" = rgb(50/255, 108/255, 100/255),
          "Jonquetella" = rgb(200/255, 200/255, 170/255),
          "Jonquetella anthropi" = rgb(199/255, 199/255, 165/255),
          "Lachnospiraceae bacterium" = rgb(83/255, 99/255, 10/255), 
          "Lachnospiraceae" = rgb(83/255, 99/255, 10/255),
          "Lachnospiraceae_unclassified" = rgb(83/255, 99/255, 10/255),
          "Lactobacillaceae_unclassified" = rgb(223/255, 162/255, 82/255),
          "Lactobacillus" = rgb(249/255, 200/255, 94/255),
          "Lactobacillus spp." = rgb(249/255, 200/255, 94/255),
          "Lactobacillus acetotolerans" = rgb(230/255, 170/255, 90/255),
          "Lactobacillus acidophilus" = "#D1AF04",
          "Lactobacillus amylolyticus" = rgb(229/255, 169/255, 91/255),
          "Lactobacillus amylovorus" = "#AA5A4A",
          "Lactobacillus apis" = rgb(227/255, 167/255, 89/255),
          "Lactobacillus bifermentans" = rgb(226/255, 163/255, 81/255),
          "Lactobacillus cacaonum" = rgb(225/255, 161/255, 80/255),
          "Lactobacillus crispatus" = "orange1",
          "Lactobacillus coleohominis" = rgb(250/255, 180/255, 100/255),
          "Lactobacillus curieae" = "#FBFB73",
          "Lactobacillus delbrueckii" = rgb(200/255, 120/255, 40/255),
          "Lactobacillus diolivorans" = rgb(220/255, 160/255, 80/255),
          "Lactobacillus equicursoris" = "#772210", 
          "Lactobacillus equi" = "#772210", 
          "Lactobacillus frumenti" = rgb(216/255, 150/255, 65/255),
          "Lactobacillus gallinarum" = rgb(215/255, 149/255, 63/255),
          "Lactobacillus gasseri" = "#FFB266",
          "Lactobacillus gigeriorum" = rgb(224/255, 149/255, 53/255),
          "Lactobacillus ginsenosidimutans" =rgb(226/255, 169/255, 91/255),
          "Lactobacillus hamsteri" = rgb(229/255, 169/255, 91/255),
          "Lactobacillus helveticus" = "#FBE673",
          "Lactobacillus hominis" = rgb(229/255, 169/255, 91/255),
          "Lactobacillus iners" = "#FB7900",
          "Lactobacillus intestinalis" = "#fec716",
          "Lactobacillus jensenii" = "#FCBC9C",
          "Lactobacillus johnsonii" = rgb(228/255, 165/255, 91/255),
          "Lactobacillus kalixensis" = rgb(229/255, 170/255, 92/255),
          "Lactobacillus kefiranofaciens" = rgb(251/255, 160/255, 44/255),
          "Lactobacillus kitasatonis" = rgb(250/255, 161/255, 44/255),
          "Lactobacillus mellifer" = rgb(253/255, 158/255, 40/255),
          "Lactobacillus mindensis" = rgb(248/255, 159/255, 42/255),
          "Lactobacillus paracasei" = rgb(249/255, 157/255, 48/255),
          "Lactobacillus pasteurii" = rgb(250/255, 159/255, 43/255),
          "Lactobacillus perolens" = rgb(249/255, 158/255, 42/255),
          "Lactobacillus psittaci" = rgb(250/255, 158/255, 42/255),
          "Lactobacillus reuteri" = rgb(206/255, 117/255, 4/255),
          "Lactobacillus rhamnosus" = rgb(205/255, 115/255, 3/255),
          "Lactobacillus rossiae" = rgb(204/255, 113/255, 4/255),
          "Lactobacillus senmaizukei" = "#EEDDB7",
          "Lactobacillus selangorensis" = rgb(221/255, 126/255, 4/255),
          "Lactobacillus sharpeae" = rgb(220/255, 124/255, 4/255),
          "Lactobacillus taiwanensis" = rgb(218/255, 123/255, 3/255),
          "Lactobacillus ultunensis" = rgb(217/255, 122/255, 3/255),
          "Lactobacillus vaginalis" = rgb(216/255, 121/255, 4/255),
          "Legionella spiritensis" = "#cbb9ff",
          "Listeria" = "cyan",
          "Massilia" = rgb(204/255, 255/255, 29/255),
          "Megasphaera" = rgb(90/255, 0/255, 7/255),
          "Megasphaera cerevisiae" = rgb(114/255, 1/255, 9/255),
          "Megasphaera elsdenii" = rgb(112/255, 0/255, 8/255),
          "Megasphaera genomosp." = rgb(155/255, 1/255, 12/255),
          "Megasphaera genomosp " = rgb(155/255, 1/255, 12/255),
          "Mesocricetibacter" = rgb(153/255, 81/255, 12/255),
          "Microbacteriaceae_unclassified" = rgb(6/255, 195/255, 209/255),
          "Mitsuokella" = rgb(11/255, 202/255, 80/255),
          "Mobiluncus" = "Red3",
          "Mobiluncus curtisii" = "Red1",
          "Mobiluncus mulieris" = "Red",
          "Moryella" = rgb(113/255, 54/255, 160/255),
          "Mycoplasma" = "orangered2",
          "Myxococcus" = "#ffb9ed",
          "Myxococcus stipitatus"= "#fcb4d5",
          "none" = "white",                              #'none' would mean 'no bacteria detected'
          "Oribacterium" = rgb(11/255, 202/255, 88/255), 
          "Other" = "darkgray",
          "Paenibacillus" = "#523a73",
          "Paenibacillus beijingensis" = "#71384c", 
          "Paenibacillus polymyxa" = "#926eae", 
          "Paeniclostridium sordellii" = rgb(100/255, 50/255, 150/255),
          "Paeniclostridium" = rgb(100/255, 50/255, 150/255),
          "Parasutterella" = rgb(80/255, 50/255, 144/255),
          "Parvimonas" = rgb(80/255, 100/255, 128/255),
          "Pasteurella" = "brown",
          "Pasteurellaceae_unclassified" = "brown2",
          "Pelomonas" = rgb(116/255, 225/255, 216/255),    
          "Peptoniphilus" = rgb(18/255, 145/255, 130/255),
          "Peptoniphilus duerdenii" = rgb(50/155, 143/255, 129/255),
          "Peptoniphilus harei" = rgb(19/155, 143/255, 129/255),
          "Peptoniphilus lacrimalis" = rgb(22/255, 142/255, 127/255),
          "Peptostreptococcus" = rgb(51/255, 201/255, 204/255), 
          "Peptostreptococcus anaerobius" = rgb(49/255, 200/255, 201/255),
          "Phyllobacterium" = rgb(188/255, 80/255, 241/255),  
          "Planococcus" = rgb(190/255, 250/255, 79/255),
          "Porphyromonadaceae_unclassified" = "royalblue4",
          "Porphyromonas" = "royalblue4",
          "Porphyromonas uenonis" = "royalblue3",
          #"Prevotella" = rgb(13/255, 89/255, 27/255), old Prevotella clr 
          "Prevotella" = "#25AF06", #new Prevotella clr (June 2019)
          "Prevotellaceae_unclassified" = rgb(18/255, 89/255, 27/255), 
          "Prevotella amnii" = rgb(13/255, 89/255, 27/255),
          "Prevotella baroniae" = rgb(14/255, 91/255, 28/255),
          "Prevotella bergensis" = rgb(9/255, 84/255, 22/255),
          "Prevotella bivia" = rgb(8/255, 82/255, 25/255),
          "Prevotella buccalis" = rgb(4/255, 73/255, 16/255),
          "Prevotella bryantii" = rgb(6/255, 73/255, 16/255),
          "Prevotella buccae" = rgb(7/255, 73/255, 16/255),
          "Prevotella copri" = rgb(5/255, 77/255, 20/255),
          "Prevotella corporis" = rgb(6/255, 78/255, 21/255),
          "Prevotella dentalis" = rgb(3/255, 70/255, 14/255),
          "Prevotella denticola" = rgb(10/255, 61/255, 19/255),
          "Prevotella disiens" = rgb(11/255, 63/255, 22/255),
          "Prevotella fusca" = rgb(10/255, 63/255, 22/255),
          "Prevotella histicola" = rgb(12/255, 61/255, 21/255),
          "Prevotella intermedia" = rgb(13/255, 61/255, 21/255),
          "Prevotella loescheii" = rgb(12/255, 62/255, 21/255),
          "Prevotella maculosa" = rgb(9/255, 58/255, 19/255),
          "Prevotella marshii" = rgb(9/255, 59/255, 19/255),
          "Prevotella melaninogenica" = rgb(1/255, 53/255, 11/255),
          "Prevotella micans" = rgb(1/255, 50/255, 8/255),
          "Prevotella multiformis" = rgb(20/255, 84/255, 32/255),
          "Prevotella multisaccharivorax" = rgb(21/255, 82/255, 30/255),
          "Prevotella oralis" = rgb(21/255, 83/255, 30/255),
          "Prevotella oris" = rgb(22/255, 87/255, 35/255),
          "Prevotella oulorum" = rgb(18/255, 82/255, 30/255),
          "Prevotella pallens" = rgb(26/255, 68/255, 34/255),
          "Prevotella pleuritidis" = rgb(26/255, 68/255, 35/255),
          "Prevotella saccharolytica" = rgb(22/255, 70/255, 32/255),
          "Prevotella scopos" = rgb(22/255, 71/255, 32/255),
          "Prevotella stercorea" = rgb(31/255, 135/255, 50/255),
          "Prevotella salivae" = "#91CD7D",
          "Prevotella timonensis" = rgb(32/255, 134/255, 49/255),
          "Prevotella_7" = rgb(13/255, 89/255, 27/255),
          "Propionibacterium" = rgb(157/255, 216/255, 30/255),
          "Pseudoalteromonas" = rgb(14/255, 14/255, 250/255),
          "Pseudoramibacter" = rgb(230/255, 220/255, 90/255),
          "Pseudoramibacter alactolyticus" = rgb(230/255, 220/255, 90/255),
          "Pseudomonas" = "yellow1",
          "Pseudomonas aeruginosa" = rgb(239/255, 239/255, 81/255),
          "Pseudomonas alcaligenes" = rgb(238/255, 238/255, 80/255),
          "Pseudomonas batumici" = rgb(232/255, 232/255, 85/255),
          "Pseudomonas bauzanensis" = rgb(237/255, 237/255, 79/255),
          "Pseudomonas cichorii" = rgb(219/255, 219/255, 63/255),
          "Pseudomonas extremaustralis" = rgb(226/255, 226/255, 40/255),
          "Pseudomonas fluorescens" = rgb(225/255, 225/255, 42/255),
          "Pseudomonas fulva" = rgb(223/255, 225/255, 42/255),
          "Pseudomonas kilonensis" = rgb(147/255, 148/255, 17/255), 
          "Pseudomonas lundensis" = rgb(196/255, 196/255, 21/255), 
          "Pseudomonas orientalis" = rgb(197/255, 197/255, 24/255),
          "Pseudomonas oryzihabitans" = rgb(198/255, 198/255, 27/255),
          "Pseudomonas plecoglossicida" = rgb(147/255, 147/255, 17/255), 
          "Pseudomonas psychrophila" = rgb(144/255, 144/255, 15/255),
          "Pseudomonas putida" = rgb(142/255, 142/255, 12/255),
          "Pseudomonas syringae" = rgb(179/255, 178/255, 12/255),
          "Pseudomonas synxantha" = rgb(178/255, 178/255, 12/255),
          "Pseudomonas trivialis" = rgb(175/255, 175/255, 14/255),
          "Pseudomonas weihenstephanensis" = rgb(239/255, 239/255, 2/255),
          "Rhizobium" = rgb(19/255, 178/255, 19/255),
          "Rhodanobacter" = rgb(29/255, 132/255, 74/255),
          "Rhodococcus" = rgb(22/255, 152/255, 57/255),
          "Roseburia" = rgb(178/255, 234/255, 232/255),
          "Roseburia inulinivorans" = rgb(177/255, 234/255, 231/255),
          "Roseburia faecis" = rgb(178/255, 233/255, 229/255),
          "Rothia" = "magenta3",
          "Ruegeria" = rgb(206/255, 39/255, 192/255),
          "Ruminococcaceae" = rgb(129/255, 150/255, 146/255),
          "Ruminococcaceae bacterium" = rgb(129/255, 150/255, 146/255),
          "Ruminococcaceae_unclassified" = rgb(129/255, 150/255, 146/255),
          "Ruminococcus" = rgb(128/255, 150/255, 147/255),
          "Ruminococcus albus" = rgb(127/255, 149/255, 145/255), 
          "Ruminococcus bicirculans" = rgb(128/255, 150/255, 147/255),
          "Ruminococcus bromii" = rgb(130/255, 152/255, 146/255),
          "Ruminococcus champanellensis" = rgb(67/255, 79/255, 77/255),
          "Ruminococcus lactaris" = rgb(30/255, 71/255, 68/255),
          "Ruminococcus obeum" = rgb(56/255, 71/255, 68/255),
          "Saccharofermentans" = rgb(66/255, 17/255, 212/255),
          "Salmonella" = "lightsalmon",
          "Selenomonas" = rgb(2/255, 2/255, 90/255),
          "Selenomonas ruminantium" = rgb(2/255, 2/255, 90/255),
          "Shigella" = rgb(163/255, 0/255, 89/255),
          "Shuttleworthia" = rgb(254/255, 255/255, 230/255),
          "Shuttleworthia satelles" = rgb(253/255, 254/255, 229/255),
          "Sneathia" = "black",
          "Sneathia amnii" = "black",
          "Sneathia sanguinegens" = "black",
          "Staphylococcus" = rgb(81/255, 99/255, 45/255),
          "Staphylococcaceae_unclassified" = rgb(81/255, 100/255, 44/255),
          "Stomatobaculum" = rgb(80/255, 90/255, 60/255),
          "Stomatobaculum longum" = rgb(80/255, 90/255, 60/255),
          "Streptococcus" = rgb(85/255, 136/255, 178/255),
          "Streptococcaceae_unclassified" = rgb(83/255, 139/255, 178/255),
          "Streptococcus agalactiae" = rgb(84/255, 135/255, 176/255),
          "Streptococcus anginosus" = rgb(83/255, 133/255, 175/255),
          "Streptococcus constellatus" = rgb(69/255, 112/255, 147/255),
          "Streptococcus iniae" = rgb(69/255, 110/255, 145/255),
          "Streptococcus macacae" = rgb(69/255, 111/255, 145/255),
          "Streptococcus mitis" = rgb(67/255, 110/255, 145/255),
          "Streptococcus oralis" = rgb(69/255, 110/255, 146/255),
          "Streptococcus parasanguinis" = rgb(67/255, 109/255, 143/255),
          "Streptococcus lutetiensis" = rgb(69/255, 108/255, 142/255),
          "Streptococcus pneumoniae" = rgb(56/255, 95/255, 127/255),
          "Streptococcus sinensis" = rgb(53/255, 93/255, 126/255),      
          "Streptococcus uberis" = rgb(54/255, 93/255, 126/255),
          "Streptococcus urinalis" = rgb(52/255, 90/255, 125/255),
          "Stenotrophomonas" = "lawngreen",
          "Succinatimonas" = rgb(70/255, 90/255, 120/255),
          "Sutterella" = rgb(0/255, 0/255, 204/255),
          "Tumebacillus" = rgb(140/255, 120/255, 220/255),
          "Veillonellaceae_unclassified" = rgb(155/255, 23/255, 125/255),
          "Veillonella" = rgb(157/255, 23/255, 127/255),
          "Vibrio" = "#00d68f", 
          "Vibrio fluvialis" = "#b2f4e1", 
          "Vibrio metschnikovii" = "#41f4b5", 
          "Vibrio navarrensis" = "#41f4a3", 
          "Vibrio parahaemolyticus" = "#41f485",
          "Vibrio ponticus" = "#41f4c7", 
          "Xanthomonadaceae_unclassified" = rgb(212/255, 110/255, 21/255),
          "Ureaplasma" = rgb(242/255, 134/255, 100/255),
          "uncultured" = rgb(212/255, 112/255, 212/255),
          "uncultured Anaerotruncus" = rgb(200/255, 110/255, 200/255),
          "uncultured Lactobacillus" = rgb(205/255, 109/255, 205/255),
          "uncultured Clostridium" = rgb(212/255, 112/255, 212/255),
          "undistinguishable" = "#b1bfbb",
          "Undistinguishable" = "#b1bfbb",
          "No bacteria detected" = "#faf7ff")

