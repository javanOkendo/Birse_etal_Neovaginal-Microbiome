####################################
#convert raw counts to proportions
####################################

getProp <- function(data){
  #for each column, 
  # each element is divided by the sum of that column
  for(i in 1:ncol(data)){
    data[ , i] <- (data[ , i] / sum(data[ , i]))
  }
  return(data)
}


###########################################################
#convert data to long format and determine hclust branches 
###########################################################

meltData <- function(data, n_branches = 1, subject_numbers = TRUE, dist = "euclidean"){
  
  #reverse each column so that bar colors are stacked in same order as legend
  data <- apply(data, 2, rev)
  
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


###########################
#Generate Painter's Plot
###########################

PPlot <- function(data, n_branches = 1, subject_numbers = TRUE, prop = FALSE, omit = "none", colors = "default", dist = "euclidean"){
  
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
  # pp <- ggplot(dl, aes(x = as.factor(Subject), y = Proportion, fill = Protein)) +
  #   geom_bar(stat = "identity") + xlab("Subject") + ylab("Microbial proportions") +
  #   scale_fill_manual(values = clrs, guide = guide_legend(title = " ")) + 
  #     theme_classic() + theme(legend.text = element_text(face = "italic"))
  
  pp <- ggplot(dl, aes(x = as.factor(Subject), y = Proportion, fill = Protein)) +
      geom_bar(stat = "identity") + xlab("Subject") + ylab("Microbial proportions") +
      scale_fill_manual(values = clrs, guide = guide_legend(title = " ")) + 
      theme_classic()

  
  # pp <- ggplot(dl, aes(x = as.factor(Subject), y = Proportion, fill = Protein)) +
  #     geom_bar(stat = "identity") + xlab("Subject") + ylab("Microbial proportions") +
  #     scale_fill_manual(values = clrs, guide = guide_legend(title = " ")) + theme_classic()
  
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


##########################################  
#view hierarchical clustering dendrograms
##########################################

plotClust <- function(data, subject_numbers = TRUE, prop = FALSE, dist = "euclidean"){
  
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

####################
#get summary by bug 
####################

getSummary <- function(data, prop = FALSE){
  
  #convert raw data into proportion
  if(prop == FALSE){
    data <- getProp(data)
  }
  
  #remove scientific notation
  options(scipen = 999)
  
  #get summary statistics
  # apply a function over the rows (1) of data
  m <- apply(data, 1, mean)       #mean
  me <- apply(data, 1, median)    #median
  mn <- apply(data, 1, min)       #min
  mx <- apply(data, 1, max)       #max
  sd <- apply(data, 1, sd)        #sd
  
  #combine all summary measures and add column names
  s <- cbind(m, me, mn, mx, sd)
  colnames(s) <- c("Mean", "Median", "Min", "Max", "SD")
  
  #return summary matrix
  return(s)
}


#########################
#get Shannon's H Bar Plot 
#########################

shanH <- function(data, n_branches = 1, prop = FALSE, omit = "none", dist = "euclidean"){
  
  #remove omitted bugs
  if(any(omit != "none")){
    rn <- rownames(data)         #data rownames
    x <- which(rn %in% omit)     #the row numbers of the omitted bugs
    data <- data[-(x), ]         #remove the rows of the omitted bugs
  }
  
  #get proportion data
  if(prop == FALSE){
    data <- getProp(data)
  }
  
  #conduct Shannon's H
  s <- diversity(t(data), "shannon")
  d <- data.frame(1:length(s), s)
  colnames(d) <- c("Subject", "S")
  
  #conduct hierarchical clustering
  if(dist == "euclidean"){
    hc <- hclust(dist(t(data)))
  }
  if(dist == "pearson"){
    hc <- hclust(dist(cor(data, method = "pearson")))
  }
  if(dist == "spearman"){
    hc <- hclust(dist(cor(data, method = "spearman")))
  }
  
  #add branch information to data
  b <- paste("Branch", as.numeric(cutree(hc, n_branches)))
  d <- cbind(d, b)
  colnames(d)[3] <- "Branch"
  
  #add clustering order to factors of data
  d$Subject <- factor(d$Subject, levels = hc$order)
  
  #use ggplot2 to make basic bar plot
  sh <- ggplot(data=d, aes(x = Subject, y = S)) + geom_bar(stat="identity", fill="black") + ylab("Shannon's H") + theme_classic()
  
  #if branches is not equal to one, this adds the different groupings
  # to the already existing plot
  if(n_branches != 1) {
    sh <- sh + facet_grid(~ Branch, scale="free_x", space = "free_x")
  }
  
  #return the bar plot
  return(sh)
}


#########################
#get HC branches 
#########################


branches <- function(data, n_branches = 1, prop = FALSE, subject_numbers = TRUE, omit = "none", dist = "euclidean"){
    
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






