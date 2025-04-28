setwd("/Volumes/Alter/LHISI/Analyses/TreeBuilding/aa_summarized/alignment/")

# Packages
library(ape)
library(phangorn)
library(seqinr)
library(ggplot2)
library(ggtree)


# List files
alignments <- list.files(pattern="*.fa")

# Get taxa table
taxa <- read.table("../../taxonomy_table.txt",header=T,stringsAsFactors = F,sep="\t")
taxa$Group <- ifelse(taxa$Order=="Phasmatodea","Ingroup","Outgroup")

# Make empty data.frame to store values
results <- data.frame(alignment=character(),
                      d_aus_z=numeric(),
                      branch_ratio=numeric(),
                      second_to_first=numeric(),
                      longest=numeric(),
                      mean=numeric(),
                      longest_to_mean=numeric(),
                      phasmid_mono=logical())

# Loop over all alignments

for(i in 1:length(alignments)){
  
# Data
test <- read.phyDat(alignments[i],
                    format="fasta",type="AA")

# Make tree, get tips
dist <- dist.ml(test,model="JC69")
tree <- NJ(dist)
tips <- tree$tip.label

# Get terminal branch lengths
## first get the node numbers of the tips
nodes<-sapply(tips,function(x,y) which(y==x),y=tree$tip.label)
## then get the edge lengths for those nodes
edge.lengths<-setNames(tree$edge.length[sapply(nodes,
                                               function(x,y) which(y==x),y=tree$edge[,2])],names(nodes))
# Now turn this into a data.frame for ggplot
# Merge it with the taxonomy table
data <- merge(data.frame(Name = names(edge.lengths),
           length = as.vector(edge.lengths)),
           taxa,by="Name")

# Now let's check a bunch of things

# Get a Z-score for the branch length of D.australis with respect to the distribution of ingroup branch lengths
lengths <- data$length[data$Group=="Ingroup"]
d_aus_z <- (data$length[data$Name=="Dryococelus_australis"] - mean(lengths)) / sd(lengths)

# How much of the tree is taken up by internal branches
# If it's lots... then paralogs maybe
ratio <- sum(tree$edge.length) / sum(data$length)

# Distance between longest branch length and second longest
lengths <- sort(tree$edge.length)
second_to_first <- max(lengths) - lengths[(length(lengths)-1)]

# Longest branch length
longest <- max(lengths)

# Mean branch length
mean_b <- mean(lengths)

# And their ratio
long_mean <- longest / mean_b

# Is phasmida monophyletic?
phasm <- is.monophyletic(phy=tree,tips=data$Name[data$Group=="Ingroup"])

results[i,1] <- alignments[i]
results[i,2] <- d_aus_z
results[i,3] <- ratio
results[i,4] <- second_to_first
results[i,5] <- longest
results[i,6] <- mean_b
results[i,7] <- long_mean
results[i,8] <- phasm


}

ggplot(results,aes(x=d_aus_z)) + geom_histogram(bins=50)
ggplot(results,aes(x=branch_ratio)) + geom_histogram(bins=50)
ggplot(results,aes(x=second_to_first)) + geom_histogram(bins=50)
ggplot(results,aes(x=longest)) + geom_histogram(bins=50)
ggplot(results,aes(x=mean)) + geom_histogram(bins=50)
ggplot(results,aes(x=longest_to_mean)) + geom_histogram(bins=50)


# Get upper and lower bounds for the branch_ratio
upper <- mean(results$branch_ratio) + sd(results$branch_ratio)
lower <- mean(results$branch_ratio) - sd(results$branch_ratio)

# Any of these conditions fulfilled is suspicious
results$suspicious <- results$d_aus_z > 2 | 
  results$d_aus_z < -2 | 
  results$branch_ratio < lower | 
  results$branch_ratio > upper | 
  results$second_to_first > 0.25 | 
  results$longest_to_mean > 25 | 
  results$longest > 0.5 | 
  results$mean > 0.625 | 
  results$phasmid_mono == F


# Write the names of the suspicious alignments
writeLines(results$alignment[results$suspicious==T],"suspiciousAlignments.txt")

# Also print a tree for all suspicious alignments
setwd("../..")
dir.create("suspicious_trees")

suspicious_trees <- results$alignment[results$suspicious==T]

for(i in 1:length(suspicious_trees)){
  
  # Data
  test <- read.phyDat(paste0("aa_summarized/alignment/",suspicious_trees[i]),
                      format="fasta",type="AA")
  
  # Make tree, get tips
  dist <- dist.ml(test,model="JC69")
  tree <- NJ(dist)

  setwd("suspicious_trees")
  png(file=paste0(suspicious_trees[i],".tree.png"),res=300,width=10,height=10,units="in")
  plot(ggtree(tree) + geom_tiplab())
  dev.off()
  setwd("..")
}
