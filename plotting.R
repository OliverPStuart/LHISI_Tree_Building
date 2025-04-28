
### File for plotting trees of phasmid + d australis

HOME_DIR="/Volumes/Alter/LHISI"
WORKING_DIR=paste0(HOME_DIR,"/Analyses/TreeBuilding")
setwd(WORKING_DIR)

library(tidytree)
library(ggplot2)
library(ggtree)
library(ape)
library(phytools)
library(phangorn)
library(treeio)
library(svglite)

### Read tree

tree <- read.nexus("aa_summarized/alignment/final_alignments/FullAlignments_AllSpecies.partition.treefile")
tree$tip.label <- gsub("_sp","_sp.",tree$tip.label)

# Subset species
include <- c("Timema_cristinae",
             "Diapherodes_gigantea",
             "Agamemnon_cornutus",
             "Tirachoidea_westwoodii",
             "Ramulus_artemis",
             "Medauroidea_extradentata",
             "Dimorphodes_sp",
             "Megacrania_phelaus",
             "Spinotectarchus_acornutus",
             "Clitarchus_hookeri",
             "Extatosoma_tiaratum",
             "Dryococelus_australis")

tree_sub <- drop.tip(tree,
                     tree$tip.label[-match(include,
                                           tree$tip.label)])

tree_sub <- root(tree_sub,
                outgroup="Timema_cristinae",
                 resolve.root = T)

tree_sub$tip.label <- gsub("_"," ",tree_sub$tip.label)

p <- ggtree(tree_sub,branch.length = "none") + 
  geom_tiplab(aes(subset=node!=1,label=label),offset = 0.1) + 
  geom_tiplab(aes(subset=node==1,label=label),fontface="bold",offset = 0.1,size=4.5) + 
  xlim(c(0,14)) +
  geom_tippoint(aes(subset=node==1),size=2) +
  geom_cladelab(node=19,label="Cladomorphinae",offset=3.8,extend=0.2,barsize=2, offset.text=0.1,barcolour="seagreen3",textcolour="seagreen3",fontface="bold") + 
  geom_cladelab(node=20,label="Clitumninae",offset=4.5,extend=0.2,barsize=2, offset.text=0.1,barcolour="palegreen3",textcolour="palegreen3",fontface="bold") + 
  geom_cladelab(node=17,label="Lanceocercata",offset=4.6,extend=0.2,barsize=2, offset.text=0.1,barcolour="palegreen4",textcolour="palegreen4",fontface="bold")

png("clade_tree.png",res=300,height=6,width=6.7,units="in")
plot(p)
dev.off()

# Another option
# Cladogram with just certain tips highlighted

# Subset species
include <- c("Timema_cristinae",
             "Dimorphodes_sp.",
             "Megacrania_phelaus",
             "Spinotectarchus_acornutus",
             "Clitarchus_hookeri",
             "Extatosoma_tiaratum",
             "Dryococelus_australis")
tree_sub <- drop.tip(tree,
                     tree$tip.label[-match(include,
                                           tree$tip.label)])
# Root tree by Timema
tree_sub <- root(tree_sub,
                 outgroup="Timema_cristinae",
                 resolve.root = T)

# Modify tip labels
tree_sub$tip.label <- gsub("_","\n",tree_sub$tip.label)

# Make tree 
p <- ggtree(tree_sub,
            branch.length = "none",
            layout="slanted",size=1,
            lineend="butt") +
  geom_tiplab(aes(subset=node!=1,label=label),
              size=2.5,
              fontface="italic",
              hjust=0.5, nudge_x=-0.4,
              lineheight=0.9) + 
  geom_tiplab(aes(subset=node==1,label=label),
              fontface="bold.italic",
              size=3,
              geom="label",
              fill="palegreen4", colour="white",
              hjust=0.5,nudge_x=-0.4,
              label.padding=unit(0.15,"lines"),
              lineheight=0.7) + 
  geom_rootedge(0.9,size=1) +
  layout_dendrogram() + 
  xlim_tree(xlim=c(-1,0.5)) +
  geom_nodepoint(aes(subset=node==8),size=2) + 
  geom_nodepoint(aes(subset=node==13),size=2)

# Add div.time data to certain nodes
p$data$div.time <- c(NA,NA,NA,NA,NA,NA,NA,
                     "121.8 mya (105.1–139.4)",
                     NA,NA,NA,NA,
                     "30.6 mya (24.3–37.3)")

# Add to plot
p1 <- p + geom_text(aes(label=div.time),
              nudge_y=0.05,nudge_x=0.2,
              size=2,hjust=0)

png("test_tree.png",res=300,width=6,height=2.5,units='in')
plot(p1)
dev.off()


# Another option
# Subset species
include <- c("Timema_cristinae",
             "Dimorphodes_sp.",
             "Megacrania_phelaus",
             "Spinotectarchus_acornutus",
             "Clitarchus_hookeri",
             "Extatosoma_tiaratum",
             "Dryococelus_australis")
tree_sub <- drop.tip(tree,
                     tree$tip.label[-match(include,
                                           tree$tip.label)])
# Root tree by Timema
tree_sub <- root(tree_sub,
                 outgroup="Timema_cristinae",
                 resolve.root = T)

# Modify tip labels
tree_sub$tip.label <- gsub("_","\n",tree_sub$tip.label)

# Make tree 
p <- ggtree(tree_sub,
            branch.length = "none",
            size=1) +
  geom_tiplab(aes(subset=node!=1,label=label),
              size=2.5,
              fontface="italic",
              hjust=0.5, nudge_x=-0.4,
              lineheight=0.9) + 
  geom_tiplab(aes(subset=node==1,label=label),
              fontface="bold.italic",
              size=3,
              geom="label",
              fill="palegreen4", colour="white",
              hjust=0.5,nudge_x=-0.4,
              label.padding=unit(0.15,"lines"),
              lineheight=0.7) + 
  geom_rootedge(0.9,size=1) +
  layout_dendrogram() + 
  xlim_tree(xlim=c(-1,0.5)) +
  geom_nodepoint(aes(subset=node==8),size=2) + 
  geom_nodepoint(aes(subset=node==13),size=2)

# Add div.time data to certain nodes
p$data$div.time <- c(NA,NA,NA,NA,NA,NA,NA,
                     "121.8 mya (105.1–139.4)",
                     NA,NA,NA,NA,
                     "30.6 mya (24.3–37.3)")

# Add to plot
p1 <- p + geom_text(aes(label=div.time),
                    nudge_y=0.05,nudge_x=0.25,
                    size=2,hjust=0)

png("test_tree2.png",res=300,width=6,height=2.5,units='in')
plot(p1)
dev.off()


svglite("test_tree2.svg", width=6,height=2.5)
plot(p1)
dev.off()
