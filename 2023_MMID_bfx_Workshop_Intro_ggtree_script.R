## Introduction to ggTree
## April 13th, 2023
## by Taylor Davedow
##********************************************************

# Data used in this workshop was modified from
# Jesser et al. 2019 
# doi: 10.3389/fpubh.2019.00066

# ggtree book resource: https://yulab-smu.top/treedata-book/ 

### Getting set up
# create a folder called "Workshop_ggtree" and within that create sub directories
# for "data", "output", "script" and "docs"
# set Workshop_ggtree as the working directory
# save script into "script"
# save workshop slides into "docs"
# "output" is where we will send tree output
# download the newick file and the one xlsx metadata file and store them in "data"


##********************************************************
### Install and load packages ----
##********************************************************
# You can use different methods to install packages
#   A) Install button in the Packages Window
#
#   B) Select from menu bar Tools > Install Packages > type "package name" > Install
#  
#   C) Remove the "#" and run the following 3 lines of code.
# 
# install.packages("readxl")
# install.packages("BiocManager")
# install.packages("treeio")
# install.packages("tidyverse")
# install.packages("phytools")

# ggtree package must be installed after installing and loading BiocManager

library(BiocManager)

BiocManager::install("ggtree")

library(readxl) # for reading in xl files
library(ggtree) # for building tree
library(treeio) # for read.newick function
library(phytools) # for midpoint.root (also has read.newick option)
library(tidyverse) 

# we also will be using ggplot2 which should automatically load in
# with ggtree


##********************************************************
##  Load in files ----
##********************************************************

# tree file
tree <- read.newick("data/msa.fasta.tree")

# metadata file 
metadata <- read_xlsx("data/metadata.xlsx")


##********************************************************
##  Create a basic tree ----
##********************************************************

# basic tree
ggtree(tree)

# or
 
# tree %>% ggtree()

# check the help page for ggtree to find further usages, arguments and
# references to explore!
?ggtree

##********************************************************
## Changing the tree layout ----
##********************************************************
##*
# LAYOUTS

ggtree(tree,
       layout = "rectangular") # did anything change? this is the default

ggtree(tree, 
       layout = "circular") 

ggtree(tree, 
       layout = "roundrect") 

# be careful with branch length
# setting branch.length = "none" will create a cladogram

ggtree(tree,
       layout = "rectangular", branch.length = "none")

# other options to consider: rooting 
# try changing the position of the root node using root.position argument

# adding a midpoint.root (part of phytools package)
ggtree(midpoint.root(tree))

# Identifying nodes
# we have to check the node numbers so we can refer to them
# when using certain functions

ggtree(tree) + 
  geom_text2(aes(subset=!isTip, 
                 label=node), 
             hjust = -.3)


# tree manipulation

# many different functions for tree manipulation
# example: scale, collapse, expand, flip or group clades
# we can also zoom in to a particular clade (viewClade function)


# zoom in to a clade and show only clade of interest
ggtree(tree) %>% 
  viewClade(node = 17)

# or show both the original and zoomed in tree
ggtree(tree) + 
  viewClade(node = 17)

# scale
ggtree(tree) %>% 
  scaleClade(node = 27, 
             scale = 5)

# collapse
ggtree(tree) %>% 
  collapse(node = 17)



##********************************************************
## Adding and customizing labels ----
##********************************************************

# basic tip labels
ggtree(tree) + 
  geom_tiplab(size = 4) + # displaying tip labels 
  coord_cartesian(clip = 'off')+ # allows us to draw outside the plot
  theme(plot.margin = margin(1,3,1,1, "cm")) # add space around the plot

# if we want to check tree tip labels
# this gives us an overview of what the tip labels will look like
# looks like the tip labels are the SRA numbers

head(tree$tip.label) 

colnames(metadata)

# since we are planning on linking the metadata with
# the tree, we need to make sure they have a variable that links to the tip
# tip label names. the biosample_id variable will create this link
# we can also use a vector to check if there are any biosample_id
# observations that are not in the tree
metadata$biosample_id[!tree$tip.label %in% metadata$biosample_id]

# character(0) means they all match match up

# what if we want to switch out the tip label for something more
# meaningful to us?
# the new operator %<+% is discussed in section 7.1 of the book

ggtree(tree) %<+% # operator used to attach annotation data to tree
  metadata + # our metadata 
  geom_tiplab(aes(label = sample_id)) + # change the tip label to strain_ID
  coord_cartesian(clip = 'off')+
  theme(plot.margin = margin(1,3,1,1, "cm"))



##********************************************************
## Example 1 ----
##********************************************************

# create a simple tree and save it as an object
gg_simple <- ggtree(tree) %<+%
  metadata + # link our metadata file here
  coord_cartesian(clip = 'off')+
  theme(plot.margin = margin(1,4,1,1, "cm"))


# we can use tree manipulation to reorient the tree
# for example, flip will flip the position of two selected branches
gg_flip <- gg_simple %>% 
  flip(25, 17)

gg_flip

# add tip labels
gg_flip +
  geom_tiplab(offset = 0.0001)


# add tip points
gg_flip +
  geom_tiplab(offset = 0.0001)+
  geom_tippoint(aes(color = `tdh/trh`, shape = matrix), 
                size = 4, 
                alpha = 0.7)

# use scale_manual to specify shape and color of the points
gg_flip +
  geom_tiplab(offset = 0.0001)+
  geom_tippoint(aes(color = `tdh/trh`, shape = matrix), 
                size = 4, 
                alpha = 0.7)+
  scale_color_manual(values = c("+/+" = "deepskyblue1",
                                "+/-" = "mediumturquoise", 
                                "-/-" = "coral1"))+
  scale_shape_manual(values = c("oyster" = 16,
                                "stool" = 17,
                                "water" = 15),
                     name = "Sample matrix") 

# add a clade label
gg_flip +
  geom_tiplab(offset = 0.0001)+
  geom_tippoint(aes(color = `tdh/trh`, shape = matrix), 
                size = 4, 
                alpha = 0.7)+
  scale_color_manual(values = c("+/+" = "deepskyblue1",
                                "+/-" = "mediumturquoise", 
                                "-/-" = "coral1"))+
  scale_shape_manual(values = c("oyster" = 16,
                                "stool" = 17,
                                "water" = 15),
                     name = "Sample matrix") +
  geom_cladelab(node = 17, label = "ST36", 
                offset = 0.0008, 
                barsize = 1.5,
                barcolor = 'grey44', 
                textcolor = 'grey44', 
                offset.text = 0.0001)

example1 <- gg_flip +
  geom_tiplab(offset = 0.0001)+
  geom_tippoint(aes(color = `tdh/trh`, shape = matrix), 
                size = 4, 
                alpha = 0.7)+
  scale_color_manual(values = c("+/+" = "deepskyblue1",
                                "+/-" = "mediumturquoise", 
                                "-/-" = "coral1"))+
  scale_shape_manual(values = c("oyster" = 16,
                                "stool" = 17,
                                "water" = 15),
                     name = "Sample matrix") +
  geom_cladelab(node = 17, label = "ST36", 
                offset = 0.0008, 
                barsize = 1.5,
                barcolor = 'grey44', 
                textcolor = 'grey44', 
                offset.text = 0.0001)+
  theme(legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.spacing.y = unit(0, "mm"), 
        panel.border = element_blank(),
        aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))

example1

# export tree

ggsave("output/example1.jpeg", example1, dpi = 300)

# can also specify height and width


# check out colors: http://sape.inf.usi.ch/quick-reference/ggplot2/colour
# scale_color_manual: https://www.rdocumentation.org/packages/ggplot2/versions/0.9.1/topics/scale_colour_manual


##********************************************************
## Exercise 1  ----
## HINT: use as.factor() around a continuous variable 
## to read as a discrete scale
##********************************************************











##********************************************************
## Adding layers to tiplab ----
##********************************************************

exercise1 +
  geom_tiplab(aes(label = matrix),
              offset = 0.0008, 
              size = 5, 
              color = "grey44")

##********************************************************
## Highlighing clades ----
##********************************************************
gg_flip +
  geom_tiplab(offset = 0.0001)+
  geom_tippoint(color = 'black', 
                size = 4, 
                alpha = 0.5)+
  geom_hilight(node = c(11,12),
               fill = "pink",
               alpha = 0.4,
               extend = 0.0005)+
  geom_hilight(node = c(13, 14, 15, 17),
               fill = "lightblue",
               alpha = 0.4,
               extend = 0.0005)+
  geom_hilight(node = 1,
               fill = "purple",
               alpha = 0.2,
               extend = 0.0005)
# or

gg_flip +
  geom_tiplab(offset = 0.0001)+
  geom_tippoint(color = 'black', 
                size = 4, 
                alpha = 0.5)+
  geom_hilight(mapping=aes(subset = wg_cluster %in% 1),
               fill = "pink",
               alpha = 0.4,
               extend = 0.0005)+
  geom_hilight(mapping=aes(subset = wg_cluster %in% 2),
               fill = "lightblue",
               alpha = 0.4,
               extend = 0.0005)+
  geom_hilight(mapping=aes(subset = wg_cluster %in% 3),
               fill = "purple",
               alpha = 0.2,
               extend = 0.0005)



##********************************************************
## END
##******************************************************** 