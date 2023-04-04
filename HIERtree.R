library(rhierbaps)
library(ggtree)
library(ggplot2)
library(phytools)
library(ape)

set.seed(1234)

setwd("C:/Users/lauri/Desktop/FINAL PROJ")

fasta.file.name <- "small_single.fa"
#fasta.file.name <- "big_single.fa"
snp.matrix <- load_fasta(fasta.file.name)
hb.results <- hierBAPS(snp.matrix, max.depth = 2, n.pops = 20, quiet = FALSE)
head(hb.results$partition.df)

  # uncomment the required line to compute correlating tree
  # small - contains 177 sequences; big - contains 213 sequences

  #FASTTREE GTR+CAT
#newick.file.name <- "small_single.tree"
#newick.file.name <- "big_single.tree"


  # RAXML/GTRGAMMA 1 - default settings, 2 - bootstrap 10, 3 - bootstrap 100
#newick.file.name <- "Best_scoring_small_1.tree"
#newick.file.name <- "Best_scoring_small_2.tree"
#newick.file.name <- "Best_scoring_small_3.tree"
#newick.file.name <- "Parsimony_small_1.tree"
#newick.file.name <- "Parsimony_small_2.tree"
#newick.file.name <- "Parsimony_small_3.tree"

  # RAXML/GTRCAT(i)
#newick.file.name <- "small_cat.tree"
#newick.file.name <- "small_cati.tree"
  #RAXML/GTRGAMMA(i)
#newick.file.name <- "small_gamma.tree"
#newick.file.name <- "small_gammai.tree"

  # bootstrapping criteria applied to "small_single.fa"
  # autoMR is read with 2nd iqtree line
#newick.file.name <- "final_bootstrap_autoMR.tree"
#newick.file.name <- "final_bootstrap_autoMRE.tree"
#newick.file.name <- "final_bootstrap_autoMREIGN.tree"
#newick.file.name <- "final_bootstrap_autoFC.tree"

#newick.file.name <- "big_final_bootstrap_autoMRE.tree"
#newick.file.name <- "small_gammai_mre.tree"
#newick.file.name <- "small_gammai_mre_B.tree"
#newick.file.name <- "big_gammai_mre.tree"

  # FASTTREE trees
#newick.file.name <- "small_jukes-cantor.tree"
#newick.file.name <- "small_NJ_jukes-cantor.tree"
#newick.file.name <- "small_NJ_positive.tree"

#newick.file.name <- "big_NJ_noNegative.tree"
#newick.file.name <- "big_NJ.tree"

iqtree <- ape::read.tree(newick.file.name)
#iqtree <- phytools::read.newick(newick.file.name)


metadata <- read.csv("metadata1.csv")
  #If metadata file does not include the first 'number' column, uncomment line 45
#metadata <- tibble::rowid_to_column(metadata, "number")

selection <- metadata[39]
 
 
gg <- ggtree(iqtree, layout="circular")
 
gg <- gg %<+% hb.results$partition.df
 
gg <- gg + geom_tippoint(aes(color = factor(`level 1`)), size=0.9) + 
    scale_color_manual(name="Population", values=c("#FF00BF", "#EADA0E", "#1AF892", "#1800FF", "#F8921A", "#7AE8FA", "#FF4848"))
 
gg <- gg %<+% metadata

gg <- gheatmap(gg, selection, offset = 0, width = 0.05, colnames=FALSE) +
   scale_fill_manual(name="Strain/Selection", values=c("#196F3D", "#A9DFBF", "#5B2C6F", "#D2B4DE"), labels = c("WT cipro -", "WT cipro +", "MutS cipro -", "MutS cipro +"))
gg

dpi = 300
# # ggsave(filename="test_big_single_org.png", width=9000/dpi, height=6000/dpi, dpi = dpi, scale = 0.4, limitsize = FALSE, units = "px")
# #ggsave(filename = "small_singe_ML3.png", plot = gg, scale = 4, limitsize = FALSE)
#ggsave(filename="small_gammai_mre.png", width=7000, height=5000, dpi = dpi, scale = 0.55, limitsize = FALSE, units = "px")
