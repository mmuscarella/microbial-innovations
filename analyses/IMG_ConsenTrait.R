################################################################################
#
# Trait Evolution Analysis of GreenGenes:
#       Ancestral State Reconstruction (ASR) of Evolved Traits
#
# Written by: Mario Muscarella
#
# Last Update: 20170515
#
# Goals:
#       1.
#       2.
#       3.
#
################################################################################


# Load Packages
require("ape")
require("phytools")
require("picante")

# Load Source Functions
source("../bin/ConsenTrait.R")

# Import Tree
IMG.um <- read.tree("../data/IMG.dated.tree")

# Import and Format Traits (Pathways)
IMG.traits <- read.delim("../data/IMG.KEGG.trimmed.txt")
dim(IMG.traits)


IMG.trait.tab <- t(IMG.traits[, -c(1, 2)])
colnames(IMG.trait.tab) <- IMG.traits[, 2]
IMG.trait.tab[1:5, 1:5]

dim(IMG.trait.tab)
sum(IMG.um$tip.label %in% colnames(IMG.trait.tab))
sum(colnames(IMG.trait.tab) %in% IMG.um$tip.label)

IMG.trait.tab <- IMG.trait.tab[rowSums(IMG.trait.tab) > 10, ]
dim(IMG.trait.tab)
IMG.trait.tab[1:5, 1:5]

# Drop Genomes Not Found in Traits and Tree
IMG.trait.tab <- IMG.trait.tab[, which(colnames(IMG.trait.tab) %in% IMG.um$tip.label)]
dim(IMG.trait.tab)

missing <- setdiff(IMG.um$tip.label, colnames(IMG.trait.tab) )
if(length(missing) > 0){
  IMG.tree <- drop.tip(IMG.um, missing)
} else {
  IMG.tree <- IMG.um
}

IMG.trait.tab[1:5, 1:5]

# Reformat Trait Table
IMG.trait.tab2 <- data.frame(Genome = colnames(IMG.trait.tab), t(IMG.trait.tab))
IMG.trait.tab2[1:5, 1:5]

# Test
# IMG.trait.tab2 <- IMG.trait.tab2[, 1:5]

# Run ConsenTrait on Pathways
print("Running ConsenTrait on Kegg Pathways", quote = F)
trait.cons <- ConsenTrait(tree = IMG.tree, traits = IMG.trait.tab2)
print("Completed ConsenTrait on Kegg Pathways")
# Save Output
write.table(trait.cons, "../data/IMG_ConsenTrait.txt", quote = F,
            sep = "\t", row.names = F, col.names = T)
