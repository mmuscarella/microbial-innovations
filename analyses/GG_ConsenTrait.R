# Load Packages
require("ape")
require("phytools")
require("picante")

# Load Source Functions
source("../bin/ConsenTrait.R")

# Import Tree
GG.um <- read.tree("../data/GG.dated.tree")


# Import and Format Traits
IMG.traits <- read.delim("../data/predicted_patyways.L3.txt", skip = 1, 
                         row.names = 1)
dim(IMG.traits)
colnames(IMG.traits) <- gsub("X", "", colnames(IMG.traits))
sum(GG.um$tip.label %in% colnames(IMG.traits))
sum(colnames(IMG.traits) %in% GG.um$tip.label)
rowSums(IMG.traits > 1)

# Drop tips not found in trait table
missing <- subset(GG.um$tip.label, 
                  !(GG.um$tip.label %in% colnames(IMG.traits)))
IMG.tree <- drop.tip(GG.um, missing)

# Reformat Trait Table
IMG.trait.tab <- data.frame(colnames(IMG.traits), t(IMG.traits))

# Run ConsenTrait
trait.cons <- ConsenTrait(tree = IMG.tree, traits = IMG.trait.tab)

# Save Output
write.table(trait.cons, "../data/IMG_ConsenTrait.txt", quote = F,
            sep = "\t", row.names = F, col.names = T)
