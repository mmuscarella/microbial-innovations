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
IMG.traits <- IMG.traits[rowSums(IMG.traits) > 10, ]

# Drop tips not found in trait table
missing <- subset(GG.um$tip.label, 
                  !(GG.um$tip.label %in% colnames(IMG.traits)))
IMG.tree <- drop.tip(GG.um, missing)

# Reformat Trait Table
IMG.trait.tab <- data.frame(colnames(IMG.traits), t(IMG.traits))

# Run ConsenTrait on Pathways
print("Running ConsenTrait on Kegg Pathways")
trait.cons <- ConsenTrait(tree = IMG.tree, traits = IMG.trait.tab)
print("Completed ConsenTrait on Kegg Pathways")
# Save Output
write.table(trait.cons, "../data/IMG_ConsenTrait.txt", quote = F,
            sep = "\t", row.names = F, col.names = T)


# Import and Format Traits
IMG.traits <- read.delim("../data/metagenome_predicted.txt", skip = 1, 
                         row.names = 1)
dim(IMG.traits)
colnames(IMG.traits) <- gsub("X", "", colnames(IMG.traits))
sum(GG.um$tip.label %in% colnames(IMG.traits))
sum(colnames(IMG.traits) %in% GG.um$tip.label)
IMG.traits <- IMG.traits[rowSums(IMG.traits) > 10, ]
dim(IMG.traits)

# Drop tips not found in trait table
missing <- subset(GG.um$tip.label, 
                  !(GG.um$tip.label %in% colnames(IMG.traits)))
IMG.tree <- drop.tip(GG.um, missing)

# Reformat Trait Table
IMG.trait.tab <- data.frame(colnames(IMG.traits), t(IMG.traits))



# Run ConsenTrait on Genes
print("Running ConsenTrait on Kegg Genes")
genetrait.cons <- ConsenTrait(tree = IMG.tree, traits = IMG.trait.tab)
print("Completed ConsenTrait on Kegg Genes")


# Save Output
write.table(genetrait.cons, "../data/IMG_Genes_ConsenTrait.txt", quote = F,
            sep = "\t", row.names = F, col.names = T)
