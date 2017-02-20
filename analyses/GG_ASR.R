# Load Packages
require("ape")
require("phytools")
require("picante")
library("methods")

# Load Source Functions
source("../bin/fitMC.R")

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
IMG.trait.tab <- data.frame(t(IMG.traits))

IMG.trait.PA <- (IMG.trait.tab > 0) * 1

out <- data.frame(matrix(NA, nrow = dim(IMG.trait.PA)[1], 
                         ncol = 4 + IMG.tree$Nnode))
colnames(out) <- c("trait", "LogL", "Rate1", "Rate2", 
                   seq(1:IMG.tree$Nnode) + length(IMG.tree$tip.label))

# Run Ancestral State Reconstruction
for(i in 1:dim(IMG.trait.PA)[2]){
  #temp <- fitMC2(phy = IMG.tree, x = IMG.trait.PA[,i], prior = c(0.999, 0.001), 
  #                  posterior = TRUE)
  temp <- fitMk(IMG.tree, IMG.trait.PA[, 1], model = "ARD", 
                pi = c(0.999, 0.001), output.liks = TRUE)
  out[i,1] <- colnames(IMG.trait.PA)[i]
  out[i,2] <- temp$logLik
  out[i,3] <- temp$rates[1]
  out[i,4] <- temp$rates[2]
  out[i,5:dim(out)[2]] <- temp$lik.anc[, 1]
  return(out)
}

# Save Output
write.table(out, "../data/IMG_ASR.txt", quote = F,
            sep = "\t", row.names = F, col.names = T)
