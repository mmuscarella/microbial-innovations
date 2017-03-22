################################################################################
#
# Trait Evolution Analysis of GreenGenes:
#       Ancestral State Reconstruction (ASR) of Evolved Traits
#
# Written by: Mario Muscarella
#
# Last Update: 20170222
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
library("methods")

# Load Source Functions
source("../bin/fitMC.R")

# Import Tree
GG.um <- read.tree("../data/GG.dated.tree")

# Import and Format Traits (Pathways)
IMG.traits <- read.delim("../data/predicted_patyways.L3.txt", skip = 1, 
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

root.dists <- as.matrix(dist.nodes(IMG.tree))[,length(IMG.tree$tip.label) + 1]

# Reformat Trait Table
IMG.trait.tab <- data.frame(t(IMG.traits))
IMG.trait.PA <- (IMG.trait.tab > 0) * 1

# Identify Traits That All Species Have
complete <- which(apply(IMG.trait.PA, 2, min) == 1)
print(paste("The following trait is found in all taxa: ", 
            names(complete), sep = ""), quote = FALSE)
print("Those triats have been removed from the analysis", quote = FALSE)

IMG.trait.PA <- IMG.trait.PA[ , -c(as.numeric(complete))]
dim(IMG.trait.PA)

out <- data.frame(matrix(NA, nrow = dim(IMG.trait.PA)[1], ncol = 8))
colnames(out) <- c("trait", "LogL", "Rate1", "Rate2", "Lik.Root", 
                   "First.Evol", "Med.Evol", "N.evol")

out.2 <- data.frame(Trait = character(), Anc  = character(), Des = character(),
                    Len = character(), S.Anc = character(), 
                    S.Des = character(),  R.Dist = character())

states <- c("Off", "On")

# Run Ancestral State Reconstruction on Pathways
print("Running ASR on Kegg Pathways", quote = FALSE)
for(i in 1:dim(IMG.trait.PA)[2]){
  print(paste("Analysing Trait ", i, " of ", dim(IMG.trait.PA)[2], ": ", 
              colnames(IMG.trait.PA)[i], sep = ""), quote = FALSE)
  #temp <- fitMC2(phy = IMG.tree, x = IMG.trait.PA[,i], prior = c(0.999, 0.001), 
  #                  posterior = TRUE)
  temp <- fitMk(IMG.tree, IMG.trait.PA[, i], model = "ARD", 
                pi = c(0.001, 0.999), output.liks = TRUE)
  temp.edge <- as.data.frame(IMG.tree$edge)
  colnames(temp.edge) <- c("Anc", "Des")
  temp.edge$Len <- IMG.tree$edge.length
  temp.edge$S.Anc <- NA
  temp.edge$S.Des <- NA
  temp.edge$R.Dist <- NA
  for(j in 1:dim(temp$lik.anc)[1]){
    temp.state <- sample(states, size = 1, prob = temp$lik.anc[j, ])
    temp.node <- rownames(temp$lik.anc)[j]
    temp.edge$S.Anc[which(temp.edge$Anc == temp.node)] <- temp.state
    temp.edge$S.Des[which(temp.edge$Des == temp.node)] <- temp.state
  }
  for(k in 1:dim(temp.edge)[1]){
    temp.edge$R.Dist[k] <- root.dists[temp.edge$Des[k]]
  }
  evol <- temp.edge[which(temp.edge$S.Anc == "Off" & temp.edge$S.Des == "On"), ]
  first <- try(min(evol$R.Dist))
  median <- try(median(evol$R.Dist))
  N.evol <- try(dim(evol)[1])
  
  out[i,1] <- colnames(IMG.trait.PA)[i]
  out[i,2] <- temp$logLik
  out[i,3] <- temp$rates[1]
  out[i,4] <- temp$rates[2]
  out[i,5] <- temp$lik.anc[1, 1]
  out[i,6] <- first
  out[i,7] <- median
  out[i,8] <- N.evol
  
  trait <- colnames(IMG.trait.PA)[i]
  evol.2 <- data.frame(Trait = rep(trait, dim(evol)[1]), evol)
  out.2 <- rbind(out.2, evol.2)
}

# Save Output
write.table(out, "../data/IMG_ASR.txt", quote = F,
            sep = "\t", row.names = F, col.names = T)

write.table(out.2, "../data/IMG_ASR_long.txt", quote = F,
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

root.dists <- as.matrix(dist.nodes(IMG.tree))[,length(IMG.tree$tip.label) + 1]

# Reformat Trait Table
IMG.trait.tab <- data.frame(t(IMG.traits))
IMG.trait.PA <- (IMG.trait.tab > 0) * 1

# Identify Traits That All Species Have
complete <- which(apply(IMG.trait.PA, 2, min) == 1)
print(paste("The following trait is found in all taxa: ", 
            names(complete), sep = ""), quote = FALSE)
print("Those triats have been removed from the analysis", quote = FALSE)

if(length(complete) > 0){
  IMG.trait.PA <- IMG.trait.PA[ , -c(as.numeric(complete))]
}
dim(IMG.trait.PA)

out <- data.frame(matrix(NA, nrow = dim(IMG.trait.PA)[1], ncol = 8))
colnames(out) <- c("trait", "LogL", "Rate1", "Rate2", "Lik.Root", 
                   "First.Evol", "Med.Evol", "N.evol")

out.2 <- data.frame(Trait = character(), Anc  = character(), Des = character(),
                    Len = character(), S.Anc = character(), 
                    S.Des = character(),  R.Dist = character())

states <- c("Off", "On")

# Run Ancestral State Reconstruction on Genes
print("Running ASR on Kegg Genes", quote = FALSE)
for(i in 1:dim(IMG.trait.PA)[2]){
  print(paste("Analysing Trait ", i, " of ", dim(IMG.trait.PA)[2], ": ", 
              colnames(IMG.trait.PA)[i], sep = ""), quote = FALSE)
  #temp <- fitMC2(phy = IMG.tree, x = IMG.trait.PA[,i], prior = c(0.999, 0.001), 
  #                  posterior = TRUE)
  temp <- fitMk(IMG.tree, IMG.trait.PA[, i], model = "ARD", 
                pi = c(0.999, 0.001), output.liks = TRUE)
  temp.edge <- as.data.frame(IMG.tree$edge)
  colnames(temp.edge) <- c("Anc", "Des")
  temp.edge$Len <- IMG.tree$edge.length
  temp.edge$S.Anc <- NA
  temp.edge$S.Des <- NA
  temp.edge$R.Dist <- NA
  for(j in 1:dim(temp$lik.anc)[1]){
    temp.state <- sample(states, size = 1, prob = temp$lik.anc[j, ])
    temp.node <- rownames(temp$lik.anc)[j]
    temp.edge$S.Anc[which(temp.edge$Anc == temp.node)] <- temp.state
    temp.edge$S.Des[which(temp.edge$Des == temp.node)] <- temp.state
  }
  for(k in 1:dim(temp.edge)[1]){
    temp.edge$R.Dist[k] <- root.dists[temp.edge$Des[k]]
  }
  evol <- temp.edge[which(temp.edge$S.Anc == "Off" & temp.edge$S.Des == "On"), ]
  first <- try(min(evol$R.Dist))
  median <- try(median(evol$R.Dist))
  N.evol <- try(dim(evol)[1])
  
  out[i,1] <- colnames(IMG.trait.PA)[i]
  out[i,2] <- temp$logLik
  out[i,3] <- temp$rates[1]
  out[i,4] <- temp$rates[2]
  out[i,5] <- temp$lik.anc[1, 1]
  out[i,6] <- first
  out[i,7] <- median
  out[i,8] <- N.evol
  
  trait <- colnames(IMG.trait.PA)[i]
  evol.2 <- data.frame(Trait = rep(trait, dim(evol)[1]), evol)
  out.2 <- rbind(out.2, evol.2)
}

# Save Output
write.table(out, "../data/IMG_genes_ASR.txt", quote = F,
            sep = "\t", row.names = F, col.names = T)

write.table(out.2, "../data/IMG_genes_ASR_long.txt", quote = F,
            sep = "\t", row.names = F, col.names = T)

