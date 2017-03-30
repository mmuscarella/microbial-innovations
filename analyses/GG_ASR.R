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
sub.trees <- subtrees(IMG.tree)

root.dists <- as.matrix(dist.nodes(IMG.tree))[,length(IMG.tree$tip.label) + 1]

# Reformat Trait Table
IMG.trait.tab <- data.frame(t(IMG.traits))
IMG.trait.PA <- (IMG.trait.tab > 0) * 1

# Identify Traits That All Species Have
complete <- which(apply(IMG.trait.PA, 2, min) == 1)
print(paste("The following trait is found in all taxa: ", 
            names(complete), sep = ""), quote = FALSE)
print("Those triats have been removed from the analysis", quote = FALSE)

write.table(complete, "../data/IMG_ASR_CompletePathways.txt", quote = F,
            sep = "\t", row.names = T, col.names = F)

IMG.trait.PA <- IMG.trait.PA[ , -c(as.numeric(complete))]
dim(IMG.trait.PA)

out <- data.frame(matrix(NA, nrow = dim(IMG.trait.PA)[1], ncol = 11))
colnames(out) <- c("trait", "Rate1", "Rate2", "LogL", "Lik.Root", 
                   "First.Evol", "Med.Evol", "N.evol", "N.origins",
                   "Origins", "Dates.Origins")

cat(colnames(out), file = "../data/IMG_ASR_PathwaysTemp.csv", sep = "; ",
    append = F, fill = F)
cat(file = "../data/IMG_ASR_PathwaysTemp.csv", append = T, fill = T)

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
                pi = c(0.999, 0.001), output.liks = TRUE)
  temp.edge <- as.data.frame(IMG.tree$edge)
  colnames(temp.edge) <- c("Anc", "Des")
  temp.edge$Len <- IMG.tree$edge.length
  temp.edge$S.Anc <- NA
  temp.edge$S.Des <- NA
  temp.edge$R.Dist <- NA
  for(j in 1:dim(temp$lik.anc)[1]){
    temp.state <- sample(states, size = 1, prob = abs(temp$lik.anc[j, ]))
    temp.node <- rownames(temp$lik.anc)[j]
    temp.edge$S.Anc[which(temp.edge$Anc == temp.node)] <- temp.state
    temp.edge$S.Des[which(temp.edge$Des == temp.node)] <- temp.state
  }
  for(k in 1:dim(temp.edge)[1]){
    temp.edge$R.Dist[k] <- root.dists[temp.edge$Des[k]]
  }
  evol <- temp.edge[which(temp.edge$S.Anc == "Off" & temp.edge$S.Des == "On"), ]
  
  evols <- evol$Des    
  origins <- vector(mode = "character", length = 0)
  positives <- vector(mode = "character", length = 0)
  
  # Loop through all subtrees and determine origins
  for (j in 1:length(sub.trees)){
    tip_names <- sub.trees[[j]]$tip.label
    tree_name <- sub.trees[[j]]$name
    if (tree_name %in% evols){
      match_test <- match(tip_names, positives)
      if (all(is.na(match_test))){
        positives <- c(positives,tip_names)
        origins <- c(origins, tree_name)
      } else {
        if (any(is.na(match_test))) {
          print("some NAs - something is weird")
        }
      }
    }
  }
  
  first <- try(min(evol$R.Dist))
  median <- try(median(evol$R.Dist))
  N.evol <- try(dim(evol)[1])
  N.origins <- try(length(origins))
  name.origins <- try(toString(origins, collapse = ", "))
  dates.origins <- try(toString(evol$R.Dist[which(evol$Des %in% origins)], collapse = ","))
  
  out[i,1] <- colnames(IMG.trait.PA)[i]
  out[i,2] <- -temp$rates[1]
  out[i,3] <- -temp$rates[2]
  out[i,4] <- temp$logLik
  out[i,5] <- temp$lik.anc[1, 1]
  out[i,6] <- first
  out[i,7] <- median
  out[i,8] <- N.evol
  out[i,9] <- N.origins
  out[i,10] <- name.origins
  out[i,11] <- dates.origins
  
  cat(unlist(out[i, ]), file = "../data/IMG_ASR_PathwaysTemp.csv", sep = "; ",
      append = T)
  cat(file = "../data/IMG_ASR_PathwaysTemp.csv", append = T, fill = T)
  
  
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

write.table(complete, "../data/IMG_ASR_CompleteGenes.txt", quote = F,
            sep = "\t", row.names = T, col.names = F)

if(length(complete) > 0){
  IMG.trait.PA <- IMG.trait.PA[ , -c(as.numeric(complete))]
}
dim(IMG.trait.PA)

out <- data.frame(matrix(NA, nrow = dim(IMG.trait.PA)[1], ncol = 11))
colnames(out) <- c("trait", "Rate1", "Rate2", "LogL", "Lik.Root", 
                   "First.Evol", "Med.Evol", "N.evol", "N.origins",
                   "Origins", "Dates.Origins")

cat(colnames(out), file = "../data/IMG_ASR_GenesTemp.csv", sep = "; ",
    append = F, fill = F)
cat(file = "../data/IMG_ASR_GenesTemp.csv", append = T, fill = T)

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
    temp.state <- sample(states, size = 1, prob = abs(temp$lik.anc[j, ]))
    temp.node <- rownames(temp$lik.anc)[j]
    temp.edge$S.Anc[which(temp.edge$Anc == temp.node)] <- temp.state
    temp.edge$S.Des[which(temp.edge$Des == temp.node)] <- temp.state
  }
  for(k in 1:dim(temp.edge)[1]){
    temp.edge$R.Dist[k] <- root.dists[temp.edge$Des[k]]
  }
  evol <- temp.edge[which(temp.edge$S.Anc == "Off" & temp.edge$S.Des == "On"), ]
  
  evols <- evol$Des    
  origins <- vector(mode = "character", length = 0)
  positives <- vector(mode = "character", length = 0)
  
  # Loop through all subtrees and determine origins
  for (j in 1:length(sub.trees)){
    tip_names <- sub.trees[[j]]$tip.label
    tree_name <- sub.trees[[j]]$name
    if (tree_name %in% evols){
      match_test <- match(tip_names, positives)
      if (all(is.na(match_test))){
        positives <- c(positives,tip_names)
        origins <- c(origins, tree_name)
      } else {
        if (any(is.na(match_test))) {
          print("some NAs - something is weird")
        }
      }
    }
  }
  
  first <- try(min(evol$R.Dist))
  median <- try(median(evol$R.Dist))
  N.evol <- try(dim(evol)[1])
  N.origins <- try(length(origins))
  name.origins <- try(toString(origins, collapse = ", "))
  dates.origins <- try(toString(evol$R.Dist[which(evol$Des %in% origins)], collapse = ","))
  
  out[i,1] <- colnames(IMG.trait.PA)[i]
  out[i,2] <- -temp$rates[1]
  out[i,3] <- -temp$rates[2]
  out[i,4] <- temp$logLik
  out[i,5] <- temp$lik.anc[1, 1]
  out[i,6] <- first
  out[i,7] <- median
  out[i,8] <- N.evol
  out[i,9] <- N.origins
  out[i,10] <- name.origins
  out[i,11] <- dates.origins
  
  cat(unlist(out[i, ]), file = "../data/IMG_ASR_GenesTemp.csv", sep = "; ",
      append = T)
  cat(file = "../data/IMG_ASR_GenesTemp.csv", append = T, fill = T)
  
  trait <- colnames(IMG.trait.PA)[i]
  evol.2 <- data.frame(Trait = rep(trait, dim(evol)[1]), evol)
  out.2 <- rbind(out.2, evol.2)
}

# Save Output
write.table(out, "../data/IMG_genes_ASR.txt", quote = F,
            sep = "\t", row.names = F, col.names = T)

write.table(out.2, "../data/IMG_genes_ASR_long.txt", quote = F,
            sep = "\t", row.names = F, col.names = T)

