################################################################################
#
# Trait Evolution Analysis of GreenGenes:
#       Ancestral State Reconstruction (ASR) of Evolved Traits
#
# Written by: Mario Muscarella
#
# Last Update: 20170525
#
# Goals:
#       1. Run Ancestral State Reconstruction on Genes in the IMG Database
#       2. Determine the Origin of Genes
#       3. Determine the rate parameters
#
################################################################################

# Load Packages
library("ape")
library("phytools")
library("picante")
library("methods")

# Load Source Function
source("../bin/fitMC.R")

# Import Tree
IMG.um <- read.tree("../data/IMG.dated.tree")

# Import and Format Traits (Pathways)
IMG.traits <- read.delim("../data/IMG.KEGG.trimmed.txt")
dim(IMG.traits)

IMG.trait.tab <- t(IMG.traits[, -c(1, 2)])
colnames(IMG.trait.tab) <- IMG.traits[, 2]

sum(IMG.um$tip.label %in% colnames(IMG.trait.tab))
sum(colnames(IMG.trait.tab) %in% IMG.um$tip.label)

# Only use genes with more than 10 entries
IMG.trait.tab <- IMG.trait.tab[rowSums(IMG.trait.tab) > 10, ]
dim(IMG.trait.tab)

# Drop Genomes Not Found in Traits and Tree
IMG.trait.tab <- IMG.trait.tab[, which(colnames(IMG.trait.tab) %in% IMG.um$tip.label)]
dim(IMG.trait.tab)

missing <- setdiff(IMG.um$tip.label, colnames(IMG.trait.tab) )
IMG.tree <- drop.tip(IMG.um, missing)
length(IMG.tree$tip.label)

sub.trees <- subtrees(IMG.tree)

root.dists <- as.matrix(dist.nodes(IMG.tree))[,length(IMG.tree$tip.label) + 1]

# Reformat Trait Table
IMG.trait.tab <- data.frame(t(IMG.trait.tab))
IMG.trait.PA <- (IMG.trait.tab > 0) * 1
IMG.trait.tab[1:5, 1:5]

# Identify Genes That All Species Have
complete <- which(apply(IMG.trait.PA, 2, min) == 1)
print(paste("The following gene is found in all taxa: ",
            names(complete), sep = ""), quote = FALSE)
print("Those genes have been removed from the analysis", quote = FALSE)

write.table(complete, "../data/IMG_ASR_CompleteGenes.txt", quote = F,
            sep = "\t", row.names = T, col.names = F)

if (length(complete) > 0){
  IMG.trait.PA <- IMG.trait.PA[ , -c(as.numeric(complete))]
}
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

#############################################
# Run Ancestral State Reconstruction on Genes
#############################################
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
