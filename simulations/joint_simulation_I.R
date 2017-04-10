#rm(list=ls())
#setwd("~/GitHub/microbial-innovations/simulations")
source("../bin/TreeSimulationFxns.R")
source("../bin/fitMC.R")
source("../bin/ConsenTrait.R")
library("methods")

# Add Basic Functions
se <- function(x, ...){sd(x, na.rm = TRUE)/sqrt(length(na.omit(x)))}
ci <- function(x, ...){1.96 * sd(x,na.rm = TRUE)}


# Simulation 1: Accuracy of ASR Predictions
test.x <- c(-0.2, -0.02, -0.002, -0.0002, -0.00002, -0.000002)
test.y <- c(-0.2, -0.02, -0.002, -0.0002, -0.00002, -0.000002)
test.comb <- data.frame(expand.grid(X = test.x, Y = test.y))
rep.comb <- data.frame(test.comb[rep(seq_len(nrow(test.comb)), each=10),])
rownames(rep.comb) <- NULL

mc.testA <- matrix(NA, ncol = 4, nrow = dim(rep.comb)[1])

tree <- Tree.sim(tips = 2000, std = 4000)

for (i in 1:dim(mc.testA)[1]){
  print(paste("Simulation", i, "out of", dim(mc.testA)[1], sep = " "))
  sim <- Trait.sim(tree = tree, x = rep.comb[i,1], y = rep.comb[i,2])
  tree <- sim$tree
  traits <- sim$traits
  
  if(sum(traits == "ON") > 0){
  # Run Ancestral State Reconstruction
  ASR <- fitMk(tree, traits, model = "ARD", 
               output.liks = TRUE, pi = c(0.001, 0.999))
  
  # Isolate Output
  rates <- -ASR$rates
  mc.testA[i, 1:2] <- c(rates[1], rates[2])
  #mc.testA[i, 3:4] <- unlist(xy_to_ab(rates[1], rates[2]))
  } else {
    mc.testA[i, 1:2] <- c(NA, NA)
  }
}

# Plot Output
#png(filename="../simulations/output/Estimates.png",
#    width = 1800, height = 900, res = 96*2)

#layout(matrix(1:2, ncol = 2, byrow = T))
#par(mar= c(3,3,1,1))
#plot(y = log10(-mc.testA[, 1]), x = log10(-rep.comb[, 1]),
#     xlab = "True X Rate (log10)", ylab = "Predicted X Rate (log10)")

#abline(0, 1)

#plot(y = log10(-mc.testA[, 2]), x = log10(-rep.comb[, 2]),
#     xlab = "True Y Rate (log10)", ylab = "Predicted Y Rate (log10)")

#abline(0, 1)

#dev.off()
#graphics.off()

#img <- readPNG("../simulations/output/Estimates.png")
#grid.raster(img)



# Simulation 2: How does 1st Evolution Change with X
test.x <- c(-0.2, -0.02, -0.002, -0.0002, -0.00002, -0.000002)
test.y <- c(-0.2, -0.02, -0.002, -0.0002, -0.00002, -0.000002)
test.comb <- data.frame(expand.grid(X = test.x, Y = test.y))
rep.comb <- data.frame(test.comb[rep(seq_len(nrow(test.comb)), each=10),])
rownames(rep.comb) <- NULL

mc.testB <- data.frame(matrix(NA, ncol = 10, nrow = dim(rep.comb)[1]))
colnames(mc.testB) <- c("X", "Y", "X.p", "Y.p", 
                       "ASR.Origins", "ASR.Mean", "ASR.SE",
                       "Con.Origins", "Con.Mean", "Con.SE")

tree <- Tree.sim(tips = 2000, std = 4000)

for (i in 1:dim(mc.testB)[1]){
  print(paste("Simulation", i, "out of", dim(mc.testB)[1], sep = " "))
  sim <- Trait.sim(tree = tree, x = rep.comb[i,1], y = rep.comb[i,2])
  tree <- sim$tree
  attributes(tree)$seed <- NULL
  traits <- sim$traits
  
  if(sum(traits == "ON") > 0){
    # Run Ancestral State Reconstruction
    ASR <- ASRTrait(tree, traits)
    ASR.r <- ASR$ASR
    ASR.o <- ASR$Origins
    
    # Run ConsenTrait
    trait.tab <- data.frame(names(traits), (traits == "On"))
    CON <- ConsenTrait(tree = tree, traits = trait.tab)
    
    # Isolate Output
    rates <- -ASR.r$rates
    mc.testB[i, 1:2] <- c(rep.comb[i,1], rep.comb[i,2])
    mc.testB[i, 3:4] <- c(rates[1], rates[2])
    mc.testB[i, 5] <- try(dim(ASR.o)[1])
    mc.testB[i, 6] <- try(mean(ASR.o$distance))
    mc.testB[i, 7] <- try(se(ASR.o$distance))
    mc.testB[i, 8] <- try(dim(CON)[1])
    mc.testB[i, 9] <- try(mean(CON$distance))
    mc.testB[i, 10] <- try(se(CON$distance))
  }
}

save.image(file = "joint_simulation_I.RData")



