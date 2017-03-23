################################################################################
#
# Trait Evolution Simulations:
#       Ancestral State Reconstruction (ASR) of Evolved Traits
#
# Written by: Mario Muscarella
#
# Last Update: 20170209
#
# Simulation Goals:
#       1. How does ASR likelihood change with trait evolution parameters
#       2. How well does ASR predict rates across evolution parameters
#       3.
#
################################################################################

# Load Dependencies
library("ape")
library("geiger")
library("expm")
library("MASS")
library("methods")

# Load Source Functions
source("../bin/TreeSimulationFxns.R")
source("../bin/fitMC.R")

# Define Trait Evolution Parameters
test.a <- seq(0.81, 0.99, 0.04)
test.b <- seq(0.81, 0.99, 0.04)
test.comb <- expand.grid(test.a, test.b)

# Define Simulation Output
simmy <- data.frame(matrix(NA, nrow = dim(test.comb)[1],
                           ncol = 6 + 99))

colnames(simmy) <- c("Alpha", "Beta", "LogL", "LogL_sem", "Mean_X", "Mean_Y",
                     seq(1:99) + 100)

simmy$Alpha <- test.comb[,1]
simmy$Beta <- test.comb[,2]

# Run Simulaton Across Parameters
for (i in 1:dim(test.comb)[1]){
  print(paste("Running Simulation", i, "of", dim(test.comb)[1]))
  sim <- TraitEvolASR.Sim(birth = 0.2, a = test.comb[i, 1], b = test.comb[i, 2],
                        nsim = 10)
  sim <- as.data.frame(sim)
  simmy$LogL[i] <- mean(sim$LogL)
  simmy$LogL_sem[i] <- sem(sim$LogL)
  simmy$Mean_X[i] <- mean(sim$X)
  simmy$Mean_Y[i] <- mean(sim$Y)
  simmy[i, 7:dim(simmy)[2]] <- apply(sim[, -c(1:3)], 2, mean)
}

# Save Simulation Output
write.table(simmy, "../simulations/output/asr_simulation.txt", quote = F,
            sep = "\t", row.names = F, col.names = T)
