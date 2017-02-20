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

# Load Source Functions
source("../bin/TreeSimulationFxns.R")
source("../bin/fitMC.R")

# Define Trait Evolution Parameters
test.a <- seq(0.51, 0.99, 0.04)
test.b <- seq(0.51, 0.99, 0.04)
test.comb <- expand.grid(test.a, test.b)

# Define Simulation Output
simmy <- data.frame(test.a = test.comb[,1], test.b = test.comb[,2],
                      logl.mean = rep(NA, dim(test.comb)[1]),
                      logl.sem = rep(NA, dim(test.comb)[1]),
                      x.mean = rep(NA, dim(test.comb)[1]),
                      y.mean = rep(NA, dim(test.comb)[1]))

# Run Simulaton Across Parameters
for (i in 1:dim(test.comb)[1]){
  print(i)
  sim <- TraitEvolASR.Sim(birth = 0.2, a = test.comb[i, 1], b = test.comb[i, 2],
                        nsim = 2)
  if (is.list(sim)){
  simmy$logl.mean[i] <- try(mean(sapply(sim[2,], '[')))
  simmy$logl.sem[i] <- try(sem(sapply(sim[2,], '[')))
  simmy$x.mean[i] <- try(mean(sapply(sim[1,], '[')[1,]))
  simmy$y.mean[i] <- try(mean(sapply(sim[1,], '[')[2,]))
  }
}

# Save Simulation Output
write.table(simmy, "../simulations/output/asr_simulation.txt", quote = F,
            sep = "\t", row.names = F, col.names = T)
