################################################################################
#
# Trait Evolution Simulations:
#       Ancestral State Reconstruction (ASR) and Root Posteriors
#
# Written by: Mario Muscarella
#
# Last Update: 20170213
#
# Simulation Goals:
#       1. How do the priors change the posterior of the root
#       2. 
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

# How does the likelihood of the root state change with A and B (no prior)
test.a <- seq(0.50, 0.95, 0.01)
test.b <- seq(0.50, 0.95, 0.01)
test.comb <- expand.grid(test.a, test.b)[-1,]

sim.lik <- data.frame(test.a = test.comb[,1], test.b = test.comb[,2],
                      sim.mean = rep(NA, dim(test.comb)[1]),
                      sim.sem = rep(NA, dim(test.comb)[1]))

for (i in 1:dim(test.comb)[1]){
  sim <- TraitEvol.sim2(birth = 0.2, a = test.comb[i, 1], b = test.comb[i, 2], 
                        nsim = 100)
  sim.lik$sim.mean[i] <- mean(sim[1, ])
  sim.lik$sim.sem[i] <- sem(sim[1, ])
}

write.table(sim.lik, "../simulations/output/TraitEvolutionPrior_Flat.txt", 
            quote = F, col.names = T, row.names = F)


# How does the likelihood of the root change with A and B (0.99, 0.01 prior)
test.a <- seq(0.54, 0.95, 0.01)
test.b <- seq(0.54, 0.95, 0.01)
test.comb <- expand.grid(test.a, test.b)[-1,]

sim.lik <- data.frame(test.a = test.comb[,1], test.b = test.comb[,2],
                      sim.mean = rep(NA, dim(test.comb)[1]),
                      sim.sem = rep(NA, dim(test.comb)[1]))

for (i in 1:dim(test.comb)[1]){
  sim <- TraitEvol.sim3(birth = 0.2, a = test.comb[i, 1], b = test.comb[i, 2], 
                        nsim = 100)
  sim.lik$sim.mean[i] <- mean(sim[1, ])
  sim.lik$sim.sem[i] <- sem(sim[1, ])
}

write.table(sim.lik, "../simulations/output/TraitEvolutionPrior_Off.txt", 
            quote = F, col.names = T, row.names = F)

