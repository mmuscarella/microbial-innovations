################################################################################
# 
# Trait Evolution Simulations: 
#       Trait Conservation of Evolved Traits
#
# Written by: Mario Muscarella
#
# Last Update: 20170215
#
# Simulation Goals:
#       1. How well Does Trait Conservation Predict Trait Origins
#       2. 
#       3. 
#
################################################################################

# Load Dependencies
require("ape")||install.packages("ape");require("ape")
require("geiger")||install.packages("geiger");require("geiger")
require("expm")||install.packages("expm");require("expm")
require("MASS")||install.packages("MASS");require("MASS")

# Load Source Functions
source("../bin/TreeSimulationFxns.R")
source("../bin/ConsenTrait.R")


# Define Trait Evolution Parameters
test.a <- seq(0.51, 0.99, 0.02)
test.b <- seq(0.51, 0.99, 0.02)
test.comb <- expand.grid(test.a, test.b)

# Define Simulation Output
simmy <- data.frame(test.a = test.comb[,1], test.b = test.comb[,2],
                    first.obs.mean = rep(NA, dim(test.comb)[1]),
                    first.pred.mean = rep(NA, dim(test.comb)[1]),
                    nevol.obs.mean = rep(NA, dim(test.comb)[1]),
                    nevol.pred.mean = rep(NA, dim(test.comb)[1]))

# Run Simulaton Across Parameters
for (i in 1:dim(test.comb)[1]){
  print(paste("ConsenTrait Simulation:", i, "of", dim(test.comb)[1]))
  sim <- TraitEvolCon.Sim(birth = 0.2, a = test.comb[i, 1], b = test.comb[i, 2], 
                   nsim = 10, level = 0.90)
  if (dim(sim)[1] > 1){
    simmy$first.obs.mean[i] <- mean(sim[, 1])
    simmy$first.pred.mean[i] <- mean(sim[, 2])
    simmy$nevol.obs.mean[i] <- mean(sim[, 3])
    simmy$nevol.pred.mean[i] <- mean(sim[, 4])
  }
}

# Save Simulation Output
write.table(simmy, "../simulations/output/ConsenTraitSimulation.txt", 
            quote = F, sep = "\t", row.names = F, col.names = T)
