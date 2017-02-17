################################################################################
# 
# Trait Evolution Simulations: 
#       First Occurence of Traits 
#
# Written by: Mario Muscarella
#
# Last Update: 20170214
#
# Simulation Goals:
#       1. How do the Markov Chain Parameters Change the Evolution of Traits
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

# Simulation 1: How does 1st Evolution Change with A
test.a <- seq(0.51, 0.99, 0.02)
test.b <- seq(0.51, 0.99, 0.02)

test1 <- data.frame(test.a = test.a, test.b = rep(0.90, length(test.a)),
                   sim.mean = rep(NA, length(test.a)),
                   sim.sem = rep(NA, length(test.a)))

for (i in 1:dim(test1)[1]){
  print(paste("Test 1:", i, "of", dim(test1)[1]), quote = F)
  sim <- TraitEvol.sim(birth = 0.2, a = test1$test.a[i], b = test1$test.b[i], 
                       nsim = 100)
  test1$sim.mean[i] <- mean(sim)
  test1$sim.sem[i] <- sem(sim)
}

write.table(test1, "../simulations/output/TraitEvolution_A.txt", 
            quote = F, col.names = T, row.names = F)

# Simulation 2: How does 1st Evolution Change with B
test2 <- data.frame(test.a = rep(0.90, length(test.b)), test.b = test.b,
                    sim.mean = rep(NA, length(test.b)),
                    sim.sem = rep(NA, length(test.b)))

for (i in 1:dim(test2)[1]){
  print(paste("Test 2:", i, "of", dim(test2)[1]))
  sim <- TraitEvol.sim(birth = 0.2, a = test2$test.a[i], b = test2$test.b[i], 
                       nsim = 100)
  test2$sim.mean[i] <- mean(sim)
  test2$sim.sem[i] <- sem(sim)
}

write.table(test2, "../simulations/output/TraitEvolution_B.txt", 
            sep = "\t", quote = F, col.names = T, row.names = F)

print("Simulation 1 & 2 Done", quote = F)

################################################################################
# Simulation 3: 1st Evolution across ranges of A and B
test.comb <- expand.grid(test.a, test.b)[-1,]

test3 <- data.frame(test.a = test.comb[,1], test.b = test.comb[,2],
                    sim.mean = rep(NA, dim(test.comb)[1]),
                    sim.sem = rep(NA, dim(test.comb)[1]))

for (i in 1:dim(test3)[1]){
  print(paste("Test 3:", i, "of", dim(test3)[1]), quote = F)
  sim <- TraitEvol.sim(birth = 0.2, a = test3$test.a[i], b = test3$test.b[i], 
                       nsim = 100)
  test3$sim.mean[i] <- mean(sim)
  test3$sim.sem[i] <- sem(sim)
}

write.table(test3, "../simulations/output/TraitEvolution_AB.txt", 
            sep = "\t", quote = F, col.names = T, row.names = F)

print("Test 3 Done", quote = F)