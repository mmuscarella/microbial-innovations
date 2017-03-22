source("../bin/TreeSimulationFxns.R")
source("../bin/fitMC.R")

mc.test <- matrix(NA, ncol = 4, nrow = 100)

for (i in 1:dim(mc.test)[1]){
  sim <- TT.sim(a = 0.8, b = 0.8, tips = 1000, std = 100)
  tree <- sim$tree
  traits <- sim$traits
  
  # Run Ancestral State Reconstruction
  #ASR <- fitMC2(phy = tree, x = traits, init.parms = init.parms,
  #              prior = prior, posterior = TRUE)
  ASR <- fitMk(tree, traits, model = "ARD", 
               output.liks = TRUE)
  
  # Isolate Output
  rates <- ASR$rates
  mc.test[i, 1:2] <- c(rates[1], rates[2])
  mc.test[i, 3:4] <- unlist(xy_to_ab(rates[1], rates[2]))
}

plot(mc.test[, 3], mc.test[, 4],
     xlab = "Alpha", ylab = "Beta")
points(0.8, 0.8, pch = 16, col = "red")

xy <- ab_to_xy(0.8, 0.8)


plot(mc.test[, 1], mc.test[, 2],
     xlab = "Alpha", ylab = "Beta")
points(xy[1], xy[2], pch = 16, col = "red")

library(MASS)
k <- kde2d(mc.test[, 3], mc.test[, 4], n=200, 
           lims = c(c(0.65, 0.9), c(0.66, 0.9)))
image(k, xlab = "Alpha", ylab = "Beta", las = 1)
points(mc.test[, 3], mc.test[, 4], pch = 16, cex = 0.5)
points(0.8, 0.8, pch = 16, col = "blue", cex = 1.5)



k <- kde2d(mc.test[, 1], mc.test[, 2], n=200, 
           lims = c(c(0.1, 0.5), c(0.1, 0.5)))
image(k, xlab = "Alpha Rate", ylab = "Beta Rate", las = 1)
points(mc.test[, 1], mc.test[, 2], pch = 16, cex = 0.5)
points(-xy[[1]], -xy[[2]], pch = 16, col = "blue", cex = 1.5)

################################################################################
## Keep Tree The Same
mc.test <- matrix(NA, ncol = 4, nrow = 100)

tree <- Tree.sim(tips = 1000, std = 100)

for (i in 1:dim(mc.test)[1]){
  sim <- Trait.sim(tree = tree, x = -0.2, y = -0.4)
  tree <- sim$tree
  traits <- sim$traits
  
  # Run Ancestral State Reconstruction
  #ASR <- fitMC2(phy = tree, x = traits, init.parms = init.parms,
  #              prior = prior, posterior = TRUE)
  ASR <- fitMk(tree, traits, model = "ARD", 
               output.liks = TRUE, pi = c(0.001, 0.999))
  
  # Isolate Output
  rates <- -ASR$rates
  mc.test[i, 1:2] <- c(rates[1], rates[2])
  #mc.test[i, 3:4] <- unlist(xy_to_ab(rates[1], rates[2]))
}


plot(mc.test[, 1], mc.test[, 2],
     xlab = "X Rate", ylab = "Y Rate")
points(-0.2, -0.4, pch = 16, col = "red")

mean(mc.test[, 1])
median(mc.test[, 1])

mean(mc.test[, 2])
median(mc.test[, 2])

library(MASS)
k <- kde2d(mc.test[, 1], mc.test[, 2], n=200, 
           lims = c(c(-0.4, 0), c(-0.6, -0.2)))
image(k, xlab = "Alpha Rate", ylab = "Beta Rate", las = 1)
points(mc.test[, 1], mc.test[, 2], pch = 16, cex = 0.5)
points(-0.2, -0.4, pch = 16, col = "blue", cex = 1.5)


hist(mc.test[, 1])
hist(mc.test[, 2])
