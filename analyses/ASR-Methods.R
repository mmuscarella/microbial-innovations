################################################################################
# 
# Ancestral State Reconstruction (ASR) Method Comparisons
#
# Written by: Mario Muscarella
#
# Last Update: 20170214
#
# Goals:
#       1. Compare Custom ASR Function to Others
#
################################################################################

# Load Source Code
source("../bin/TreeSimulationFxns.R")
source("../bin/fitMC.R")
library("png")
library("grid")

# Simulate Tree and Traits: Yule Tree and Markov Chain Trait Evolution
a <- 0.8
b <- 0.6
sim <- TT.sim(birth = 0.2, a = a, b = b, tips = 100)
Hx <- ((a - 1)/(a + b - 2)) * log(a + b - 1)
Hy <- ((b - 1)/(a + b - 2)) * log(a + b - 1)
phy <- sim$tree
x <- sim$traits

# Run Ancestral State Reconstruction Methods
fitmc <- fitMC2(phy = phy, x = x, pp = c(0.5, 0.5), pi = c(1, 1) , post = TRUE)
ace1 <- ace(as.factor(x), phy, type = "discrete", model = "ARD", marginal = TRUE)
ace2 <- ace(as.factor(x), phy, type = "discrete", model = "ARD", marginal = FALSE)
rr1 <- rerootingMethod(phy, x, model = "ARD")
fMK <- fitMk(phy, x, model = "ARD", output.liks = TRUE)

# Compare Parameter Estimates
fitmc$fit$par
ace1$rates
rr1$Q
fMK$rates

# Plot Comparisions of Custom Method with Others
png(filename="../figures/ASR_Methods.png",
    width = 1800, height = 600, res = 96*2)

par(mar = c(4, 2, 1, 1) + 0.5, oma = c(1, 4, 1, 1))
layout(matrix(1:4, 1, 4, byrow = T))


plot(x = ace1$lik.anc[, 1], y = fitmc$liks$liks[, 1],
     xlab = "ACE (marginal = T)", ylab = "", las = 1)
abline(0, 1)
mtext(side = 2, "Custom (no prior)", line = 3)

plot(x = ace2$lik.anc[, 1], y = fitmc$liks$liks[, 1],
     xlab = "ACE (marginal = F)", ylab = "", las = 1)
abline(0, 1)

plot(x = rr1$marginal.anc[, 1], y = fitmc$liks$liks[, 1],
     xlab = "RerootingMethod", ylab = "", las = 1)
abline(0, 1)

plot(x = fMK$lik.anc[, 1], y = fitmc$liks$liks[, 1],
     xlab = "fit MK", ylab = "", las = 1)
abline(0, 1)


# Close Plot Device
dev.off()
graphics.off()

## Show Plot
img <- readPNG("../figures/ASR_Methods.png")
grid.raster(img)

######################
# SYM Transition Rates
######################

ace3 <- ace(as.factor(x), phy, type = "discrete", model = "SYM", marginal = FALSE)
ace4 <- ace(as.factor(x), phy, type = "discrete", model = "SYM", marginal = TRUE)
rr2 <- rerootingMethod(phy, x, model = "SYM")



