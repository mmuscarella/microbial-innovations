################################################################################
# 
# Ancestral State Reconstruction (ASR) Method Comparisons
#
# Written by: Mario Muscarella
#
# Last Update: 20170213
#
# Goals:
#       1. Compare Custom ASR Function to Others
#
################################################################################

# Simulate Tree and Traits: Yule Tree and Markov Chain Trait Evolution
a <- 0.58
b <- 0.90
Hx <- ((a - 1)/(a + b - 2)) * log(a + b - 1)
Hy <- ((b - 1)/(a + b - 2)) * log(a + b - 1)

sim <- TT.sim(birth = 0.2, a = a, b = b, tips = 100)
phy <- sim$tree
x <- sim$traits

mc1 <- fitMC(phy = phy, x = x, ip = 0.5)
mc2 <- fitMC2(phy = phy, x = x, pp = c(0.5, 0.5), pi = c(1, 1) , post = TRUE)
ace1 <- ace(as.factor(x), phy, type = "discrete", model = "ARD", marginal = TRUE)
ace2 <- ace(as.factor(x), phy, type = "discrete", model = "ARD", marginal = FALSE)
rr1 <- rerootingMethod(phy, x, model = "ARD")
fMK <- fitMk(phy, x, model = "ARD", output.liks = TRUE)

a <- 0.8
b <- 0.6
sim <- TT.sim(birth = 0.2, a = a, b = b, tips = 100)
Hx <- ((a - 1)/(a + b - 2)) * log(a + b - 1)
Hy <- ((b - 1)/(a + b - 2)) * log(a + b - 1)
phy <- sim$tree
x <- sim$traits



mc1$par
mc2$fit$par
ace1$rates
rr1$Q
fMK$rates

dim(ace1$lik.anc)
dim(mc2$liks$liks)

par(mar = c(4, 2, 1, 1) + 0.5, oma = c(1, 4, 1, 1))
layout(matrix(1:4, 1, 4, byrow = T))
# ARD
plot(x = ace1$lik.anc[, 1], y = mc2$liks$liks[, 1],
     xlab = "ACE (marginal = T)", ylab = "", las = 1)
abline(0, 1,)
mtext(side = 2, "Custom (no prior)", line = 3)

plot(x = ace2$lik.anc[, 1], y = mc2$liks$liks[, 1],
     xlab = "ACE (marginal = F)", ylab = "", las = 1)
abline(0, 1,)

plot(x = rr1$marginal.anc[, 1], y = mc2$liks$liks[, 1],
     xlab = "RerootingMethod", ylab = "", las = 1)
abline(0, 1,)

plot(x = fMK$lik.anc[, 1], y = mc2$liks$liks[, 1],
     xlab = "fit MK", ylab = "", las = 1)
abline(0, 1,)


plot(x = rr1$marginal.anc[, 1], y = ace1$lik.anc[, 1],
     xlab = "RerootingMethod", ylab = "ACE (marginal = T)",
     main = "All Rates Different Model")

plot(x = rr1$marginal.anc[, 1], y = ace2$lik.anc[, 1],
     xlab = "RerootingMethod", ylab = "ACE (marginal = F)",
     main = "All Rates Different Model")

plot(x = fMK$lik.anc[, 1], y = ace2$lik.anc[, 1],
     xlab = "fitMK", ylab = "ACE (marginal = F)",
     main = "All Rates Different Model")

plot(x = fMK$lik.anc[, 1], y = mc2$liks$liks[, 1],
     xlab = "fitMK", ylab = "Custom (no prior)",
     main = "All Rates Different Model")

plot(x = fMK$lik.anc[, 1], y = rr1$marginal.anc[, 1],
     xlab = "fitMK", ylab = "RerootingMethod",
     main = "All Rates Different Model")



# SYM

ace3 <- ace(as.factor(x), phy, type = "discrete", model = "SYM", marginal = FALSE)
ace4 <- ace(as.factor(x), phy, type = "discrete", model = "SYM", marginal = TRUE)
rr2 <- rerootingMethod(phy, x, model = "SYM")
plot(x = rr2$marginal.anc[, 1], y = ace3$lik.anc[, 1],
     xlab = "RerootingMethod", ylab = "ACE (marginal = F)",
     main = "Symmetric Rates Model")

plot(x = rr2$marginal.anc[, 1], y = ace4$lik.anc[, 1],
     xlab = "RerootingMethod", ylab = "ACE (marginal = T)",
     main = "Symmetric Rates Model")


