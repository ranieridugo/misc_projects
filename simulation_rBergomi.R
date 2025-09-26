# Simulation of Volterra process and BM

source("fCovariancem_rBergomi.R")
iM = 1e5
v = 1; iSteps = 30; u = 1/iSteps; H = 0.07; rho = - 0.9; 
ksi0 = 0.235 ^ 2; eta = 1.9; S0 = 100

iN = round(abs(v - u) * iSteps + 1, 0)
mIn = fMCov_rBer(v = v, u = u, iSteps = iSteps, H = H, rho = rho)
mLT = t(chol(mIn))
mGN = matrix(rnorm(n = 2 * iN * iM), nrow = 2 * iN, ncol = iM)
mSim = mLT %*% mGN
mVolterra = rbind(
  rep(0, iM),
  mSim[1 : iN, ]
)
mMB = rbind(
  rep(0, iM),
  mSim[(iN + 1) : (2 * iN), ]
)
vT = seq(from = 0, to = v, by = 1/iSteps)

# Verification
k = sample(1:iM, 1)
plot(mMB[, k], type = "l")

# Verify distribution mB
plot(apply(mMB, 1, mean), type = "l")
abline(h = 0)
plot(apply(mMB, 1, var), type = "l")
lines(vT, col = "red")

# Verify variance Volterra process
plot(x = vT,
     y = mVolterra[, k], type = "l")

plot(apply(mVolterra, 1, var))
lines(vT ^ (2 * H), col = "red")

# Verifying corellation
k = sample(1:iN, 1)
cor(mVolterra[k,], mMB[k, ])

# Simulation of the rough Bergomi

# variance
mVar = matrix(NA, nrow = nrow(mVolterra),
              ncol = ncol(mVolterra))
for(i in 1 : ncol(mVar)){
  mVar[, i] = ksi0 * exp(eta * mVolterra[, i] - 0.5 * eta ^ 2 * vT ^ (2 * H))
}

par(mfrow = c(2, 2), mar = c(0, 0, 0, 0))
for(i in runif(4, 1, iM)){
  plot(x = vT, y = mVar[, i], type = "l")
  abline(h = 0, col = "red")
}
par(mfrow = c(1, 1), mar = c(2, 2, 2, 2))

plot(apply(mVar, 1, mean), type = "l")
abline(h = ksi0, col = "red")

# price
mBM_diff = diff(mMB)
dt = 1/iSteps
mIncrements = mVar[- nrow(mVar), ] ^ 0.5 * mBM_diff - 0.5 * mVar[- nrow(mVar), ] * dt
mIntegral = apply(mIncrements, 2, cumsum)
mPrice = S0 * rbind(rep(1, iM), exp(mIntegral))

par(mfrow = c(2, 2))
for(i in runif(4, 1, iM)){
  plot(x = vT, y = mPrice[, i], type = "l")
  abline(h = 0, col = "red")
}
par(mfrow = c(1, 1))

plot(apply(mPrice, 1, mean))
abline(h = S0, col = "red")

plot(x = vT, y = - 2 * apply(log(mPrice), 1, mean))
lines(x = vT, y = ksi0 * vT, col = "red")

# Plot resulting price density at final time T, their density, and writing 
# the simulated spot trajectories to csv

hist(mPrice[iN + 1, ], breaks = 200, freq = FALSE)
curve(dnorm(x, mean = 100, sd=sd(mPrice[iN + 1, ])),
      col="red4", lwd=2, add=TRUE, yaxt="n")

matplot(mPrice[, 1 : 10] * S0, type = "l")

write.csv(mPrice, "sim_spot_m1e5_rbergomi.csv", row.names = FALSE)

