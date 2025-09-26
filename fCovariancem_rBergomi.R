if (!requireNamespace("hypergeo", quietly = TRUE)) {
  install.packages("hypergeo")
}

# covariance matrix of the Volterra process driving the volatility

fG = function(x, H){
  gam = 1/2 - H
  value = 
    ((1 - 2 * gam) / (1 - gam)) * x ^ gam *
    hypergeo::hypergeo(1, gam, 2 - gam, x)
  return(Re(value))
}

fCov_WvWu = 
  function(v, u, H) {
    u ^ (2 * H) * fG(x = v/u, H = H)
  }

fMCov_W = function(v, u, iSteps, H){
  iN = round(abs(v - u) * iSteps + 1, 0)
  M = matrix(NA, nrow = iN, ncol = iN)
  for(i in 0 : (iN - 2)){
    u1 = u + i / iSteps
    M[i + 1, (i + 1) : iN ] = 
      sapply(X = seq(from = u1, to = v, by = 1/iSteps),
             FUN = fCov_WvWu, u = u1, H = H)
  } 
  M[iN, iN] = fCov_WvWu(v = v, u = v, H = H)
  M[lower.tri(M)] = t(M)[lower.tri(M)]
  return(M)
}

# covariance between the Volterra process driving the volatility
# and the Brownian motion driving the spot prices

fCov_WvZu = function(v, u, rho, H){ # v > u 
  Dh = sqrt(2 * H)/(H + 1/2)
  value = rho * Dh * (v ^ (H + 1/2) - (v - u) ^ (H + 1/2))
  return(value)
  
}

fCov_ZvWu = # v > u 
  function(rho, H, u) {
    Dh = sqrt(2 * H)/(H + 1/2)
    value = rho * Dh * u ^ (H + 1/2)
    return(value)
  }


fMCov_WZ = function(v, u, iSteps, H, rho) {
  iN = round(abs(v - u) * iSteps + 1, 0)
  M = matrix(NA, nrow = iN, ncol = iN)
  for(i in 0 : (iN - 2)){
    u1 = u + i / iSteps
    M[(i + 1) : iN, i + 1] = 
      sapply(X = seq(from = u1, to = v, by = 1/iSteps),
             FUN = fCov_WvZu, rho = rho , H = H, u = u1)
  }
  M[iN, iN] = fCov_WvZu(v, v, rho, H)
  M_ret = t(apply(M, 1L, zoo::na.locf))
  return(M_ret)
}

eigen(fMCov_WZ(1, 1/2, 2, 0.1, -0.9))

fAutocov_fBM = 
  
  function(s, t, H, sigma = 1){ 
    
    gamma = sigma ^ 2 / 2 * (abs(s) ^ (2 * H) + abs(t) ^ (2 * H) - abs(t - s) ^ (2 * H))
    
    return(gamma)
  }


fCov_BvBu = function(v, u, sigma){
  return(sigma ^ 2 * u)
}

fBM_CovMat = 
  function(v, u, sigma, iSteps, corr = FALSE){
    
    iN = round(abs(v - u) * iSteps + 1, 0)
    M = matrix(NA, nrow = iN, ncol = iN)
    for(i in 0 : (iN - 2)){
      u1 = u + i / iSteps
      M[i + 1, (i + 1) : iN ] = 
        sapply(X = seq(from = u1, to = v, by = 1/iSteps),
               FUN = fCov_BvBu, sigma = sigma, u = u1)
    } 
    M[iN, iN] = fCov_BvBu(u = v, v = v, sigma = sigma)
    M[lower.tri(M)] = t(M)[lower.tri(M)]
    return(M)
    
  }


# joint overall covariance matrix

fMCov_rBer = 
  function(v, u, iSteps, H, rho, sigma = 1){
    M1 = fMCov_W(v = v, u = u, iSteps = iSteps, H = H)
    M2 = fMCov_WZ(v = v, u = u, iSteps = iSteps, H = H, rho = rho)
    M3 = fBM_CovMat(v = v, u = u, iSteps = iSteps, sigma = sigma)
    M = 
      rbind(
        cbind(M1, M2),
        cbind(t(M2), M3)
      )
    return(M)
  }
