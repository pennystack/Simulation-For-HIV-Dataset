# LIKELIHOOD ESTIMATION FUNCTIONS ----------------------------------------------

make_freq <- function(df) {
  
  freqij <- matrix(0, nrow = nstates, ncol = nstates)
  
  for (i in 1:nstates) {
    for(j in 1:nstates) {
      if(i != j){
        freqij[i, j] <- df[state == i & state_next == j, .N]
      }
      else{
        freqij[i, j] <- 0
      }
    }
  }
  
  return (freqij / rowSums(freqij))
}



loglikelihood_fast <- function(params, obstimes) {
  
  vij <- abs(matrix(params[(1):(nstates**2)], nstates, nstates)) / parscale
  sij <- abs(matrix(params[(nstates**2 + 1):(2 * nstates**2)], nstates, nstates)) / parscale
  aij <- matrix(params[(2 * nstates**2 + 1):(3 * nstates ** 2)], nstates, nstates)
  diag(aij) = 0
  aij[nstates, nstates-1] <- - sum(aij[nstates, 1:(nstates-2)])
  diag(aij) = 0
  aij[, nstates] = - rowSums(aij[, 1:(nstates-1)])
  diag(aij) = 0
  bij <- matrix(params[(3 * nstates**2 + 1):(4 * nstates ** 2)], nstates, nstates)
  diag(bij) = 0
  bij[nstates, nstates-1] <- 1 - sum(bij[nstates, 1:(nstates-2)])
  bij[, nstates] = 1 - rowSums(bij[, 1:(nstates-1)])
  
  
  return(cpp_loglikelihood(obstimes, aij = c(aij), bij = c(bij), vij = c(vij), sij = c(sij), nstates = nstates))
}

my_gradient <- function(y, obstimes) numDeriv::grad(func = function(z) {loglikelihood_fast(z, obstimes)}, x = y)


# -----------------------------------------------------------------------------------------------
est_vij <- function(params, nstates, parscale){
  vij <- abs(matrix(params[(1):(nstates**2)], nstates, nstates)) / parscale
  diag(vij) = 0
  return(vij)
}

est_sij <- function(params, nstates, parscale){
  sij <- abs(matrix(params[(nstates**2 + 1):(2 * nstates**2)], nstates, nstates)) / parscale
  diag(sij) = 0
  return(sij)
}

est_aij <- function(params, nstates, parscale){
  aij <- matrix(params[(2 * nstates**2 + 1):(3 * nstates ** 2)], nstates, nstates)
  diag(aij) = 0
  aij[nstates, nstates-1] <- - sum(aij[nstates, 1:(nstates-2)])
  diag(aij) = 0
  aij[, nstates] = - rowSums(aij[, 1:(nstates-1)])
  diag(aij) = 0
  return(aij)
}

est_bij <- function(params, nstates, parscale){
  bij <- matrix(params[(3 * nstates**2 + 1):(4 * nstates ** 2)], nstates, nstates)
  diag(bij) = 0
  bij[nstates, nstates-1] <- 1 - sum(bij[nstates, 1:(nstates-2)])
  bij[, nstates] = 1 - rowSums(bij[, 1:(nstates-1)])
  return(bij)
}
