library(rcarbon)
library(nimbleCarbon)
library(truncnorm)
library(coda)
library(parallel)
library(MCMCvis)

binsize = 50
ncores = 3
nsim = 1000
runm = 50
method = "uncalsample"

load("perm_prep.RData")
load("mcmc_prep.RData")

# Permutation test

pt1 <- permTest(x = calDates.ap, marks = as.factor(ap$Affiliation.num), timeRange = c(2000,100),
                nsim = 1000, bins = bins.ap, spdnormalised = FALSE, datenormalised = FALSE)

pt2 <- permTest(x = calDates.bp, marks = as.factor(bp$Affiliation.num), timeRange = c(2000,100),
                nsim = 1000, bins = bins.bp, spdnormalised = FALSE, datenormalised = FALSE)

save(pt1, pt2,
     file = "perm_res.RData")

# Bayesian MCMC

cl <- makeCluster(ncores)

runFun <- function(seed, data, constants, m.dates) {
  library(nimbleCarbon)
  initsFunction = function() list(r1=rnorm(1,sd=0.0004), #0.0004 original value
                                  r2=rnorm(1,sd=0.0004),
                                  chp=round(truncnorm::rtruncnorm(1,mean=600,a=150,b=1600)))
  
  m.dblexp <- nimbleCode({
    for (i in 1:N){
      theta[i] ~ dDoubleExponentialGrowth(a=1600,b=150,r1=r1,r2=r2,mu=changept);
      mu[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]);
      sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[]);
      sd[i] <- (sigma[i]^2+sigmaCurve[i]^2)^(1/2);
      X[i] ~ dnorm(mean=mu[i],sd=sd[i]);
    }
    r1 ~ dnorm(0,sd=0.004);
    r2 ~ dnorm(0,sd=0.004); 
    chp ~ T(dnorm(600,sd=100),150,1600);
    changept <- round(chp);
  })  
  
  results <- nimbleMCMC(code = m.dblexp,constants = constants,data = data,niter = 50000, nchains = 3, thin=2, nburnin = 5000, 
                        progressBar = FALSE, monitors=c('r1', 'r2', 'chp'), inits=initsFunction, 
                        samplesAsCodaMCMC=TRUE,setSeed=seed) # lower the iterations and burnins to speed up, at own risk
  
  return(results)
}

seeds = c(123, 456, 789)

chain_output = parLapply(cl = cl, X = seeds, fun = runFun, data = data, constants=constants,m.dates=medDates) # run chains
stopCluster(cl)

chain_coda = coda::mcmc.list(c(chain_output[[1]], chain_output[[2]], chain_output[[3]]))

rhat = gelman.diag(chain_output[[1]])
ess = effectiveSize(chain_output[[1]])
save(ess, rhat, file="mcmc_diags.RData")

params = list(r1 = c(chain_coda$chain1[,'r1'], chain_coda$chain2[,'r1'], chain_coda$chain3[,'r1']),
              r2 = c(chain_coda$chain1[,'r2'], chain_coda$chain2[,'r2'], chain_coda$chain3[,'r2']),
              mu = c(chain_coda$chain1[,'chp'], chain_coda$chain2[,'chp'], chain_coda$chain3[,'chp']))
save(params, file='mcmc_params.RData')

pp.check=postPredSPD(obs.data$CRA,obs.data$Error, # run posterior check
                     model = dDoubleExponentialGrowth,a = 1600,b=150,params=params, calCurve = "intcal20",
                     nsim = 500,ncores = 5,verbose=TRUE,method=method, spdnormalised = FALSE,datenormalised = TRUE)

save(chain_output, chain_coda, pp.check, file="mcmc_res.RData")
