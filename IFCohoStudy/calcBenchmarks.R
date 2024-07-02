
require('dplyr')

# Functions
sGenOptimum <- function ( S, theta ) {
  # Function called from sGenSolver
  loga <- theta[1]
  b <- theta[2]
  prt <- S * exp( loga - b * S)
  sMSY <- ( 1 - gsl::lambert_W0 (exp ( 1 - loga) ) ) / b
  epsilon <- log(sMSY) - log(prt)
  nLogLike <- - sum( dnorm ( epsilon, 0, 1, log = T))
  return( nLogLike )
}


sGenSolver <- function (loga, b) {
  # Function to estimate Sgen from loga and b Ricker parameters
  theta <- c(loga, b)
  sMSY <- (1 - gsl::lambert_W0(exp(1 - loga))) / b
  fit <- optimize(f = sGenOptimum, interval = c(0, sMSY),
                  theta = theta)
  return(fit$minimum)
}

sGenSolver(adjProd,mbeta)

SMSY <- function(loga, beta){
  out <- (1 - gsl::lambert_W0(exp(1 - loga))) / beta
  # out <- (loga/beta)*(0.5-0.07*loga)
  return(out)
}

alpha_adj <- function(loga, gamma){
  ave_m <- 0.011584
  out <- loga + gamma*log(ave_m)
  return(out)
} 




# Read in data
mc1 <- read.csv(here::here("SamSimInputs","SR_IndivRicker_Surv_mcmc.csv"))
mc2 <- read.csv(here::here("SamSimInputs","SR_IndivRicker_SurvCap_mcmc.csv"))
# Combine posteriors
mc <- rbind(mc1,mc2)


mc <- mc %>% 
  mutate (logalpha = alpha) %>% 
  mutate(logalpha_prime = alpha_adj(logalpha, gamma))

gamma <- mc %>% group_by(stk) %>% summarize(m.gamma = mean(gamma))

# on  line 277 of fitmcmc.R, logA is called post_longlapha, so I think this is called alpha at output not logA
# on lline 315 adn 319, adjProd is calculated as exp(logalpha-prime)

# mc <- mc %>% 
  # mutate(Sgen=sGenSolver(logalpha_prime, beta)) %>% #sGEnSOlver does not vectorize well.... see below
  # mutate(Smsy = SMSY(logalpha_prime, beta))


Sgen <- NA
Smsy <- NA
for(j in 1:length(mc$alpha)) {
  Sgen[j] <- sGenSolver(mc$logalpha_prime[j], mc$beta[j])
  Smsy[j] <- SMSY(mc$logalpha_prime[j], mc$beta[j])
}
mc <- mc %>% mutate(Sgen=Sgen, Smsy=Smsy, Smsy2=Smsy2)

alpha <- mc %>% group_by(stk) %>% summarize(m.alpha=mean(alpha))

logalpha_prime <- mc %>% group_by(stk) %>% summarize(m.logalpha_prime=mean(logalpha_prime))
# adjProd in the paper-- is not logged

# adjProd <- (exp(logalpha_prime$m.logalpha_prime))

## These adjProd values differ from Table C1, maybe because adjProd is estimated within TMB code using age-specific  mortality rates??

beta <- mc %>% group_by(stk) %>% summarize(m.beta=mean(beta))

Sgen <- mc %>% group_by(stk) %>% summarize(m.Sgen=mean(Sgen))
# This is the same as Michael Arbeider's email Feb 14, 10:31am

Smsy <- mc %>% group_by(stk) %>% summarize(m.Smsy=mean(Smsy))
# These values are not the same as MA's



