testmcmc <- read.csv("samSimInputs/Ricker_mcmc.csv")
# testmcmc <- read.csv("samSimInputs/Even_mcmc.csv")
# testmcmc <- read.csv("samSimInputs/Same_mcmc.csv")
# testmcmc <- read.csv("samSimInputs/SameProd_mcmc.csv")
# testmcmc <- read.csv("samSimInputs/SameBeta_mcmc.csv")
testmcmc <- read.csv("samSimInputs/SameSREP_mcmc.csv")

s<-(1:1000)*100
r <- exp(testmcmc$alpha[2])*s*exp(-testmcmc$beta[2]*s)
plot(x=s, y=r, ylim=c(0,15000), xlim=c(0,100000), type="l")
for (i in 1:10) lines(x=s, y=exp(testmcmc$alpha[i])*s*exp(-testmcmc$beta[i]*s))
title("Same SREP")


# When Productivity is the same but SREP varies (and hence beta), then LRP is not sensitivity to ER
# Plot Sgens for each curve in SameProd, add replacement line, and replacement line with high ER
# Each curve with different SREPs should be equally sensitive to high ER as population status they are scaled to SREP via Sgen, ???


# When SREP (or beta) is the same and productivity varies, then LRP is slightly sensitive to ER

# However when Productivity  and SREP change, then LRP is very sensitive to ER

