#===============================================================================
# Code to derive uncertainty intevals for UMSY given uncertainty in productivity
# For Brown et al. (in revisions, 2024)
# Using Hilborn and Walters (1992) approximation
#===============================================================================


# Estimate of productivity from life-cycle model for WCVI CK
lna <- 1
# Uncertainty in productivity from life-cycle model for WCVI CK
sig.lna <- 0.5

U_MSY  <- function (lna){0.5 * lna - 0.07 *lna^2 }

lna.vec <- rnorm(10000,lna, sig.lna)

umsy <- U_MSY(lna.vec)
quantile(umsy, p=c(0.025, 0.5,0.975))

