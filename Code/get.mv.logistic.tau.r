# Author: Brigitte Dorner
# Last Modified: 12 May 2021 by Carrie Holt

# A quick-and-dirty grid search for finding a good tau parameter for the multivariate logistic distribution
# used to represent relative proportions of recruits in different age groups from
# J.T. Schnute and L.J. Richards 1995.
# The Influence of error on population estimates from catch-age models
# CJFAS 52:2063-2077  (see equations S.9 and S.10 on p. 2067)
#
# parameters:
# age.struct.data: a matrix or data frame containing age structure data for one stock
# columns in the age structure data should be headed "RecruitX.Y",
# where X is the number of winters in freshwater
# and Y is the number of winters in the ocean.
# Note: recruit numbers should be the actual counts of individual fish, i.e., NOT scaled to
# 100s or 1000s of fish. This is because the program adds 1 to all recruit counts to avoid
# problems with having to take the log of 0.
#
# return value:
# the tau value that creates the best fit to the age structure data
#
# side effects:
# plots are created of
# (1) value of the objective function (SSQ of differences in standard error between actual and simulated data)
# 		vs. tau
# (2) simulated vs. actual pdf for proportions in each age group
#THese are by age-structured data sets (by BY)
#dum<-as.data.frame(read.table("H:/FromUSB/DFO/Benchmarks/Modelling/Simulation/ESages.dat",header=TRUE))

#    "get.mv.logistic.tau"(dum)

"get.mv.logistic.tau" <- function(age.struct.data)
{
	# # get rid of everything except age structure data
	# if (length(grep("Recruit?.?", names(age.struct.data))) > 0)
	# 	df <- age.struct.data[,grep("Recruit?.?", names(age.struct.data))]
	# else {if (length(grep("R?.?", names(age.struct.data))) > 0)
	# 	df <- age.struct.data[,grep("R?.?", names(age.struct.data))]
	# else stop("could not find age structure info in input data!")
	# df <- df+1
	# R <- apply(df, 1, sum, na.rm=T)											# get total number of recruits for each brood year
	# props <- sweep(df, 1, R, FUN="/") 											# convert recruit numbers to proportions
	#

  # If input data = proportions
  props <- age.struct.data


	target.fun <- function(tau, expected.log.props, target.stderr, n = 10000) # objective function: minimize difference between ssq of observed vs. simulated standard deviations for proportions
	{
		cat(". ")
		Tau <- rep(tau, n)
		simulated.props <- t(sapply(Tau, mv.logistic.log.props, log.props= expected.log.props)) 	# n draws from the multivariate logistic with parameter tau
		sim.stderr <- apply(simulated.props, 2, sd)													# get standard deviation for simulated proportions
		return(sum((sim.stderr - target.stderr)^2, na.rm=T))
	}

	"plot.recruit.props" <- function(props, tau, n=10000)
	{
   		tau <- rep(tau, n)
   		sim.props <- t(sapply(tau, mv.logistic.log.props, log.props=apply(log(props), 2, mean)))
		actual <- apply(props, 2, density, n=100, from=0, to=1, na.rm=T)
		sim <- apply(sim.props, 2, density, n=100, from=0, to=1, na.rm=T)
		windows()#graphsheet()
		for (p in which(!is.na(props[1, ])))
		{
			yrange <- c(0,1)#c(0, max(c(as.numeric(actual[[p]]), as.numeric(sim[[p]])), na.rm=T))
			plot(x=0, y=0, type="n", xlim=c(0, 1), ylim=yrange, ylab="density", xlab="")
			title(main=paste("Proportion of ", dimnames(props)[[2]][p], "\nempirical (dots) vs. modeled distribution", sep="") )
			points(actual[[p]], col=3, type="p", pch=18)
			points(sim[[p]], col=4, type="l")
		}
	}

	# grid search for best tau:
	tau <- seq(from = 0, to=3, by=0.1)
	result <- data.frame(tau=tau, objective=rep(NA, length(tau)))
	cat("finding scale parameter for variability in relative recruit proportions:\n")
	for (t in tau)
	{
		cat(". ")
		result[t*10 + 1, "objective"] <- target.fun(t, apply(log(props), 2, mean), apply(props, 2, sd))
	}
	cat("\n")
	best.tau <- tau[result[, "objective"] == min(result[, "objective"])]
	windows()#graphsheet()
	plot(tau, result$objective, xlab="tau", ylab="objective function")
	plot.recruit.props(props, best.tau)
	return(list(obj.fun=result, best.tau=best.tau))
}


#------------------------------------------------------------------------------------------
# This is the multivariate logistic distribution used to model stochastic fluctuations
# in age group proportions in
#
# J.T. Schnute and L.J. Richards 1995.
# The Influence of error on population estimates from catch-age models
# CJFAS 52:2063-2077  (see equations S.9 and S.10 on p. 2067)
#
# The function takes a vector of expected or predicted proportions and adds errors
# (observation errors in Schnute & Richards), thus generating a vector of actual (or observed)
# proportions.  The parameter tau controls the amount of error added. For tau = 0 predicted
# and actual proportions are identical. For tau -> Inf actual proportions become random and
# unrelated to the predicted values.
#

"mv.logistic.log.props" <- function(tau, log.props)
{
	p <- log.props[!is.na(log.props)]
	N <- length(p)
	rnum<-rnorm(N, sd=1)
	eps <- tau *rnum
	x <- p + eps - (1/N) * sum(p + eps)
	result <- rep(NA, N)
	result[!is.na(log.props)] <- exp(x)/sum(exp(x))
	result
}

# This function simulates some ppns from means and a tau
"get.traj" <- function(means,tau)
{
	log.props<-log(means)
	sim<-"mv.logistic.log.props"(tau,log.props)
  sim
	}

#"get.traj"(c(0.268,0.727, 0.007),0.5)	# Togiak mean values for age 1.2, 1.3, 1.4
#"get.traj"(c(0.012,0.937, 0.051),0.9)	# Chilko mean values for age 1.1, 1.2, 1.3


