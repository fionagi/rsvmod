library(deSolve)
library(conmat)
library(dplyr)
library(pracma)
library(ggplot2)

load("data/data.greaterPerth.2011.rda")
load("data/data.contact.monthly.rda")
source("R/functions.r")
source("R/parameters.r")
source("R/models.r")
source("R/diagnostic.r")

#Plot age-to-risk function using fixed and fitted values
#Figure 3 in paper
plot_ageToRisk(A = A_fit, B = B_fit, C = C_fit, D = D_fit,
               E = E_fit, maxAge = 12)

#Set up model
age <- age_structure(year = "2011")

mixing <- data.contact.monthly #contact matrix

#create vectors
omega_vect <- as.vector(rep(1, age$nAges))
omega_vect[!(age$age_years < 10)] <- omega #reduced infectiousness
                                            #for older age groups
sigma_vect <- as.vector(c(sigma, rep(1, (age$nAges-6)))) #reduced susceptibility
                                                    #due to maternal protection

# output time steps
n_times <- 400
times     <- seq(0, n_times, 0.25)

#Run base model (without preterm risk-group separation)
#Initialise y - state variables
y0_seir <- matrix(0, age$nAges, 8)
y0_seir[,1] <- 0.99 * age$pop_groups #S0 term
y0_seir[,3] <- 0.01 * age$pop_groups #I0 term
y0 <- cbind(y0_seir, matrix(0, age$nAges, 6))

#set parameter values
pars_base <- list(
  b0 = b0_fit,
  b1 = b1_fit,
  phi = phi_fit,
  delta = delta,
  gamma0 = gamma0,
  gamma1 = gamma1,
  nu = nu,
  omega_vect = omega_vect,
  omegaE = omegaE,
  A = A_fit,
  B = B_fit,
  C = C_fit,
  D = D_fit,
  age_in = age$age_in,
  age_out = age$age_out,
  age_months = age$age_years*12,
  sigma_vect = sigma_vect,
  sigmaE = sigmaE,
  mixing = mixing,
  nAges = age$nAges)

#Solve ODEs
deSolve_out <- deSolve::ode(y0, times, deSolve_base, pars_base)

#Aggregate output into age groups for plots
modelOut_base <- aggregate_output(deSolveOut = deSolve_out,
                             times = times,
                             nAges = age$nAges,
                             model = "base",
                             old = 1)
#Find index range matching up to complete calendar years
#spanning observational data and assuming peak month is July
iRange <- findIndexRange(counts = modelOut_base$combined, nYears = 13)

#Plot to compare model output to observed. Observed data can not be provided
#due to confidentiality but if available would be in the same format as
#model output. Default inputs for plot function assumes obs1 is observed data
#and obs2 is model output. Here we have put in model output twice
comp_data(obs1 = modelOut_base$combined[iRange$start:iRange$end,],
          obs2 = modelOut_base$combined[iRange$start:iRange$end,])

#Run Risk model (with preterm stratification)
#Initialise y - state variables - pre-terms
y0_seir <- matrix(0, age$nAges, 16)
y0_seir[,1] <- (1 - alpha) * 0.99 * age$pop_groups #S0 term
y0_seir[,3] <- (1 - alpha) * 0.01 * age$pop_groups #I0 term
y0_seir[,9] <- alpha * 0.99 * age$pop_groups #S0_bar preterm
y0_seir[,11] <- alpha * 0.01 * age$pop_groups #I0_bar preterm
y0 <- cbind(y0_seir, matrix(0, age$nAges, 14))

#set parameter values
pars_pT <- list(
  b0 = b0_fit,
  b1 = b1_fit,
  phi = phi_fit,
  delta = delta,
  gamma0 = gamma0,
  gamma1 = gamma1,
  nu = nu,
  omega_vect = omega_vect,
  omegaE = omegaE,
  alpha = alpha,
  A = A_fit,
  B = B_fit,
  C = C_fit,
  D = D_fit,
  E = E_fit,
  age_in = age$age_in,
  age_out = age$age_out,
  age_months = age$age_years*12,
  sigma_vect = sigma_vect,
  sigmaE = sigmaE,
  mixing = mixing,
  nAges = age$nAges)

deSolve_out <- deSolve::ode(y0, times, deSolve_preTerm, pars_pT)

modelOut_pT <- aggregate_output(deSolveOut = deSolve_out,
                             times = times,
                             nAges = age$nAges,
                             model = "risk",
                             old = 1)

mrange <- findIndexRange(counts = modelOut_pT$combined, month = 7, nYears = 13)

#Plot model predicted output. Usually would model against observed data but
#as above, we put the modelled output in the format for the observed data and
#hence plot model output twice
dummyObs <- cbind(modelOut_pT$term[mrange$start:mrange$end,1],
                  modelOut_pT$pre[mrange$start:mrange$end,1] )
for(i in 2:5) dummyObs <- cbind(dummyObs,
                                modelOut_pT$term[mrange$start:mrange$end,i],
                                modelOut_pT$pre[mrange$start:mrange$end,i])
all_plots <- plot_riskFit(term = modelOut_pT$term[mrange$start:mrange$end,],
             pTerm = modelOut_pT$pre[mrange$start:mrange$end,], obs = dummyObs,
                years = 2000:2013, old = 1)
#First plot gives preterm and term predictions on same plot
all_plots[[1]]
#Second, just term
all_plots[[2]]
#Third, just preterm
all_plots[[3]]
