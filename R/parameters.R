#Set parameter values

#Age group structure
final_age <- 80
numMonthly <- 5*12 #number of age groups (from age 0), in 1 month groups
                   #Assume the first 5 years are in monthly groups

#Fixed values
omega <- 1 #reduced infectiousness in older age groups
omegaE <- 0.7 #reduced infectiousness from prior exposure
nu <- 0.132 #immunity period 230 days 1/(230/(365.25/12))
delta <- 7.610 #latent period 4 days 1/(4/(365.25/12))
gamma0 <- 3.044 #infectious period after first exposure 10 days 1/(10/(365.25/12))
gamma1 <- 4.348 #infectious period after subsequent exposures 7 days 1/(7/(365.25/12))
sigma <- as.vector(c(0.08, 0.45, 0.45, 1, 1, 1)) #maternal protection for first 6 months, scaled susceptibility
sigmaE <- 0.77 #reduced susceptibility from prior exposure
alpha <- 0.0849 #proportion of births born pre-term (<37 weeks)

#Fitted values
b0_fit <- 0.02040069
b1_fit <- 0.33962388
phi_fit <- 0.98456 #phase shift

#Risk of Hosp. = (Ae^(-Bi)+C)D*E
A_fit <- 0.51442670 #max. risk of hospitalisation (at age 0)
B_fit <- 0.37756892 #decay constant
C_fit <- 0.015 #min. risk across all ages
D_fit <- 0.2 #scaling factor for previous hospitalisations
E_fit <- 2.63292393 #scaling of risk for pre-terms
