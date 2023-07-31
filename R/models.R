#'ODE model with multiple exposures and pre-term risk group
#'
#'This function is used with R package deSolve to solve the ODEs
#'
#' @param t a vector of times at which model output is wanted
#' @param y a matrix of the initial values of the state variables.
#'          Each column corresponds to a state variable.
#' @param parms the required model parameters and their values. For this model,
#'              parms <- pars_pT <- list(b0 = b0,
#'                                       b1 = b1,
#'                                       phi = phi,
#'                                       delta = delta,
#'                                       gamma0 = gamma0,
#'                                       gamma1 = gamma1,
#'                                       nu = nu,
#'                                       omega_vect = omega_vect,
#'                                       omegaE = omegaE,
#'                                       alpha = alpha,
#'                                       A = A,
#'                                       B = B,
#'                                       C = C,
#'                                       D = D,
#'                                       E = E,
#'                                       age_in = age_in,
#'                                       age_out = age_out,
#'                                       age_months = age_months,
#'                                       sigma_vect = sigma_vect,
#'                                       sigmaE = sigmaE,
#'                                       mixing = mixing,
#'                                       nAges = nAges)
#' @return list
#' @export
deSolve_preTerm <- function(t, y, parms) {

  y <- matrix(y, nrow = 75, ncol = 30)


  with(as.list(c(y, parms)), {

    #term population
    S0 <- y[, 1]
    E0 <- y[, 2]
    I0 <- y[, 3]
    R0 <- y[, 4]

    S1 <- y[, 5]
    E1 <- y[, 6]
    I1 <- y[, 7]
    R1 <- y[, 8]

    #pre-term population
    S0_bar <- y[, 9]
    E0_bar <- y[, 10]
    I0_bar <- y[, 11]
    R0_bar <- y[, 12]

    S1_bar <- y[, 13]
    E1_bar <- y[, 14]
    I1_bar <- y[, 15]
    R1_bar <- y[, 16]

    N <- S0 + E0 + I0 + R0 + S1 + E1 + I1 + R1 +
      S0_bar + E0_bar + I0_bar + R0_bar + S1_bar + E1_bar + I1_bar + R1_bar

    tau_shift <- c((1 - alpha) * age_in[1], age_in[-1])
    tau <- c(age_out[-nAges], (1 - alpha) * age_out[nAges])
    tauBar_shift <- c(alpha * age_in[1], age_in[-1])
    tauBar <- c(age_out[-nAges], alpha * age_out[nAges])

    temp <- omega_vect * (I0 + omegaE*I1 + I0_bar + omegaE*I1_bar) / N
    s <- sweep(mixing, MARGIN = 2, temp, FUN = "*")
    lambda <- b0 * (1 + b1 * cos(2* pi *t / 12 + phi)) * rowSums(s)

    infect0 <- lambda * sigma_vect * S0
    infect1 <- lambda * sigma_vect * sigmaE * S1
    infect0_bar <- lambda * sigma_vect * S0_bar
    infect1_bar <- lambda * sigma_vect * sigmaE * S1_bar

    S0_shift <- c(1, S0[-nAges])
    E0_shift <- c(0, E0[-nAges])
    I0_shift <- c(0, I0[-nAges])
    R0_shift <- c(0, R0[-nAges])

    S1_shift <- c(0, S1[-nAges])
    E1_shift <- c(0, E1[-nAges])
    I1_shift <- c(0, I1[-nAges])
    R1_shift <- c(0, R1[-nAges])

    S0bar_shift <- c(1, S0_bar[-nAges])
    E0bar_shift <- c(0, E0_bar[-nAges])
    I0bar_shift <- c(0, I0_bar[-nAges])
    R0bar_shift <- c(0, R0_bar[-nAges])

    S1bar_shift <- c(0, S1_bar[-nAges])
    E1bar_shift <- c(0, E1_bar[-nAges])
    I1bar_shift <- c(0, I1_bar[-nAges])
    R1bar_shift <- c(0, R1_bar[-nAges])

    dS0 <- tau_shift * S0_shift - infect0 - tau * S0
    dE0 <- tau_shift * E0_shift + infect0 - delta * E0 - tau * E0
    dI0 <- tau_shift * I0_shift + delta*E0 - gamma0 * I0 - tau * I0
    dR0 <- tau_shift * R0_shift + gamma0 * I0 - nu * R0 - tau * R0

    dS1 <- tau_shift * S1_shift - infect1 + nu * (R1+ R0) - tau * S1
    dE1 <- tau_shift * E1_shift + infect1 - delta * E1 - tau * E1
    dI1 <- tau_shift * I1_shift + delta*E1 - gamma1 * I1 - tau * I1
    dR1 <- tau_shift * R1_shift + gamma1 * I1 - nu * R1 - tau * R1

    dS0_bar <- tauBar_shift * S0bar_shift - infect0_bar - tauBar * S0_bar
    dE0_bar <- tauBar_shift * E0bar_shift + infect0_bar - delta * E0_bar - tauBar * E0_bar
    dI0_bar <- tauBar_shift * I0bar_shift + delta*E0_bar - gamma0 * I0_bar - tauBar * I0_bar
    dR0_bar <- tauBar_shift * R0bar_shift + gamma0 * I0_bar - nu * R0_bar - tauBar * R0_bar

    dS1_bar <- tauBar_shift * S1bar_shift - infect1_bar + nu * (R1_bar + R0_bar) - tauBar * S1_bar
    dE1_bar <- tauBar_shift * E1bar_shift + infect1_bar - delta * E1_bar - tauBar * E1_bar
    dI1_bar <- tauBar_shift * I1bar_shift + delta*E1_bar - gamma1 * I1_bar - tauBar * I1_bar
    dR1_bar <- tauBar_shift * R1bar_shift + gamma1 * I1_bar - nu * R1_bar - tauBar * R1_bar

    #Term: First exposure - Incidence and detected incidence (hospitalisations)
    dinc0 <- infect0
    ddetinc0 <- (A * exp(-B * age_months) + C) * infect0

    #Term: Subsequent exposures
    dinc1 <- infect1
    ddetinc1 <- (A * exp(-B * age_months) + C) * D * infect1

    #Term: combined
    dincT <- dinc0 + dinc1
    ddetincT <- ddetinc0 + ddetinc1

    #Preterm: First exposure - Incidence and detected incidence (hospitalisations)
    dinc0_bar <- infect0_bar
    ageRiskF <- (A * exp(-B * age_months) + C) * E
    ageRiskF[which(ageRiskF >1)] <- 1
    ddetinc0_bar <- ageRiskF * infect0_bar

    #Preterm: Subsequent exposures
    dinc1_bar <- infect1_bar
    ageRiskS <- (A * exp(-B * age_months) + C) * E * D
    ageRiskS[which(ageRiskS >1)] <- 1
    ddetinc1_bar <- ageRiskS * infect1_bar

    #Preterm: combined
    dincPT <- dinc0_bar + dinc1_bar
    ddetincPT <- ddetinc0_bar + ddetinc1_bar

    #Term + Preterm: combined
    dinc <- dincT + dincPT
    ddetinc <- ddetincT + ddetincPT

    return(list(c(dS0,
                  dE0,
                  dI0,
                  dR0,
                  dS1,
                  dE1,
                  dI1,
                  dR1,
                  dS0_bar,
                  dE0_bar,
                  dI0_bar,
                  dR0_bar,
                  dS1_bar,
                  dE1_bar,
                  dI1_bar,
                  dR1_bar,
                  dinc0,
                  ddetinc0,
                  dinc1,
                  ddetinc1,
                  dincT,
                  ddetincT,
                  dinc0_bar,
                  ddetinc0_bar,
                  dinc1_bar,
                  ddetinc1_bar,
                  dincPT,
                  ddetincPT,
                  dinc,
                  ddetinc)))
  })
}

#'Base ODE model with multiple exposures
#'
#'This function is used with R package deSolve to solve the ODEs
#'
#' @param t a vector of times at which model output is wanted
#' @param y a matrix of the initial values of the state variables.
#'          Each column corresponds to a state variable.
#' @param parms the required model parameters and their values. For this model,
#'              parms <- pars_pT <- list(b0 = b0,
#'                                       b1 = b1,
#'                                       phi = phi,
#'                                       delta = delta,
#'                                       gamma0 = gamma0,
#'                                       gamma1 = gamma1,
#'                                       nu = nu,
#'                                       omega_vect = omega_vect,
#'                                       omegaE = omegaE,
#'                                       A = A,
#'                                       B = B,
#'                                       C = C,
#'                                       D = D,
#'                                       age_in = age_in,
#'                                       age_out = age_out,
#'                                       age_months = age_months,
#'                                       sigma_vect = sigma_vect,
#'                                       sigmaE = sigmaE,
#'                                       mixing = mixing,
#'                                       nAges = nAges)
#' @return list
#' @export
deSolve_base <- function(t, y, parms) {

  y <- matrix(y, nrow = 75, ncol = 14)


  with(as.list(c(y, parms)), {

    #term population
    S0 <- y[, 1]
    E0 <- y[, 2]
    I0 <- y[, 3]
    R0 <- y[, 4]

    S1 <- y[, 5]
    E1 <- y[, 6]
    I1 <- y[, 7]
    R1 <- y[, 8]

    N <- S0 + E0 + I0 + R0 + S1 + E1 + I1 + R1

    tau_shift <- c(age_in[1], age_in[-1])
    tau <- c(age_out[-nAges], age_out[nAges])

    temp <- omega_vect * (I0 + omegaE*I1) / N
    s <- sweep(mixing, MARGIN = 2, temp, FUN = "*")
    lambda <- b0 * (1 + b1 * cos(2* pi *t / 12 + phi)) * rowSums(s)

    infect0 <- lambda * sigma_vect * S0
    infect1 <- lambda * sigma_vect * sigmaE * S1

    S0_shift <- c(1, S0[-nAges])
    E0_shift <- c(0, E0[-nAges])
    I0_shift <- c(0, I0[-nAges])
    R0_shift <- c(0, R0[-nAges])

    S1_shift <- c(0, S1[-nAges])
    E1_shift <- c(0, E1[-nAges])
    I1_shift <- c(0, I1[-nAges])
    R1_shift <- c(0, R1[-nAges])

    dS0 <- tau_shift * S0_shift - infect0 - tau * S0
    dE0 <- tau_shift * E0_shift + infect0 - delta * E0 - tau * E0
    dI0 <- tau_shift * I0_shift + delta*E0 - gamma0 * I0 - tau * I0
    dR0 <- tau_shift * R0_shift + gamma0 * I0 - nu * R0 - tau * R0

    dS1 <- tau_shift * S1_shift - infect1 + nu * (R1+ R0) - tau * S1
    dE1 <- tau_shift * E1_shift + infect1 - delta * E1 - tau * E1
    dI1 <- tau_shift * I1_shift + delta*E1 - gamma1 * I1 - tau * I1
    dR1 <- tau_shift * R1_shift + gamma1 * I1 - nu * R1 - tau * R1

    #First exposure - Incidence and detected incidence (hospitalisations)
    dinc0 <- infect0
    ddetinc0 <- (A * exp(-B * age_months) + C) * infect0

    #Subsequent exposures
    dinc1 <- infect1
    ddetinc1 <- (A * exp(-B * age_months) + C) * D * infect1

    #Combined
    dinc <- infect0 + infect1
    ddetinc <- (A * exp(-B * age_months) + C) * (infect0 + D * infect1)

    return(list(c(dS0,
                  dE0,
                  dI0,
                  dR0,
                  dS1,
                  dE1,
                  dI1,
                  dR1,
                  dinc0,
                  ddetinc0,
                  dinc1,
                  ddetinc1,
                  dinc,
                  ddetinc)))
  })
}

#'Base ODE model with multiple exposures but NO continuous ageing
#'
#'This function is used with R package deSolve to solve the ODEs
#'
#' @param t a vector of times at which model output is wanted
#' @param y a matrix of the initial values of the state variables.
#'          Each column corresponds to a state variable.
#' @param parms the required model parameters and their values. For this model,
#'              parms <- pars_pT <- list(b0 = b0,
#'                                       b1 = b1,
#'                                       phi = phi,
#'                                       delta = delta,
#'                                       gamma0 = gamma0,
#'                                       gamma1 = gamma1,
#'                                       nu = nu,
#'                                       omega_vect = omega_vect,
#'                                       omegaE = omegaE,
#'                                       A = A,
#'                                       B = B,
#'                                       C = C,
#'                                       D = D,
#'                                       age_months = age_months,
#'                                       sigma_vect = sigma_vect,
#'                                       sigmaE = sigmaE,
#'                                       mixing = mixing,
#'                                       nAges = nAges)
#' @return list
#' @export
deSolve_cohortAgeing <- function(t, y, parms) {

  y <- matrix(y, nrow = 75, ncol = 14)


  with(as.list(c(y, parms)), {

    #term population
    S0 <- y[, 1]
    E0 <- y[, 2]
    I0 <- y[, 3]
    R0 <- y[, 4]

    S1 <- y[, 5]
    E1 <- y[, 6]
    I1 <- y[, 7]
    R1 <- y[, 8]

    N <- S0 + E0 + I0 + R0 + S1 + E1 + I1 + R1

    temp <- omega_vect * (I0 + omegaE*I1) / N
    s <- sweep(mixing, MARGIN = 2, temp, FUN = "*")
    lambda <- b0 * (1 + b1 * cos(2* pi *t / 12 + phi)) * rowSums(s)
    infect0 <- lambda * sigma_vect * S0
    infect1 <- lambda * sigma_vect * sigmaE * S1

    dS0 <- -infect0
    dE0 <- infect0 - delta * E0
    dI0 <- delta*E0 - gamma0 * I0
    dR0 <- gamma0 * I0 - nu * R0

    dS1 <- -infect1 + nu * (R1+ R0)
    dE1 <- infect1 - delta * E1
    dI1 <- delta * E1 - gamma1 * I1
    dR1 <- gamma1 * I1 - nu * R1

    #First exposure - Incidence and detected incidence (hospitalisations)
    dinc0 <- infect0
    ddetinc0 <- (A * exp(-B * age_months) + C) * infect0

    #Subsequent exposures
    dinc1 <- infect1
    ddetinc1 <- (A * exp(-B * age_months) + C) * D * infect1

    #Combined
    dinc <- infect0 + infect1
    ddetinc <- (A * exp(-B * age_months) + C) * (infect0 + D * infect1)

    return(list(c(dS0,
                  dE0,
                  dI0,
                  dR0,
                  dS1,
                  dE1,
                  dI1,
                  dR1,
                  dinc0,
                  ddetinc0,
                  dinc1,
                  ddetinc1,
                  dinc,
                  ddetinc)))
  })
}

