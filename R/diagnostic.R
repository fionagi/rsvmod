############################################################################
# DIAGNOSTIC FUNCTIONS
############################################################################
#' Comparison plot of 2 data sets over same time range
#'
#' @param obs1 table of monthly hospitalisation data in age groups (colored dots)
#' @param obs2 table of monthly hospitalisation data in age groups (dashed line)
#' @param ylabel default is Hospitalisations
#' @param years year labels for x-axis
#' @param text plot text, default is RSV hospitalisations: observed (colour .) and model (black --)
#' @return ggplot
#'
comp_data <- function(obs1, obs2, ylabel = "Hospitalisations", years = 2000:2013,
                      text = "RSV hospitalisations: observed (colour .) and model (black --)")
{
  attr <- c("Month", "Age", ylabel)

  melt_obs1 <- reshape2::melt(as.matrix(obs1))
  colnames(melt_obs1) <- attr

  melt_obs2 <- reshape2::melt(as.matrix(obs2))
  colnames(melt_obs2) <- attr

  g1 <- ggplot() +
    geom_point(data = melt_obs1, aes_string(x = attr[1], y = attr[3], col = attr[2])) +
    geom_line(data = melt_obs2, aes_string(x = attr[1], y = attr[3]), col = "black", linetype = "dashed") +
    scale_x_continuous(name = "Year", breaks = seq(0, nrow(obs1), 12), labels = years)+
    ggplot2::theme_bw()+
    theme(legend.position="none")+
    ggtitle(text)+
    theme(plot.title = element_text(size = 10))+
    facet_wrap(~Age, ncol = 1)

  return(g1)
}

############################################################################
#' FOR OLD MODEL - POSSIBLY DELETE/ARCHIVE
#' Compare model(s) predicted hospitalisations to observed hospitalisation data
#'
#' @param pars1 list of parameter values from model 1
#' @param pars2 list of compartment counts from model 2 if comparing
#' @param obs table of monthly hospitalisation data for 4 age groups
#' @param years year labels for x-axis
#' @param max_t the "run-in" time used in fitting the model
#' @return ggplot
#'
plot_fitToData <- function(pars1, pars2 = NA, obs, years = 2000:2013, max_t = 200)
{
  nMonths <- nrow(obs)
  parSet <- list(pars1, pars2)

  modelList = NA
  for(i in 1:length(parSet))
  {
    if(is.na(parSet[[i]])) break

    counts <- run_rsv_model(b0 = parSet[[i]]$b0,
                            b1 = parSet[[i]]$b1,
                            phi = parSet[[i]]$phi,
                            h_k1 = parSet[[i]]$h_k1,
                            h_k2 = parSet[[i]]$h_k2,
                            h_k3 = parSet[[i]]$h_k3,
                            h_k4 = parSet[[i]]$h_k4,
                            max_t = max_t,
                            init_conds_from_file = 0,
                            save_init_conds = 0)

    mrange <- findIndexRange(counts = counts, peakMethod = 1, month = 7,
                             ageGroup = 5, obs = obs)

    mod_DetIncidence <- as.data.frame(counts$DetIncidence[mrange$start:mrange$end,])
    DetInc_1_3 <- apply(mod_DetIncidence[, 1:3], 1, sum)
    DetInc_4_6 <- apply(mod_DetIncidence[, 4:6], 1, sum)
    DetInc_7_12 <- apply(mod_DetIncidence[, 7:12], 1, sum)
    DetInc_13_24 <- apply(mod_DetIncidence[, 13:24], 1, sum)
    model <- cbind(DetInc_1_3, DetInc_4_6, DetInc_7_12, DetInc_13_24)
    colnames(model) <- c("<3months", "3-<6months", "6-<12months", "12-<24months")
    rownames(model) <- 1:nMonths

    # calculate the log-likelihood
    loglikeli <- sum(obs$`<3months` * log(DetInc_1_3) - DetInc_1_3) +
      sum(obs$`3-<6months` * log(DetInc_4_6) - DetInc_4_6) +
      sum(obs$`6-<12months` * log(DetInc_7_12) - DetInc_7_12) +
      sum(obs$`12-<24months` * log(DetInc_13_24) - DetInc_13_24)

    if(is.na(modelList))
    {
      modelList <- list(model, loglikeli)
    }else{
      modelList <- c(modelList, list(model, loglikeli))
    }
  }

  # compare model output and observations
  melt_model1 <- reshape2::melt(modelList[[1]])
  colnames(melt_model1) <- c("Month", "Age", "Hospitalisations")

  if(!is.na(pars2))
  {
    melt_model2 <- reshape2::melt(modelList[[3]])
    colnames(melt_model2) <- c("Month", "Age", "Hospitalisations")
  }

  melt_obs <- reshape2::melt(as.matrix(obs))
  colnames(melt_obs) <- c("Month", "Age", "Hospitalisations")

  if(!is.na(pars2))
  {
    g1 <- ggplot2::ggplot() +
      ggplot2::geom_point(data = melt_obs, ggplot2::aes(x = Month, y = Hospitalisations, col = Age)) +
      ggplot2::geom_line(data = melt_model1, ggplot2::aes(x = Month, y = Hospitalisations, col = Age)) +
      ggplot2::geom_line(data = melt_model2, ggplot2::aes(x = Month, y = Hospitalisations), col = "black", linetype = "dashed") +
      ggplot2::theme_bw()+
      ggplot2::scale_x_continuous(name = "Year", breaks = seq(0, nrow(obs), 12), labels = years)+
      ggplot2::theme(legend.position="none")+
      ggplot2::ggtitle("RSV hospitalisations: observed (colour .), old model fit (colour -) and new model fit (black --)")+
      ggplot2::theme(plot.title = element_text(size = 10))+
      ggplot2::facet_wrap(~Age, ncol = 1)
    return(c(list("g1" = g1, "likelihood_1" = modelList[[2]], "likelihood_2" = modelList[[4]])))
  }

  g1 <- ggplot2::ggplot() +
    ggplot2::geom_point(data = melt_obs, ggplot2::aes(x = Month, y = Hospitalisations, col = Age)) +
    ggplot2::geom_line(data = melt_model1, ggplot2::aes(x = Month, y = Hospitalisations, col = Age)) +
    ggplot2::theme_bw()+
    ggplot2::theme(legend.position="none")+
    ggtitle("RSV hospitalisations: observed (colour .) and model fit (colour -)")+
    ggplot2::theme(plot.title = element_text(size = 10))+
    ggplot2::facet_wrap(~Age, ncol = 1)
  return(c(list("g1" = g1, "likelihood" = modelList[[2]])))
}

#################################################################################
#' CHANGE FOR NEW MODEL
#' Compare model predicted hospitalisations to observed hospitalisation data
#' with 95% confidence interval from mcmc derived distributions
#'
#' @param parTab table of values, with following headings:
#'               names, bestPars, lower2.5 and upper97.5
#' @param obs table of monthly hospitalisation data for 4 age groups
#' @param max_t the "run-in" time used in fitting the model
#' @return ggplot
#'
plot_fitToData95 <- function(parTab, obs, max_t = 200)
{
  nMonths <- nrow(obs)
  parSet <- c("bestPars", "lower2.5", "upper97.5")

  modelList <- NA
  for(pars in parSet)
  {
    counts <- run_rsv_model(b0 = as.numeric(parTab[which(parTab$names == "b0"), pars]),
                            b1 = as.numeric(parTab[which(parTab$names == "b1"), pars]),
                            phi = as.numeric(parTab[which(parTab$names == "phi"), pars]),
                            h_k1 = as.numeric(parTab[which(parTab$names == "h_k1"), pars]),
                            h_k2 = as.numeric(parTab[which(parTab$names == "h_k2"), pars]),
                            h_k3 = as.numeric(parTab[which(parTab$names == "h_k3"), pars]),
                            h_k4 = as.numeric(parTab[which(parTab$names == "h_k4"), pars]),
                            max_t = max_t,
                            init_conds_from_file = 0,
                            save_init_conds = 0)

    mrange <- findIndexRange(counts = counts, peakMethod = 1, month = 7,
                             ageGroup = 5, obs = obs)

    mod_DetIncidence <- as.data.frame(counts$DetIncidence[mrange$start:mrange$end,])
    DetInc_1_3 <- apply(mod_DetIncidence[, 1:3], 1, sum)
    DetInc_4_6 <- apply(mod_DetIncidence[, 4:6], 1, sum)
    DetInc_7_12 <- apply(mod_DetIncidence[, 7:12], 1, sum)
    DetInc_13_24 <- apply(mod_DetIncidence[, 13:24], 1, sum)
    model <- cbind(DetInc_1_3, DetInc_4_6, DetInc_7_12, DetInc_13_24)
    colnames(model) <- c("<3months", "3-<6months", "6-<12months", "12-<24months")
    rownames(model) <- 1:nMonths

    # calculate the log-likelihood
    loglikeli <- sum(obs$`<3months` * log(DetInc_1_3) - DetInc_1_3) +
      sum(obs$`3-<6months` * log(DetInc_4_6) - DetInc_4_6) +
      sum(obs$`6-<12months` * log(DetInc_7_12) - DetInc_7_12) +
      sum(obs$`12-<24months` * log(DetInc_13_24) - DetInc_13_24)

    if(is.na(modelList))
      modelList <- list(model, loglikeli)
    else
      modelList <- c(modelList, list(model, loglikeli))
  }

  # compare model output and observations
  melt_best <- reshape2::melt(modelList[[1]])
  colnames(melt_best) <- c("Month", "Age", "Hospitalisations")

  melt_lower <- reshape2::melt(modelList[[3]])
  colnames(melt_lower) <- c("Month", "Age", "Hospitalisations")

  melt_upper <- reshape2::melt(modelList[[5]])
  colnames(melt_upper) <- c("Month", "Age", "Hospitalisations")

  melt_model <- cbind(melt_best, melt_lower$Hospitalisations, melt_upper$Hospitalisations)
  colnames(melt_model) <- c("Month", "Age", "Best", "Lower", "Upper")

  melt_obs <- reshape2::melt(as.matrix(obs))
  colnames(melt_obs) <- c("Month", "Age", "Hospitalisations")

  g1 <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = melt_model, ggplot2::aes(x = Month, ymin = Lower, ymax = Upper), fill = "grey70") +
    ggplot2::geom_point(data = melt_obs, ggplot2::aes(x = Month, y = Hospitalisations, col = Age)) +
    ggplot2::geom_line(data = melt_model, ggplot2::aes(x = Month, y = Best), col = "black",) +
    ggplot2::theme_bw()+
    ggplot2::theme(legend.position="none")+
    ggplot2::ggtitle("RSV hospitalisations: observed (colour .), best pars model fit (black -) and 95% CI (grey shading)")+
    ggplot2::theme(plot.title = element_text(size = 10))+
    ggplot2::facet_wrap(~Age, ncol = 1)

  return(c(list("g1" = g1), modelList))
}

###############################################################################################
#' CHANGE TO YEARS
#' Plot model output - no comparison to observed
#'
#' @param outT table of monthly hospitalisation/infection data for 4 age groups
#' @param hos if T, y axis label is "Hospitalisation" else "Infections". Default is T
#' @param plotText title text. Default is no title
#'
#' @return ggplot
#'
plot_mod <- function(outT, hos = T, plotText = '')
{
  ylabel <- ifelse(hos, "Hospitalisations", "Infections")

  melt_outT <- reshape2::melt(as.matrix(outT))
  colnames(melt_outT) <- c("Month", "Age", "Incidence")

  g1 <- ggplot() +
    geom_line(data = melt_outT, aes(x = Month, y = Incidence , col = Age), linetype = "dotted", size = 1.1) +
    ggplot2::theme_bw()+
    theme(legend.position="none")+
    ylab(ylabel)+
    ggtitle(plotText)+
    theme(plot.title = element_text(size = 10))+
    facet_wrap(~Age, ncol = 1)

  return(g1)
}
##################################################################################################
#' Takes output from deSolve and pulls out hospitalisations aggregated in the
#' age groups used for fitting for term first exposure, term subsequent exposures,
#' term combined, preterm first exposure, preterm subsequent exposures,
#' preterm combined, combined first exposure, combined second exposure, combined
#'
#' @param deSolveOut output matrix
#' @param times time vector
#' @param nAges number of age groups
#' @param model either base or risk
#' @param acc if continuous ageing, need to make detInc not accumulative
#' @param old 1 = include older age group 24-<60 months,
#'            2 = include 26-<36 months, 36->48 months, 48-<60 months
#' @return list
#'
aggregate_output <- function(deSolveOut, times, nAges, model, acc = 1, old = 0)
{
  if(model == "base")
  {
    nState <- 8
    nRet <- 14
  }else{
    #model == "risk"
    nState <- 16
    nRet <- 30
  }
  deSolveOut <- array(deSolveOut[,-1], dim = c(length(times), nAges, nRet))
  deSolveOut <- deSolveOut[which(times == round(times)),,]

  if(model == "base")
  {
    listNames <- c("deSolve", "first", "second", "combined")

  }else{
    listNames <- c("deSolve", "term0", "term1", "term",
                 "pre0", "pre1", "pre",
                 "combined")
  }
  modelOut <- sapply(listNames,function(x) NULL)

  colNames <- c("<3months", "3-<6months", "6-<12months", "12-<24months")
  if(old == 1) colNames <- c(colNames, "24-<60months")
  if(old == 2) colNames <- c(colNames, "24-<36months", "36-<48months", "48-<60months")

  if(acc)
  {
    for(i in (nState+1):nRet)
      deSolveOut[, , i] <- apply(as.data.frame(deSolveOut[, , i]), 2, function(x) diff(c(0, x)))
  }
  modelOut[[listNames[1]]] <- deSolveOut

  for(i in 1:length(listNames[-1]))
  {
    detInc <- as.data.frame(deSolveOut[, , nState + i*2])
    DetInc_1_3 <- apply(detInc[, 1:3], 1, sum)
    DetInc_4_6 <- apply(detInc[, 4:6], 1, sum)
    DetInc_7_12 <- apply(detInc[, 7:12], 1, sum)
    DetInc_13_24 <- apply(detInc[, 13:24], 1, sum)
    if(old == 1) DetInc_25_60 <- apply(detInc[, 25:60], 1, sum)
    if(old == 2){
      DetInc_25_36 <- apply(detInc[, 25:36], 1, sum)
      DetInc_37_48 <- apply(detInc[, 37:48], 1, sum)
      DetInc_49_60 <- apply(detInc[, 49:60], 1, sum)
    }

    detIC <- cbind(DetInc_1_3, DetInc_4_6, DetInc_7_12, DetInc_13_24)
    if(old == 1) detIC <- cbind(detIC, DetInc_25_60)
    if(old == 2) detIC <- cbind(detIC, DetInc_25_36, DetInc_37_48, DetInc_49_60)
    colnames(detIC) <- colNames

    modelOut[[listNames[i+1]]] <- detIC
  }

  return(modelOut)
}
##################################################################################################
#' Takes output from deSolve and pulls out hospitalisations aggregated in the
#' age groups used for fitting for term first exposure, term subsequent exposures,
#' term combined, preterm first exposure, preterm subsequent exposures,
#' preterm combined, combined first exposure, combined second exposure, combined
#'
#' @param deSolveOut output matrix
#' @param times time vector
#' @param nAges number of age groups
#' @param numImm number of immunisation categories
#' @param old 1 = include older age group 24-<60 months,
#'            2 = include 26-<36 months, 36->48 months, 48-<60 months
#' @return list
#'
aggregate_output_immi <- function(deSolveOut, times, nAges, numImm, old = 0)
{browser()
  nState <- numImm * 8
  nRet <- nState + numImm * 6

  deSolveOut <- array(deSolveOut[,-1], dim = c(length(times), nAges, nRet))
  deSolveOut <- deSolveOut[which(times == round(times)),,]

  listNames <- c("deSolve", "unP_first", "unP_second", "unP_combined",
                            "M_first", "M_second", "M_combined",
                            "Prior_first", "Prior_second", "Prior_combined",
                            "all_first", "all_second", "all_combined")
  modelOut <- sapply(listNames,function(x) NULL)

  colNames <- c("<3months", "3-<6months", "6-<12months", "12-<24months")
  if(old == 1) colNames <- c(colNames, "24-<60months")
  if(old == 2) colNames <- c(colNames, "24-<36months", "36-<48months", "48-<60months")


  for(i in (nState+1):nRet)
    deSolveOut[, , i] <- apply(as.data.frame(deSolveOut[, , i]), 2, function(x) diff(c(0, x)))
  modelOut[[listNames[1]]] <- deSolveOut

  for(i in 1:length(listNames[-1]))
  {
    detInc <- as.data.frame(deSolveOut[, , nState + i*2])
    DetInc_1_3 <- apply(detInc[, 1:3], 1, sum)
    DetInc_4_6 <- apply(detInc[, 4:6], 1, sum)
    DetInc_7_12 <- apply(detInc[, 7:12], 1, sum)
    DetInc_13_24 <- apply(detInc[, 13:24], 1, sum)
    if(old == 1) DetInc_25_60 <- apply(detInc[, 25:60], 1, sum)
    if(old == 2){
      DetInc_25_36 <- apply(detInc[, 25:36], 1, sum)
      DetInc_37_48 <- apply(detInc[, 37:48], 1, sum)
      DetInc_49_60 <- apply(detInc[, 49:60], 1, sum)
    }

    detIC <- cbind(DetInc_1_3, DetInc_4_6, DetInc_7_12, DetInc_13_24)
    if(old == 1) detIC <- cbind(detIC, DetInc_25_60)
    if(old == 2) detIC <- cbind(detIC, DetInc_25_36, DetInc_37_48, DetInc_49_60)
    colnames(detIC) <- colNames

    modelOut[[listNames[i+1]]] <- detIC
  }

  return(modelOut)
}

###############################################################################################
#' Plot age-to-risk function
#'
#' @param A
#' @param B
#' @param C
#' @param D
#' @param E
#' @param maxAge age in months to plot to
#' @return ggplot
#'
plot_ageToRisk <- function(A, B, C, D, E, maxAge = 12)
{

  age <- seq(0, maxAge, by = 0.25) #months
  Y_T0 <- data.frame(Age = age, Risk = A * exp(-B * age) + C)
  Y_T1 <- data.frame(Age = age, Risk = (A * exp(-B * age) + C) * D)
  Y_PT0 <- data.frame(Age = age, Risk = (A * exp(-B * age) + C) * E)
  Y_PT1 <- data.frame(Age = age, Risk = (A * exp(-B * age) + C) * D * E)

  if(any(Y_T0$Risk > 1)) Y_T0$Risk[which(Y_T0$Risk > 1)] <- 1
  if(any(Y_T1$Risk > 1)) Y_T1$Risk[which(Y_T1$Risk > 1)] <- 1
  if(any(Y_PT0$Risk > 1)) Y_PT0$Risk[which(Y_PT0$Risk > 1)] <- 1
  if(any(Y_PT1$Risk > 1)) Y_PT1$Risk[which(Y_PT1$Risk > 1)] <- 1

  risk_df <- cbind(Y_T0, Y_T1$Risk, Y_PT0$Risk, Y_PT1$Risk)
  colnames(risk_df) <- c("Age", "T0", "T1", "PT0", "PT1")

  maxLimY <- plyr::round_any(((A + C)*E*1.05), 0.1, f = ceiling)
  if(maxLimY > 1) maxLimY <- 1.01

  g1 <- ggplot() +
    geom_line(data = Y_PT0, aes(x = Age, y = Risk, col = "red"), size = 1.1) +
    geom_line(data = Y_PT1, aes(x = Age, y = Risk, col = "red"), size = 1.1, linetype = "dashed") +
    geom_line(data = Y_T0, aes(x = Age, y = Risk, col = "blue"), size = 1.1) +
    geom_line(data = Y_T1, aes(x = Age, y = Risk, col = "blue"), size = 1.1, linetype = "dashed") +
    ggplot2::theme_bw()+
    theme(legend.position="none")+
    xlab("Age (months)")+
    ylab("Risk of hospitalisation")+
    ggtitle(label =   "Term (red), Preterm (blue), First exposure (line), Subsequent (dashed)",
            subtitle =  paste("A = ", A, " B = ", B, " C = ", C, " D = ", D, " E = ", E))+
    theme(plot.title = element_text(size = 12))+
    scale_x_continuous(limits = c(0, maxAge + 0.1), expand = c(0, 0), breaks = 0:maxAge)+
    scale_y_continuous(limits = c(0, maxLimY), expand = c(0,0),
                       breaks = seq(0, maxLimY, 0.1), labels = c(0, seq(0.1, maxLimY, 0.1)))


  return(list(g1, risk_df))
}

###############################################################################################
#' Plot fit to data for risk model (with preterms)
#'
#' @param term table of model predicted monthly hospitalisations
#'             for term infants
#' @param pTerm table of model predicted monthly hospitalisations
#'              for preterm infants
#' @param obs table of observed term and preterm population monthly
#'            hospitalisation data
#' @param years year labels for x-axis
#' @param old 0 = 4 age groups - <3months, 3-<6months, 6-<12months, 12-<24months
#'            1 = include older age group 24-<60months,
#'            2 = include older ages in yearly groups
#' @param lineCol default line colour is black
#' @param lineType default line type is dashed
#' @param lineSize default line size is 1
#' @param plotText if plotText is NA (as default), plotText is "RSV
#'                 hospitalisations observed (dots) and model (dashed)"
#' @param fileNames if fileName is NA (as default), resulting plots are not
#'                  saved. Otherwise fileName is base for 3 resulting jpeg files
#' @return ggplot
#'
plot_riskFit <- function(term, pTerm, obs, years = 2000:2013, old = 0,
                         lineCol = "black", lineType = "dashed", lineSize = 1,
                         plotText = NA, fileName = NA)
{
  browser()
  if(is.na(plotText)) plotText <- "RSV hospitalisations observed (dots) and model (dashed)"

  colNames <- c("<3months","3-<6months", "6-<12months", "12-<24months")
  if(old == 1) colNames <- c(colNames, "24-<60months")
  if(old == 2) colNames <- c(colNames, "24-36months", "36-<48months", "48-<60months")
  colnames(term) <- colNames
  colnames(pTerm) <- colNames

  melt_term <- reshape2::melt(as.matrix(term))
  colnames(melt_term) <- c("Month", "Age", "Hospitalisations")

  obs_term <- obs[, seq(1, ncol(obs), 2)]
  colnames(obs_term) <- colNames
  melt_obs_term <- reshape2::melt(as.matrix(obs_term))
  colnames(melt_obs_term) <- c("Month", "Age", "Hospitalisations")

  melt_pTerm <- reshape2::melt(as.matrix(pTerm))
  colnames(melt_pTerm) <- c("Month", "Age", "Hospitalisations")

  obs_pTerm <- obs[, seq(2, ncol(obs), 2)]
  colnames(obs_pTerm) <- colNames
  melt_obs_pTerm <- reshape2::melt(as.matrix(obs_pTerm))
  colnames(melt_obs_pTerm) <- c("Month", "Age", "Hospitalisations")

  g1 <- ggplot() +
    geom_point(data = melt_obs_term, aes(x = Month, y = Hospitalisations)) +
    geom_point(data = melt_obs_pTerm, aes(x = Month, y = Hospitalisations, col = Age)) +
    geom_line(data = melt_pTerm, aes(x = Month, y = Hospitalisations, col = Age), linetype = lineType) +
    geom_line(data = melt_term, aes(x = Month, y = Hospitalisations), col = lineCol, linetype = lineType) +
    scale_x_continuous(name = "Year", breaks = seq(0, nrow(obs), 12), labels = years)+
    ggplot2::theme_bw()+
    theme(legend.position="none")+
    ggtitle(label = plotText,
            subtitle = "Term (black) and preterm (colour)")+
    theme(plot.title = element_text(size = 12))+
    facet_wrap(~Age, ncol = 1)

  g2 <- ggplot() +
    geom_point(data = melt_obs_term, aes(x = Month, y = Hospitalisations, col = Age)) +
    geom_line(data = melt_term, aes(x = Month, y = Hospitalisations, col = Age), linetype = lineType) +
    scale_x_continuous(name = "Year", breaks = seq(0, nrow(obs), 12), labels = years)+
    ggplot2::theme_bw()+
    theme(legend.position="none")+
    ggtitle(label = plotText,
            subtitle = "Term")+
    theme(plot.title = element_text(size = 12))+
    facet_wrap(~Age, ncol = 1)

  g3 <- ggplot() +
    geom_point(data = melt_obs_pTerm, aes(x = Month, y = Hospitalisations, col = Age)) +
    geom_line(data = melt_pTerm, aes(x = Month, y = Hospitalisations, col = Age),  linetype = lineType) +
    scale_x_continuous(name = "Year", breaks = seq(0, nrow(obs), 12), labels = years)+
    ggplot2::theme_bw()+
    theme(legend.position="none")+
    ggtitle(label = plotText,
            subtitle =  "Preterm")+
    theme(plot.title = element_text(size = 12))+
    facet_wrap(~Age, ncol = 1)

  if(!is.na(fileName))
  {
    print(g1)
    ggsave(filename = paste(fileName, ".jpeg", sep = ''), dpi = 900)
    ggsave(filename = paste(fileName, ".pdf", sep = ''), dpi = 900)
    print(g2)
    ggsave(filename = paste(fileName, "_term.pdf", sep = ''), dpi = 900)
    print(g3)
    ggsave(filename = paste(fileName, "_preterm.pdf", sep = ''), dpi = 900)
  }

  return(list(g1, g2, g3))
}

############################################################################
#' Plot of 2 time series of hospitalisation or infection data for 1 age group
#'
#' @param obs1 vector of monthly hospitalisation data (colored dots)
#' @param obs2 vector of monthly hospitalisation data (dashed line)
#' @param ylabel default is Hospitalisations
#' @param years year labels for x-axis
#' @param text plot text, default is RSV hospitalisations: observed (colour .) and model (black --)
#' @return ggplot
#'
comp_series <- function(obs1, obs2, ylabel = "Hospitalisations children < 5yo", years = 2000:2013,
                      text = "")
{
  plotData <- cbind(1:length(obs1), obs1, obs2)
  colnames(plotData) <- c("Month", "obs1", "obs2")

  #Change 0 average age to NA so just plots gap in line graph
  #plotData[5, 2] <- NA

  g1 <- ggplot(as.data.frame(plotData), aes(x=Month)) +
    geom_line(aes(y = obs1), color = "darkred", size = 1.1) +
    geom_line(aes(y = obs2), color="steelblue", linetype="twodash", size = 1.1) +
    scale_x_continuous(name = "Year", breaks = seq(0, nrow(plotData), 12), labels = years)+
    ylab(ylabel)+
    ggplot2::theme_bw()+
    ggtitle(text)+
    theme(plot.title = element_text(size = 10))

  return(g1)
}

