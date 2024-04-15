# rsvmod
A mathematical compartmental transmission model for RSV

We developed a dynamic compartmental model of Respiratory Syncytial Virus (RSV) fitted to individually-linked population-based laboratory, perinatal and hospitalisation data for 2000-2012 from metropolitan Western Australia (WA), stratified by age and prior exposure. We account for the differential risk of RSV-hospitalisation in full-term and preterm infants (defined as <37 weeks gestation). We formulated a function relating age, RSV exposure history, and preterm status to the risk of RSV-hospitalisation given infection. Further details of this work can be found xx.  The wrapper file, **main.r** provides examples of how to run the model and produce some of the key plots.

A basic summary of the R code files follows. The functions are described in further detail in the comments above the relevant code.

## data.r ##
This file contains descriptions of the two data files used. Note that the timeseries of monthly WA RSV-hospitalisations used in model fitting is not provided due to confidentiality issues.

## parameters.r ##
Key parameter values are set here, including the estimated fitted values of the seasonal forcing function and exponential function that relates infections (by age) to hospitalisations (the age-to-risk function).

## models.r ##
Contains the functions that define the model ODEs. These functions are used by the R package _deSolve_ to solve the ODEs. The deSolve_base model does not stratify the population by gestational age whereas deSolve_preTerm stratifies the population into preterm (born <37 weeks gestational age) and full-term. 

## functions.r ##
This file contains functions to set up the age structure used in the model, a function to create the contact matrix (using the R package _conmat_ and with resulting matrix saved as a .rda data file), functions used to set up the likelihood functions for model fitting using an mcmc process as implemented in the R package _lazymcmc_, and functions to find average age values.

## diagnostics.r ##
This file contains functions to produce a range of different plots.

## Linked publications ##

If you wish to use any part of this code, please reference

xx
