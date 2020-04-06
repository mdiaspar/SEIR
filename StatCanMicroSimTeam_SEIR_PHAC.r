rm(list = ls())
#################################################################################
#
# Modified PHAC S(E)IR model script with n age classes to estimate the number of
# COVID-19 cases (infected, hospitalized, etc).
#
#==================================================================================
#
#
# Author of the modifications: Maikol Diasparra, Claude Nadeau and Joel Barnes
# Date: 25 MAR 2020
# Modification objective: to extend the model to get the estimates by age-groups
# Original script provided by PHAC Analysts: Antoinette Ludwig & Erin Rees.
#
#
# Ex.The input file contains 5 age classes (age groups: 1=under 20, 2=20-59, 3=60-69, 4=70-79, 5=80+)
# But the script can be used for n age classes.
#
#
# We assume that the population is closed
# (no immigration, births or natural deaths are considered, although we include
# disease induced death) and that there is only one disease in operation.
#
#     
# 
#################################################################################

print("SPECIAL NOTE: c is a primitive function in R, so it was replace by c_")

#==========================================================================
#  Define the parameter values
#==========================================================================
# agegrp = age group
# cagegrp = age group by row (to write beta)
# ragegrp = age group by column (to write beta)
# c_ = contact rate (age-groups x age-groups).          
# cr = contact rate for refractory people (age-groups x age-groups). Number contacts per day individuals in age-group A make with ndividuals in age-group B
# cq =  contact rate for quarantined people (age-groups x age-groups)
# beta = transmission probability  (age-groups x age-groups)
# sigma = 1 / latency period
# lambda = proportion of exposed (incubant) individuals are detected and placed in quarantine (contact tracing) (by age-groups)
# rho = quarantine complience rate (by age-groups)
# epsilon = 1/pre-symptomatic infectious period
# epsilonq = 1/(pre-symptomatic infectious period +duration between onset of symptoms and diagnostic)
# alpha = percentage of infectious (symptomatic) that develop mild symptoms (by age-groups)
# delta = percentage of infectious pre-symptomatic who will develop symptoms (by age-groups)
# upsilon = 1/duration of the asymptomatic period between pre-symptomatic and recovery
# num = 1/duration of symptomatic period for mild cases before recover
# nus = 1/duration of symptomatic period for severe cases before recover
# nud = 1/duration of symptomatic period for severe cases before dying
# kappa =  1/duration between onset of symptoms and diagnostic
# feim = percentage of mild cases who go in isolation (by age-groups)
# feisi =  percentage of severe cases who go in isolation (by age-groups)
# feish = percentage of severe cases who go in hospital (by age-groups)
# feimq = percentage of mild cases in quarantine who go in isolation (by age-groups)
# feimr = percentage of mild cases refractory who go in isolation (by age-groups)
# mu = percentage of severe cases dying (by age-groups)
# tau = Tau = complience rate with social distancing measures (by age-groups)
# phi = modulator (percentage) for the effect of quarantine infetious on transmission

#==========================================================================
#  Define the components
#==========================================================================

# S = susceptible
# L = latent(incubant)
# L_q = latent (incubant) in quarantine
# L_r = latent (icubant) refractory to social distancing
# I_a = pre-symptomatic infectious
# I_ar = pre-symptomatic infectious refractory to social distancing
# I_aq = pre-symptomatic infectious in quarantine and complient
# I_aqn = pre-symptomatic infectious in quarantine and not complient
# I_sm = symptomtic with mild symptoms before diagnostic
# I_ss = symptomatic with severe symptoms before diagnostic
# I_smr = initial symptomatic until testing with mild symptoms refractory to social distancing
# I_ssr = initial symptomatic until testing with severe symptoms refractory to social distancing
# I_smis = symptomatic with mild symptoms after diagnostic in isolation
# I_smisn = symptomatic with mild symptoms after diagnostic not in isolation
# I_ssis = symptomatic with severe symptoms after diagnostic in isolation
# I_ssisn = symptomatic with severe symptoms after diagnostic not in isolation
# I_ssh =  symptomatic with severe symptoms after diagnostc hospitalized
# I_smrisn = initial symptomatic not isolated with mild symptoms refractory to social distancing
# I_ssrisn = initial symptomatic not isolated with severe symptoms refractory to social distancing
# I_smqisn = symptomatic in quarantine with mild symptoms who don't continue in isolation
#
# R = recovered
# D= Dead from infection

#==========================================================================
#  packages
#==========================================================================


package_names <- c("janitor","readxl","dplyr","deSolve","tidyr","ggplot2", "ggpubr", "tidyverse", "viridis") #"matrixStats","pryr"
load_packages <- lapply(package_names, require, character.only = TRUE)

#==========================================================================
#   Data
#==========================================================================


### User input parameters
WDir <- "C:/Users/maiko/Downloads/SEIR_Model/"    # working directory 
setwd(WDir)

source("UtilitiesChunks.R") 

#file_name <- "input_sheet_5agegrp.xls"
file_name <- "input_sheet_1agegrp.xls"

time_stuff   <- as.data.frame.from.tbl( readxl::read_excel(paste(WDir,file_name,sep=""),sheet = "time") )   # other parameters
time_stuff_m <- as.data.frame.from.tbl( readxl::read_excel(paste(WDir,file_name,sep=""), sheet = "time2") ) # c_, cr, cq, and beta
input_stuff  <- as.data.frame.from.tbl( readxl::read_excel(paste(WDir,file_name,sep=""), sheet = "input") ) # initial values

nagegrp <- length(unique(time_stuff$agegrp))        # number of age groups

nrow_   <- dim(time_stuff)[1]/nagegrp

time_stuff <- dplyr::arrange(time_stuff, tmin, agegrp) # sort by tmin agegrp
time_stuff <- time_stuff %>%
  mutate(isim = rep(1:nrow_, each=nagegrp)) 

time_stuff_m <- arrange(time_stuff_m, tmin, cagegrp, ragegrp) # sort by tmin cagegrp  ragegrp
time_stuff_m <- time_stuff_m %>%
  mutate(isim = rep(1:nrow_, each = nagegrp*nagegrp))
#===================================================================
# Initial values (components)
#===================================================================
#initial values
input_stuff_age_columns = setdiff(colnames(input_stuff), "NAME")
init_list <- list()
for(k in input_stuff$NAME){
  init_list[[k]] <- as.matrix( subset(input_stuff, NAME == k)[,input_stuff_age_columns] )
}
#==========================================================================
#  Main routine
#==========================================================================
# The SEIR model with N age classes
#
SEIR.n.Age.Classes <- function( time=NULL, age.parms = NULL, age.age.parms = NULL,list.inits = NULL, not.parms=  c("tmin", "tmax", "agegrp", "cagegrp", "ragegrp", "isim"))
{
  nage = nrow(age.parms)
  
  if (is.null(age.parms))
    stop("undefined 'age.parms'")
  
  if (is.null(age.age.parms))
    stop("undefined 'age.age.parms'")
  
  if (is.null(time))
    stop("undefined 'time'")
  
  if (is.null(list.inits))
    stop("undefined 'list.inits'")
  
  
  list.parms <- list()
  for(k in setdiff(names(age.parms), not.parms )){
    list.parms[[k]] <- age.parms[,k]
  }
  for(k in setdiff(names(age.age.parms), not.parms ))
  {
    temp<- array(NA, c(nage,nage))
    temp[cbind(age.age.parms$cagegrp, age.age.parms$ragegrp)] <- age.age.parms[,k]
    list.parms[[k]] <- temp
    if(any(is.na(temp)))
      stop(paste0(k," matrix has some missing entries"))
  }
  
  ### to write the system of differential equations
  calculate_derivatives <- function(time, vec.inits, list.parms, names.inits) {
    
    deriv.char <- paste0("d", paste(names.inits, collapse =", d")) 
    iota <- seq(length(vec.inits)/length(names.inits)) 
    list.inits <- list()
    for(k in names.inits){
      if (length(iota) > 1) {
        list.inits[[k]] <- vec.inits[paste0(k, iota)] 
      }else{list.inits[[k]] <- vec.inits[k] }
    }
    
    with(as.list(c(list.inits, list.parms)),{
      
      I_sum <- I_a + I_aqn + I_sm + I_ss + I_smisn + I_ssisn + I_ar + I_smr + I_ssr + I_smrisn + I_ssrisn + phi * I_aq
      BetaCrOneMinusLambda <- BetaC_OneMinusLambda <- BetaCqLambda <- C_OneMinusLambda <- CqLambda <-  CrOneMinusLambda <- CqLambda <- one <- as.matrix(rep(1,nage))
      OneMinusLambda <- (one - lambda)
      C_OneMinusLambda     <- c_ %*% OneMinusLambda          
      CrOneMinusLambda     <- cr %*% OneMinusLambda          
      CqLambda             <- cq %*% OneMinusLambda          
      BetaCrOneMinusLambda <- beta %*% CrOneMinusLambda  
      BetaC_OneMinusLambda <- beta %*% C_OneMinusLambda  
      BetaCqLambda         <- beta %*% CqLambda          
      
      
      # rates of change that depends on matrices
      #--------------------------------------------
      dS   <- -(BetaCrOneMinusLambda + BetaC_OneMinusLambda + BetaCqLambda) * S * I_sum
      dL   <- BetaC_OneMinusLambda * tau * S * I_sum - sigma * L
      dL_q <- BetaCqLambda* S * I_sum - sigma * L_q
      dL_r <- BetaCrOneMinusLambda * (one -  tau) * S * I_sum - sigma * L_r
      
      # rates of change that depends on vectors
      #--------------------------------------------
      dI_a      <- sigma * L - I_a * delta * epsilon - I_a * (one -  delta) * upsilon
      dI_aq     <- sigma * rho * L_q - I_aq * delta * epsilon - I_aq * (one -  delta) * upsilon
      dI_ar     <- sigma * L_r - I_ar * delta * epsilon - I_ar * (one -  delta) * upsilon
      dI_aqn    <- sigma * (one -  rho) * L_q - I_aqn * delta * epsilon - I_aqn * (one -  delta) * upsilon
      dI_sm     <- (I_a + I_aqn) * delta * epsilon * alpha - kappa * I_sm # updated
      dI_ss     <- (I_a + I_aqn) * delta * epsilon * OneMinusLambda - kappa * I_ss # updated
      dI_smr    <- I_ar * delta * epsilon * alpha - kappa * I_smr
      dI_ssr    <- I_ar * delta * epsilon * OneMinusLambda - kappa * I_ssr
      dI_smis   <- kappa * feim * I_sm + kappa * feimr * I_smr + delta * alpha * epsilonq * feimq * I_aq - num * I_smis
      dI_smisn  <- kappa * (one -  feim) * I_sm - num * I_smisn
      dI_ssis   <- kappa * feisi * (I_ss + I_ssr) - I_ssis * ((one -  mu) * nus + mu * nud)
      dI_ssisn  <- kappa * ((one - feisi-feish) * I_ss)  - I_ssisn * ((one -  mu) * nus + mu * nud) # updated
      dI_ssh    <- kappa * feish * (I_ss + I_ssr) + delta * OneMinusLambda * epsilonq * I_aq - I_ssh * ((one - mu) * nus + mu * nud)
      dI_smrisn <- kappa * (one - feimr) * I_smr - num * I_smrisn
      dI_ssrisn <- kappa * (one - feisi-feish) * (I_ssr) - I_ssrisn * ((one -  mu) * nus + mu * nud)
      dI_smqisn <- I_aq * delta * alpha * epsilonq * (one - feimq) - num * I_smqisn
      dR        <- (I_a + I_aq + I_aqn + I_ar) * (one - delta) * upsilon + num * (I_smis + I_smisn + I_smqisn + I_smrisn) + (I_ssis + I_ssisn + I_ssh + I_ssrisn) * (one - mu) * nus
      dD        <- mu * nud * (I_ssis + I_ssisn + I_ssh + I_ssrisn)
      
      out <- eval(parse(text=paste0("c(", paste0("d", paste(names.inits, collapse=", d")),")"))) # out = c(dS, dL, .... , dR, dD)
      
      out.print <- out
      names(out.print) <- c()
      #print(out.print)
      names(out)<-names(vec.inits)
      list(out)
    })
    
  } #end of function calculate_derivatives
  
  ###---------------------------------------------------------------------------------
  ### Solver for Ordinary Differential Equations (ODE), Switching Automatically
  ### Between Stiff and Non-stiff Methods
  ###---------------------------------------------------------------------------------
  # let's determine the values of the components at time (t)
  
  output <-  lsoda(y = unlist(list.inits),
                   times = time,
                   func =  calculate_derivatives,
                   parms = list.parms,
                   names.inits = names(list.inits))
  
  return(output)
}  # END  of function SEIR.n.Age.Classes
################################################################
#        To run the example
################################################################
# the current example includes 5 age-groups (Ex. 1=under 20 years, 2=20-59, 3=60-69, 4=70-79, 5=80 years or more)
excluded_names <- c("tmin", "tmax","agegrp","cagegrp","ragegrp","isim")
sprintf("S(E)IR model script to estimate the number of COVID-19 cases")
sprintf("Number of age groups considered: %s", nagegrp)
sprintf("Components:")
sprintf( names(init_list) )
sprintf("Parameters that change with age (age-groups):")
sprintf( setdiff(colnames(time_stuff  ), excluded_names) )
sprintf("Parameters that change with age and contact with others (age-groups x age-groups):")
sprintf( setdiff(colnames(time_stuff_m), excluded_names) )
sprintf("...Computing ... ")

nSim <- max(time_stuff$isim)
listOut <- list()
previous.tmax <- 0

for(i in seq(1, nSim, 1)){
  parameter.by.age     <- subset(time_stuff  , isim == i)
  parameter.by.age.age <- subset(time_stuff_m, isim == i)
  
  tmin <- unique(c(parameter.by.age$tmin, parameter.by.age.age$tmin))
  tmax <- unique(c(parameter.by.age$tmax, parameter.by.age.age$tmax)) 
  
  if(length(tmin)>1 || length(tmax)>1 || tmin>=tmax )
    stop(paste0("Unexpected pattern in tmin, tmax for interval ", i))
  
  tt <- seq(0, tmax - tmin, by = 1)
  #tt <- seq(tmin, tmax, by = 1)
  
  if(tmin != previous.tmax)
    stop(paste(interval.label , "\n  Interval lower bound not equal to previous interval upper bound"))
  
  previous.tmax <- tmax
  out <- SEIR.n.Age.Classes( time=tt,
                             age.parms = parameter.by.age,
                             age.age.parms = parameter.by.age.age,
                             list.inits = init_list)
  
  out <- as.data.frame(out)
  # ode/lsoda Output diagnostic #######################
  #diagn <- diagnostics.deSolve(out)
  
  out$time <- seq(tmin,tmax,1) #tt
  out_for_init <- out %>%
    slice(nrow(out)) %>%
    pivot_longer(-time)
  init <- out_for_init$value
  names(init) <- out_for_init$name      #new initial values
  
  rowns <- names(select(out,-c(time)))
  out <- out %>%
    mutate(N_tot = rowSums(.[rowns]))  # Total number of individuals 
  
  #updating the initial values  
  #====================================================================
  for(k in 1:length(init_list)){
    init_list[[k]][1:nagegrp] <- init[seq(nagegrp*(k-1)+1,nagegrp*k)] 
  }
  
  # Add outputs to the list
  listOut[[i]] <- out
}

# Merge the data
big_out <- bind_rows(listOut, .id = "column_label") %>% distinct(time, .keep_all= TRUE)
xx <- yy <- df <- df2 <- NULL
for (p in 1: nagegrp){
  if (nagegrp>1){varsc<-names(big_out)[grepl(p,names(big_out))]}else{varsc<-names(big_out)}
  df <-big_out %>% 
    select(one_of(varsc))
  xx <- df %>%
    select_at(vars(starts_with("L"))) %>% 
    rowSums()
  yy <-df %>% 
    select_at(vars(starts_with("I"))) %>% 
    rowSums()
  df2 <- cbind(xx,yy)
  big_out <- cbind(big_out,df2)
  names(big_out)[c(dim(big_out)[2]-1,dim(big_out)[2])]<-c(paste0(c("L_tot","I_tot"),p))
  varsc<-names(big_out)[grepl(p,names(big_out))]
}

write.csv(big_out,paste0(WDir,paste0("SEIR_results_",paste0(nagegrp,"agegrp.csv"))), row.names = FALSE)

#==========================================================================
#  Graphics 
#==========================================================================
# Variables of interest for plots
if (nagegrp > 1){
  variables_of_interest <- as.vector(sapply(c("S","L_tot","I_tot","R","D"), function(x) paste0(x, 1:nagegrp)))
  timelimit <- 365.25
}else{
  variables_of_interest <- c("S","L_tot1","I_tot1","R","D")
  timelimit <- 1500
}
big_out_graphs <- big_out %>%
  select(c("time", variables_of_interest)) %>%
  filter(time < timelimit) # I set a limit of days for the graphics


# Add lookup table for age groups (if nagegrp!=5 or nagegrp=1, please write here the labels)
if (nagegrp==1){lookup0 <- c("all age groups")}
if (nagegrp==5){lookup0 <- c("< 20 year-olds", "20- to 59-year-olds", "60- to 69-year-olds", "70-79-year-olds", "80+ year-olds")}

# if needs nagegrp and lookup0
get_plot <- function(data, age_group,lookup=lookup0) {
  # Subset the data frame to include only the vectors of interest
  if (nagegrp>1){
    data_subset <- filter(data, meta_key %in% paste0(c("S","L_tot","I_tot","R","D"), age_group))
  }else{data_subset <- filter(data, meta_key %in% c("S","L_tot1","I_tot1","R","D"))}
  
  # Refactor the meta_key vector so that levels no longer represented in the vector are removed
  data_subset$meta_key <- factor(data_subset$meta_key)
  
  # Add labels to the factor, which will also appear in the legend
  #data_subset$meta_key <- factor(data_subset$meta_key, levels = rev(levels(data_subset$meta_key)), labels = c("Susceptible", "Latent", "Infected", "Recovered", "Dead"))
  data_subset$meta_key <- factor(data_subset$meta_key, levels = levels(data_subset$meta_key), labels = c("Susceptible", "Latent", "Infected", "Recovered", "Dead"))
    # Output the plot
  ggplot(data_subset, aes(x = time, y = meta_value)) + 
    geom_line(aes(color = meta_key), size = 0.55) +
    ggtitle(paste0("SEIR model, age group ", lookup[age_group])) +
    xlab("Time (days)") +
    ylab("N (individuals)") +
    scale_color_viridis(discrete = TRUE) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      legend.title = element_blank()
    )
}

# Reshape SEIR model output from wide to long format
big_out_long <- gather(big_out_graphs, key = meta_key, value = meta_value, 2:ncol(big_out_graphs), factor_key = TRUE)
# Create a time series plot for each age group and add it to a list object
plots <- lapply(1:nagegrp, FUN = get_plot, data = big_out_long)
# Output the plots in a panel
ggarrange(plotlist = plots, ncol = 2, nrow = ceiling(length(plots)/ 2), common.legend = TRUE)

#################################################################################
#      
#     ^^
#    (oo)     WORK IN PROGRESS......
#   _(  )_                          
#     ""
#
#################################################################################
