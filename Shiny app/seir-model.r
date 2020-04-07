time_stuff <- as.data.frame.from.tbl( readxl::read_excel(file_to_read$datapath, sheet = "time") )   # other parameters
            time_stuff_m <- as.data.frame.from.tbl( readxl::read_excel(file_to_read$datapath, sheet = "time2") ) # c_, cr, cq, and beta
            input_stuff  <- as.data.frame.from.tbl( readxl::read_excel(file_to_read$datapath, sheet = "input") ) # initial values
            
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