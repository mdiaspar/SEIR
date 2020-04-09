package_names <- c("janitor","readxl","dplyr","deSolve","tidyr","ggplot2", "ggpubr", "tidyverse", "viridis", "shinycssloaders", "DT", "scales", "plotly") 
load_packages <- lapply(package_names, require, character.only = TRUE)

# Define UI
ui <- fluidPage(
  # Application title
  tags$head(
    tags$style(HTML("
      // CSS goes here
    "))
  ),
  titlePanel("SEIR model app"),
  fluidRow(
    column(6, fileInput("file", "Please upload your model parameters (Excel file)")),
    column(6, uiOutput("age_group", width = "500px"))
  ),
  tabsetPanel(
    tabPanel("Plot", br(), plotlyOutput("plot") %>% withSpinner(color = "#337ab7")),
    tabPanel("Model output", br(), DTOutput("table") %>% withSpinner(color = "#337ab7"))
  )
)

# Define server logic 
server <- function(input, output) {
  output$age_group <- renderUI({
    if(is.null(run_model()$nagegrp)) { 
      return() 
    } else {
      age_groups <- 1:run_model()$nagegrp
      radioButtons('age_group', paste0("Please select one of the ", run_model()$nagegrp, " age groups detected."), choiceNames = paste0("Age group ", age_groups), choiceValues = age_groups, selected = character(0))
    }
  })
  
  run_model <- reactive({
    # generate bins based on input$bins from ui.R
    file_to_read <- input$file
    if(is.null(file_to_read)) {
      return(list(big_out = data.frame(), nagegrp = NULL, columns = NULL))
    } else {
      source("UtilitiesChunks.R") 
      time_stuff   <- as.data.frame.from.tbl( readxl::read_excel(file_to_read$datapath, sheet = "time") )  # other parameters
      time_stuff_m <- as.data.frame.from.tbl( readxl::read_excel(file_to_read$datapath, sheet = "time2") ) # c_, cr, cq, and beta
      input_stuff  <- as.data.frame.from.tbl( readxl::read_excel(file_to_read$datapath, sheet = "input") ) # initial values
      
      nagegrp <- length(unique(time_stuff$agegrp))        # number of age groups
      
      nrow_   <- dim(time_stuff)[1]/nagegrp
      
      time_stuff <- dplyr::arrange(time_stuff, tmin, agegrp) 
      time_stuff <- time_stuff %>%
        mutate(isim = rep(1:nrow_, each=nagegrp)) 
      
      time_stuff_m <- arrange(time_stuff_m, tmin, cagegrp, ragegrp) 
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
            
            I_sum <- I_a + I_aqn + I_sm + I_ss + I_smisn + I_ssisn + I_ar + I_smr + I_ssr + I_smrisn + I_ssrisn + phi*I_aq
            I_ssss <- I_ssis + I_ssisn + I_ssh + I_ssrisn
            BetaCrOneMinusLambda <- BetaC_OneMinusLambda <- BetaCqLambda <- C_OneMinusLambda <- CqLambda <-  CrOneMinusLambda <- CqLambda <- one <- as.matrix(rep(1,nage))
            
            OneMinusLambda       <-    (one - lambda)
            C_OneMinusLambda     <- c_ %*% OneMinusLambda          
            CrOneMinusLambda     <- cr %*% OneMinusLambda          
            CqLambda             <- cq %*% OneMinusLambda          
            BetaCrOneMinusLambda <- beta %*% CrOneMinusLambda  
            BetaC_OneMinusLambda <- beta %*% C_OneMinusLambda  
            BetaCqLambda         <- beta %*% CqLambda          
            
            # rates of change that depends on matrices
            #--------------------------------------------
            
            dS   <- -(BetaC_OneMinusLambda*tau + BetaCqLambda + BetaCrOneMinusLambda*(one - tau))*S*I_sum # updated
            dL_r <- ((BetaCrOneMinusLambda *(one-tau))*S*I_sum) - (sigma*L_r)                                       
            dL_q <- (BetaCqLambda*S*I_sum) - (((sigma*(one - rho)) + (sigma*rho))*L_q)                      
            dL   <- ((BetaC_OneMinusLambda*tau) *S*I_sum) - (sigma*L)                                          
            
            # rates of change that depends on vectors
            #--------------------------------------------
            dI_a      <- sigma*L - I_a*delta*epsilon - I_a*(one -  delta)*upsilon
            dI_aq     <- sigma*rho*L_q - I_aq*delta*epsilonq - I_aq*(one -  delta)*upsilon # updated
            dI_ar     <- sigma*L_r - I_ar*delta*epsilon - I_ar*(one -  delta)*upsilon
            dI_aqn    <- sigma*(one -  rho)*L_q - I_aqn*delta*epsilon - I_aqn*(one -  delta)*upsilon
            dI_sm     <- (I_a + I_aqn)*delta*epsilon*alpha - kappa*I_sm 
            dI_ss     <- (I_a + I_aqn)*delta*epsilon*(one - alpha)    - kappa*I_ss # updated
            dI_smr    <- I_ar*delta*epsilon*alpha - kappa*I_smr
            dI_ssr    <- I_ar*delta*epsilon*(one - alpha)    - kappa*I_ssr # updated
            dI_smis   <- kappa*feim*I_sm + kappa*feimr*I_smr + delta*alpha*epsilonq*feimq*I_aq - num*I_smis
            dI_smisn  <- kappa*(one -  feim)*I_sm - num*I_smisn
            dI_ssis   <- kappa*feisi*(I_ss + I_ssr) - I_ssis*((one -  mu)*nus + mu*nud)
            dI_ssisn  <- kappa*((one - feisi-feish)*I_ss)  - I_ssisn*((one -  mu)*nus + mu*nud) 
            dI_ssh    <- kappa*feish*(I_ss + I_ssr) + delta*(one-alpha)*epsilonq*I_aq - I_ssh*((one - mu)*nus + mu*nud) # updated
            dI_smrisn <- kappa*(one - feimr)*I_smr - num*I_smrisn
            dI_ssrisn <- kappa*(one - feisi-feish)*(I_ssr) - I_ssrisn*((one -  mu)*nus + mu*nud)
            dI_smqisn <- I_aq*delta*alpha*epsilonq*(one - feimq) - num*I_smqisn
            dR        <- (I_a + I_aq + I_aqn + I_ar)*(one - delta)*upsilon + num*(I_smis + I_smisn + I_smqisn + I_smrisn) + (one - mu)*nus*I_ssss
            dD        <- mu*nud*I_ssss
            
            out <- eval(parse(text=paste0("c(", paste0("d", paste(names.inits, collapse=", d")),")"))) # out = c(dS, dL, .... , dR, dD)
            
            out.print <- out
            names(out.print) <- c()
            names(out)<-names(vec.inits)
            list(out)
          })
          
        } #end of function calculate_derivatives
        
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
        
        if(tmin != previous.tmax)
          stop(paste(interval.label , "\n  Interval lower bound not equal to previous interval upper bound"))
        
        previous.tmax <- tmax
        out <- SEIR.n.Age.Classes( time=tt,
                                   age.parms = parameter.by.age,
                                   age.age.parms = parameter.by.age.age,
                                   list.inits = init_list)
        
        out <- as.data.frame(out)
        
        out$time <- seq(tmin,tmax,1) 
        out_for_init <- out %>%
          slice(nrow(out)) %>%
          pivot_longer(-time)
        init <- out_for_init$value
        names(init) <- out_for_init$name     
        
        rowns <- names(select(out,-c(time)))
        out <- out %>%
          mutate(N_tot = rowSums(.[rowns]))  # Total number of individuals 
        
        #updating the initial values  
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
      return(list(big_out = big_out[,-1], nagegrp = nagegrp))
    }
  })
  
  output$plot <- renderPlotly({
    # generate bins based on input$bins from ui.R
    if(nrow(run_model()$big_out) < 1 | is.null(input$age_group)) {
      return()
    } else {
      big_out <- run_model()$big_out
      nagegrp <- run_model()$nagegrp
      if (nagegrp > 1){
        variables_of_interest <- as.vector(sapply(c("S","L_tot","I_tot","R","D"), function(x) paste0(x, input$age_group)))
        timelimit <- 365.25
      }else{
        variables_of_interest <- c("S","L_tot1","I_tot1","R","D")
        timelimit <- 1500
      }
      big_out_graphs <- big_out %>%
        select(c("time", variables_of_interest)) %>%
        filter(time < timelimit) # I set a limit of days for the graphics
      
      
      # Add lookup table for age groups (if nagegrp!=5 or nagegrp=1, please write here the labels)
      #if (nagegrp==1){lookup0 <- c("all age groups")}
      #if (nagegrp==5){lookup0 <- c("< 20 year-olds", "20- to 59-year-olds", "60- to 69-year-olds", "70-79-year-olds", "80+ year-olds")}
      
      # if needs nagegrp and lookup0
      get_plot <- function(data, age_group) {
        # Subset the data frame to include only the vectors of interest
        if (nagegrp>1){
          data_subset <- filter(data, meta_key %in% paste0(c("S","L_tot","I_tot","R","D"), age_group))
        }else{data_subset <- filter(data, meta_key %in% c("S","L_tot1","I_tot1","R","D"))}
        
        # Refactor the meta_key vector so that levels no longer represented in the vector are removed
        data_subset$meta_key <- factor(data_subset$meta_key)
        
        # Add labels to the factor, which will also appear in the legend
        #data_subset$meta_key <- factor(data_subset$meta_key, levels = rev(levels(data_subset$meta_key)), labels = c("Susceptible", "Latent", "Infected", "Recovered", "Dead"))
        data_subset$meta_key <- factor(data_subset$meta_key, levels = levels(data_subset$meta_key), labels = c("Susceptible", "Latent", "Infected", "Recovered", "Dead"))
        
        names(data_subset) <- c("Day", "Compartment", "Count")
        data_subset$Count <- as.integer(data_subset$Count)
        
        todays_date <- as.Date("03/13/20", "%m/%d/%y")
        data_subset$Date <- todays_date + data_subset$Day
        #data_subset$Date <- format(data_subset$Date, "%b %d")
        
        # Output the plot
          p <- ggplot(data_subset, aes(x = Day, y = Count)) +
          geom_line(aes(color = Compartment), size = 0.85) +
          ggtitle(paste0("SEIR model, age group ", age_group)) +
          xlab("Time (days)") +
          ylab("N (individuals)") +
          scale_y_continuous(labels = comma) +
          theme_minimal() +
          theme(
            plot.title = element_text(size = 12),
            axis.title.x = element_text(size = 12),
            axis.title.y = element_text(size = 12),
            legend.text = element_text(size = 12),
            legend.title = element_blank()
          )
        p <- ggplotly(p)
      }
      
      # Reshape SEIR model output from wide to long format
      big_out_long <- gather(big_out_graphs, key = meta_key, value = meta_value, 2:ncol(big_out_graphs), factor_key = TRUE)
      # Create a time series plot for each age group and add it to a list object
      #plots <- lapply(1:nagegrp, FUN = get_plot, data = big_out_long)
      
      # Output the plots in a panel
      #print(ggarrange(plotlist = ggplotly(plots), ncol = 2, nrow = ceiling(length(plots)/ 2), common.legend = TRUE))
      get_plot(big_out_long, input$age_group)
    }
  })
  
  output$table = renderDT(
    run_model()$big_out %>% round(),
    extensions = c("Buttons", "Scroller"), 
    rownames = FALSE,
    options = list(
      columnDefs = list(list(visible = FALSE, targets = c())),
      pageLength = 50, 
      dom = "Bfrtip", 
      buttons = c("colvis", "copy", "csv", "excel", "pdf"), 
      deferRender = TRUE, 
      searchDelay = 500,
      initComplete = JS(
        "function(settings, json) {",
        "$(this.api().table().header()).css({'background-color': '#fff', 'color': '#111'});",
        "}"
      )
    )
  )
}

# Run the application 
shinyApp(ui = ui, server = server)