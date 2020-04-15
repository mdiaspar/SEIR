# Load packages to extend base R


required_Packages <-c("janitor","readxl","dplyr","deSolve","tidyr","ggplot2", "ggpubr", "tidyverse", "viridis", "shinycssloaders", "DT", "scales", "plotly","matrixcalc")
for(p in required_Packages){
  if(!require(p,character.only = TRUE)) {install.packages(p)}
}
load_packages <- lapply(required_Packages, require, character.only = TRUE)

# Define UI
ui <- fluidPage(
  # Add CSS here
  tags$head(
    tags$style(HTML("
      // CSS goes here
    "))
  ),
  # Application title
  titlePanel("COVID-19 model app"),
  
  # Layout
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload model parameters (Excel file)"),
      uiOutput("age_group"),
      uiOutput("compartment"),
      uiOutput("start_date"),
      width = 3
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Plot", br(), plotlyOutput("plot") %>% withSpinner(color = "#337ab7")),
        tabPanel("Summary statistics", br(), DTOutput("summary_statistics") %>% withSpinner(color = "#337ab7")),
        tabPanel("Model output", br(), DTOutput("model_output") %>% withSpinner(color = "#337ab7")),
        tabPanel("Model input (time)", br(), DTOutput("model_inputs_time1") %>% withSpinner(color = "#337ab7")),
        tabPanel("Model input (time2)", br(), DTOutput("model_inputs_time2") %>% withSpinner(color = "#337ab7")),
        tabPanel("Model input (inputs)", br(), DTOutput("model_inputs") %>% withSpinner(color = "#337ab7"))
      ),
      width = 9
    )
  )
)

# Define server logic 
server <- function(input, output) {
  
  # Build age group menu based on the number of age groups detected in the uploaded Excel file
  output$age_group <- renderUI({
    if(is.null(run_model()$nagegrp)) { 
      return() 
    } else {
      age_groups <- 1:run_model()$nagegrp
      names(age_groups) <- paste0("Age group ", age_groups)
      selectInput("age_group", "Select age group", choices = age_groups)
    }
  })
  
  # Build compartment menu
  output$compartment <- renderUI({
    if(is.null(run_model()$nagegrp)) { 
      return() 
    } else {
      checkboxGroupInput("compartment", label = "Select compartments", choices = list("Susceptible" = "S", "Latent" = "L_tot", "Infected" = "I_tot", "Recovered" = "R", "Dead" = "D"), selected = c("S", "L_tot", "I_tot", "R", "D"))
    }
  })
  
  # Build start date field
  output$start_date <- renderUI({
    if(is.null(run_model()$nagegrp)) { 
      return() 
    } else {
      dateInput("start_date", "Select start date (for Time = 0)", value="")
    }
  })
  
  # Read sheets from uploaded Excel file
  get_inputs <- reactive({
    file_to_read <- input$file
    if(is.null(file_to_read)) {
      return(list(time1 = data.frame(), time2 = data.frame(), inputs = data.frame(), columns = NULL))
    } else {
      time_stuff   <- as.data.frame.from.tbl( readxl::read_excel(file_to_read$datapath, sheet = "time") )  # other parameters
      time_stuff_m <- as.data.frame.from.tbl( readxl::read_excel(file_to_read$datapath, sheet = "time2") ) # c_, cr, cq, and beta
      input_stuff  <- as.data.frame.from.tbl( readxl::read_excel(file_to_read$datapath, sheet = "input") ) # initial values
      return(list(time1 = time_stuff, time2 = time_stuff_m, inputs = input_stuff))
    }
  })
  
  # Compute summary statistics based on the age group and compartment(s) selected with Time converted to YYYY-mm-dd if the start date field is populated
  get_statistics <- reactive({
    if(nrow(run_model()$big_out) < 1) {
      return()
    } else {
      # Function to compute min counts 
      get_min_count <- function(df, variable_of_interest) {
        df <- df %>% filter(Time < 365.25) %>% filter(!!rlang::sym(all_of(variable_of_interest)) == min(!!rlang::sym(all_of(variable_of_interest)))) %>% select(Time, all_of(variable_of_interest))
        start_date <- paste0("|", input$start_date, collapse = "")
        if(nchar(gsub("[|]", "", start_date)) == 10) {
          start_date <- as.Date(gsub("[|]", "", start_date), format = "%Y-%m-%d")
          df$Time <- as.character(start_date + df$Time)
        }
        df <- tibble(Compartment = names(df)[2], Min = as.integer(df[,2]), Day = df$Time)
      }
      
      # Function to compute max counts
      get_max_count <- function(df, variable_of_interest) {
        df <- df %>% filter(Time < 365.25) %>% filter(!!rlang::sym(variable_of_interest) == max(!!rlang::sym(variable_of_interest))) %>% select(Time, variable_of_interest)
        start_date <- paste0("|", input$start_date, collapse = "")
        if(nchar(gsub("[|]", "", start_date)) == 10) {
          start_date <- as.Date(gsub("[|]", "", start_date), format = "%Y-%m-%d")
          df$Time <- as.character(start_date + df$Time)
        }
        df <- tibble(Compartment = names(df)[2], Max = as.integer(df[,2]), Day = df$Time)
      }
      
      # Data frame with min stats
      df1 <- as.data.frame(t(sapply(paste0(input$compartment, input$age_group), function(x) get_min_count(run_model()$big_out, x))))
      for(i in 1:ncol(df1)) {
        df1[,i] <- unname(unlist(df1[,i]))
      }
      
      # Data frame with max stats
      df2 <- as.data.frame(t(sapply(paste0(input$compartment, input$age_group), function(x) get_max_count(run_model()$big_out, x))))
      for(i in 1:ncol(df2)) {
        df2[,i] <- unname(unlist(df2[,i]))
      }
      
      # Merge min and max data frames
      df <- cbind(df1, df2 %>% select(-Compartment))
      
      # Convert variable names to long form
      lookup <- tibble(short = c("S", "L_tot", "I_tot", "R", "D"), long = c("Susceptible", "Latent", "Infected", "Recovered", "Dead"))
      df$Compartment <- lookup$long[match(gsub('[0-9]+', '', df$Compartment), lookup$short)]
      
      return(list(df = df))
    }
  })
  
  # Run SEIR model by age groups
  run_model <- reactive({
    # generate bins based on input$bins from ui.R
    file_to_read <- input$file
    if(is.null(file_to_read)) {
      return(list(big_out = data.frame(), nagegrp = NULL, columns = NULL))
    } else {
      source("UtilitiesChunks.R") 
      time_stuff   <- get_inputs()$time1
      time_stuff_m <- get_inputs()$time2
      input_stuff  <- get_inputs()$inputs
      
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
            # check the symmetry condition
            if(isSymmetric(beta,check.attributes=FALSE)==FALSE){stop("beta is not a symmetric matrix")}
            if(isSymmetric(c_,check.attributes=FALSE)==FALSE){stop("c_ is not a symmetric matrix")}
            if(isSymmetric(cq,check.attributes=FALSE)==FALSE){stop("cq is not a symmetric matrix")}
            if(isSymmetric(cr,check.attributes=FALSE)==FALSE){stop("cr is not a symmetric matrix")}
            
            I_sum <- I_a + I_aqn + I_sm + I_ss + I_smisn + I_ssisn + I_ar + I_smr + I_ssr + I_smrisn + I_ssrisn + phi*I_aq
            I_ssss <- I_ssis + I_ssisn + I_ssh + I_ssrisn
            c_BetaIprop <- CrBetaIprop <- CqBetaIprop <- C_OneMinusLambda <- one <- as.matrix(rep(1,nage))
            
            OneMinusLambda      <- (one - lambda)
            c_BetaIprop         <- hadamard.prod(c_,beta) %*% I_sum          
            CqBetaIprop         <- hadamard.prod(cq,beta) %*% I_sum         
            CrBetaIprop         <- hadamard.prod(cr,beta) %*% I_sum          
            
            
            # rates of change that depends on matrices 
            #--------------------------------------------
            dS   <- -((OneMinusLambda * tau * c_BetaIprop) + (lambda * CqBetaIprop) + (OneMinusLambda * (one-tau) * CrBetaIprop)) * S #new formula
            dL   <- (OneMinusLambda * tau * c_BetaIprop) * S  - sigma * L                                                             #new formula
            dL_q <- (lambda * CqBetaIprop) * S  - sigma * L_q                                                                         #new formula
            dL_r <- (OneMinusLambda * (one-tau) * CrBetaIprop) * S  - sigma * L_r                                                     #new formula
            
            # rates of change that depends on vectors
            #--------------------------------------------
            dI_a      <- sigma * L - I_a * delta * epsilon - I_a * (one -  delta) * upsilon
            dI_aq     <- sigma * rho * L_q - I_aq * delta * epsilonq - I_aq * (one -  delta) * upsilon # updated 
            dI_ar     <- sigma * L_r - I_ar * delta * epsilon - I_ar * (one -  delta) * upsilon
            dI_aqn    <- sigma * (one -  rho) * L_q - I_aqn * delta * epsilon - I_aqn * (one -  delta) * upsilon
            dI_sm     <- (I_a + I_aqn) * delta * epsilon * alpha - kappa * I_sm 
            dI_ss     <- (I_a + I_aqn) * delta * epsilon * (one - alpha)    - kappa * I_ss # updated 
            dI_smr    <- I_ar * delta * epsilon * alpha - kappa * I_smr
            dI_ssr    <- I_ar * delta * epsilon * (one - alpha)    - kappa * I_ssr # updated 
            dI_smis   <- kappa * feim * I_sm + kappa * feimr * I_smr + delta * alpha * epsilonq * feimq * I_aq - num * I_smis
            dI_smisn  <- kappa * (one -  feim) * I_sm - num * I_smisn
            dI_ssis   <- kappa * feisi * (I_ss + I_ssr) - I_ssis * ((one -  mu) * nus + mu * nud)
            dI_ssisn  <- kappa * ((one - feisi-feish) * I_ss)  - I_ssisn * ((one -  mu) * nus + mu * nud) 
            dI_ssh    <- kappa * feish * (I_ss + I_ssr) + delta * (one-alpha) * epsilonq * I_aq - I_ssh * ((one - mu) * nus + mu * nud) # updated 
            dI_smrisn <- kappa * (one - feimr) * I_smr - num * I_smrisn
            dI_ssrisn <- kappa * (one - feisi-feish) * (I_ssr) - I_ssrisn * ((one -  mu) * nus + mu * nud)
            dI_smqisn <- I_aq * delta * alpha * epsilonq * (one - feimq) - num * I_smqisn
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
      out<-NULL
      
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
      big_out <- big_out %>% select(order(colnames(big_out))) %>% select(time, everything())
      colnames(big_out)[1] <- "Time"
      return(list(big_out = big_out, nagegrp = nagegrp))
    }
  })
  
  # Build line plot based on the age group and compartment(s) selected with Time converted to YYYY-mm-dd if the start date field is populated
  output$plot <- renderPlotly({
    # generate bins based on input$bins from ui.R
    if(nrow(run_model()$big_out) < 1 | is.null(input$age_group) | is.null(input$compartment)) {
      return()
    } else {
      big_out <- run_model()$big_out
      nagegrp <- run_model()$nagegrp
      if (nagegrp > 1){
        #variables_of_interest <- as.vector(sapply(c("S","L_tot","I_tot","R","D"), function(x) paste0(x, input$age_group)))
        variables_of_interest <- as.vector(sapply(input$compartment, function(x) paste0(x, input$age_group)))
        timelimit <- 365.25
      } else {
        #variables_of_interest <- c("S","L_tot1","I_tot1","R","D")
        variables_of_interest <- input$compartment
        variables_of_interest <- gsub("L_tot", "L_tot1", variables_of_interest)
        variables_of_interest <- gsub("I_tot", "I_tot1", variables_of_interest)
        timelimit <- 1500
      }
      big_out_graphs <- big_out %>%
        select(c("Time", all_of(variables_of_interest))) %>%
        filter(Time < timelimit) # I set a limit of days for the graphics
      
      # if needs nagegrp and lookup0
      get_plot <- function(data, age_group) {
        # Subset the data frame to include only the vectors of interest
        if (nagegrp > 1) {
          data_subset <- filter(data, meta_key %in% paste0(c("S","L_tot","I_tot","R","D"), age_group))
        } else {
          data_subset <- filter(data, meta_key %in% c("S","L_tot1","I_tot1","R","D"))
        }
        
        # Refactor the meta_key vector so that levels no longer represented in the vector are removed
        data_subset$meta_key <- factor(data_subset$meta_key)
        
        # Add labels to the factor, which will also appear in the legend
        lookup <- tibble(short = c("S", "L_tot", "I_tot", "R", "D"), long = c("Susceptible", "Latent", "Infected", "Recovered", "Dead"))
        data_subset$meta_key <- factor(data_subset$meta_key, levels = levels(data_subset$meta_key), labels = lookup$long[match(gsub('[0-9]+', '', variables_of_interest), lookup$short)])
        
        names(data_subset) <- c("Day", "Compartment", "Individuals")
        data_subset$Individuals <- as.integer(data_subset$Individuals)
        
        start_date <- paste0("|", input$start_date, collapse = "")
        if(nchar(gsub("[|]", "", start_date)) == 10) {
          start_date <- as.Date(gsub("[|]", "", start_date), format = "%Y-%m-%d")
          data_subset$Day <- start_date + data_subset$Day
          x_lab_label <- paste0("Time (since ", format(start_date, format = "%B %d, %Y"), ")")
        } else {
          x_lab_label = "Time (days)"
        }
        
        # Output the plot
        p <- ggplot(data_subset, aes(x = Day, y = Individuals)) +
          geom_line(aes(color = Compartment), size = 0.85) +
          ggtitle(paste0("SEIR model, age group ", age_group)) +
          xlab(x_lab_label) +
          ylab("Count (individuals)") +
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
      
      # Output the plots in a panel
      get_plot(big_out_long, input$age_group)
    }
  })
  
  # Render the summary statistics in searchable/sortable table
  output$summary_statistics <- renderDT(
    get_statistics()$df,
    extensions = c("Buttons", "Scroller"), 
    rownames = FALSE,
    options = list(
      columnDefs = list(list(visible = FALSE, targets = c())),
      pageLength = 10, 
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
  
  # Render uploaded Excel file (the "time" sheet) in searchable/sortable table
  output$model_inputs_time1 <- renderDT(
    get_inputs()$time1,
    extensions = c("Buttons", "Scroller"), 
    rownames = FALSE,
    options = list(
      columnDefs = list(list(visible = FALSE, targets = c())),
      pageLength = 10, 
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
  
  # Render uploaded Excel file (the "time2" sheet) in searchable/sortable table
  output$model_inputs_time2 <- renderDT(
    get_inputs()$time2,
    extensions = c("Buttons", "Scroller"), 
    rownames = FALSE,
    options = list(
      columnDefs = list(list(visible = FALSE, targets = c())),
      pageLength = 10, 
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
  
  # Render uploaded Excel file (the "inputs" sheet) in searchable/sortable table
  output$model_inputs <- renderDT(
    get_inputs()$inputs,
    extensions = c("Buttons", "Scroller"), 
    rownames = FALSE,
    options = list(
      columnDefs = list(list(visible = FALSE, targets = c())),
      pageLength = 10, 
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
  
  # Render the SEIR model output in searchable/sortable table
  output$model_output = renderDT(
    run_model()$big_out[,grepl(paste0("Time|", input$age_group, collapse = ""), names(run_model()$big_out))] %>% round(),
    extensions = c("Buttons", "Scroller"), 
    rownames = FALSE,
    options = list(
      columnDefs = list(list(visible = FALSE, targets = c())),
      pageLength = 10, 
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