# Define packages that will be used to extend base R
package_names <- c("janitor","readxl","dplyr","deSolve","tidyr","ggplot2", "ggpubr", "tidyverse", "shiny", "shinycssloaders", "DT", "scales", "plotly", "matrixcalc") 

# Install any packages that do not exist
install_packages <- lapply(package_names, FUN = function(x) if(! require(x, character.only = TRUE)) install.packages(x))

# Load the packages
load_packages <- lapply(package_names, require, character.only = TRUE)

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
      uiOutput("other_outcome"),
      uiOutput("start_date"),
      uiOutput("max_time"),
      width = 3
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Summary", htmlOutput("plot_title"), plotlyOutput("compartment_plot") %>% withSpinner(color = "#337ab7"), htmlOutput("summary_statistics_title"), DTOutput("summary_statistics"), br()),
        tabPanel("Model output", br(), DTOutput("model_output") %>% withSpinner(color = "#337ab7")),
        tabPanel("Initial conditions", br(), DTOutput("initial_conditions") %>% withSpinner(color = "#337ab7")),
        tabPanel("Parameters by age", br(), DTOutput("parameters_by_age") %>% withSpinner(color = "#337ab7")),
        tabPanel("Parameters by age x age", br(), DTOutput("parameters_by_age_x_age") %>% withSpinner(color = "#337ab7"))
      ),
      width = 9
    )
  )
)

# Define server logic 
server <- function(input, output) {
  # Read sheets from uploaded Excel file
  get_inputs <- reactive({
    file_to_read <- input$file
    if(is.null(file_to_read)) {
      return(list(parameters_by_age = data.frame(), parameters_by_age_x_age = data.frame(), initial_conditions = data.frame(), columns = NULL))
    } else {
      time_stuff   <- as.data.frame.from.tbl(readxl::read_excel(file_to_read$datapath, sheet = run_model()$sheet_names$parms.1d))  # other parameters
      time_stuff_m <- as.data.frame.from.tbl(readxl::read_excel(file_to_read$datapath, sheet = run_model()$sheet_names$parms.2d)) # c_, cr, cq, and beta
      input_stuff  <- as.data.frame.from.tbl(readxl::read_excel(file_to_read$datapath, sheet = run_model()$sheet_names$initial.conditions)) # initial values
      return(list(parameters_by_age = time_stuff, parameters_by_age_x_age = time_stuff_m, initial_conditions = input_stuff))
    }
  })
  
  # Function to process user input from the start date field
  get_start_date <- function() {
    start_date <- paste0("|", input$start_date, collapse = "")
    if(nchar(gsub("[|]", "", start_date)) == 10) {
      start_date <- as.Date(gsub("[|]", "", start_date), format = "%Y-%m-%d")
    } else {
      start_date <- NA
    }
    return(start_date)
  }
  
  # Compute summary statistics based on the age group and compartment(s) selected with Time converted to YYYY-mm-dd if the start date field is populated
  get_statistics <- reactive({
    if(nrow(run_model()$big_out) < 1 | (is.null(input$compartment) & is.null(input$other_outcome)) | is.null(input$max_time)) {
      return()
    } else {
      # Function to compute min counts 
      get_min_count <- function(df, variable_of_interest) {
        df <- df %>% filter(Time <= input$max_time) %>% filter(!!rlang::sym(variable_of_interest) == min(!!rlang::sym(variable_of_interest))) %>% select(c(Time, all_of(variable_of_interest)))
        start_date <- get_start_date()
        if(! is.na(start_date)) {
          df$Time <- as.character(start_date + df$Time)
        } 
        df <- tibble(Description = names(df)[2], Min = as.integer(df[,2]), Day1 = df$Time)
      }
      
      # Function to compute max counts
      get_max_count <- function(df, variable_of_interest) {
        df <- df %>% filter(Time <= input$max_time) %>% filter(!!rlang::sym(variable_of_interest) == max(!!rlang::sym(variable_of_interest))) %>% select(c(Time, all_of(variable_of_interest)))
        start_date <- get_start_date()
        if(! is.na(start_date)) {
          df$Time <- as.character(start_date + df$Time)
        } 
        df <- tibble(Description = names(df)[2], Max = as.integer(df[,2]), Day2 = df$Time)
      }
      
      if (run_model()$nagegrp > 1){
        variables_of_interest <- as.vector(sapply(c(input$compartment, input$other_outcome), function(x) paste0(x, input$age_group)))
      } else {
        variables_of_interest <- c(input$compartment, input$other_outcome)
        variables_of_interest <- gsub("L_tot", "L_tot1", variables_of_interest)
        variables_of_interest <- gsub("I_tot", "I_tot1", variables_of_interest)
      }
      
      # Data frame with min stats
      df1 <- as.data.frame(t(sapply(variables_of_interest, function(x) get_min_count(run_model()$big_out, x))))
      for(i in 1:ncol(df1)) {
        df1[,i] <- unname(unlist(df1[,i]))
      }
      
      # Data frame with max stats
      df2 <- as.data.frame(t(sapply(variables_of_interest, function(x) get_max_count(run_model()$big_out, x))))
      for(i in 1:ncol(df2)) {
        df2[,i] <- unname(unlist(df2[,i]))
      }
      
      # Merge min and max data frames
      df <- cbind(Description = df1$Description, Rate = rep(NA, nrow(df1)), df1 %>% select(-Description), df2 %>% select(-Description))
      
      # Convert variable names to long form
      df$Description <- run_model()$lookup$long[match(gsub('[0-9]+', '', df$Description), run_model()$lookup$short)]
      
      # Compute attack rate
      attack_rate <- df %>% filter(Description == "Susceptible compartment")
      if(nrow(attack_rate) == 1) {
        df <- df %>% add_row(Description = "Attack rate", "Rate" = as.integer((attack_rate$Max - attack_rate$Min) / attack_rate$Max * 100), "Min" = NA, "Day1" = NA, "Max" = NA, "Day2" = NA)
      }
      
      # Order data frame alphabetically by description
      df <- df[order(df$Description),]
      
      names(df)[c(4,6)] <- rep("Day", 2)
      
      return(list(df = df))
    }
  })
  
  # Get variable names for model output tab, depending on how many total age groups are in the model
  get_variable_names <- function() {
    if (run_model()$nagegrp > 1){
      output <- run_model()$big_out[,grepl(paste0("Time|", input$age_group, collapse = ""), names(run_model()$big_out))]
    } else {
      output <- run_model()$big_out %>% select(-column_label)
    }
  }
  
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
  
  # Build line plot based on the age group and compartment(s) selected with Time converted to YYYY-mm-dd if the start date field is populated
  output$compartment_plot <- renderPlotly({
    # generate bins based on input$bins from ui.R
    if(nrow(run_model()$big_out) < 1 | is.null(input$age_group) | (is.null(input$compartment) & is.null(input$other_outcome)) | is.null(input$max_time)) {
      return()
    } else {
      big_out <- run_model()$big_out
      nagegrp <- run_model()$nagegrp
      if (nagegrp > 1){
        #variables_of_interest <- as.vector(sapply(c("S","L_tot","I_tot","R","D"), function(x) paste0(x, input$age_group)))
        variables_of_interest <- as.vector(sapply(c(input$compartment, input$other_outcome), function(x) paste0(x, input$age_group)))
        timelimit <- input$max_time
      } else {
        #variables_of_interest <- c("S","L_tot1","I_tot1","R","D")
        variables_of_interest <- c(input$compartment, input$other_outcome)
        variables_of_interest <- gsub("L_tot", "L_tot1", variables_of_interest)
        variables_of_interest <- gsub("I_tot", "I_tot1", variables_of_interest)
        timelimit <- input$max_time
      }
      big_out_graphs <- big_out %>%
        select(c(Time, all_of(variables_of_interest))) %>%
        filter(Time <= timelimit) # I set a limit of days for the graphics
      
      get_plot <- function(data, age_group) {
        # Subset the data frame to include only the vectors of interest
        if (nagegrp > 1) {
          data_subset <- filter(data, meta_key %in% paste0(c("S","L_tot","I_tot","R","D", input$other_outcome), age_group))
        } else {
          data_subset <- filter(data, meta_key %in% c("S","L_tot1","I_tot1","R","D", input$other_outcome))
        }
        
        # Refactor the meta_key vector so that levels no longer represented in the vector are removed
        data_subset$meta_key <- factor(data_subset$meta_key)
        
        # Add labels to the factor, which will also appear in the legend
        data_subset$meta_key <- factor(data_subset$meta_key, levels = levels(data_subset$meta_key), labels = gsub("compartment", "", run_model()$lookup$long[match(gsub('[0-9]+', '', variables_of_interest), run_model()$lookup$short)]))
        
        names(data_subset) <- c("Day", "Compartment", "Individuals")
        data_subset$Individuals <- as.integer(data_subset$Individuals)
        
        start_date <- get_start_date()
        if(! is.na(start_date)) {
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
  
  # Render uploaded Excel file (the "inputs" sheet) in searchable/sortable table
  output$initial_conditions <- renderDT(
    get_inputs()$initial_conditions,
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
  
  # Build other_outcome menu
  output$other_outcome <- renderUI({
    if(is.null(run_model()$nagegrp)) { 
      return() 
    } else {
      checkboxGroupInput("other_outcome", label = "Select other outcomes", choices = list("Incidence (new cases per day)" = "IncI", "Cumulative incidence" = "cumI", "Hospitalized" = "Iss_hosp", "Quarantined" = "I_aq", "Isolated" = "Isolat"), selected = c(""))
    }
  })
  
  # Build max_time menu based on the number of age groups detected in the uploaded Excel file
  output$max_time <- renderUI({
    if(is.null(run_model()$time_min) | is.null(run_model()$time_max)) { 
      return() 
    } else {
      sliderInput("max_time", "Select maximum time", min = run_model()$time_min + 1, max = run_model()$time_max, value = run_model()$time_max, step = 1)
    }
  })
  
  # Render the SEIR model output in searchable/sortable table
  output$model_output = renderDT(
    get_variable_names() %>% round(3),
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
  output$parameters_by_age <- renderDT(
    get_inputs()$parameters_by_age,
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
  output$parameters_by_age_x_age <- renderDT(
    get_inputs()$parameters_by_age_x_age,
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
  
  # Render HTML title for the plot
  output$plot_title <- renderText({
    if(is.null(run_model()$nagegrp)) { 
      return() 
    } else {
      return(paste("<br>", "<h4>Plot</h4>", "<br>"))
    }
  })
  
  
  # Build start date field
  output$start_date <- renderUI({
    if(is.null(run_model()$nagegrp)) { 
      return() 
    } else {
      dateInput("start_date", "Select start date (for Time = 0)", value = "", format = "yyyy-mm-dd")
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
  
  # Render HTML title for the summary statistics
  output$summary_statistics_title <- renderText({
    if(is.null(run_model()$nagegrp)) { 
      return() 
    } else {
      return(paste("<br>", "<br>", "<br>", "<h4>Descriptives</h4>", "<br>"))
    }
  })
  
  # Run SEIR model by age groups
  run_model <- reactive({
    # geneother_outcome bins based on input$bins from ui.R
    file_to_read <- input$file
    if(is.null(file_to_read)) {
      return(list(big_out = data.frame(), nagegrp = NULL, time_min = NULL, time_max = NULL, lookup = NULL, columns = NULL))
    } else {
      source("UtilitiesChunks.R")
      source("SEIR.n.Age.Classes.R")
      
      sheet_names <- list(initial.conditions = "Initial conditions", parms.1d = "Parameters by Age" , parms.2d = "Parameters by Age x Age", model.flow = "Model Specs (lazy)", auxiliary.vars = "Intermediate calculations")
      results <- SEIR.n.Age.Classes(file_to_read$datapath, sheet_names.lazy)
      
      # continue on below with listOut as before but should consider using results$solution
      listOut = results$solution 
      nagegrp = ncol(results$input.info$initial.conditions) - 1
      
      # Merge the data
      big_out <- bind_rows(listOut, .id = "column_label") %>% distinct(time, .keep_all= TRUE)
      xx <- yy <- df <- df2 <- NULL
      for (p in 1:nagegrp){
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
      
      parameters_by_age <- as.data.frame.from.tbl(readxl::read_excel(file_to_read$datapath, sheet = sheet_names$parms.1d))
      
      # Compute other outcomes
      for(i in 1:nagegrp) {	
        v <- c()	
        I_tot <- unname(unlist(big_out %>% select(all_of(paste0("I_tot", i)))))	
        L_tot <- unname(unlist(big_out %>% select(all_of(paste0("L_tot", i)))))	
        Lq <- unname(unlist(big_out %>% select(all_of(paste0("Lq", ifelse(nagegrp > 1, i, ""))))))	
        Iq_pres <- unname(unlist(big_out %>% select(all_of(paste0("Iq_pres", ifelse(nagegrp > 1, i, ""))))))	
        Iaq_r <- unname(unlist(big_out %>% select(all_of(paste0("Iaq_r", ifelse(nagegrp > 1, i, ""))))))	
        Ism_iso <- unname(unlist(big_out %>% select(all_of(paste0("Ism_iso", ifelse(nagegrp > 1, i, ""))))))
        Iss_isohome <- unname(unlist(big_out %>% select(all_of(paste0("Iss_isohome", ifelse(nagegrp > 1, i, ""))))))
        
        for(j in 1:nrow(big_out)) {	
          if(j == 1) {	
            v[j] <- I_tot[1]	
          } else {	
            sigma <- 1 / as.numeric(parameters_by_age %>% filter(agegrp == i & big_out$time[j] > tmin & big_out$time[j] <= tmax ) %>% select(t_latency))
            v[j] <- L_tot[j - 1] * sigma
          }	
        }	
        
        # Incidence	(number of new cases each day, over time [incidence])
        big_out[[paste0("IncI", ifelse(nagegrp > 1, i, ""))]] <- v
        
        # Cumulative incidence
        big_out[[paste0("cumI", ifelse(nagegrp > 1, i, ""))]] <- cumsum(v)
        
        # Quarantined
        big_out[[paste0("I_aq", ifelse(nagegrp > 1, i, ""))]] <- Lq + Iq_pres + Iaq_r
        
        # Isolated
        big_out[[paste0("Isolat", ifelse(nagegrp > 1, i, ""))]] <- Ism_iso + Iss_isohome
      } 
      
      # Organize model output vectors alphabetically
      big_out <- big_out %>% select(order(colnames(big_out))) %>% select(time, everything())
      
      # Rename the time vector
      colnames(big_out)[1] <- "Time"
      
      # Return the model output
      return(list(big_out = big_out, nagegrp = nagegrp, sheet_names = sheet_names, time_min = min(parameters_by_age$tmin), time_max = max(parameters_by_age$tmax), lookup = tibble(short = c("S", "L_tot", "I_tot", "R", "D", "IncI", "cumI", "Iss_hosp", "I_aq", "Isolat"), long = c("Susceptible compartment", "Latent compartment", "Infected compartment", "Recovered compartment", "Dead compartment", "Incidence", "Cumulative incidence", "Hospitalized", "Quarantined", "Isolated"))))
    }
  })
}

# Run the application 
shinyApp(ui = ui, server = server)