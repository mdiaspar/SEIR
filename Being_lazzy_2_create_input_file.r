#===========================================
# being lazzy to do it manually
#===========================================

#install.packages("xlsx")
#library(xlsx)
library(readxl)
library(tidyr)
library(dplyr)
WDir<-"C:/Users/maiko/Downloads/seir_model-master/seir_model-master/" # working directory 
setwd(WDir)

file_name<-"input_sheet.xlsx" #original PHAC file
nagegrp<-5 # user-defined 


df.template <- read_excel(paste(WDir,file_name,sep=""), sheet = "time")

# parameters by age-groups 
excluded_vars <- c("c_", "cr","cq","beta")
time_stuff <- df.template %>% 
  select(-one_of(excluded_vars)) %>% 
  mutate(agegrp = 1)

for (i in 2:nagegrp){
  time_stuff_B<- df.template%>% 
    select(-one_of(excluded_vars)) %>% 
    mutate(agegrp = i) 
  
  time_stuff <- bind_rows(list(time_stuff,time_stuff_B))
}

# parameters more complex (for age-group x age-group matrices): 
# parameters by age-groups 
exclusive_vars <- c("tmin","tmax","c_", "cr","cq","beta")


for (i in 1:nagegrp){
  for (j in 1:nagegrp){
    if(j==1 & i==1){
      time_stuff_m <- df.template %>% 
        select(one_of(exclusive_vars)) %>% 
        mutate(cagegrp = 1) %>% 
        mutate(ragegrp = 1)
    }else{
      time_stuff_B<- df.template%>% 
        select(one_of(exclusive_vars)) %>% 
        mutate(cagegrp = i) %>% 
        mutate(ragegrp = j)
      
      time_stuff_m <- bind_rows(list(time_stuff_m,time_stuff_B))
    }
  }
}




new_file<-"input_sheet_by_agegrp.xlsx"
#write.xlsx(time_stuff, paste(WDir,file_name,sep=""), sheetName = "time", 
#           col.names = TRUE, row.names = TRUE, append = FALSE)
#write.xlsx(time_stuff_m, paste(WDir,file_name,sep=""), sheetName = "time2", 
#          col.names = TRUE, row.names = TRUE, append = FALSE)

# a little issue with java :-(
library(writexl)
#writexl::write_xlsx(time_stuff_m, new_file,col_names=TRUE)
writexl::write_xlsx(time_stuff, new_file,col_names=TRUE)