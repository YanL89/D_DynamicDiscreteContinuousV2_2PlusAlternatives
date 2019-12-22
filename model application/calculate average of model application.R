rm(list=ls())
library(dcmodels)
temp_hh = 456
temp_time = 10
temp_look = 4
temp_sim = 3

source("C:/Users/Yan Liu/Documents/UMD research/recursive probit/3dynamic multivariate D-C model - series covariances/simulate_data/dynData.R")
source("C:/Users/Yan Liu/Documents/UMD research/recursive probit/3dynamic multivariate D-C model - series covariances/simulate_data/create_X.R")
source("C:/Users/Yan Liu/Documents/UMD research/recursive probit/3dynamic multivariate D-C model - series covariances/simulate_data/dynUtils.R")

#basic information---------------------------------------------------------------------
nobs_sim = temp_hh
nalt_sim = 3
nvarD_disc = 1
nvarG_disc = 7
nvarD_reg = 1
time_sim = temp_time
n_sim = 1000

apply_choice_0 = matrix(0, temp_sim, temp_time-temp_look)
apply_choice_1 = matrix(0, temp_sim, temp_time-temp_look)
apply_choice_2 = matrix(0, temp_sim, temp_time-temp_look)
apply_miles = matrix(0, temp_sim, temp_time-temp_look)
true_choice_0 = matrix(0, temp_sim, temp_time-temp_look)
true_choice_1 = matrix(0, temp_sim, temp_time-temp_look)
true_choice_2 = matrix(0, temp_sim, temp_time-temp_look)
true_miles = matrix(0, temp_sim, temp_time-temp_look)

for(sim in 1:temp_sim){
  
  dir_path= paste("sim_",temp_hh,"_",temp_time,"_",temp_look,"_", sim, sep="")
  dir_path= paste("C:/Users/Yan Liu/Documents/UMD research/recursive probit/3dynamic multivariate D-C model - series covariances/simulate_data/",dir_path, sep="")
  setwd(dir_path)
  
  #read data
  choice_apply = read.table("choice_application.txt", sep = " ", head = T)
  miles_apply = read.table("miles_application.txt", sep = " ", head = T)
  choice_true = read.table("choice.txt", sep = " ", head = T)
  miles_true = read.table("miles.txt", sep = " ", head = T)
  
  #save average value-------------------------------------------------------
  for(t in 1:(temp_time-temp_look)){
      for(i in 1:temp_hh){
        if(choice_apply[i, t]==0)
          apply_choice_0[sim, t] = apply_choice_0[sim, t] + 1/temp_hh
        else if(choice_apply[i, t]==1)
          apply_choice_1[sim, t] = apply_choice_1[sim, t] + 1/temp_hh
        else apply_choice_2[sim, t] = apply_choice_2[sim, t] + 1/temp_hh
        
        apply_miles[sim, t] = apply_miles[sim, t] + miles_apply[i, t]/temp_hh
        
        if(choice_true[i, t]==0)
          true_choice_0[sim, t] = true_choice_0[sim, t] + 1/temp_hh
        else if(choice_true[i, t]==1)
          true_choice_1[sim, t] = true_choice_1[sim, t] + 1/temp_hh
        else true_choice_2[sim, t] = true_choice_2[sim, t] + 1/temp_hh
        
        true_miles[sim, t] = true_miles[sim, t] + miles_true[i, t]/temp_hh
      }  
  }
  
}

#save choice and miles for each simulation
choice_name = paste("t", seq(1, temp_time-temp_look, 1), sep="")
sim_choice = rbind(choice_name, apply_choice_0, choice_name, true_choice_0, choice_name, choice_name, apply_choice_1, choice_name, true_choice_1, choice_name, choice_name, apply_choice_2, choice_name, true_choice_2)
sim_miles = rbind(choice_name, apply_miles, choice_name, true_miles)
write.table(sim_choice, file = "sim_application_choice_for_plot.txt", row.names=FALSE, col.names=FALSE)
write.table(sim_miles, file = "sim_application_miles_for_plot.txt", row.names=FALSE, col.names=FALSE)
