# install package "dcmodels_1.0.01" if you have not installed it before
#install.packages("dcmodels_1.0.01.tar.gz", repos=NULL, type="source")

rm(list=ls())
library(dcmodels)
temp_hh = 456
temp_time = 18
temp_look = 2
temp_sim = 10

setwd("C:/Users/Yan Liu/Documents/UMD research/recursive probit/3dynamic multivariate D-C model - series covariances")

source("dcmodels/R/probitUtils.R")
source("dcmodels/R/models.R")
source("dcmodels/R/dynUtils.R")
source("dcmodels/R/dyn.R")
source("dcmodels/R/dyn_dc.R")
source("dcmodels/R/dynData.R")

Sys.time()
for(sim in 4:temp_sim){   
  
dir_path= paste("sim_",temp_hh,"_",temp_time,"_",temp_look,"_",sim, sep="")

spec = list(
  D = paste("simulate_data/",dir_path,"/dyn",sep=""),
  Global = paste("simulate_data/",dir_path,"/global.txt",sep=""),
  miles = paste("simulate_data/",dir_path,"/miles.txt",sep=""),
  choices = paste("simulate_data/",dir_path,"/choice.txt",sep=""),
  
  generic = list(c("type.0", "type.1", "type.2")),
  specific = list(c(), c("asc", "sex", "edu", "inc", "kid", "res"), c("asc", "sex", "edu", "inc", "kid", "res")),
  reg = c("asc", "sex", "age", "inc", "res", "gas"),
  
  nvar_altspec = 1,   #number of alternative specific variables
  modifyD = function(D,Global,t){ 
    Dyan = cbind(D[[t]], Global)
    return(Dyan)
  },
  SD = "hessian",
  ASC = FALSE, # specify this now
  nTime = temp_time - temp_look,  #actual time, not total time
  nLook = temp_look,  #look ahead period
  transition = matrix(1, 3, 3),
  stopAlt = c(), # HERE, if you have more than one alternative 
  # that halts the decision process, put that in 
  # a vector (like stopAlt = c(1,2,10) if 
  # alt 1, 2 and 10 are stopping alternatives)
  reltol = 1e-5,
  delta = 1e-1,
  method = "pmvnorm",
  SD = "hessian",    #bootstrap OR none
  nboot = 100,
  verbose = TRUE
)

spec = checkFillDynSpec(spec)
modelFns = dyn_dc
D = NULL

Sys.time()
Mdyn = model_Yan(modelFns, spec, D)
Sys.time()
Mdyn

#calculate LL at 0
spec = checkFillDynSpec(spec)
Args = dyn_dc$computeArgs(spec,NULL)
#zero = rep(0, length(dyn_dc$computeStart(spec,NULL)))
b = rep(0, 19)
for(t in 1:Args$spec$nTime){
  b = c(b, c(0, 1, 0, 0, 1))
}
LL_0 = sum(log(dyn_dc$LLVec(b,Args)))
cat("LL at 0: " ,LL_0,"\n")

R_sq = 1-Mdyn$maxLL/LL_0
ad_R_sq = 1-Mdyn$maxLL*(temp_hh*(temp_time-temp_look)-1)/LL_0/(temp_hh*(temp_time-temp_look)-length(Mdyn$results$beta_hat))
true_beta = c(0.5, -1, -0.4, -0.1, 0.3, 0.2, -0.1, -2, -0.3, -0.2, 0.5, 0.3, -0.2, 3.0, -0.3, -0.2, 0.1, -0.1, -6.0)
for(tt in 1:(temp_time - temp_look)){
  true_beta = c(true_beta, 0.05*tt, (1+0.05*tt))
}

RMSE = sqrt(sum((true_beta-Mdyn$results$beta_hat)^2)/length(true_beta))
RMSE_19 = sqrt(sum((true_beta[1:19]-Mdyn$results$beta_hat[1:19])^2)/19)

result_save = cbind(Mdyn$results, Mdyn$maxLL, LL_0, R_sq, ad_R_sq, RMSE, RMSE_19)

#save results
path_save = paste("C:/Users/Yan Liu/Documents/UMD research/recursive probit/3dynamic multivariate D-C model - series covariances/simulate_data/",dir_path, "/output_save_", sim, ".csv", sep="")
write.csv(result_save, file = path_save, row.names=FALSE, col.names=FALSE)

}
#END of model estimation 
Sys.time()




#re-organize data (please ignore the following part if you don't need it)
output_beta_hat = matrix(NA, 19+5*(temp_time-temp_look), temp_sim)
output_SD = matrix(NA, 19+5*(temp_time-temp_look), temp_sim)
output_t = matrix(NA, 19+5*(temp_time-temp_look), temp_sim)
output_p = matrix(NA, 19+5*(temp_time-temp_look), temp_sim)
output_maxLL = matrix(NA, 1, temp_sim)
output_LL0 = matrix(NA, 1, temp_sim)
output_R_sq = matrix(NA, 1, temp_sim)
output_ad_R_sq = matrix(NA, 1, temp_sim)
output_RMSE = matrix(NA, 1, temp_sim)
output_RMSE_19 = matrix(NA, 1, temp_sim)

for(jj in 1:temp_sim){
  dir_path= paste("sim_",temp_hh,"_",temp_time,"_",temp_look,"_",jj, sep="")
  dir_output = paste("C:/Users/Yan Liu/Documents/UMD research/recursive probit/3dynamic multivariate D-C model - series covariances/simulate_data/",dir_path, "/output_save_", jj, ".csv", sep="")
  data = read.csv(dir_output, head = T, sep = ",")
  output_beta_hat[,jj] = data$beta_hat
  output_SD[,jj] = data$SD
  output_t[,jj] = data$t
  output_p[,jj] = data$p
  output_maxLL[,jj] = unique(data$Mdyn.maxLL)
  output_LL0[,jj] = unique(data$LL_0)
  output_R_sq[,jj] = unique(data$R_sq)
  output_ad_R_sq[,jj] = unique(data$ad_R_sq)
  output_RMSE[,jj] = unique(data$RMSE)
  output_RMSE_19[,jj] = unique(data$RMSE_19)
}
dir_path= paste("sim_",temp_hh,"_",temp_time,"_",temp_look,"_",temp_sim, sep="")
#save beta_hat
dir_beta_hat = paste("C:/Users/Yan Liu/Documents/UMD research/recursive probit/3dynamic multivariate D-C model - series covariances/simulate_data/",dir_path, "/output_beta_hat.csv", sep="")
write.csv(output_beta_hat, file = dir_beta_hat, row.names=FALSE, col.names=FALSE)
#save SD
dir_SD = paste("C:/Users/Yan Liu/Documents/UMD research/recursive probit/3dynamic multivariate D-C model - series covariances/simulate_data/",dir_path, "/output_SD.csv", sep="")
write.csv(output_SD, file = dir_SD, row.names=FALSE, col.names=FALSE)
#save t
dir_t = paste("C:/Users/Yan Liu/Documents/UMD research/recursive probit/3dynamic multivariate D-C model - series covariances/simulate_data/",dir_path, "/output_t.csv", sep="")
write.csv(output_t, file = dir_t, row.names=FALSE, col.names=FALSE)
#save p
dir_p = paste("C:/Users/Yan Liu/Documents/UMD research/recursive probit/3dynamic multivariate D-C model - series covariances/simulate_data/",dir_path, "/output_p.csv", sep="")
write.csv(output_p, file = dir_p, row.names=FALSE, col.names=FALSE)
#save maxLL
dir_maxLL = paste("C:/Users/Yan Liu/Documents/UMD research/recursive probit/3dynamic multivariate D-C model - series covariances/simulate_data/",dir_path, "/output_maxLL.csv", sep="")
write.csv(output_maxLL, file = dir_maxLL, row.names=FALSE, col.names=FALSE)
#save SD
dir_LL0 = paste("C:/Users/Yan Liu/Documents/UMD research/recursive probit/3dynamic multivariate D-C model - series covariances/simulate_data/",dir_path, "/output_LL0.csv", sep="")
write.csv(output_LL0, file = dir_LL0, row.names=FALSE, col.names=FALSE)
#save t
dir_R_sq = paste("C:/Users/Yan Liu/Documents/UMD research/recursive probit/3dynamic multivariate D-C model - series covariances/simulate_data/",dir_path, "/output_R_sq.csv", sep="")
write.csv(output_R_sq, file = dir_R_sq, row.names=FALSE, col.names=FALSE)
#save p
dir_ad_R_sq = paste("C:/Users/Yan Liu/Documents/UMD research/recursive probit/3dynamic multivariate D-C model - series covariances/simulate_data/",dir_path, "/output_ad_R_sq.csv", sep="")
write.csv(output_ad_R_sq, file = dir_ad_R_sq, row.names=FALSE, col.names=FALSE)
#save RMSE
dir_RMSE = paste("C:/Users/Yan Liu/Documents/UMD research/recursive probit/3dynamic multivariate D-C model - series covariances/simulate_data/",dir_path, "/output_RMSE.csv", sep="")
write.csv(output_RMSE, file = dir_RMSE, row.names=FALSE, col.names=FALSE)
#save RMSE_19
dir_RMSE_19 = paste("C:/Users/Yan Liu/Documents/UMD research/recursive probit/3dynamic multivariate D-C model - series covariances/simulate_data/",dir_path, "/output_RMSE_19.csv", sep="")
write.csv(output_RMSE_19, file = dir_RMSE_19, row.names=FALSE, col.names=FALSE)