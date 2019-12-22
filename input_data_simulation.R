  rm(list=ls())
  library(dcmodels)
  library("truncnorm")
  temp_hh = 456
  temp_time = 18
  temp_look = 2
  temp_sim = 10
  
  source("C:/Users/Yan Liu/Documents/UMD research/recursive probit/3dynamic multivariate D-C model - series covariances/simulate_data/dynData.R")
  source("C:/Users/Yan Liu/Documents/UMD research/recursive probit/3dynamic multivariate D-C model - series covariances/simulate_data/create_X.R")
  source("C:/Users/Yan Liu/Documents/UMD research/recursive probit/3dynamic multivariate D-C model - series covariances/simulate_data/dynUtils.R")
  
  for(sim in 4:temp_sim){
    
    dir_path= paste("sim_",temp_hh,"_",temp_time,"_",temp_look,"_", sim, sep="")
    dir_path= paste("C:/Users/Yan Liu/Documents/UMD research/recursive probit/3dynamic multivariate D-C model - series covariances/simulate_data/",dir_path, sep="")
    dir.create(dir_path)
    setwd(dir_path)
  
  #model description
  #U1 = type1*beta_typ
  #U2 = type2*beta_typ + asc*beta_asc2 + gender2*beta_sex2 + education2*beta_edu2 
  #     + income2*beta_inc2 + kids2*beta_kid2 + residential_density2*beta_res2
  #U3 = type3*beta_typ + asc*beta_asc3 + gender3*beta_sex3 + education3*beta_edu3 
  #     + income3*beta_inc3 + kids3*beta_kid3 + residential_density3*beta_res3
  #Y = asc*beta_asc4 + gender4*beta_sex4 + age4*beta_age4 + income4*beta_inc4 
  #     + residential_density4*beta_res4 + driving_cost4*beta_gas4 
  
  #data generation--------------------------------------------------------------------- 
  nobs_sim = temp_hh   #number of households 456 
  nalt_sim = 3    #number of alternatives
  nvarD_disc = 1   #number of dynamic variables in discrete part (typ)
  nvarG_disc = 7   #number of global variables (asc, sex, edu, kid, inc, age, res)
  nvarD_reg = 1   #number of dynamic variables in the regression (gas)
  time_sim = temp_time   #number of study time periods 
  n_sim = 1000   #number of simulations to calculate decision variables  
  
  #simulate "global" data
  asc_sim = rep(1, nobs_sim)
  sex_sim = rbinom (nobs_sim, 1, 0.5)
  edu_sim = rep(0, nobs_sim)
  kid_sim = rep(0, nobs_sim)
  inc_sim = rep(0, nobs_sim)
  age_sim = rep(0, nobs_sim)
  for(i in 1:nobs_sim){
    edu_sim[i] = which(rmultinom(1, 1, c(0.0130, 0.1771, 0.2354, 0.1166, 0.2549, 0.2030))==1)
    kid_sim[i] = which(rmultinom(1, 1, c(0.7149, 0.1620, 0.0929, 0.0238, 0.0022, 0.0042))==1)
    inc_sim[i] = which(rmultinom(1, 1, c(0.1123, 0.2376, 0.2441, 0.1533, 0.1620, 0.0670, 0.0237))==1)
    age_sim[i] = which(rmultinom(1, 1, c(0.0994, 0.2225, 0.1447, 0.2030, 0.3283, 0.0021))==1)
  }
  res_sim = rtruncnorm (nobs_sim, 0.04, 17.50, 2.49, 2.75)
  
  global = matrix(data = NA, nrow = nobs_sim, ncol = nvarG_disc, byrow = FALSE)
  global[, 1] = asc_sim
  global[, 2] = sex_sim
  global[, 3] = edu_sim
  global[, 4] = kid_sim
  global[, 5] = inc_sim
  global[, 6] = age_sim
  global[, 7] = res_sim
  global = rbind(c("asc", "sex", "edu", "kid", "inc", "age", "res"), global)
  write.table(global, file = "global.txt", row.names=FALSE, col.names=FALSE)
#global = read.table("global.txt", sep = " ", head = T)

#simulate "dyn" data
dyn = list()
for(t in 1:time_sim){
  dyn[[t]] = list()
  dyn[[t]] = matrix(data = NA, nrow = nobs_sim, ncol = nvarD_disc * nalt_sim + nvarD_reg, byrow = FALSE)
  typ1_sim = runif(nobs_sim, min = 0 + 0.02 * t, max = 2 + 0.02 * t)
  typ2_sim = runif(nobs_sim, min = 1 + 0.02 * t, max = 3 + 0.02 * t)
  typ3_sim = runif(nobs_sim, min = 2 + 0.02 * t, max = 4 + 0.02 * t)
  gas = rtruncnorm (nobs_sim, 0.03, 0.43, 0.13 + 0.005 * t, 0.06)
  dyn[[t]] = cbind(typ1_sim, typ2_sim, typ3_sim, gas)
  dyn[[t]] = rbind(c("type.0", "type.1", "type.2", "gas"), dyn[[t]])
  data_name = paste("dyn", t, sep = "")
  write.table(dyn[[t]], file = paste(data_name, "txt", sep = "."), row.names=FALSE, col.names=FALSE)
  #dyn[[t]] = read.table(paste(data_name, "txt", sep = "."), sep = " ", head = T)
}
  
#define true coefficients - beta-------------------------------------------------------
#beta = c(beta_typ, beta_asc2, beta_sex2, beta_edu2, beta_inc2, beta_kid2, beta_res2, 
#  beta_asc3, beta_sex3, beta_edu3, beta_inc3, beta_kid3, beta_res3, 
#  beta_asc4, beta_sex4, beta_age4, beta_inc4, beta_res4, beta_gas4) 
beta = c(0.5, -1, -0.4, -0.1, 0.3, 0.2, -0.1, -2, -0.3, -0.2, 0.5, 0.3, -0.2, 3.0, -0.3, -0.2, 0.1, -0.1, -6.0)
beta_disc = beta[1:13]
beta_reg = beta[14:19]

#define D-C covariance matrix
S_disc = diag(nalt_sim)
M = matrix(c(-1, -1, 1, 0, 0 ,1), nalt_sim-1, nalt_sim)
Diff_S_disc = M %*% S_disc %*% t(M)
M1 = matrix(c(-1, -1, -1, 1, 0, 0, 0, 1, 0, 0, 0, 1), nalt_sim, nalt_sim+1)
S = list()
Diff_S = list()
Diff_S_1 = list()
scale = list()
TrueDiff_S = list()

for(t in 1:time_sim){
  S_reg = diag(1)*(1 + 0.05*t)
  S12 = as.matrix(c(0, 0.03*t, 0.05*t), nalt_sim, 1)
  S[[t]] = diag(nalt_sim+1)
  S[[t]][1:nalt_sim, 1:nalt_sim] = S_disc
  S[[t]][1:nalt_sim, nalt_sim+1] = S12
  S[[t]][nalt_sim+1, 1:nalt_sim] = t(S12)
  S[[t]][nalt_sim+1, nalt_sim+1] = S_reg
  Diff_S12 = as.matrix(c(0.03*t, 0.05*t), nalt_sim-1, 1)
  Diff_S[[t]] = diag(nalt_sim+1-1)
  Diff_S[[t]][1:nalt_sim-1, 1:nalt_sim-1] = Diff_S_disc
  Diff_S[[t]][1:nalt_sim-1, nalt_sim] = Diff_S12
  Diff_S[[t]][nalt_sim, 1:nalt_sim-1] = t(Diff_S12)
  Diff_S[[t]][nalt_sim+1-1, nalt_sim+1-1] = S_reg
  Diff_S_1[[t]] = Diff_S[[t]][1, 1]
  scale[[t]] = sqrt(Diff_S_1[[t]])
  TrueDiff_S[[t]]= Diff_S[[t]]/Diff_S_1[[t]]
}
#Diff_S[[t]] is the true covariance in difference at time t


#calculate utilities-------------------------------------------------------------------
#define model specification
spec = list(
  D = "dyn",
  Global = "global.txt",
  nvarD_disc = 1,
  nvarG_disc = 7,
  nvarD_reg = 1,
  
  generic = list(c("type.0", "type.1", "type.2")),
  specific = list(c(), c("asc", "sex", "edu", "inc", "kid", "res"), c("asc", "sex", "edu", "inc", "kid", "res")),
  reg = c("asc", "sex", "age", "inc", "res", "gas"),
  
  modifyD = function(D,Global,t){ 
    Dyan = cbind(D[[t]], Global)
    return(Dyan)
  },
  nTime = temp_time - temp_look,  #actual time, not total time
  nLook = temp_look,  #look ahead period
  transition = matrix(1, nalt_sim, nalt_sim),
  stopAlt = c() , #if you have more than one alternative 
  # that halts the decision process, put that in 
  # a vector (like stopAlt = c(1,2,10) if 
  # alt 1, 2 and 10 are stopping alternatives)
  method = "pmvnorm"
)

spec = checkFillDynSpec(spec)

#define function computeArgs
computeArgs = function(spec, D){
  X = list()
  RegX = list()
  
  for(t in 1:spec$nTime){
    Dreg = spec$modifyD(spec$D, spec$Global, t)
    RegX[[t]] = as.matrix(Dreg[, spec$reg])
    X[[t]] = list()
    for(l in 1:(spec$nLook+1)){
      Dt = spec$modifyD(spec$D, spec$Global, t+l-1)   #get all attributes at t
      X[[t]][[l]] = create_X(spec$generic, spec$specific, Dt)
    }
  }
  
  # check if we reached stop point already
  stopMat = matrix(0, spec$nObs, spec$nTime) 
  
  list(spec = spec, X = X, RegX = RegX, stop = stopMat, method = spec$method)    
}

args = computeArgs(spec, D)  
args = c(args, spec)


#calculate utility and choice
U = list()
choice = matrix(data = NA, nrow = nobs_sim, ncol = spec$nTime, byrow = FALSE)
Y = matrix(data = NA, nrow = nobs_sim, ncol = spec$nTime, byrow = FALSE)


for(t in 1:args$spec$nTime){
  #calculate utility for discrte part
  for(l in 1:(args$spec$nLook+1)){
    U[[l]] = matrix(args$X[[t]][[l]] %*% beta_disc, args$spec$nObs, args$spec$nAlt, byrow = T)
    #U[[l]] = U[[l]]/scale
  }
  
  if(nalt_sim > 2){
    fU = getfU(U, args$spec$transition)  #logsum of recursive logit
  } else if(nalt_sim == 2) {
    fU = getProbitU(U, Diff_S_1, args$spec$transition)  #obtain instant utility plus downstream utility at  time t
  }
  
  #calculate Y for continuous part
  Y[, t] = args$RegX[[t]] %*% beta_reg
  
  #simulate error components AND obatin final utilities
  Err = matrix(data = NA, nrow = nobs_sim, ncol = nalt_sim + 1, byrow = FALSE)
  Utility = matrix(data = NA, nrow = nobs_sim, ncol = nalt_sim, byrow = FALSE)
  Prob = matrix(0, nobs_sim, nalt_sim)
  #RegErr = rep(0, nobs_sim)
  for(j in 1:n_sim){
  Err = rmvnorm(n = nobs_sim, mean = rep(0, nalt_sim + 1), sigma = S[[t]])
  Utility = fU + Err[, 1:nalt_sim]
    for(h in 1:nobs_sim){
      if(Utility[h, 1] == max(Utility[h,])){
        Prob[h, 1] = Prob[h, 1] + 1/n_sim
      } else if(Utility[h, 2] == max(Utility[h,])){
        Prob[h, 2] = Prob[h, 2] + 1/n_sim
      }else
      Prob[h, 3] = 1-Prob[h, 1]-Prob[h, 2]
    }
  }
  Y[, t] = Y[, t] + Err[, nalt_sim + 1]
  
  for(h in 1:nobs_sim){
    temp = runif(1, min=0, max=1) #general a random  probability
    if(temp <= Prob[h, 1]){
      choice[h, t] = 0
    } else if(temp <= Prob[h, 1]+Prob[h, 2]){choice[h, t] = 1
    } else {choice[h, t] = 2
    }
  }
}

#save discrete dependent variable - "choice"
choice_name = paste("t", seq(1, spec$nTime, 1), sep="")
choice = rbind(choice_name, choice)
write.table(choice, file = "choice.txt", row.names=FALSE, col.names=FALSE)
#save continuous dependent variable - "miles"
miles_name = paste("vmt", seq(1, spec$nTime, 1), sep="")
miles = rbind(miles_name, Y)
write.table(miles, file = "miles.txt", row.names=FALSE, col.names=FALSE)
#mile = read.table(paste("miles", "txt", sep = "."), sep = " ", head = T)

}