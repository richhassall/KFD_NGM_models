#Authors: Richard Hassall and Nienke Hartemink 
#Date: 14/09/21

# NGM using habitat as a scaling factor

# The script extracts values from a csv file outlining the range of parameter values and uses 
# these to run multiple NGM models.
# Sensitivity and elasticity values are also calculated for individual parameters. 
# There is an option to exclude the egg type-at-infection 'remove_TOT', which will 
# remove the first row and column of the matrix 

####______________________________________####
#### 1. Load required libraries and files ####
####______________________________________####

library(Deriv)
library(popbio) 
library(stringr) 
library(reshape) 
library(tidyr)
library(ggplot2)
library(forcats)
library(dplyr)
library(VGAM)


## wd and csv file ##
  
params <- read.csv("./parameters/Parameters_scenarios_habitat.csv",stringsAsFactors = FALSE) ### load in csv with parameter ranges 


####____________________________####
####2. Inputs to control script ####
####____________________________####

scenarios <- rep(c("EVR","DCF","PLA","PAD","HAB"), times=2) ## column names

trans <- rep(c("TOT","NO_TT"),times=rep(length(scenarios)/2,times=2)) # tag to add to file name and control transovarial transmission

runs =100 ## number of runs 

for(sc in 1:length(scenarios)){ ## loop through scenarios

remove_TOT <- ifelse(trans[sc]=="TOT",FALSE,TRUE)
## should transovarial transmission be excluded

run_name <- paste0(scenarios[sc],"_",trans[sc]) ## name for saving results files

param_set <- params[,scenarios[sc]] ### select the column with range of parameters to use in models


####________________________####
#### 3.  Storage of outputs ####
####________________________####

Results=array(0,dim=c(runs,7)) #this sets up a matrix to store the results 

colnames(Results) <- c("TOT","NST","SYS","H1","H2","H3","R0") ### name the columns in results matrix

param_sens <- matrix(nrow=nrow(params),ncol=runs) ### df to store sensitivity and elasticity values

param_elas <- matrix(nrow=nrow(params),ncol=runs) ### df to store sensitivity and elasticity values

cat(paste0("Starting scenario ",run_name,"\n"))

pb = txtProgressBar(min = 0, max = runs, initial = 0,style=3) ## set up progress bar for loop

####___________________####
####4. Simulations loop####
####___________________####

for (r in 1:runs) {
  
  expres_df <- as.data.frame(matrix(nrow=7,ncol=7)) ## data frame to store formula for each parameter for sensitivty and elasiticity if individual parameters
  ####______________________________________####  
  #### 4.1 Extract range for all parameters ####
  ####______________________________________####
  ## This extracts range given in csv file (format: 2 ; 10) for each parameter
  ## Should change to extract values based on parameter name rather than row in future##
  
  ### extracts all parameters listed in csv and samples within range defined under scenario
  for(pr in 1:length(param_set)){
    assign(params$Parameter[pr],
           runif(1, min=as.numeric(str_trim(sub("\\;.*", "", param_set[pr]))), max=as.numeric(str_trim(sub('.*;', '', param_set[pr])))))
  }

  ####__________________________________####
  #### 4.2 Formulae for matrix elements ####
  ####__________________________________####
  ## mat element formula followed by code to store formula for sensitivity analysis ##
  
  ## Transovarial transmission ##
  
  #no eggs infected by one tick that was infected as an egg
  k11=Np*E*Sl*Sn*Sa*Ra	
  
  expres_df[1,1] <- expression(Np*E*Sl*Sn*Sa*Ra) 
  
  #no ticks infected as an egg (transovarial) produced by tick infected as a larva
  k12=Np*E*Sn*Sa*Ra  
  
  expres_df[1,2] <- expression(Np*E*Sn*Sa*Ra)
  
  #no eggs infected by  tick that was infected while feeding as a nymph
  k13=Np*E*Sa*Ra 
  
  expres_df[1,3] <- expression(Np*E*Sa*Ra)
  
  #no eggs infected by tick that was infected while feeding as an adult
  k14=Np*E*Ra 
  
  expres_df[1,4] <- expression(Np*E*Ra)
  
  k15=0 #0 because vertebrate host cannot pass infection to tick eggs
  
  expres_df[1,5] <- NA #expression()
  
  k16=0 #0 because vertebrate host cannot pass infection to tick eggs
  
  expres_df[1,6] <- NA # expression()
  
  k17=0 #0 because vertebrate host cannot pass infection to tick eggs
  
  expres_df[1,7] <- NA # expression()
  
  ##elements arising from non-systemic (cofeeding transmission)##
  
  #no larvae infected by one tick that was infected as an egg
  k21=Sl*Hc*(PbH4*TT*(Np*CllH4*Cs) + PbH3*TT*(Np*CllH3*Cs) + PbH2*TT*(Np*CllH2*Cs) + PbH1*TT*(Np*CllH1*Cs)) +  
    Sl*Sn*Hc*(PbH4*TT*(Np*ClnH4*Cs) + PbH3*TT*(Np*ClnH3*Cs) + PbH2*TT*(Np*ClnH2*Cs) + PbH1*TT*(Np*ClnH1*Cs)) +
    Sl*Sn*Sa*Hc*(PbH4*TT*(Np*ClaH4*Cs) + PbH3*TT*(Np*ClaH3*Cs) + PbH2*TT*(Np*ClaH2*Cs) + PbH1*TT*(Np*ClaH1*Cs)) 
  
  expres_df[2,1] <- expression(Sl*Hc*(PbH4*TT*(Np*CllH4*Cs) + PbH3*TT*(Np*CllH3*Cs) + PbH2*TT*(Np*CllH2*Cs) + PbH1*TT*(Np*CllH1*Cs)) +  
                                 Sl*Sn*Hc*(PbH4*TT*(Np*ClnH4*Cs) + PbH3*TT*(Np*ClnH3*Cs) + PbH2*TT*(Np*ClnH2*Cs) + PbH1*TT*(Np*ClnH1*Cs)) +
                                 Sl*Sn*Sa*Hc*(PbH4*TT*(Np*ClaH4*Cs) + PbH3*TT*(Np*ClaH3*Cs) + PbH2*TT*(Np*ClaH2*Cs) + PbH1*TT*(Np*ClaH1*Cs)))
  
  # number of larva infected by tick-infected-as-larva
  k22=Sn*Hc*(PbH4*TT*(Np*ClnH4*Cs)+PbH3*TT*(Np*ClnH3*Cs) + PbH2*TT*(Np*ClnH2*Cs) + PbH1*TT*(Np*ClnH1*Cs))+
    Sn*Sa*Hc*(PbH4*TT*(Np*ClaH4*Cs)+PbH3*TT*(Np*ClaH3*Cs) + PbH2*TT*(Np*ClaH2*Cs) + PbH1*TT*(Np*ClaH1*Cs))
  
  expres_df[2,2] <- expression(Sn*Hc*(PbH4*TT*(Np*ClnH4*Cs)+PbH3*TT*(Np*ClnH3*Cs) + PbH2*TT*(Np*ClnH2*Cs) + PbH1*TT*(Np*ClnH1*Cs))+
                                 Sn*Sa*Hc*(PbH4*TT*(Np*ClaH4*Cs)+PbH3*TT*(Np*ClaH3*Cs) + PbH2*TT*(Np*ClaH2*Cs) + PbH1*TT*(Np*ClaH1*Cs)))
  
  # number of larva infected by tick-infected-as-nymph
  k23=Sa*Hc*(PbH4*TT*(Np*ClaH4*Cs)+PbH3*TT*(Np*ClaH3*Cs) + PbH2*TT*(Np*ClaH2*Cs) + PbH1*TT*(Np*ClaH1*Cs)) 
  
  expres_df[2,3] <- expression(Sa*Hc*(PbH4*TT*(Np*ClaH4*Cs)+PbH3*TT*(Np*ClaH3*Cs) + PbH2*TT*(Np*ClaH2*Cs) + PbH1*TT*(Np*ClaH1*Cs)))
  
  k24=0
  expres_df[2,4] <- NA#expression() 
  
  # number of nymphs infected by tick-infected-as-egg
  k31=Sl*Hc*(PbH4*TT*(Np*CnlH4*Cs)+PbH3*TT*(Np*CnlH3*Cs) + PbH2*TT*(Np*CnlH2*Cs) + PbH1*TT*(Np*CnlH1*Cs)) + 
    Sl*Sn*Hc*(PbH4*TT*(Np*CnnH4*Cs)+PbH3*TT*(Np*CnnH3*Cs) + PbH2*TT*(Np*CnnH2*Cs) + PbH1*TT*(Np*CnnH1*Cs)) +
    Sl*Sn*Sa*Hc*(PbH4*TT*(Np*CnaH4*Cs)+PbH3*TT*(Np*CnaH3*Cs) + PbH2*TT*(Np*CnaH2*Cs) + PbH1*TT*(Np*CnaH2*Cs))
  
  expres_df[3,1] <- expression(Sl*Hc*(PbH4*TT*(Np*CnlH4*Cs)+PbH3*TT*(Np*CnlH3*Cs) + PbH2*TT*(Np*CnlH2*Cs) + PbH1*TT*(Np*CnlH1*Cs)) + 
                                 Sl*Sn*Hc*(PbH4*TT*(Np*CnnH4*Cs)+PbH3*TT*(Np*CnnH3*Cs) + PbH2*TT*(Np*CnnH2*Cs) + PbH1*TT*(Np*CnnH1*Cs)) +
                                 Sl*Sn*Sa*Hc*(PbH4*TT*(Np*CnaH4*Cs)+PbH3*TT*(Np*CnaH3*Cs) + PbH2*TT*(Np*CnaH2*Cs) + PbH1*TT*(Np*CnaH2*Cs)))
  
  # number of nymphs infected by tick-infected-as-larvae
  k32=Sn*Hc*(PbH4*TT*(Np*CnnH4*Cs)+PbH3*TT*(Np*CnnH3*Cs) + PbH2*TT*(Np*CnnH2*Cs) + PbH1*TT*(Np*CnnH1*Cs))+
    Sn*Sa*Hc*(PbH4*TT*(Np*CnaH4*Cs)+PbH3*TT*(Np*CnaH3*Cs) + PbH2*TT*(Np*CnaH2*Cs) + PbH1*TT*(Np*CnaH1*Cs))
  
  expres_df[3,2] <- expression(Sn*Hc*(PbH4*TT*(Np*CnnH4*Cs)+PbH3*TT*(Np*CnnH3*Cs) + PbH2*TT*(Np*CnnH2*Cs) + PbH1*TT*(Np*CnnH1*Cs))+
                                 Sn*Sa*Hc*(PbH4*TT*(Np*CnaH4*Cs)+PbH3*TT*(Np*CnaH3*Cs) + PbH2*TT*(Np*CnaH2*Cs) + PbH1*TT*(Np*CnaH1*Cs)))
  
  # number of nymphs infected by tick-infected-as-nymph
  k33=Sa*Hc*(PbH4*TT*(Np*CnaH4*Cs)+PbH3*TT*(Np*CnaH3*Cs) + PbH2*TT*(Np*CnaH2*Cs) + PbH1*TT*(Np*CnaH1*Cs)) 
  
  expres_df[3,3] <- expression(Sa*Hc*(PbH4*TT*(Np*CnaH4*Cs)+PbH3*TT*(Np*CnaH3*Cs) + PbH2*TT*(Np*CnaH2*Cs) + PbH1*TT*(Np*CnaH1*Cs)))
  
  
  k34=0 #0 because tick infected as an adult (during final blood meal) can only pass infection to her eggs
  
  expres_df[3,4] <- NA #expression()
  
  #no adults infected by one tick that was infected as an egg
  k41=Sl*Hc*(PbH4*TT*(Np*CalH4*Cs)+PbH3*TT*(Np*CalH3*Cs) + PbH2*TT*(Np*CalH2*Cs) + PbH1*TT*(Np*CalH1*Cs)) + 
    Sl*Sn*Hc*(PbH4*TT*(Np*CanH4*Cs)+PbH3*TT*(Np*CanH3*Cs) + PbH2*TT*(Np*CanH2*Cs) + PbH1*TT*(Np*CanH1*Cs)) +
    Sl*Sn*Sa*Hc*(PbH4*TT*(Np*CaaH4*Cs)+PbH3*TT*(Np*CaaH3*Cs) + PbH2*TT*(Np*CaaH2*Cs) + PbH1*TT*(Np*CaaH1*Cs))
  
  expres_df[4,1] <- expression(Sl*Hc*(PbH4*TT*(Np*CalH4*Cs)+PbH3*TT*(Np*CalH3*Cs) + PbH2*TT*(Np*CalH2*Cs) + PbH1*TT*(Np*CalH1*Cs)) + 
                                 Sl*Sn*Hc*(PbH4*TT*(Np*CanH4*Cs)+PbH3*TT*(Np*CanH3*Cs) + PbH2*TT*(Np*CanH2*Cs) + PbH1*TT*(Np*CanH1*Cs)) +
                                 Sl*Sn*Sa*Hc*(PbH4*TT*(Np*CaaH4*Cs)+PbH3*TT*(Np*CaaH3*Cs) + PbH2*TT*(Np*CaaH2*Cs) + PbH1*TT*(Np*CaaH1*Cs)))
  
  #no adults infected by one tick that was infected while feeding as a larva
  k42=Sn*Hc*(PbH4*TT*(Np*CanH4*Cs)+PbH3*TT*(Np*CanH3*Cs) + PbH2*TT*(Np*CanH2*Cs) + PbH1*TT*(Np*CanH1*Cs)) +
    Sn*Sa*Hc*(PbH4*TT*(Np*CaaH4*Cs)+PbH3*TT*(Np*CaaH3*Cs) + PbH2*TT*(Np*CaaH2*Cs) + PbH1*TT*(Np*CaaH1*Cs))
  
  expres_df[4,2] <- expression(Sn*Hc*(PbH4*TT*(Np*CanH4*Cs)+PbH3*TT*(Np*CanH3*Cs) + PbH2*TT*(Np*CanH2*Cs) + PbH1*TT*(Np*CanH1*Cs)) +
                                 Sn*Sa*Hc*(PbH4*TT*(Np*CaaH4*Cs)+PbH3*TT*(Np*CaaH3*Cs) + PbH2*TT*(Np*CaaH2*Cs) + PbH1*TT*(Np*CaaH1*Cs)))
  
  #no adults infected by one tick that was infected while feeding as a nymph
  k43=Sa*Hc*(PbH4*TT*(Np*CaaH4*Cs)+PbH3*TT*(Np*CaaH3*Cs) + PbH2*TT*(Np*CaaH2*Cs) + PbH1*TT*(Np*CaaH1*Cs))
  expres_df[4,3] <- expression(Sa*Hc*(PbH4*TT*(Np*CaaH4*Cs)+PbH3*TT*(Np*CaaH3*Cs) + PbH2*TT*(Np*CaaH2*Cs) + PbH1*TT*(Np*CaaH1*Cs)))
  
  k44=0  #0 because tick infected as an adult (during final blood meal) can only pass infection to her eggs
  expres_df[4,4] <- NA # expression()
  
  ### Systemic transmission ###
  
  #number of larvae infected by one infected H1
  k25=(Np*NlhH1*Pl*IH1)/Dl 
  expres_df[2,5] <- expression((Np*NlhH1*Pl*IH1)/Dl)
  
  #number of larvae infected by one infected H2
  k26=(Np*NlhH2*Pl*IH2)/Dl 
  expres_df[2,6] <- expression((Np*NlhH2*Pl*IH2)/Dl)

  #number of larvae infected by one infected H3
  k27=(Np*NlhH3*Pl*IH3)/Dl 
  expres_df[2,7] <- expression((Np*NlhH3*Pl*IH3)/Dl)
  
  #number of nymphs infected by one infected H1
  k35=(Np*NnhH1*Pn*IH1)/Dn 
  expres_df[3,5] <- expression((Np*NnhH1*Pn*IH1)/Dn)

  #number of nymphs infected by one infected H2
  k36=(Np*NnhH2*Pn*IH2)/Dn 
  expres_df[3,6] <- expression((Np*NnhH2*Pn*IH2)/Dn)
 
  #number of nymphs infected by one infected H3
  k37=(Np*NnhH3*Pn*IH3)/Dn 
  expres_df[3,7] <- expression((Np*NnhH3*Pn*IH3)/Dn )

  
  #number of adults infected by one infected H1
  
  k45=(Np*NahH1*Pa*IH1)/Da 
  expres_df[4,5] <- expression((Np*NahH1*Pa*IH1)/Da)

  #number of adults infected by one infected H2
  k46=(Np*NahH2*Pa*IH2)/Da 
  expres_df[4,6] <- expression((Np*NahH2*Pa*IH2)/Da)

  #number of adults infected by one infected H3
  k47=(Np*NahH3*Pa*IH3)/Da
  expres_df[4,7] <- expression((Np*NahH3*Pa*IH3)/Da)

  
  # numbers of Host 1 infected by ......
  # tick infected as egg 
  k51=Sl*QlH1*Hc*PbLH1 + Sl*Sn*QnH1*Hc*PbNH1 + Sl*Sn*Sa*QaH1*Hc*PbAH1
  
  expres_df[5,1] <- expression(Sl*QlH1*Hc*PbLH1 + Sl*Sn*QnH1*Hc*PbNH1 + Sl*Sn*Sa*QaH1*Hc*PbAH1)
  # tick infected as larva
  k52=Sl*Sn*QnH1*Hc*PbNH1 + Sl*Sn*Sa*QaH1*Hc*PbAH1
  
  expres_df[5,2] <- expression(Sl*Sn*QnH1*Hc*PbNH1 + Sl*Sn*Sa*QaH1*Hc*PbAH1)
  # tick infected as nymph
  k53=Sl*Sn*Sa*QaH1*Hc*PbAH1
  
  expres_df[5,3] <- expression(Sl*Sn*Sa*QaH1*Hc*PbAH1)
  
  k54=0 # 0 as tick infected as adult can only infect eggs
  
  expres_df[5,4] <- NA #expression()
  
  k55=0 # 0 as hosts can't infect each other
  
  expres_df[5,5] <- NA #expression()
  
  k56=0 # 0 as hosts can't infect each other
  
  expres_df[5,6] <- NA #expression()
  
  k57=0 # 0 as hosts can't infect each other
  
  expres_df[5,7] <- NA #expression()
  
  # numbers of Host 2 infected by ......
  # tick infected as egg 
  # tick infected as egg 
  k61=Sl*QlH2*Hc*PbLH2 + Sl*Sn*QnH2*Hc*PbNH2 + Sl*Sn*Sa*QaH2*Hc*PbAH2
  
  expres_df[6,1] <- expression(Sl*QlH2*Hc*PbLH2 + Sl*Sn*QnH2*Hc*PbNH2 + Sl*Sn*Sa*QaH2*Hc*PbAH2)
  # tick infected as larva
  k62=Sl*Sn*QnH2*Hc*PbNH2 + Sl*Sn*Sa*QaH2*Hc*PbAH2
  expres_df[6,2] <- expression(Sl*Sn*QnH2*Hc*PbNH2 + Sl*Sn*Sa*QaH2*Hc*PbAH2)
  
  # tick infected as nymph
  k63=Sl*Sn*Sa*QaH2*Hc*PbAH2
  expres_df[6,3] <- expression(Sl*Sn*Sa*QaH2*Hc*PbAH2)
  
  k64=0 # 0 as tick infected as adult can only infect eggs
  expres_df[6,4] <- NA #expression()
  
  k65=0 # 0 as hosts can't infect each other
  expres_df[6,5] <- NA #expression()
  
  k66=0 # 0 as hosts can't infect each other
  expres_df[6,6] <- NA #expression()
  
  k67=0 # 0 as hosts can't infect each other
  expres_df[6,7] <- NA #expression()
  
  
  # numbers of Host 3 infected by ......
  # tick infected as egg
  k71=Sl*QlH3*Hc*PbLH3 + Sl*Sn*QnH3*Hc*PbNH3 + Sl*Sn*Sa*QaH3*Hc*PbAH3
  expres_df[7,1] <- expression(Sl*QlH3*Hc*PbLH3 + Sl*Sn*QnH3*Hc*PbNH3 + Sl*Sn*Sa*QaH3*Hc*PbAH3)#
  
  # tick infected as larva
  k72=Sl*Sn*QnH3*Hc*PbNH3 + Sl*Sn*Sa*QaH3*Hc*PbAH3
  expres_df[7,2] <- expression(Sl*Sn*QnH3*Hc*PbNH3 + Sl*Sn*Sa*QaH3*Hc*PbAH3)
  
  # tick infected as nymph
  k73=Sl*Sn*Sa*QaH3*Hc*PbAH3
  expres_df[7,3] <- expression(Sl*Sn*Sa*QaH3*Hc*PbAH3)
  
  k74=0 # 0 as tick infected as adult can only infect eggs
  expres_df[7,4] <-NA #expression()
  
  k75=0 # 0 as hosts can't infect each other
  expres_df[7,5] <-NA # expression()
  
  k76=0 # 0 as hosts can't infect each other
  expres_df[7,6] <-NA # expression()
  
  k77=0 # 0 as hosts can't infect each other
  expres_df[7,7] <-NA # expression()
  
  
  
  
  ####__________________####    
  #### 4.3 Build matrix ####
  ####__________________####
  
  F=c(
    k11,k12,k13,k14,k15,k16,k17,
    k21,k22,k23,k24,k25,k26,k27,
    k31,k32,k33,k34,k35,k36,k37,
    k41,k42,k43,k44,k45,k46,k47,
    k51,k52,k53,k54,k55,k56,k57,
    k61,k62,k63,k64,k65,k66,k67,
    k71,k72,k73,k74,k75,k76,k77)
  A=matrix(F, nrow=7, byrow=TRUE) ## puts in matrix form with dimensions 7x7
  
  ## if TOT is to be excluded get rid of the first row and column of the matrix and expression data frame
  if(remove_TOT==TRUE){
    A <- A[-1,]
    A <- A[,-1]
    
    expres_df <- expres_df[-1,]
    expres_df <- expres_df[,-1]
  }
  
  ####____________________________________________________________________####
  #### 4.4 Get R0, elasticity and sensitivity values for matrix elements ####
  ####____________________________________________________________________####
  
  ##R0 is largest eigenvalue of the NGM matrix A
  R0=Re(lambda(A))  # largest eigenvalue in the matrix 
  
  MES=eigen(A) 
  REV=MES$vec[,1] 
  TM=aperm(A,c(2,1)) 
  TMES=eigen(TM)
  LEV=TMES$vec[,1] 
  Sens=outer(LEV,REV,"*")/c(LEV%*%REV) ## Sensitivity values 
  Elas=1/MES$values[1]*Sens*A # Elasticity values
  
  if(r==1){
    ### build array for storing elasticity matrices
    if(remove_TOT==FALSE){
    Elas_arr <- array(NA,dim=c(7,7,runs))
    }else{
      Elas_arr <- array(NA,dim=c(6,6,runs)) 
    }
    }
  Elas_arr[,,r] <- Elas
  ## if TOT is not exluded assume 7 x 7 matrix ##
  if(remove_TOT==FALSE){
    ## Contribution of TOT to R0  
    TOT=Elas[1,1]+Elas[1,2]+Elas[1,3]+Elas[1,4] 
    
    ## Contribution of Systemic transmission to R0
    ST=Elas[5,1]+Elas[5,2]+Elas[5,3]+ 
      Elas[6,1]+Elas[6,2]+Elas[6,3]+
      Elas[7,1]+Elas[7,2]+Elas[7,3]+ 
      Elas[2,5]+Elas[2,6]+Elas[2,7]+
      Elas[3,5]+Elas[3,6]+Elas[3,7]+
      Elas[4,5]+Elas[4,6]+Elas[4,7]
    
    ## Contribution of non systemic transmission to R0
    NST=Elas[2,1]+Elas[2,2]+Elas[2,3]+
      Elas[3,1]+Elas[3,2]+Elas[3,3]+
      Elas[4,1]+Elas[4,2]+Elas[4,3] 
    
    ## Contribution of H1 to systemic transmission
    ElasH1=Elas[2,5]+Elas[3,5]+Elas[4,5]+
      Elas[5,1]+Elas[5,2]+Elas[5,3]
    
    ## Contribution of H2 to systemic transmission
    ElasH2=Elas[2,6]+Elas[3,6]+Elas[4,6]+
      Elas[6,1]+Elas[6,2]+Elas[6,3] 
    
    ## Contribution of H3 to systemic transmission
    ElasH3=Elas[2,7]+Elas[3,7]+Elas[4,7]+
      Elas[7,1]+Elas[7,2]+Elas[7,3] 
    
    Results[r,]=c(Re(TOT),Re(NST),Re(ST), Re(ElasH1), Re(ElasH2),Re(ElasH3), R0) ## store results
    
  }else{ ## else TOT is excluded assume 6 x 6 matrix ##
    TOT=0
    
    ## Contribution of Systemic transmission to R0
    ST=Elas[4,1]+Elas[4,2]+ 
      Elas[5,1]+Elas[5,2]+
      Elas[6,1]+Elas[6,2]+ 
      Elas[1,4]+Elas[1,5]+Elas[1,6]+
      Elas[2,4]+Elas[2,5]+Elas[2,6]+
      Elas[3,4]+Elas[3,5]+Elas[3,6]
    
    ## Contribution of non systemic transmission to R0
    NST=Elas[1,1]+Elas[1,2]+
      Elas[2,1]+Elas[2,2]+
      Elas[3,1]+Elas[3,2] 
    
    ## relative contribution of H1 to systemic transmission
    ElasH1=Elas[1,4]+Elas[2,4]+Elas[3,4]+
      +Elas[4,1]+Elas[4,2]
    
    ## relative contribution of H2 to systemic transmission
    ElasH2=Elas[1,5]+Elas[2,5]+Elas[3,5]+
      +Elas[5,1]+Elas[5,2] 
    
    ## relative contribution of H3 to systemic transmission
    ElasH3=Elas[1,6]+Elas[2,6]+Elas[3,6]+
      Elas[6,1]+Elas[6,2] 
    
    Results[r,]=c(TOT,NST,ST, ElasH1, ElasH2,ElasH3, R0) ## store results
    
  }
  
  #### 4.5 Get sensitivity and elasticity values for individual parameters ####
  
  for(e in 1:nrow(param_sens)){ ### loop to cycle though each parameter
    
    focal.parameter <- params[e,"Parameter"] ## store name of focal parameter
    
    par_mat_elements <- character(0) ## character vector to store matrix elements that parameter contributes to
    
    ## find the matrix elements that focal parameter contributes to 
    
    for(expi in 1:nrow(expres_df)){
      for(expj in 1:ncol(expres_df)){
        if(grepl(focal.parameter,as.character(expres_df[expi,expj]),fixed=TRUE)){ ## if parameter is in expression
          par_mat_elements[length(par_mat_elements)+1] <- paste0(expi,",",expj) ## store index for matrix
        }
      }
    }
    
    if(length(par_mat_elements)>0){## if parameter contributes to at least one  matrix element calculate sensitivity (included so code doesn't break when TOT is excluded) 
      param_sens_vec <- rep(NA,length(par_mat_elements)) ### create vector to store calculation for each relevant matrix element
      
      for(me in 1:length(par_mat_elements)){ ### for each matrix element that parameter contributes to
        
        i <- as.numeric(str_trim(sub("\\,.*", "", par_mat_elements[me]))) ## extract row of focal matrix element
        j <- as.numeric(str_trim(sub('.*,', '', par_mat_elements[me]))) ## extract column of focal matrix element
        
        dir <- as.expression(Deriv(expres_df[i,j],focal.parameter)) ### extract expression from expression matrix and get formula for partial derivative 
        
        sol <- eval(dir) ## Solve formula for partial derivative (dkij/da)
        
        param_sens_vec[me] <- sol*Re(Sens[i,j]) ## calculate sensitivity of R0 to individual parameter (dkij/da * dR0/kij)
      }## end me loop
    }else{param_sens_vec <- NA} ## if parameter is excluded due to TOT being excluded give NA
    
    param_sens[e,r] <- abs(sum(param_sens_vec)) ### sum sensitivity values across all matrix elements focal parameter contributes to
    
    param_elas[e,r] <-  abs(sum(param_sens_vec))*(get(focal.parameter)/R0) ### calculate elasticity for individual parameters
  } ## end e loop
  
  setTxtProgressBar(pb,r) ## update progress bar 
} ## end run loop



#### 5. Store results from sensitivity and elasiticty of individual parameters ####
sens_elas_results <- data.frame(Parameter=params$Parameter,Sens=NA,Sens_lower=NA,Sens_upper=NA,Elas=NA,Elas_lower=NA,Elas_upper=NA)

## calculate mean and quantiles for sensitivity and elasticity values

for(se in 1:nrow(sens_elas_results)){  
  sens_elas_results$Sens[se] <- mean(param_sens[se,],na.rm=TRUE)
  sens_elas_results$Sens_lower[se] <- as.numeric(quantile(param_sens[se,],0.025,na.rm=TRUE))
  sens_elas_results$Sens_upper[se] <- as.numeric(quantile(param_sens[se,],0.975,na.rm=TRUE))

  sens_elas_results$Elas[se] <- mean(param_elas[se,],na.rm=TRUE)
  sens_elas_results$Elas_lower[se] <- as.numeric(quantile(param_elas[se,],0.025,na.rm=TRUE))
  sens_elas_results$Elas_upper[se] <- as.numeric(quantile(param_elas[se,],0.975,na.rm=TRUE))
}

## print results to check 
Results.df <- as.data.frame(Re(Results)) ## create data frame from results array

cat(paste0("\n","R0=",mean(Results.df[,"R0"])," (",as.numeric(quantile(Results.df[,"R0"],c(0.025),na.rm=TRUE)),";",as.numeric(quantile(Results.df[,"R0"],c(0.975),na.rm=TRUE)),")","\n",  
"TOT=",mean(Results.df[,"TOT"])," (",as.numeric(quantile(Results.df[,"TOT"],c(0.025),na.rm=TRUE)),";",as.numeric(quantile(Results.df[,"TOT"],c(0.975),na.rm=TRUE)),")","\n",
"NST=",mean(Results.df[,"NST"])," (",as.numeric(quantile(Results.df[,"NST"],c(0.025),na.rm=TRUE)),";",as.numeric(quantile(Results.df[,"NST"],c(0.975),na.rm=TRUE)),")","\n",
"SYS=",mean(Results.df[,"SYS"])," (",as.numeric(quantile(Results.df[,"SYS"],c(0.025),na.rm=TRUE)),";",as.numeric(quantile(Results.df[,"SYS"],c(0.975),na.rm=TRUE)),")","\n",
"H1=",mean(Results.df[,"H1"])," (",as.numeric(quantile(Results.df[,"H1"],c(0.025),na.rm=TRUE)),";",as.numeric(quantile(Results.df[,"H1"],c(0.975),na.rm=TRUE)),")","\n",
"H2=",mean(Results.df[,"H2"])," (",as.numeric(quantile(Results.df[,"H2"],c(0.025),na.rm=TRUE)),";",as.numeric(quantile(Results.df[,"H2"],c(0.975),na.rm=TRUE)),")","\n",
"H3=",mean(Results.df[,"H3"])," (",as.numeric(quantile(Results.df[,"H3"],c(0.025),na.rm=TRUE)),";",as.numeric(quantile(Results.df[,"H3"],c(0.975),na.rm=TRUE)),")","\n",
"All trans=",mean(Results.df[,"TOT"],na.rm=TRUE)+mean(Results.df[,"NST"],na.rm=TRUE)+mean(Results.df[,"SYS"],na.rm=TRUE),"\n",
"All hosts=",mean(Results.df[,"H1"],na.rm=TRUE)+mean(Results.df[,"H2"],na.rm=TRUE)+mean(Results.df[,"H3"],na.rm=TRUE),"\n"
))




####6.  Save files ####

 write.csv(sens_elas_results,paste0(run_name,"_sens_elas.csv"),row.names = FALSE) ## save sensitivity and elasticity results for individual parameters
 
 write.csv(Results.df,paste0(run_name,"_results.csv"),row.names = FALSE) ## save results 
 
 saveRDS(Elas_arr, paste0(run_name,"_elas_array.rds"))
 
} # end sc loop
