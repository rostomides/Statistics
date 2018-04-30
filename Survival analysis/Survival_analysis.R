#Author : Larbi Bedrani
#Date : 30 April 2018
#version: 0.6

#The survival package provides no robust C-index and hence no confdence interval for the latter.
#This code allow to generate a robust C-index with 95% CI using a bootstrapping approach.

library(boot)
library(survival)

#############################
#Function performing CoxPH Model
#############################
coxPH_fit = function(data, TimeCol, eventsCol, mainIndependantVariableCol
                        , covariatesCols=NULL){
  #Build a table
  if(is.null(covariatesCols)){
    tab = data[,c(TimeCol, eventsCol, mainIndependantVariableCol)]
  }else{
    tab = data[,c(TimeCol, eventsCol, mainIndependantVariableCol, covariatesCols)]
    }
  names(tab)[1:2]=c("Time", "events")
  #Fit Cox PH model
  cox.surv = coxph(Surv(Time, events, type ="right") ~., data=tab)
  return(cox.surv)
}
##########################################################################################################
#Function for bootstrapping that calculates the Cox PH unbiased C-index as well as its confidence interval
##########################################################################################################
bs <- function(formula, data, indices) {
  d <- data[indices,] 
  fit <- coxph(formula, data=d)
  return(summary(fit)$concordance[1])
}

c_index = function(coxphObject, data, TimeCol, eventsCol, mainIndependantVariableCol
                   , nboot=10000, covariatesCols=NULL){
  #Create the sub table with relevant independant variable and covariates
  if(is.null(covariatesCols)){
    tab = data[,c(TimeCol, eventsCol, mainIndependantVariableCol)]
  }else{
    tab = data[,c(TimeCol, eventsCol, mainIndependantVariableCol, covariatesCols)]
  }
  names(tab)[1:2]=c("Time", "events")
  #Perform bootstrapping 
  booResults <- boot(data=tab, statistic=bs,
                     R=nboot, 
                     formula=coxphObject$formula)
  CI = boot.ci(booResults, type="bca", index = 1)
  return(CI)
}

########################################################
#Costom cox PH function with C-index confidence interval
########################################################
coxPH_costom = function(data, TimeCol, eventsCol, mainIndependantVariableCol
                        , covariatesCols=NULL){
  #Fit the cox PH model
  coxphObject =coxPH_fit(data = data, TimeCol = TimeCol, eventsCol=eventsCol,
                         mainIndependantVariableCol=mainIndependantVariableCol, 
                         covariatesCols=covariatesCols)
  #Calculate the unbiased c-index and its confidence interval
  cox_CI =c_index(coxphObject=coxphObject, data = data, TimeCol = TimeCol, eventsCol=eventsCol,
                  mainIndependantVariableCol=mainIndependantVariableCol, covariatesCols=covariatesCols)
  #Return the final object with the reslts
  final_list = list("CoxPH_Object" = coxphObject, "Unbiased_C_index" = cox_CI$t0, 
                    "C_index_confidence_interval"= cox_CI$bca[1,c(4,5)])
  return(final_list)
}


#Testing the code

data("lung")
names(lung)


model = coxPH_costom(data=lung, TimeCol="time", eventsCol="status", mainIndependantVariableCol="ph.karno",
                 covariatesCols=c("age","sex"))

#C-index
model$Unbiased_C_index

#C-index 95% CI
model$C_index_confidence_interval

#P value
summary(model$CoxPH_Object)$coefficients[1,5]

#Hazard ratio
summary(model$CoxPH_Object)$coefficients[1,2]

#Wald test p pvalue
summary(model$CoxPH_Object)$waldtest[3]

