library(Hmisc)
library(parallel)
library(MASS)
library(plyr)
library(dplyr)
library(ggplot2)
library(pomp)
library(reshape2)
library(tidyverse)
library(tidyr)
library(matrixStats)
#library(DescTools)
library(extdplyr)



cmv_sim = function(params){
  
  setwd("/home/sgazea01/projects/def-craigm/sgazea01/NewCon1")
  #setwd("C:/Users/soniu/OneDrive/Dokumenty/BrazilTorun7")
  #setwd("C:/Users/c-byr.LAPTOP-UG0F9R6R/OneDrive - math.ubc.ca/cmv/cmv transmission model/POLYMOD/original contact matrix, no daycare/CMV_transmission_model40")
  
  RoundTo = function (x, multiple = 1, FUN = round) 
  {
    if (is.function(FUN)) {
      fct <- FUN
      FUN <- "fct"
      FUN <- gettextf("%s", FUN)
    }
    return(eval(parse(text = gettextf("%s(x/multiple) * multiple", 
                                      FUN))))
  }
  
  t = 0
  delta_t = 1/(12) #set delta_t = 2 weeks
  #params = c(-2.4016959, -1.4045884, -0.4425165, -0.8195270,  1.4725439,  0.3689964,  1.0000000,  0.5000000)
  
  params = tail(params,8)
  params[1:6] = 10^params[1:6]
  
  child = params[1]
  adult = params[1]*params[2]
  secondary = params[3]
  chronic = params[1]*params[2]*params[3]*params[4]
  breast = params[5]
  diaper = params[6]
  immune_waning = params[7]
  immune_const = params[8]
  duration = 1.5 #how long until children hit chronic level of infection
  
  params = c(adult,child,chronic,breast,diaper,secondary,immune_waning, immune_const,duration)
  params = as.data.frame(t(params))
  colnames(params) = c("adult_infectivity_val",
                       "child_infectivity_val",
                       "chronic_infectivity_val",
                       "breast_factor_val",
                       "diaper_factor_val",
                       "secondary_factor_val",
                       "immune_waning",
                       "immune_const",
                       "duration_child_val")
  
  
  
  rand = round(runif(1,1,100000),2)
  print(rand)
  while(length(list.files(path = ".", pattern = as.character(rand)))>0){
    rand = round(runif(1,1,100000),2)
  }
  write.table(data.frame(setup = 0), paste0("sim_final",rand,".csv"), sep = ",", col.names = T,row.names = F, append = F)
  


  #BIRTH RATE
  a = 0.1259 #rate at which women aged 15-19 have babies per year
  b = 0.2448 #rate at which women aged 20-24 have babies per year
  c = 0.2449 # rate at which women aggood_initial_fit_data_10000_newed 25-29 have babies per year
  d = 0.2065 #rate at which women aged 29-34 have babies per year 
  e = 0.1308 #rate at which women age 35-39 have babies per year 
  f = 0.0364 #rate at which women aged 40-44 have babies per year
  #rates at which adults give birth for each age 
  birth =c(rep(0,15),rep(c(a,b,c,d,e,f),each = 5),rep(0,5))*2/4 #scales rates so that women have an ave of 2 children each
  birth = data.frame(age = c(0:49),gender = 1,rate_birth = birth)
  birth =rbind(birth,data.frame(age = c(0:49),gender = 0,rate_birth = 0))
  
  #CHRONIC INFECTION
  #taken from hogea et al 2015 cmv math model
  #0-5 y.o. = 24 months = 2 years
  #10 y.o.+ = 2 months = 1/6 years
  #exponential decay inbetween
  r = -log((1/6)/duration)/(10-5)
  chronic = data.frame(age_of_infection = seq(0,5,delta_t),time_to_chronic = duration)
  for(i in seq((5+delta_t),10,delta_t)){
    chronic = rbind(chronic,data.frame(age_of_infection = i, time_to_chronic = duration*exp(-r*(i-5))))
  }
  chronic = rbind(chronic,data.frame(age_of_infection = seq((10+delta_t),(50-delta_t),delta_t), time_to_chronic = (1/6)))
  
  
  #TRANSMISSIBILITY PARAMS
  #max infectivity for every age
  #age 0-5 have infectivity of child_infectivity_val
  #age 10+ have infectivity of adult_infectivity_val
  #age 5-10 have exponential decay of infectivity 
  r = -log(params$adult_infectivity_val/params$child_infectivity_val)/(10-5)
  max_infectivity = data.frame(age_of_infection = seq(0,5,delta_t),max.infectivity =rep(params$child_infectivity_val))
  for(i in seq((5+delta_t),10,delta_t)){
    max_infectivity = rbind(max_infectivity,data.frame(age_of_infection = i, max.infectivity = params$child_infectivity_val*exp(-r*(i-5))))
  }
  max_infectivity = rbind(max_infectivity,data.frame(age_of_infection =  seq((10+delta_t),(50-delta_t),delta_t),max.infectivity = params$adult_infectivity_val))
  max_infectivity$max.infectivity[is.nan(max_infectivity$max.infectivity)]=0
  infectivity_params = data.frame(chronic,max.infectivity=max_infectivity$max.infectivity)
  infectivity_params$age_of_infection = round(infectivity_params$age_of_infection,2)
  
  #IMMUNE WANING
  #want immune waning to happen once individuals reach the same transmissibility as that seen during a reinfection
  
  max_infectivity$reinfection.infectivity = infectivity_params$max.infectivity*params$secondary_factor_val
  max_infectivity$time_to_wane = log(max_infectivity$reinfection.infectivity/max_infectivity$max.infectivity)/(log(params$chronic_infectivity_val/infectivity_params$max.infectivity)/
                                                                                                                 chronic$time_to_chronic)
  #let immune protection be the probability of reinfection with time_since_last_infection=0 being when waning starts 
  r = -params$immune_waning
  immune_protection = data.frame(time_since_last_infection = rep(seq(0,50,delta_t),2),infection_status=rep(c(1,2),each = length(seq(0,50,delta_t))),immune_protection = rep(c(0,1),each = length(seq(0,50,delta_t))))
  for (i in seq(0,50,delta_t)){
    immune_protection = rbind(immune_protection,data.frame(time_since_last_infection = i,infection_status=3, immune_protection = exp(r*i)))
  }
  immune_protection$immune_protection[which(immune_protection$immune_protection<params$immune_const&immune_protection$infection_status==3)]=params$immune_const
  
  immune_protection$time_since_last_infection = round(immune_protection$time_since_last_infection,3)
  immune_protection[which(immune_protection$immune_protection<0.000001),"immune_protection"]=0
  #CONGENITAL PARAMS
  #probability of healthy vs infected birth given infection status of mother
  #know from Kenneson and Cannon 2007 that 32% of primary infected mothers will pass infection to baby
  #know from Hogea et al that 33.4% of primary infected mothers and 8.5% of recurrently infected mothers will pass infection to baby
  #thus all adults start with a 33.4% chance of congenital transmission when they acquire primary CMV
  #this probability decays over the 2 months it takes them to establish chronic infection
  #at chronic infection, probability of infection should be 0.334*chronic_param/adult_param
  prob_congen_chronic = 0.334*params$chronic_infectivity_val/params$adult_infectivity_val
  r = -6*log(params$chronic_infectivity_val/params$adult_infectivity_val)
  if(params$adult_infectivity_val==0)
    r=0
  prob_inf = NULL
  for (i in seq(0,50,delta_t)){
    prob_inf = rbind(prob_inf,data.frame(infection_status = c(2,3), time_since_last_infection = i, prob_congenital = c((0.334*exp(-r*i)),(0.334*exp(-r*i)*params$secondary_factor_val))))
  }
  if(params$adult_infectivity_val==0)
    prob_inf$prob_congenital = 0
  prob_inf[which(prob_inf$prob_congenital<=prob_congen_chronic),c("prob_congenital")] = prob_congen_chronic
  prob_inf[which(prob_inf$prob_congenital<0.000001),c("prob_congenital")] = 0
  prob_inf$time_since_last_infection = round(prob_inf$time_since_last_infection,2)

  
  #family interactions
  ave_mom_child = as.numeric(read.table("ave_mom_child.txt"))*365
  
  #toilet training
  fit_vals_toilet = as.numeric(read.table("toilet_training_gompertz.txt"))
  
  gompertz <- function(time, a, mu, lambda){
    y <- a*exp(-exp(mu*exp(1)/a*(lambda-time)+1))
    return(data.frame(age=time, frac=y))
  }
  diapers <- gompertz(seq(0,(5*12), delta_t*12), fit_vals_toilet[1], fit_vals_toilet[2], fit_vals_toilet[3])
  
  diapers[(diapers$frac>1),"frac"]=1
  prob_of_convert = 1-(1-diapers$frac[2:nrow(diapers)])/(1-diapers$frac[1:(nrow(diapers)-1)])
  prob_of_convert[is.na(prob_of_convert)]=1
  diapers = data.frame(diapers[1:(nrow(diapers)-1),], prob_of_convert = c(prob_of_convert))
  diapers$age = diapers$age/12
  diapers$age = round(diapers$age,2)
  
  #breastfeeding
  fit_vals_breast = as.numeric(read.table("breast_feeding_exponential.txt"))
  
  exponential <- function(time, a,b){
    y <- b*exp(-a*time)
    return(data.frame(age=time, frac=y))
  }
  
  breast = exponential(seq(0,(5*12), delta_t*12), fit_vals_breast[1], fit_vals_breast[2])
  
  breast[(breast$frac>1),"frac"]=1
  prob_of_convert = 1-(breast$frac[2:nrow(breast)])/(breast$frac[1:(nrow(breast)-1)])
  prob_of_convert[is.na(prob_of_convert)]=1
  breast = data.frame(breast[1:(nrow(breast)-1),], prob_of_convert = c(prob_of_convert))
  breast$age = breast$age/12
  breast_prob_to_convert = breast$prob_of_convert[1]
  
  #infectivity of each individual to an uninfected individual
  infectivity_fun = function(data){
    infectivity_info = data.frame(data, age_of_infection = (data$age-data$time_since_last_infection))
    infectivity_info$age_of_infection = round(infectivity_info$age_of_infection,2)
    infectivity_info = merge(infectivity_info,infectivity_params,by = "age_of_infection",all.x = TRUE, all.y = FALSE)
    infectivity_info = data.frame(infectivity_info,
                                  infectivity = infectivity_info$max.infectivity*exp(log(params$chronic_infectivity_val/infectivity_info$max.infectivity)/
                                                                                       infectivity_info$time_to_chronic*
                                                                                       (infectivity_info$time_since_last_infection)))
    infectivity_info[is.nan(infectivity_info$infectivity),"infectivity"]=0
    infectivity_info$time_since_last_infection = round(infectivity_info$time_since_last_infection,3)
    infectivity_info = merge(infectivity_info,immune_protection, by = c("time_since_last_infection","infection_status"),all.x = TRUE, all.y = FALSE)
    #individuals who are experiencing 2ndary infection aren't as infectious - scale 
    infectivity_info$infectivity[infectivity_info$infection_status==3] = infectivity_info$infectivity[infectivity_info$infection_status==3]*params$secondary_factor_val
    #now correct infectivity so that infected individuals have infectivity no less than the chronic amount
    #and so that uninfected individuals have no infectivity 
    infectivity_info$infectivity = rowMaxs(as.matrix(data.frame(infectivity = infectivity_info$infectivity,min = params$chronic_infectivity_val)))
    infectivity_info$infectivity[infectivity_info$infection_status==1]=0
    #add 5 year age categories
    infectivity_info = data.frame(infectivity_info,age_round = RoundTo(infectivity_info$age,multiple = 5, FUN = floor))
    infectivity_info = infectivity_info[,c("person_id","family_id","age_round","age_class","gender","infection_status","breast_status","diaper_status","infectivity","immune_protection","time_since_last_infection")]
    return(infectivity_info)
  }
  
  #Use init data to figure out contacts
  contact = read.table("contact_edited3BRA.txt",header  = TRUE)
  contact = unique(contact[,c(1,2,5,6,7)])
  contact[c("male_interact","female_interact")] = contact[c("male_interact","female_interact")]*365
  
  baby_stat_total = matrix(0,ncol = 5,nrow = 3000)
  baby_stat_total = as.data.frame(baby_stat_total)
  colnames(baby_stat_total) = c("time","mother_status","mother_time_since_last_infection","prob_congenital","baby_status")
  baby_count = 0
  birth_binary=0
  
  #read in initial data
  init_data = read.csv("good_initial_fit_data_10000_new.csv",header = TRUE)
  data = init_data

  #total number of people
  person_id = max(data$person_id)
  family_id = max(data$family_id)

  infectivity_info = infectivity_fun(data)
  
  system.time({
    
    for (n in 1:(70*(1/delta_t))){
      
      new_birth_data <- NULL
      
      #POPULATION-LEVEL INFECTION
      #Use infectivity info to calculate the rate at which population-level infection happens
      #calculate average infectivity for each age group and each gender
      contact_info = as.data.frame(infectivity_info%>%group_by(age_round,gender)%>%dplyr::summarize(average_infectivity = mean(infectivity)))
      #have seperate column for male and female infectivity
      contact_info = dcast(contact_info, age_round~gender, value.var="average_infectivity")
      colnames(contact_info) = c("cont.age.cat","ave_infectivity_male_cont","ave_infectivity_female_cont")
      #merge these infectivities with the contact matrix describing how many contacts there are 
      contact_info = merge(contact_info,contact,by = "cont.age.cat",all.x = TRUE, all.y = FALSE)
      #get rate of infection by men and women by multiplying rate of contact by infectivity of each contact
      contact_info$male_infectivity = contact_info$male_interact*contact_info$ave_infectivity_male_cont
      contact_info$female_infectivity = contact_info$female_interact*contact_info$ave_infectivity_female_cont
      
      #sum rates of infection for each type of individual 
      contact_info = as.data.frame(contact_info%>%group_by(part.age.cat,part.gender)%>%dplyr::summarize(pop_inf_rate = sum(male_infectivity+female_infectivity)))
      #combine this with data - get rate of pop infection for each individual
      pop_inf_rate = infectivity_info[,c("person_id","family_id","age_round","gender","infection_status","immune_protection","time_since_last_infection")]
      colnames(pop_inf_rate) = c("person_id","family_id","part.age.cat","part.gender","part.inf.status","part.immune.protection","time_since_last_infection")
      pop_inf_rate = merge(pop_inf_rate,contact_info,by = c("part.age.cat","part.gender"),all.x = TRUE)
      #now adjust for participant's infection status so that people who are currently experiencing primary infectio cannot get infected 
      #and so that people who are chronically/reinfected/waning infected have resistance
      pop_inf_rate$pop_inf_rate[pop_inf_rate$part.inf.status%in%c(2)]=0
      pop_inf_rate$pop_inf_rate[pop_inf_rate$part.inf.status==3]= pop_inf_rate$pop_inf_rate[pop_inf_rate$part.inf.status==3]*(1-pop_inf_rate$part.immune.protection[pop_inf_rate$part.inf.status==3])
      
      transmission_rates = pop_inf_rate
      
      
      #WITHIN-FAMILY INFECTION
      #all families need to be updated every time step since infectivity changes with age 
      #update rates of events for families that have changed
      person_pairs_data = subset(data,family_id!=1)
      if(nrow(person_pairs_data)>0){
        #only want to include families that have more than one person in them
        person_pairs_info = as.data.frame(person_pairs_data%>%dplyr::group_by(family_id)%>%tally())
        person_pairs_info = subset(person_pairs_info,n>1)
        if(nrow(person_pairs_info)>0){
          person_pairs_data = subset(person_pairs_data,family_id%in%person_pairs_info$family_id)
          person_pairs = as.data.frame(person_pairs_data %>% dplyr::group_by(family_id) %>% dplyr::do(data.frame(t(combn(.$person_id, 2)))))
          colnames(person_pairs) = c("family_id","person_id_V1","person_id_V2")
          #V1 is the infection receiver and V2 is the infection transmitter
          person_pairs = rbind(person_pairs,data.frame(family_id=person_pairs$family_id,
                                                       person_id_V1 = person_pairs$person_id_V2,
                                                       person_id_V2 = person_pairs$person_id_V1))
          merging_info_1 = infectivity_info[,c("person_id","age_class","infection_status","breast_status")]
          colnames(merging_info_1) = paste0(colnames(merging_info_1),"_V1")
          merging_info_2 = infectivity_info[,c("person_id","infection_status","infectivity","age_class","diaper_status")]
          colnames(merging_info_2) = paste0(colnames(merging_info_2),"_V2")
          person_pairs = merge(person_pairs,merging_info_2,by = "person_id_V2")
          person_pairs = merge(person_pairs,merging_info_1, by = "person_id_V1")
          #assume siblings are not infectious to eachother - remove sibling pairs
          person_pairs = subset(person_pairs,!(age_class_V2==1&age_class_V1==1))
          person_pairs$infectivity_V2 = person_pairs$infectivity_V2*ave_mom_child
          #a mother's infectivity to her child will depend on whether the child is still breastfeeding
          person_pairs$infectivity_V2[person_pairs$age_class_V2==3&person_pairs$breast_status_V1==1]= 
            person_pairs$infectivity_V2[person_pairs$age_class_V2==3&person_pairs$breast_status_V1==1]*params$breast_factor_val
          person_pairs$infectivity_V2[person_pairs$age_class_V2==3&person_pairs$breast_status_V1==0]= 0
          #a child's infectivity to their mother will depend on whether the child is still in diapers 
          person_pairs$infectivity_V2[person_pairs$age_class_V2==1&person_pairs$diaper_status_V2==1&person_pairs$age_class_V1==3] =
            person_pairs$infectivity_V2[person_pairs$age_class_V2==1&person_pairs$diaper_status_V2==1&person_pairs$age_class_V1==3]*params$diaper_factor_val
          person_pairs$infectivity_V2[person_pairs$age_class_V2==1&person_pairs$diaper_status_V2==0&person_pairs$age_class_V1==3] =0
          #add all infectivities from each contact together
          person_pairs = as.data.frame(person_pairs%>%group_by(person_id_V1,infection_status_V1)%>%
                                         dplyr::summarize(fam_inf_rate = sum(infectivity_V2)))
          colnames(person_pairs)[1:2] = c("person_id","infection_status")
          transmission_rates = merge(transmission_rates,person_pairs[,c("person_id","fam_inf_rate")],by = "person_id",all = TRUE)
          transmission_rates$fam_inf_rate[is.na(transmission_rates$fam_inf_rate)] = 0
          transmission_rates$fam_inf_rate[transmission_rates$part.inf.status==2]=0
          transmission_rates$fam_inf_rate[transmission_rates$part.inf.status==3]=transmission_rates$fam_inf_rate[transmission_rates$part.inf.status==3]*(1-transmission_rates$part.immune.protection[transmission_rates$part.inf.status==3])
          
          
        }else{
          transmission_rates = data.frame(transmission_rates,fam_inf_rate=0)
        }
      }else{
        transmission_rates = data.frame(transmission_rates,fam_inf_rate=0)
      }
      
      # #CHOOSE INFECTIONS TO OCCUR
      event_long = as.data.frame(transmission_rates[,c("person_id","family_id","pop_inf_rate","fam_inf_rate")])
      
      event_long = arrange(event_long,
                           pop_inf_rate,
                           fam_inf_rate)
      
      event_stat = ddply(event_long,colnames(event_long[,3:ncol(event_long)]),nrow)
      num_people = event_stat$V1
      event_stat = event_stat[,1:(ncol(event_stat)-1)]
      
      event_data = NULL
      for(i in 1:nrow(event_stat)){
        ra_sub = as.numeric(event_stat[i,])
        (choice = reulermultinom(num_people[i],1,ra_sub,delta_t))
        chosen = which(t(choice)==1, arr.ind=TRUE)
        if(nrow(chosen)>0){
          index_chosen = sum(num_people[1:i])-num_people[i]+chosen[,1]
          family_person_chosen = event_long[index_chosen,c("family_id","person_id")]
          event_chosen = colnames(event_stat)[chosen[,2]]
          event_data = rbind(event_data,data.frame(family_person_chosen, event = event_chosen))
        }
      }
      
      #CHOOSE BIRTHS TO OCCUR
      birth_long =data[,c("person_id","family_id","gender","age","infection_status")]
      birth_long$age = RoundTo(birth_long$age,multiple = 1, FUN = floor)
      birth_long = merge(birth_long,birth,by = c("age","gender"),all.x = TRUE,all.y = FALSE)
      
      birth_long = arrange(birth_long,rate_birth)
      
      birth_stat = ddply(birth_long,c("rate_birth"),nrow)
      event_stat = ddply(event_long,colnames(event_long[,3:ncol(event_long)]),nrow)
      num_people = birth_stat$V1
      birth_stat = birth_stat[,1:(ncol(birth_stat)-1)]
      for(i in 1:length(birth_stat)){
        ra_sub = as.numeric(birth_stat[i])
        choice = reulermultinom(num_people[i],1,ra_sub,delta_t)
        chosen = which(t(choice)==1, arr.ind=TRUE)
        if(nrow(chosen)>0){
          index_chosen = sum(num_people[1:i])-num_people[i]+chosen[,1]
          family_person_chosen = birth_long[index_chosen,c("family_id","person_id")]
          event_data = rbind(event_data,data.frame(family_person_chosen,event = "birth"))
        }
      }
      
      #CHOOSE STOP BREASTFEEDING
      breast_children = subset(data,breast_status==1)
      if(nrow(breast_children)>0){
        breast_children = data.frame(breast_children,prob = breast_prob_to_convert)
        data[data$person_id%in%breast_children$person_id,"breast_status"] = rbinom(nrow(breast_children),1,(1-breast_prob_to_convert))
      }
      
      #CHOOSE POTTY TRAINING
      diaper_children = subset(data,diaper_status==1)
      diaper_children$age = round(diaper_children$age,2)
      if(nrow(diaper_children)>0){
        diaper_children = merge(diaper_children,diapers,by = "age",all.x = TRUE,all.y = FALSE)
        diaper_children = diaper_children[,c("person_id","age","prob_of_convert")]
        diaper_children = arrange(diaper_children,prob_of_convert)
        event_stat = ddply(diaper_children,"prob_of_convert",nrow)
        num_people = event_stat$V1
        event_stat = event_stat[,1:(ncol(event_stat)-1)]
        diaper_person = NULL
        for(i in 1:length(event_stat)){
          ra_sub = as.numeric(event_stat[i])
          choice = rbinom(num_people[i],1,ra_sub)
          chosen = which((choice)==1)
          if(length(chosen)>0){
            index_chosen = sum(num_people[1:i])-num_people[i]+chosen
            diaper_person = c(diaper_person,diaper_children[index_chosen,c("person_id")])
          }
        }
        if(length(diaper_person)>0){
          data[data$person_id%in%diaper_person,"diaper_status"]=0
        }
      }
      
      
      if(!is.null(event_data)){
        
        infection_status_increase = subset(event_data,event%in%c("pop_inf_rate","fam_inf_rate"))$person
        infection_status_increase = subset(data,person_id%in%infection_status_increase)
        index_1 = data$person_id%in%infection_status_increase$person_id&data$infection_status==1
        data[index_1,c("infection_status","last_infection","time_since_last_infection","time_since_first_infection")] = data.frame(infection_status=2,last_infection=1,time_since_last_infection=-delta_t,time_since_first_infection = -delta_t)
        index_3 = data$person_id%in%infection_status_increase$person_id&data$infection_status==3
        data[index_3,c("infection_status","last_infection","time_since_last_infection")] = data.frame(infection_status=3,last_infection=3,time_since_last_infection=-delta_t)
        
        birth_change = subset(event_data,event=="birth")
        if(nrow(birth_change)>0){
          #make sure all people who were chosen to give birth are in age class 3 (mothers)
          data[data$person_id%in%birth_change$person,c("age_class")]=3
          #look at the data of these people
          current_people_data = data[data$person_id%in%birth_change$person,]
          #some women might be new mothers - update family number of these people
          family_update = subset(current_people_data,family_id==1)$person_id
          if(length(family_update)>0){
            data[data$person_id%in%family_update,c("family_id")] = (family_id+1):(family_id+length(family_update)) #women have now been assigned a new family number
            family_id = family_id+length(family_update) #update family_id count
            current_people_data = data[data$person_id%in%birth_change$person,] #update our info on the women having children to include their new family ids
          }
          data_baby = NULL
          #now choose infection status of new babies
          #first have uninfected mothers - these will all have no risk of transmission
          uninfected_mothers= subset(current_people_data,infection_status==1)
          if (nrow(uninfected_mothers)>0){
            data_baby = rbind(data_baby,
                              data.frame(family_id =uninfected_mothers$family_id,
                                         person_id = (person_id+1):(person_id+nrow(uninfected_mothers)),
                                         age_class= 1,
                                         age = -delta_t,
                                         infection_status = 1,
                                         last_infection = 1,
                                         time_since_last_infection = -delta_t,
                                         time_since_first_infection = NA))
            person_id = person_id+nrow(uninfected_mothers)
            #update info on births to be used in summary statistics
            new_birth_data = rbind(new_birth_data,data.frame(time = rep((t+delta_t),nrow(uninfected_mothers)),
                                                             mother_status =1,
                                                             mother_time_since_last_infection = NA,
                                                             prob_congenital = 0,
                                                             baby_status = 1))
          }
          #want to figure out the maximum probability of congenital infection mothers have had in the last 9 months
          infected_mothers = subset(current_people_data,infection_status%in%c(2,3))
          if(nrow(infected_mothers)>0){
            infected_mothers$infection_status=2
            infected_mothers$time_since_last_infection = infected_mothers$time_since_first_infection-9/12
            infected_mothers$time_since_last_infection = ifelse(infected_mothers$time_since_last_infection<0,0,infected_mothers$time_since_last_infection)
            
            reinfected_mothers = subset(current_people_data,infection_status==3)
            if(nrow(reinfected_mothers)>0){
              reinfected_mothers$time_since_last_infection = reinfected_mothers$time_since_last_infection-9/12
              reinfected_mothers$time_since_last_infection = ifelse(reinfected_mothers$time_since_last_infection<0,0,reinfected_mothers$time_since_last_infection)
              infected_mothers = rbind(infected_mothers,reinfected_mothers)[,c("person_id","family_id","time_since_last_infection","infection_status")]
              infected_mothers$time_since_last_infection = round(infected_mothers$time_since_last_infection,2)
              infected_mothers = merge(infected_mothers,prob_inf,by = c("time_since_last_infection","infection_status"),all.x = TRUE,all.y = FALSE)
              infected_mothers2 = as.data.frame(infected_mothers%>%dplyr::group_by(person_id,family_id)%>%dplyr::summarize(prob_congenital = max(prob_congenital)))
              infected_mothers = merge(infected_mothers2,infected_mothers,by = colnames(infected_mothers2),all.x = TRUE,all.y = FALSE)
              infected_mothers3 = as.data.frame(infected_mothers%>%dplyr::group_by(person_id,family_id,prob_congenital)%>%dplyr::summarize(infection_status = max(infection_status)))
              infected_mothers = merge(infected_mothers3,infected_mothers,all.x = TRUE,all.y = FALSE)
            }else{
              infected_mothers = infected_mothers[,c("person_id","family_id","time_since_last_infection","infection_status")]
              infected_mothers$time_since_last_infection = round(infected_mothers$time_since_last_infection,2)
              infected_mothers = merge(infected_mothers,prob_inf,by = c("time_since_last_infection","infection_status"),all.x = TRUE,all.y = FALSE)
              infected_mothers = infected_mothers[,c("person_id","family_id","prob_congenital","infection_status","time_since_last_infection")]
            }
            current_people_data = infected_mothers
            if(nrow(current_people_data)>0){
            for (i in 1:nrow(current_people_data)){
              rate = current_people_data$prob_congenital[i]
              if(!is.numeric(rate))
                write.table(data.frame(column = i, current_people_data),"birth_issue.csv",sep = ",",col.names = FALSE,append = TRUE)
              baby_inf = as.numeric(rbinom(1,1,rate))
              data_baby = rbind(data_baby,
                                data.frame(family_id =current_people_data$family_id[i],
                                           person_id = (person_id+1),
                                           age_class= 1,
                                           age = -delta_t,
                                           infection_status = baby_inf+1,
                                           last_infection = baby_inf+1,
                                           time_since_last_infection = -delta_t,
                                           time_since_first_infection = c(NA,-delta_t)[baby_inf+1]))
              person_id = person_id+1
              #update info on births to be used in summary statistics
              
            }
            new_birth_data = rbind(new_birth_data,data.frame(time = t+delta_t,
                                                             mother_status = current_people_data$infection_status,
                                                             mother_time_since_last_infection = current_people_data$time_since_last_infection,
                                                             prob_congenital = current_people_data$prob_congenital,
                                                             baby_status = tail(data_baby,nrow(current_people_data))$infection_status))
            }
          }
          
          
          data = rbind(data,data.frame(data_baby[,1:5],
                                       breast_status = rbinom(nrow(data_baby),1,fit_vals_breast[2]),
                                       diaper_status=1,
                                       gender = sample(c(0,1),nrow(data_baby),TRUE),data_baby[,6:8]))
        }
      }
      
      
      data[,c("age","time_since_last_infection","time_since_first_infection")] = data[,c("age","time_since_last_infection","time_since_first_infection")]+delta_t
      
      
      #age children into adults
      #if this makes mothers childless, put mother back in adult class
      children = subset(data,age_class==1&round(age,2)>=5)
      if (nrow(children)>0){
        num_children_family = as.data.frame(subset(data,family_id%in%children$family_id&age_class==1)%>%
                                              group_by(family_id)%>%dplyr::tally()) #count current number of children - we are about to lose 1
        data[data$person_id%in%children$person_id,c("family_id")]=1 #put child into family_id=1
        data[data$person_id%in%children$person_id,c("age_class")]=2 #age child into age class 2 (adult)
        data[data$person_id%in%children$person_id,c("breast_status")]=0 #make sure child is no longer breastfeeding
        data[data$person_id%in%children$person_id,c("diaper_status")]=0 #make sure child is no longer is diapers
        no_children = subset(num_children_family,n==1) #families that will no longer have children
        if(nrow(no_children)>0){
          data[data$family_id%in%no_children$family_id,c("age_class")]=2 #put mother back into age class 2
          data[data$family_id%in%no_children$family_id,c("family_id")]=1 #put mother back into family id 1
        }
      }
      
      #Age adults into death
      adult_age = subset(data,age>50)$person_id
      if(length(adult_age)>0){
        data[data$person_id%in%adult_age,c("age_class")]=4
      }
      
      data = subset(data,age_class<=3)
      
      #put individuals in reinf/immune waning class based on their time_to_wane
      
      waning_people = subset(data,infection_status==2)
      if (nrow(waning_people)>0){
        waning_people = data.frame(person_id = waning_people$person_id,
                                   infection_status=2,
                                   time_since_last_infection = waning_people$time_since_last_infection, 
                                   age_of_infection = waning_people$age-waning_people$time_since_last_infection)
        if(nrow(waning_people)>0){
          waning_people$age_of_infection =round(waning_people$age_of_infection,3)
          waning_info = max_infectivity[,c(1,4)]
          waning_info$age_of_infection = round(waning_info$age_of_infection,3)
          waning_people = merge(waning_people,waning_info,by = "age_of_infection",all.x = TRUE,all.y = FALSE)
          waning_people = waning_people[waning_people$time_since_last_infection>=(waning_people$time_to_wane),]
          data[data$person_id%in%waning_people$person_id,"infection_status"]=3
          data[data$person_id%in%waning_people$person_id,"time_since_last_infection"]=delta_t
          data[data$person_id%in%waning_people$person_id,"last_infection"]=2
        }
      }
      
      
      #infectivity of each individual to an uninfected individual
      infectivity_info = infectivity_fun(data)
      
      if (sum(data$infection_status<1)>0){
        print("infection status wrong")
        stop()
      }
      
      
      
      t = t+delta_t
      #print(t)
      
      if(sum(data$infection_status==6)>0){
        print("stop")
        stop()
      }
      
      if(sum(data$time_since_last_infection>data$age)>0){
        print("stop")
        stop()
      }
      
      if(round(t,2)>=50){
        if(round(t,2)==50){
          write.table(data.frame(data,time = t), paste0("sim_final",rand,".csv"), sep = ",",col.names = TRUE,row.names = FALSE, append = FALSE)
        }else{
          write.table(data.frame(data,time = t), paste0("sim_final",rand,".csv"), sep = ",", col.names = FALSE,row.names = FALSE, append = TRUE)
        }
        
        
        if(!is.null(nrow(new_birth_data))){
          baby_stat_total[(baby_count+1):(baby_count+nrow(new_birth_data)),] = data.frame(new_birth_data)
          baby_count = baby_count+nrow(new_birth_data)
        }
        
      }
    }
  })
  
  #Summary Statistics
  baby_stat_total = subset(baby_stat_total,mother_status!=0)
  data_total = read.csv(paste0("sim_final",rand,".csv"))
  
  #Want to know the rate at which people seroconvert as a function of age (percent_uninfected)
  data_total$age_rounded = floor(data_total$age)
  
  percent_uninfected = subset(as.data.frame(data_total%>%group_by(age_rounded,infection_status)%>%pct_routine()),infection_status==1)
  percent_all = data.frame(age_rounded =0:49, infection_status = 1)
  percent_uninfected = as.matrix(merge(percent_uninfected, percent_all,by = c("age_rounded", "infection_status"), all = TRUE))
  percent_uninfected[is.na(percent_uninfected)] = 0
  percent_uninfected = as.data.frame(percent_uninfected)
  percent_uninfected = arrange(percent_uninfected,age_rounded)

  #Want to know how many infants were born with congenital CMV
  baby_stat_total = subset(baby_stat_total, baby_status!=0)
  congenital_percent  = mean(baby_stat_total$prob_congenital)

 
  #Important sim info
  all_residuals = rbind(data.frame(val = 1-percent_uninfected$pct[2:50], test =paste0("infected_",1:49,"_val")),
                        data.frame(val = congenital_percent,test = "congenital_val")
  )
  
  predicted = c(congenital_percent,1-percent_uninfected$pct[2:50])
  actual = c(0.01,NHANES$percent_weighted[1:49])
  
  perc_error = 1/length(actual)*sum(abs(actual-predicted)/actual)

  
  write.table(rbind(all_residuals,data.frame(val = as.numeric(params[1,]), test= colnames(params))),
              paste0("residuals_",rand,".csv"), sep = ",",col.names =TRUE,row.names = FALSE)
  
  #file.remove(paste0("birth_data",rand,".csv"))
  write.table(subset(data_total,round(time,2)==70)[,1:12], paste0("sim_final",rand,".csv"), sep = ",", col.names = T,row.names = F, append = F)
  
  print(paste0(rand," all good"))
  return(c(perc_error,rand))
  
}













