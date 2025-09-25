library(Hmisc)
library(parallel)
library(MASS)
library(plyr)
library(dplyr)
library(ggplot2)
library(pomp)
library(reshape2)
#library(tidyverse)
library(tidyr)
library(matrixStats)
#library(DescTools)
library(extdplyr)
library(data.table)

cmv_sim = function(input){
  input = as.numeric(input)
  params = input[1:8]
  sim_num = input[9]
  vacc_waning = input[10]
  vacc_const = input[11]
  vacc_prob = input[12]
  age_vacc = input[13]
  booster_times = input[14]
  rand = input[15]
  
  
  setwd("/home/sgazea01/projects/def-craigm/sgazea01/LowSero1")
  #setwd("C:/Users/c-byr.LAPTOP-UG0F9R6R/OneDrive - math.ubc.ca/cmv/cmv transmission model/POLYMOD/original contact matrix, no daycare/CMV_transmission_model40_2")
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
  
  t = -5
  delta_t = 1/(12) #set delta_t = 2 weeks
  
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
  params = as.data.table(t(params))
  colnames(params) = c("adult_infectivity_val",
                       "child_infectivity_val",
                       "chronic_infectivity_val",
                       "breast_factor_val",
                       "diaper_factor_val",
                       "secondary_factor_val",
                       "immune_waning",
                       "immune_const",
                       "duration_child_val")
  
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
  #0-5 y.o. = 24 months = 1.5 years
  #10 y.o.+ = 2 months = 1/6 years
  #exponential decay inbetween
  r = -log((1/6)/duration)/(10-5)
  chronic = data.table(age_of_infection = seq(0,5,delta_t),time_to_chronic = duration)
  for(i in seq((5+delta_t),10,delta_t)){
    chronic = rbind(chronic,data.table(age_of_infection = i, time_to_chronic = duration*exp(-r*(i-5))))
  }
  chronic = rbind(chronic,data.table(age_of_infection = seq((10+delta_t),(50-delta_t),delta_t), time_to_chronic = (1/6)))
  
  
  #TRANSMISSIBILITY PARAMS
  #max infectivity for every age
  #age 0-5 have infectivity of child_infectivity_val
  #age 10+ have infectivity of adult_infectivity_val
  #age 5-10 have exponential decay of infectivity 
  r = -log(params$adult_infectivity_val/params$child_infectivity_val)/(10-5)
  max_infectivity = data.table(age_of_infection = seq(0,5,delta_t),max.infectivity =rep(params$child_infectivity_val))
  for(i in seq((5+delta_t),10,delta_t)){
    max_infectivity = rbind(max_infectivity,data.table(age_of_infection = i, max.infectivity = params$child_infectivity_val*exp(-r*(i-5))))
  }
  max_infectivity = rbind(max_infectivity,data.table(age_of_infection =  seq((10+delta_t),(50-delta_t),delta_t),max.infectivity = params$adult_infectivity_val))
  max_infectivity$max.infectivity[is.nan(max_infectivity$max.infectivity)]=0
  infectivity_params = data.table(chronic,max.infectivity=max_infectivity$max.infectivity)
  infectivity_params$age_of_infection = round(infectivity_params$age_of_infection,2)
  
  #IMMUNE WANING
  #want immune waning to happen once individuals reach the same transmissibility as that seen during a reinfection
  
  max_infectivity$reinfection.infectivity = infectivity_params$max.infectivity*params$secondary_factor_val
  max_infectivity$time_to_wane = log(max_infectivity$reinfection.infectivity/max_infectivity$max.infectivity)/(log(params$chronic_infectivity_val/infectivity_params$max.infectivity)/
                                                                                                                 chronic$time_to_chronic)
  #let immune protection be the probability of reinfection with time_since_last_infection=0 being when waning starts 
  r_regular = -params$immune_waning
  r_vacc = -vacc_waning
  immune_protection = data.table(time_since_last_infection = rep(seq(0,50,delta_t),2),infection_status=rep(c(1,2),each = length(seq(0,50,delta_t))),immune_protection = rep(c(0,1),each = length(seq(0,50,delta_t))))
  for (i in seq(0,50,delta_t)){
    immune_protection = rbind(immune_protection,data.table(time_since_last_infection = i,infection_status=c(3,4,5), immune_protection = c(exp(r_regular*i),exp(r_vacc*i),exp(r_vacc*i))))
  }
  immune_protection$immune_protection[which(immune_protection$immune_protection<params$immune_const&immune_protection$infection_status==3)]=params$immune_const
  immune_protection$immune_protection[which(immune_protection$immune_protection<vacc_const&immune_protection$infection_status%in%c(4,5))]=params$immune_const
  
  immune_protection$time_since_last_infection = round(immune_protection$time_since_last_infection,3)
  
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
  for (i in seq(0,100,delta_t)){
    prob_inf = rbind(prob_inf,data.table(infection_status = c(2,3), time_since_last_infection = i, prob_congenital = c((0.334*exp(-r*i)),(0.334*exp(-r*i)*params$secondary_factor_val))))
  }
  if(params$adult_infectivity_val==0)
    prob_inf$prob_congenital = 0
  prob_inf[which(prob_inf$prob_congenital<=prob_congen_chronic),c("prob_congenital")] = prob_congen_chronic
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
  diapers = data.table(diapers[1:(nrow(diapers)-1),], prob_of_convert = c(prob_of_convert))
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
  breast = data.table(breast[1:(nrow(breast)-1),], prob_of_convert = c(prob_of_convert))
  breast$age = breast$age/12
  breast_prob_to_convert = breast$prob_of_convert[1]
  
  #infectivity of each individual to an uninfected individual
  infectivity_fun = function(vaccinated_infected_people,data2){
    #calculate immunity using new vaccinated states
    immunity_info  = data.table(data2, age_of_infection = (data2$age-data2$time_since_last_infection))
    immunity_info$time_since_last_infection = round(immunity_info$time_since_last_infection,3)
    immunity_info = merge(immunity_info,immune_protection, by = c("time_since_last_infection","infection_status"),all.x = TRUE, all.y = FALSE)
    immunity_info = immunity_info[,c("person_id","immune_protection")]
    
    #calculate infectivity from statuses before vaccination. - assume vaccination does not impact trajectory of infectivity
    data = data2
    if(!is.null(vaccinated_infected_people)){
      chosen_vaccinated = vaccinated_infected_people[vaccinated_infected_people$person_id%in%data$person_id,]
      chosen_vaccinated = chosen_vaccinated[order(-person_id,decreasing=TRUE),]
      data[data$person_id%in%vaccinated_infected_people$person_id,c("infection_status","last_infection","time_since_last_infection","time_since_first_infection")]= chosen_vaccinated[,c("infection_status","last_infection","time_since_last_infection","time_since_first_infection")]
    }
    infectivity_info = data.table(data, age_of_infection = (data$age-data$time_since_last_infection))
    infectivity_info$age_of_infection = round(infectivity_info$age_of_infection,2)
    infectivity_info = merge(infectivity_info,infectivity_params,by = "age_of_infection",all.x = TRUE, all.y = FALSE)
    infectivity_info = data.table(infectivity_info,
                                  infectivity = infectivity_info$max.infectivity*exp(log(params$chronic_infectivity_val/infectivity_info$max.infectivity)/
                                                                                       infectivity_info$time_to_chronic*
                                                                                       (infectivity_info$time_since_last_infection)))
    infectivity_info[is.nan(infectivity_info$infectivity),"infectivity"]=0
    #individuals who are experiencing 2ndary infection aren't as infectious - scale 
    infectivity_info$infectivity[infectivity_info$infection_status==3] = infectivity_info$infectivity[infectivity_info$infection_status==3]*params$secondary_factor_val
    #now correct infectivity so that infected individuals have infectivity no less than the chronic amount
    infectivity_info$infectivity = rowMaxs(as.matrix(data.table(infectivity = infectivity_info$infectivity,min = params$chronic_infectivity_val)))
    #uninfected individuals or vaccinated-uninfected individuals have no infectivity 
    infectivity_info$infectivity[infectivity_info$infection_status%in%c(1,4)]=0
    #add 5 year age categories
    infectivity_info = data.table(infectivity_info,age_round = RoundTo(infectivity_info$age,multiple = 5, FUN = floor))
    infectivity_info = infectivity_info[,c("person_id","family_id","age_round","age_class","gender","infection_status","breast_status","diaper_status","infectivity","time_since_last_infection")]
    infectivity_info = merge(infectivity_info,immunity_info,all = TRUE)
    infectivity_info = infectivity_info[order(person_id)]
    infectivity_inf =infectivity_info[,c(1:9,11,10)]
    
    return(infectivity_info)
  }
  
  #Use init data to figure out contacts
  contact = read.table("contact_edited3BRA.txt",header  = TRUE)
  contact = unique(contact[,c(1,2,5,6,7)])
  contact[c("male_interact","female_interact")] = contact[c("male_interact","female_interact")]*365
  contact = as.data.table(contact)
  
  #determine ages that people are to be vaccinated
  if (booster_times!=0&r_vacc!=0)
    for (i in 1:booster_times)
      age_vacc = c(age_vacc,(age_vacc[1]-log(2)/r_vacc*i))
  age_vacc=age_vacc*(1/delta_t) #want age_vacc to occur at a time that matches with time steps
  age_vacc = round(age_vacc)
  age_vacc = age_vacc*delta_t
  age_vacc = round(age_vacc,2)
  age_vacc = unique(age_vacc)
  
  #read in initial data
  
  data = fread(paste0("sim_final",sim_num,".csv"),header =TRUE)[,c(1:11)]
  data =data[order(-person_id,decreasing=TRUE),]
  data$vaccinated=0
  Reff = NULL
  new_infections = NULL
  
  vaccinated_list = NULL
  vaccinated_infected_people = NULL
  
  #total number of people
  person_id = max(data$person_id)
  family_id = max(data$family_id)
  
  infectivity_info = infectivity_fun(vaccinated_infected_people,data)
  
  
  all_mothers = NULL
  source_infection = NULL
  child_data = NULL
  
  system.time({
    
    for (n in 1:(60*(1/delta_t)-1)){
      
      #POPULATION-LEVEL INFECTION
      #Use infectivity info to calculate the rate at which population-level infection happens
      #calculate average infectivity for each age group and each gender
      contact_info = infectivity_info[,.(average_infectivity = mean(infectivity)),by = list(age_round,gender)]
      #have seperate column for male and female infectivity
      contact_info = dcast(contact_info, age_round~gender, value.var="average_infectivity")
      colnames(contact_info) = c("cont.age.cat","ave_infectivity_male_cont","ave_infectivity_female_cont")
      #merge these infectivities with the contact matrix describing how many contacts there are 
      contact_info = merge(contact_info,contact,by = "cont.age.cat",all.x = TRUE, all.y = FALSE)
      #get rate of infection by men and women by multiplying rate of contact by infectivity of each contact
      contact_info$male_infectivity = contact_info$male_interact*contact_info$ave_infectivity_male_cont
      contact_info$female_infectivity = contact_info$female_interact*contact_info$ave_infectivity_female_cont
      #sum rates of infection for each type of individual 
      contact_info = contact_info[,.(pop_inf_rate = sum(male_infectivity+female_infectivity)),by = list(cont.age.cat,part.age.cat,part.gender)]
      #combine this with data - get rate of pop infection for each individual
      pop_inf_rate = infectivity_info[,c("person_id","family_id","age_round","gender","infection_status","immune_protection","time_since_last_infection")]
      colnames(pop_inf_rate) = c("person_id","family_id","part.age.cat","part.gender","part.inf.status","part.immune.protection","time_since_last_infection")
      pop_inf_rate = merge(pop_inf_rate,contact_info,by = c("part.age.cat","part.gender"),all.x = TRUE,allow.cartesian=TRUE)
      #now adjust for participant's infection status so that people who are currently experiencing primary infectio cannot get infected 
      #and so that people who are chronically/reinfected/waning infected have resistance
      pop_inf_rate$pop_inf_rate[pop_inf_rate$part.inf.status%in%c(2)]=0
      pop_inf_rate$pop_inf_rate[pop_inf_rate$part.inf.status%in%c(3,4,5)]= pop_inf_rate$pop_inf_rate[pop_inf_rate$part.inf.status%in%c(3,4,5)]*(1-pop_inf_rate$part.immune.protection[pop_inf_rate$part.inf.status%in%c(3,4,5)])
      child_inf_rate = pop_inf_rate[cont.age.cat==0,c(1:7,9)]
      child_inf_rate = child_inf_rate[order(person_id)]
      colnames(child_inf_rate)[8] = "child_inf_rate"
      adult_inf_rate = pop_inf_rate[cont.age.cat!=0,.(adult_inf_rate = sum(pop_inf_rate)),by =list(part.age.cat,part.gender,person_id,family_id,part.inf.status,part.immune.protection,time_since_last_infection)]
      adult_inf_rate = adult_inf_rate[order(person_id)]
      transmission_rates = data.table(adult_inf_rate,child_inf_rate=child_inf_rate$child_inf_rate)

      #WITHIN-FAMILY INFECTION
      #all families need to be updated every time step since infectivity changes with age 
      #update rates of events for families that have changed
      person_pairs_data = data[family_id!=1]
      mothers_id = person_pairs_data[age>5]$person_id
      if(nrow(person_pairs_data)>0){
        #only want to include families that have more than one person in them
        person_pairs_info = as.data.table(person_pairs_data%>%dplyr::group_by(family_id)%>%dplyr::tally()%>%dplyr::filter(n>1))
        if(nrow(person_pairs_info)>0){
          person_pairs = as.data.table(person_pairs_data%>%dplyr::filter(family_id%in%person_pairs_info$family_id) %>% dplyr::group_by(family_id) %>% dplyr::do(data.table(t(combn(.$person_id, 2)))))
          colnames(person_pairs) = c("family_id","person_id_V1","person_id_V2")
          #V1 is the infection receiver and V2 is the infection transmitter
          person_pairs = rbind(person_pairs,data.table(family_id=person_pairs$family_id,
                                                       person_id_V1 = person_pairs$person_id_V2,
                                                       person_id_V2 = person_pairs$person_id_V1))
          merging_info_1 = infectivity_info[,c("person_id","age_class","infection_status","breast_status")] #v1 is the person who might be infected (breast status makes more susceptible to v2)
          colnames(merging_info_1) = paste0(colnames(merging_info_1),"_V1")
          merging_info_2 = infectivity_info[,c("person_id","infection_status","infectivity","age_class","diaper_status")] #v2 is the person who is doing the infecting (diaper status makes more infectious to v1)
          colnames(merging_info_2) = paste0(colnames(merging_info_2),"_V2")
          person_pairs = merge(person_pairs,merging_info_2,by = "person_id_V2")
          person_pairs = merge(person_pairs,merging_info_1, by = "person_id_V1")
          #assume siblings are not infectious to eachother - remove sibling pairs
          person_pairs = person_pairs[!(age_class_V2==1&age_class_V1==1)]
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
          person_pairs = person_pairs[,.(fam_inf_rate = sum(infectivity_V2)),list(person_id_V1,infection_status_V1)]
          # person_pairs = as.data.table(person_pairs%>%group_by(person_id_V1,infection_status_V1)%>%
          #                                dplyr::summarise(fam_inf_rate = sum(infectivity_V2)))
          colnames(person_pairs)[1:2] = c("person_id","infection_status")
          transmission_rates = merge(transmission_rates,person_pairs[,c("person_id","fam_inf_rate")],by = "person_id",all = TRUE)
          transmission_rates$fam_inf_rate[is.na(transmission_rates$fam_inf_rate)] = 0
          transmission_rates$fam_inf_rate[transmission_rates$part.inf.status==2]=0
          transmission_rates$fam_inf_rate[transmission_rates$part.inf.status%in%c(3,4,5)]=transmission_rates$fam_inf_rate[transmission_rates$part.inf.status%in%c(3,4,5)]*(1-transmission_rates$part.immune.protection[transmission_rates$part.inf.status%in%c(3,4,5)])
          
        }else{
          transmission_rates = data.table(transmission_rates,fam_inf_rate=0)
        }
      }else{
        transmission_rates = data.table(transmission_rates,fam_inf_rate=0)
      }
      
      age_of_infection = data$age-data$time_since_first_infection
      age_of_infection = age_of_infection[!is.na(age_of_infection)]
      age_of_infection_child =summary(age_of_infection[age_of_infection<5])
      age_of_infection_adult = summary(age_of_infection[age_of_infection>=5])
      age_of_infection = summary(age_of_infection)
      
      #determine Reff info - calculate rate of primary infection and reinfection
      rate_infs = transmission_rates
      rate_infs$family_id[rate_infs$family_id!=1]=2
      rate_infs$part.age.cat[rate_infs$part.age.cat>0]=5
      rate_infs$part.inf.status[rate_infs$part.inf.status%in%c(4)]=1
      rate_infs$part.inf.status[rate_infs$part.inf.status%in%c(3,5)]=2
      uninfs = subset(rate_infs,part.inf.status==1)
      uninfs$inf_rate = uninfs$adult_inf_rate+uninfs$child_inf_rate+uninfs$fam_inf_rate
      uninfs = uninfs[,c("person_id","inf_rate")]
      uninfs = merge(uninfs,data[,c("person_id","age")])
      weighted_age_acquisition = sum(uninfs$inf_rate/sum(uninfs$inf_rate)*uninfs$age)
      num_infected = sum(rate_infs$part.inf.status==2)
      sum_infected = sum(uninfs$inf_rate)
      Reff = rbind(Reff, data.frame(new_per_year = sum_infected, num_infected =num_infected, age_acq=weighted_age_acquisition, time = t))

      
      #CHOOSE INFECTIONS TO OCCUR
      event_long = transmission_rates[,c("person_id","family_id","adult_inf_rate","child_inf_rate","fam_inf_rate")]
      
      event_long = event_long[order(
        adult_inf_rate,
        child_inf_rate,
        fam_inf_rate)]
      
      event_stat = ddply(event_long,colnames(event_long[,3:ncol(event_long)]),nrow)
      num_people = event_stat$V1
      event_stat = event_stat[,1:(ncol(event_stat)-1)]
      
      event_data = NULL
      for(i in 1:nrow(event_stat)){
        ra_sub = as.numeric(event_stat[i,])
        choice = reulermultinom(num_people[i],1,ra_sub,delta_t)
        chosen = which(t(choice)==1, arr.ind=TRUE)
        if(nrow(chosen)>0){
          index_chosen = sum(num_people[1:i])-num_people[i]+chosen[,1]
          family_person_chosen = event_long[index_chosen,c("family_id","person_id")]
          event_chosen = colnames(event_stat)[chosen[,2]]
          event_data = rbind(event_data,data.table(family_person_chosen, event = event_chosen))
        }
      }
      
      #CHOOSE BIRTHS TO OCCUR
      birth_long =data#[,c("person_id","family_id","gender","age","infection_status")]
      birth_long$age = RoundTo(birth_long$age,multiple = 1, FUN = floor)
      birth_long = merge(birth_long,birth,by = c("age","gender"),all.x = TRUE,all.y = FALSE)
      birth_long = birth_long[order(rate_birth)]
      
      birth_stat = ddply(birth_long,c("rate_birth"),nrow)
      event_stat = ddply(event_long,colnames(event_long[,3:ncol(event_long)]),nrow)
      num_people = birth_stat$V1
      birth_stat = birth_stat[,1:(ncol(birth_stat)-1)]
      for(i in 1:length(birth_stat)){
        ra_sub = as.numeric(birth_stat[i])
        (choice = reulermultinom(num_people[i],1,ra_sub,delta_t))
        chosen = which(t(choice)==1, arr.ind=TRUE)
        if(nrow(chosen)>0){
          index_chosen = sum(num_people[1:i])-num_people[i]+chosen[,1]
          family_person_chosen = birth_long[index_chosen,c("family_id","person_id")]
          event_data = rbind(event_data,data.table(family_person_chosen,event = "birth"))
        }
      }
      event_data
      
      #record how likely each woman of child-bearing age is for having a congenitally-infected infant
      mothers = NULL
      possible_birth = birth_long[rate_birth!=0]
      
      if(sum(possible_birth$person_id%in%vaccinated_infected_people$person_id)>0){
        chosen_vaccinated = vaccinated_infected_people[vaccinated_infected_people$person_id%in%possible_birth$person_id,]
        chosen_vaccinated = chosen_vaccinated[order(-person_id,decreasing=TRUE),]
        possible_birth[possible_birth$person_id%in%vaccinated_infected_people$person_id,c("person_id","infection_status","last_infection","time_since_last_infection","time_since_first_infection")]=
          chosen_vaccinated[,c("person_id","infection_status","last_infection","time_since_last_infection","time_since_first_infection")]
      }
      if(nrow(possible_birth)>0){
        uninfected_mothers = possible_birth[infection_status%in%c(1,4)]
        if(nrow(uninfected_mothers)>0){
          uninfected_mothers$prob_congenital = 0
          mothers = uninfected_mothers[,c("person_id","family_id","prob_congenital","infection_status","time_since_last_infection","rate_birth")]
        }
        infected_mothers = possible_birth[infection_status%in%c(2,3)]
        if(nrow(infected_mothers)>0){
          infected_mothers_new = NULL
          infected_mothers = infected_mothers[!is.na(time_since_first_infection)]
          if(nrow(infected_mothers)>0){
            infected_mothers$infection_status=2
            infected_mothers$time_since_last_infection = infected_mothers$time_since_first_infection-9/12
            infected_mothers$time_since_last_infection = ifelse(infected_mothers$time_since_last_infection<0,0,infected_mothers$time_since_last_infection)
            infected_mothers$time_since_last_infection = round(infected_mothers$time_since_last_infection,2)
            infected_mothers_new = infected_mothers[,c("person_id","family_id","time_since_last_infection","infection_status","rate_birth")]
          }
          reinfected_mothers = possible_birth[infection_status%in%c(3)]
          if(nrow(reinfected_mothers)>0){
            reinfected_mothers$time_since_last_infection = reinfected_mothers$time_since_last_infection-9/12
            reinfected_mothers$time_since_last_infection = ifelse(reinfected_mothers$time_since_last_infection<0,0,reinfected_mothers$time_since_last_infection)
            reinfected_mothers$time_since_last_infection = round(reinfected_mothers$time_since_last_infection,2)
            reinfected_mothers = reinfected_mothers[,c("person_id","family_id","time_since_last_infection","infection_status","rate_birth")]
            infected_mothers_new = rbind(reinfected_mothers,infected_mothers_new)
          }
          infected_mothers_new = merge(infected_mothers_new,prob_inf,by = c("time_since_last_infection","infection_status"),all.x = TRUE,all.y = FALSE)
          infected_mothers_new$infection_status[infected_mothers_new$prob_congenital>as.numeric(tail(prob_inf,1)[1,3])&infected_mothers_new$infection_status==2]=20
          infected_mothers_new$infection_status[infected_mothers_new$prob_congenital>as.numeric(tail(prob_inf,1)[1,3])&infected_mothers_new$infection_status==3]=30
          infected_mothers_new$infection_status[infected_mothers_new$infection_status==2]=3
          infected_mothers_new = as.data.table(infected_mothers_new%>%dplyr::group_by(person_id,family_id)%>%dplyr::slice(which.max(prob_congenital)))
          infected_mothers_new = infected_mothers_new[,c("person_id","family_id","prob_congenital","infection_status","time_since_last_infection","rate_birth")]
          if (sum(duplicated(infected_mothers_new$person_id))>0)
            stop()
          mothers = rbind(mothers,infected_mothers_new)
          mothers_chronic = subset(mothers,infection_status==3)
          mothers_not_chronic = subset(mothers,infection_status!=3)
          mothers_chronic = merge(data[,c("person_id","infection_status")],mothers_chronic[,c(1:3,5:6)],all.x = FALSE,all.y = TRUE)
          mothers = rbind(mothers_chronic,mothers_not_chronic)
        }
        if(t>=-5)
          all_mothers = rbind(all_mothers,data.table(mothers,time_of_birth = round(t,2)))
      }
      
      #CHOOSE STOP BREASTFEEDING
      breast_children = data[breast_status==1]
      if(nrow(breast_children)>0){
        breast_children = data.table(breast_children,prob = breast_prob_to_convert)
        data[data$person_id%in%breast_children$person_id,"breast_status"] = rbinom(nrow(breast_children),1,(1-breast_prob_to_convert))
      }
      
      #CHOOSE POTTY TRAINING
      diaper_children = data[diaper_status==1]
      diaper_children$age = round(diaper_children$age,2)
      if(nrow(diaper_children)>0){
        diaper_children = merge(diaper_children,diapers,by = "age",all.x = TRUE,all.y = FALSE)
        diaper_children = diaper_children[,c("person_id","age","prob_of_convert")]
        diaper_children = diaper_children[order(prob_of_convert)]
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
        
        infection_status_increase = event_data[event%in%c("adult_inf_rate","child_inf_rate","fam_inf_rate")]
        addition = infection_status_increase[,c("person_id","event")]
        if (nrow(addition)>0)
          source_infection = rbind(source_infection,data.table(time_of_infection = round(t,2),addition))
        infection_status_increase = infection_status_increase$person
        infection_status_increase = data[person_id%in%infection_status_increase]
        index_1 = data$person_id%in%infection_status_increase$person_id&data$infection_status==1
        data[index_1,c("infection_status","last_infection","time_since_last_infection","time_since_first_infection")] = data.table(infection_status=2,last_infection=1,time_since_last_infection=-delta_t,time_since_first_infection = -delta_t)
        if(sum(index_1)>0)
          new_infections = rbind(new_infections,data.frame(data[index_1,"age"],time=t))
        index_3 = data$person_id%in%infection_status_increase$person_id&data$infection_status==3
        data[index_3,c("infection_status","last_infection","time_since_last_infection")] = data.table(infection_status=3,last_infection=3,time_since_last_infection=-delta_t)
        index_4 = data$person_id%in%infection_status_increase$person_id&data$infection_status==4
        data[index_4,c("infection_status","last_infection","time_since_last_infection")] = data.table(infection_status=3,last_infection=1,time_since_last_infection=-delta_t)
        if(sum(index_4)>0)
          new_infections = rbind(new_infections,data.frame(data[index_4,"age"],time=t))
        index_5 = data$person_id%in%infection_status_increase$person_id&data$infection_status==5
        data[index_5,c("infection_status","last_infection","time_since_last_infection")] = data.table(infection_status=3,last_infection=3,time_since_last_infection=-delta_t)
        person_id_5 = data$person_id[index_5]
        vaccinated_infected_people = vaccinated_infected_people[!(person_id%in%person_id_5)]
        
        birth_change = event_data[event=="birth"]
        if(nrow(birth_change)>0){
          #make sure all people who were chosen to give birth are in age class 3 (mothers)
          data[data$person_id%in%birth_change$person,c("age_class")]=3
          #look at the data of these people
          current_people_data = data[data$person_id%in%birth_change$person,]
          #some women might be new mothers - update family number of these people
          family_update = current_people_data[family_id==1]$person_id
          if(length(family_update)>0){
            data[data$person_id%in%family_update,c("family_id")] = (family_id+1):(family_id+length(family_update)) #women have now been assigned a new family number
            family_id = family_id+length(family_update) #update family_id count
            current_people_data = data[data$person_id%in%birth_change$person,] #update our info on the women having children to include their new family ids
          }
          data_baby = NULL
          
          #now choose infection status of new babies
          #first have uninfected mothers or vaccinated mothers - these will all have no risk of transmission
          #if mothers are of infection status 5, make congenital based on infection status before vaccination
          uninfected_mothers= current_people_data[infection_status%in%c(1,4)]
          if (nrow(uninfected_mothers)>0){
            data_baby = rbind(data_baby,
                              data.table(family_id =uninfected_mothers$family_id,
                                         person_id = (person_id+1):(person_id+nrow(uninfected_mothers)),
                                         age_class= 1,
                                         age = -delta_t,
                                         infection_status = 1,
                                         last_infection = 1,
                                         time_since_last_infection = -delta_t,
                                         time_since_first_infection = NA))
            person_id = person_id+nrow(uninfected_mothers)
            #update info on births to be used in summary statistics
          }
          
          #want to figure out the maximum probability of congenital infection mothers have had in the last 9 months
          infected_mothers =current_people_data[infection_status%in%c(2,3,5)]
          if(nrow(infected_mothers)>0){
            infected_mothers = merge(infected_mothers[,c("person_id","family_id")],mothers[,c("person_id","infection_status","prob_congenital","time_since_last_infection")],by = "person_id",all.x = TRUE,all.y = FALSE)
            current_people_data = infected_mothers
            if(nrow(current_people_data)>0){
              for (i in 1:nrow(current_people_data)){
                rate = current_people_data$prob_congenital[i]
                baby_inf = as.numeric(rbinom(1,1,rate))
                data_baby = rbind(data_baby,
                                  data.table(family_id =current_people_data$family_id[i],
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
            }
          }
          
          data_baby = data_baby[order(-person_id,decreasing=TRUE),]
          data = rbind(data,data.table(data_baby[,1:5],
                                       breast_status = rbinom(nrow(data_baby),1,fit_vals_breast[2]),
                                       diaper_status=1,
                                       gender = sample(c(0,1),nrow(data_baby),TRUE),
                                       data_baby[,6:8],
                                       vaccinated = 0))
          congen_birth = data_baby[infection_status==2]
          if(nrow(congen_birth)>0)
            source_infection = rbind(source_infection,data.table(time_of_infection = t,person_id = congen_birth$person_id,event = "mother"))
        }
      }
      
      data[,c("age","time_since_last_infection","time_since_first_infection")] = data[,c("age","time_since_last_infection","time_since_first_infection")]+delta_t
      if(!is.null(vaccinated_infected_people))
        vaccinated_infected_people[,c("age","time_since_last_infection","time_since_first_infection")] = vaccinated_infected_people[,c("age","time_since_last_infection","time_since_first_infection")]+delta_t
      
      #age children into adults
      #if this makes mothers childless, put mother back in adult class
      children = data[age_class==1&round(age,2)>=5]
      if (nrow(children)>0){
        num_children_family = as.data.table(data%>%dplyr::filter(family_id%in%children$family_id&age_class==1)%>%
                                              dplyr::group_by(family_id)%>%dplyr::tally()) #count current number of children - we are about to lose 1
        data[data$person_id%in%children$person_id,c("family_id")]=1 #put child into family_id=1
        data[data$person_id%in%children$person_id,c("age_class")]=2 #age child into age class 2 (adult)
        data[data$person_id%in%children$person_id,c("breast_status")]=0 #make sure child is no longer breastfeeding
        data[data$person_id%in%children$person_id,c("diaper_status")]=0 #make sure child is no longer is diapers
        no_children = num_children_family[n==1] #families that will no longer have children
        if(nrow(no_children)>0){
          data[data$family_id%in%no_children$family_id,c("age_class")]=2 #put mother back into age class 2
          data[data$family_id%in%no_children$family_id,c("family_id")]=1 #put mother back into family id 1
        }
      }
      
      #Age adults into death
      adult_age = data[age>50]$person_id
      if(length(adult_age)>0){
        data[data$person_id%in%adult_age,c("age_class")]=4
      }
      
      data = data[age_class<=3]
      
      
      #put individuals in reinf/immune waning class based on their time_to_wane
      waning_people = data[infection_status==2]
      if (nrow(waning_people)>0){
        waning_people = data.table(person_id = waning_people$person_id,
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
      
      #repeat for those who were vaccinated when infected
      if(!is.null(vaccinated_infected_people)){
        #children into adults
        children = vaccinated_infected_people[age_class==1&round(age,2)>=5]
        if (nrow(children)>0){
          num_children_family = as.data.table(vaccinated_infected_people%>%dplyr::filter(family_id%in%children$family_id&age_class==1)%>%
                                                dplyr::group_by(family_id)%>%dplyr::tally()) #count current number of children - we are about to lose 1
          vaccinated_infected_people[vaccinated_infected_people$person_id%in%children$person_id,c("family_id")]=1 #put child into family_id=1
          vaccinated_infected_people[vaccinated_infected_people$person_id%in%children$person_id,c("age_class")]=2 #age child into age class 2 (adult)
          vaccinated_infected_people[vaccinated_infected_people$person_id%in%children$person_id,c("breast_status")]=0 #make sure child is no longer breastfeeding
          vaccinated_infected_people[vaccinated_infected_people$person_id%in%children$person_id,c("diaper_status")]=0 #make sure child is no longer is diapers
        }
        
        
        #age to adults
        adult_age = vaccinated_infected_people[age>50]$person_id
        if(length(adult_age)>0){
          vaccinated_infected_people[vaccinated_infected_people$person_id%in%adult_age,c("age_class")]=4
        }
        vaccinated_infected_people = vaccinated_infected_people[age_class<=3]
        #move into infection class 3
        waning_people = vaccinated_infected_people[infection_status==2]
        if (nrow(waning_people)>0){
          waning_people = data.table(person_id = waning_people$person_id,
                                     infection_status=2,
                                     time_since_last_infection = waning_people$time_since_last_infection, 
                                     age_of_infection = waning_people$age-waning_people$time_since_last_infection)
          if(nrow(waning_people)>0){
            waning_people$age_of_infection =round(waning_people$age_of_infection,3)
            waning_info = max_infectivity[,c(1,4)]
            waning_info$age_of_infection = round(waning_info$age_of_infection,3)
            waning_people = merge(waning_people,waning_info,by = "age_of_infection",all.x = TRUE,all.y = FALSE)
            waning_people = waning_people[waning_people$time_since_last_infection>=(waning_people$time_to_wane),]
            vaccinated_infected_people[vaccinated_infected_people$person_id%in%waning_people$person_id,"infection_status"]=3
            vaccinated_infected_people[vaccinated_infected_people$person_id%in%waning_people$person_id,"time_since_last_infection"]=delta_t
            vaccinated_infected_people[vaccinated_infected_people$person_id%in%waning_people$person_id,"last_infection"]=2
          }
        }
      }
      if(t>=5){
        #vaccinate those who are correct age at the correct fraction
        age_rounded = round(data$age,2)
        vacc_potential = which(age_rounded%in%age_vacc[1]&data$infection_status%in%c(1)) #let infection_status=4 represent those who are vaccinated
        if (length(vacc_potential)>0){
          chosen = rbinom(length(vacc_potential),1,vacc_prob)
          vacc_person =  vacc_potential[which(chosen==1)]
          data[vacc_person,"last_infection"]=data[vacc_person,"infection_status"]
          data[vacc_person,"infection_status"]=4
          data[vacc_person,"time_since_last_infection"]=0
          data[vacc_person,"vaccinated"]=1
        }
        vacc_potential = which(age_rounded%in%age_vacc[1]&data$infection_status%in%c(2,3)) #let infection_status=4 represent those who are vaccinated
        #if already infected when vaccinated, you get added to a special list - previous infection will determine infectivity - vaccination will determine immune protection
        if (length(vacc_potential)>0){
          chosen = rbinom(length(vacc_potential),1,vacc_prob)
          vacc_person =  vacc_potential[which(chosen==1)]
          addition_vacc = data[vacc_person,]
          vaccinated_infected_people = rbind(vaccinated_infected_people,addition_vacc)
          #vaccinated_infected_people = vaccinated_infected_people[order(-person_id,decreasing=TRUE),]
          data[vacc_person,"last_infection"]=data[vacc_person,"infection_status"]
          data[vacc_person,"infection_status"]=5
          data[vacc_person,"time_since_last_infection"]=0
          data[vacc_person,"vaccinated"]=1
        } 
        
        #vaccinated older 
        if(length(age_vacc)>1){
          vacc_potential = which(age_rounded%in%age_vacc[2:length(age_vacc)]&data$infection_status%in%c(1,4)&data$vaccinated>=1) #let infection_status=4 represent those who are vaccinated
          if (length(vacc_potential)>0){
            chosen = rbinom(length(vacc_potential),1,1)
            vacc_person =  vacc_potential[which(chosen==1)]
            data[vacc_person,"last_infection"]=data[vacc_person,"infection_status"]
            data[vacc_person,"infection_status"]=4
            data[vacc_person,"time_since_last_infection"]=0
            data[vacc_person,"vaccinated"]=data[vacc_person,"vaccinated"]+1
          }
          vacc_potential = which(age_rounded%in%age_vacc[2:length(age_vacc)]&data$infection_status%in%c(5)&data$vaccinated>=1)
          #if already vaccinated when vaccinated again, you've already been added to a special list - previous infection will determine infectivity - vaccination will determine immune protection
          if (length(vacc_potential)>0){
            chosen = rbinom(length(vacc_potential),1,1)
            vacc_person =  vacc_potential[which(chosen==1)]
            data[vacc_person,"last_infection"]=data[vacc_person,"infection_status"]
            data[vacc_person,"infection_status"]=5
            data[vacc_person,"time_since_last_infection"]=0
            data[vacc_person,"vaccinated"]=data[vacc_person,"vaccinated"]+1
          }
          vacc_potential = which(age_rounded%in%age_vacc[2:length(age_vacc)]&data$infection_status%in%c(2,3)&data$vaccinated>=1) #let infection_status=4 represent those who are vaccinated
          #if already infected when vaccinated, you get added to a special list - previous infection will determine infectivity - vaccination will determine immune protection
          if (length(vacc_potential)>0){
            chosen = rbinom(length(vacc_potential),1,1)
            vacc_person =  vacc_potential[which(chosen==1)]
            addition_vacc = data[vacc_person,]
            vaccinated_infected_people = rbind(vaccinated_infected_people,addition_vacc)
            #vaccinated_infected_people = vaccinated_infected_people[order(-person_id,decreasing=TRUE),]
            data[vacc_person,"last_infection"]=data[vacc_person,"infection_status"]
            data[vacc_person,"infection_status"]=5
            data[vacc_person,"time_since_last_infection"]=0
            data[vacc_person,"vaccinated"]=data[vacc_person,"vaccinated"]+1
          } 
        } 
        
      }
      #infectivity of each individual to an uninfected individual
      infectivity_info = infectivity_fun(vaccinated_infected_people,data)
      
      chosen_children = which(data$age<5)
      chosen_children_dat = data.frame(data[chosen_children,],infectivity_info[chosen_children,c("infectivity","immune_protection")])
      child_data = rbind(child_data,data.frame(chosen_children_dat,time = t))
      
      t = t+delta_t
      print(t)
      
      
    }
  })
  
  all_mothers_all = all_mothers
  #all_mothers_new = merge(all_mothers,source_infection,by = "person_id",all.x = TRUE,all.y = FALSE,allow.cartesian = TRUE)
  all_mothers_new1 = subset(all_mothers,infection_status==1)
  all_mothers_new1$event="uninfected"
  all_mothers_new2 = subset(all_mothers,infection_status==4)
  all_mothers_new2$event = "uninfected vaccinated"
  all_mothers_new3 = subset(all_mothers,infection_status%in%c(2,3))
  all_mothers_new3$event="chronic"
  all_mothers_new4 = subset(all_mothers,infection_status%in%c(20,30))
  all_mothers_new4 = merge(all_mothers_new4,source_infection,by = "person_id",all.x = TRUE,all.y = FALSE,allow.cartesian = TRUE)
  all_mothers_new4  = subset(all_mothers_new4,time_of_birth>time_of_infection)
  all_mothers_new4$event = revalue(all_mothers_new4$event,c("adult_inf_rate"="adult",
                                                            "fam_inf_rate"="child",
                                                            "child_inf_rate"="child"))
  all_mothers_new4  = as.data.table(all_mothers_new4%>%
                                      dplyr::group_by(person_id,time_of_birth,rate_birth) %>% 
                                      dplyr::slice(which.min(time_of_birth-time_of_infection)))
  all_mothers_new4 = all_mothers_new4[,c(1:7,9)]
  all_mothers_new5 = subset(all_mothers,infection_status==5)
  all_mothers_new5$event = "chronic vaccinated"
  all_mothers_new5 = rbind(all_mothers_new1,all_mothers_new2,all_mothers_new3,all_mothers_new4,all_mothers_new5)
  all_mothers_new5$time_of_birth = floor(all_mothers_new5$time_of_birth/5)*5
  all_mothers_new5 = as.data.table(all_mothers_new5%>%dplyr::group_by(infection_status,event,time_of_birth)%>%dplyr::summarize(n = sum(rate_birth*delta_t),n_congen = sum(prob_congenital*rate_birth*delta_t))%>%dplyr::group_by(time_of_birth)%>%dplyr::mutate(freq = n/sum(n),freq_congen=n_congen/sum(n_congen)))
  all_mothers_new5 = data.frame(rand = rand, sim_num = sim_num,vacc_waning = vacc_waning,vacc_const = vacc_const,vacc_prob = vacc_prob, age_vacc = age_vacc[1],booster_times=booster_times,all_mothers_new5)
  
  Reff2 = as.data.frame(Reff%>%dplyr::group_by(time = (floor(time/5)*5))%>%dplyr::summarize(new_per_year = mean(new_per_year),num_infected = mean(num_infected),age_acq=mean(age_acq)))
  Reff2 = data.frame(rand = rand, sim_num = sim_num, vacc_waning = vacc_waning, vacc_const = vacc_const, vacc_prob = vacc_prob, age_vacc = age_vacc[1],booster_times = booster_times,Reff2)
  #Reff = gather(Reff,source_inf,rate,adult_inf_rate:child_inf_rate, factor_key=TRUE)
  
  child_data$time_round = floor(child_data$time/5)*5
  
  #want ave infectivity of infected children
  ave_infectivity = as.data.frame(child_data%>%dplyr::filter(!infection_status%in%c(1,4))%>%
                                    dplyr::group_by(time_round)%>%
                                    dplyr::summarize(infectivity = mean(infectivity),immune_protection = mean(immune_protection)))
  
  #ggplot(ave_infectivity,aes(x = time_round,y = infectivity))+geom_line()
  
  
  #delaying infection?
  new_children = subset(child_data,age==0)
  
  new_children = subset(child_data,person_id%in%new_children$person_id)
  
  new_children = subset(new_children,!infection_status%in%c(1,4,5))
  
  new_children = arrange(new_children,person_id,age)
  
  new_children = new_children[!duplicated(new_children[,c("person_id")]),]
  
  age_inf = as.data.frame(new_children%>%dplyr::group_by(time_round)%>%dplyr::summarize(age_inf = mean(age)))
  
  #ggplot(age_inf,aes(x = time_round,y = age_inf))+geom_line()
  
  
  #preventing infection?
  total_children = as.data.frame(child_data%>%group_by(time_round)%>%tally())
  colnames(total_children)[2]="total_n"
  
  children_infected = child_data
  children_infected$infection_status = as.character(children_infected$infection_status)
  children_infected$infection_status = revalue(children_infected$infection_status,c("1"="uninfected","2"="infected","3"="infected","4"="uninfected","5"="infected"))
  children_infected = as.data.frame(children_infected%>%group_by(time_round,infection_status,vaccinated)%>%tally())
  children_infected = merge(children_infected,total_children)
  children_infected$freq = children_infected$n/children_infected$total_n
  children_infected$status =paste0(children_infected$vaccinated,"_",children_infected$infection_status)
  
  #ggplot(children_infected,aes(x = time_round,y = freq,fill = paste0(infection_status,"_",vaccinated)))+geom_bar(stat="identity")         
  
  one = merge(age_inf,ave_infectivity)
  
  two = merge(children_infected,one)
  two = data.frame(rand = rand, sim_num = sim_num, vacc_waning = vacc_waning, vacc_const = vacc_const, vacc_prob = vacc_prob, age_vacc = age_vacc[1],booster_times = booster_times,two)
  
  fwrite(all_mothers_new5,paste0("/home/sgazea01/projects/def-craigm/sgazea01/LowSero1/summarized_mother",rand,".csv"),row.names = FALSE) 
  fwrite(Reff2,paste0("/home/sgazea01/projects/def-craigm/sgazea01/LowSero1/Reff",rand,".csv"),row.names = FALSE) 
  fwrite(two,paste0("/home/sgazea01/projects/def-craigm/sgazea01/LowSero1/children",rand,".csv"),row.names = FALSE) 
  
  output = paste0(sim_num," done")
  return(output)
  
}

library(parallel)

setwd("/home/sgazea01/projects/def-craigm/sgazea01/LowSero1")
#setwd("C:/Users/c-byr.LAPTOP-UG0F9R6R/OneDrive - math.ubc.ca/cmv/cmv transmission model/POLYMOD/original contact matrix, no daycare/CMV_transmission_model40_2")
NHANES = subset(fread("NHANES_CMV_antibodies_combined_weighted.csv",header = TRUE),age<=50)
#init_data = fread("good_fit_initial_data_adjusted.csv",header = TRUE)
nodeslist <- unlist(strsplit(Sys.getenv("NODESLIST"), split=" "))

#clusterExport(cl,"init_data")


pars =read.table("ABC_Lenormand_params.txt",header = FALSE)
stats = read.table("ABC_Lenormand_stats.txt",header = FALSE)
pars = cbind(pars[,1:8],stats[,2])


vacc_prob = c(1,2/3,1/3)
age_vacc = c(2/12,2,12,25)
booster_times = c(0)



input = NULL
for (k in 1:length(vacc_prob)){
  for (l in 1:length(age_vacc)){
    for (m in 1:length(booster_times)){
      input = rbind(input,cbind(pars,0,pars[,8],vacc_prob[k],age_vacc[l],booster_times[m]))
    }
  }
}



colnames(input) = c(paste0("p",1:8),"sim_num","vacc_waning","vacc_const","vacc_prob","age_vacc","booster_times")
input$rand = 1:nrow(input)
count = nrow(input)


for (i in 1:(count/600)){
  cl <- makeCluster(nodeslist, type = "PSOCK",outfile = "/home/sgazea01/projects/def-craigm/sgazea01/LowSero1/errors.txt")

  clusterEvalQ(cl, {
    library(plyr)
    library(dplyr)
    library(ggplot2)
    library(pomp)
    library(reshape2)
    library(tidyr)
    #library(DescTools)
    library(matrixStats)
    library(extdplyr)
    library(data.table)
  })

  clusterExport(cl,"cmv_sim")
  input_chosen = split(input[((i-1)*(600)+1):(i*600),], seq(600))
  output = clusterApplyLB(cl,input_chosen, cmv_sim)
  output = do.call(rbind.data.frame,output)
  fwrite(output,paste0("/home/sgazea01/projects/def-craigm/sgazea01/LowSero1/vaccine_output",i,".csv"),row.names = FALSE)
  stopCluster(cl)
}

dir.create("/home/sgazea01/projects/def-craigm/sgazea01/LowSero1/combined_child_vaccine_2")
setwd("/home/sgazea01/projects/def-craigm/sgazea01/LowSero1")

nrow_input = nrow(input)
for (i in 1:(count/600)){
  all_mother = NULL
  all_reff = NULL
  all_children = NULL
  for (j in ((i-1)*600+1):(i*600)){
    new =read.csv(paste0("summarized_mother",j,".csv"))
    all_mother = rbind(all_mother,new)
    new2 =read.csv(paste0("Reff",j,".csv"))
    all_reff = rbind(all_reff,new2)
    new3 = read.csv(paste0("children",j,".csv"))
    all_children = rbind(all_children,new3)
  }
  write.csv(all_mother,paste0("/home/sgazea01/projects/def-craigm/sgazea01/LowSero1/combined_child_vaccine_2/summarized_mother_combined",i,".csv"),row.names=FALSE)
  write.csv(all_reff,paste0("/home/sgazea01/projects/def-craigm/sgazea01/LowSero1/combined_child_vaccine_2/Reff_combined",i,".csv"),row.names=FALSE)
  write.csv(all_children,paste0("/home/sgazea01/projects/def-craigm/sgazea01/LowSero1/combined_child_vaccine_2/children_combined",i,".csv"),row.names=FALSE)
}

file.copy("sterilizing_vaccine.r", paste0("/home/sgazea01/projects/def-craigm/sgazea01/LowSero1/combined_child_vaccine_2"), overwrite = TRUE )
file.copy("vaccine_cluster.sh", paste0("/home/sgazea01/projects/def-craigm/sgazea01/LowSero1/combined_child_vaccine_2"), overwrite = TRUE )



