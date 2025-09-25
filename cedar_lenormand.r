setwd("/home/sgazea01/projects/def-craigm/sgazea01/LowSero")
  
source("contact_matrix_use_continuous_progress.r")
source("EasyABC_internal_modified.r")
source("ABC_sequential.r")
library("lgarch")
library("mnormt")
library("lhs")


#library(EasyABC)

prior_new <- list(c("unif",-3,-1.5), #child primary  
                  c("unif",-1.5,-0.5), #adult primary factor
                  c("unif",-1,0), #secondary factor
                  c("unif",-2.5,-0.5), #chronic factor
                  c("unif",0,3.5),   #breast factor
                  c("unif",0,1),  #diaper factor
                  c("unif",0.5,1.5),  #immune waning
                  c("unif",0.1,0.6)  #immune const
)

NHANES = subset(read.csv("NHANES_CMV_antibodies_combined_weighted.csv",header = TRUE),age<=50)
init_data = read.csv("good_initial_fit_data_10000_new.csv",header = TRUE)



ABC_Lenormand<-ABC_sequential(method="Lenormand", model=cmv_sim,alpha =1000/1200,
                              inside_prior = TRUE,use_seed = TRUE,n_cluster = 100,
                              prior=prior_new, nb_simul=120, summary_stat_target=c(0,0),
                              p_acc_min=0.1,verbose=TRUE,progress_bar = TRUE)

write.table(ABC_Lenormand$stats, file = paste0("ABC_Lenormand_stats.txt"),
            sep = " ", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
write.table(ABC_Lenormand$nsim, file = paste0("ABC_Lenormand_nsim.txt"),
            sep = " ", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
write.table(ABC_Lenormand$epsilon, file = paste0("ABC_Lenormand_epsilon.txt"),
            sep = " ", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
write.table(ABC_Lenormand$weights, file = paste0("ABC_Lenormand_weights.txt"),
            sep = " ", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
write.table(ABC_Lenormand$stats_normalization, file = paste0("ABC_Lenormand_normalization.txt"),
            sep = " ", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
write.table(ABC_Lenormand$computime, file = paste0("ABC_Lenormand_computime.txt"),
            sep = " ", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
write.table(ABC_Lenormand$param, file = paste0("ABC_Lenormand_params.txt"),
            sep = " ", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
