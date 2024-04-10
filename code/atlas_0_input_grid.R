
# Prepare metabolomics dataset for the atlas
params=list(
  rank_normalization_yes=T, # determine if rank normalize the data
  rescale_yes=F, # determine after rank normalize if rescale the data with the original std
  imp_file_code="metacov",
  output_path="output/atlas/imputed",
  b1_imp_path_0="data/Imputed_metab_datasets/b1_zero.csv", # the original file (without imputed values)
  b1_imp_path_1="data/Imputed_metab_datasets/b1_metacov_1.csv",
  b1_imp_path_2="data/Imputed_metab_datasets/b1_metacov_2.csv",
  b1_imp_path_3="data/Imputed_metab_datasets/b1_metacov_3.csv",
  b1_imp_path_4="data/Imputed_metab_datasets/b1_metacov_4.csv",
  b1_imp_path_5="data/Imputed_metab_datasets/b1_metacov_5.csv",
  b2_imp_path_0="data/Imputed_metab_datasets/b2_zero.csv", # the original file (without imputed values)
  b2_imp_path_1="data/Imputed_metab_datasets/b2_metacov_1.csv",
  b2_imp_path_2="data/Imputed_metab_datasets/b2_metacov_2.csv",
  b2_imp_path_3="data/Imputed_metab_datasets/b2_metacov_3.csv",
  b2_imp_path_4="data/Imputed_metab_datasets/b2_metacov_4.csv",
  b2_imp_path_5="data/Imputed_metab_datasets/b2_metacov_5.csv",
  mdl1_covars="age gender center background",
  mdl2_covars="age gender center background bmi",
  mdl3_covars="age gender center background bmi alcohol_use cigarette_use gpaq_total_metG ahei2010"
)
# Source functions
source("sol_metabolomics_functions.R")
# System parameters
options(scipen = 999, stringsAsFactors = F)
options(survey.lonely.psu="remove")

library_list<-c("renv","tidyverse","hms","labelled","purrr","tableone","stringr","plyr","dplyr","ggrepel","ggpubr","htmlwidgets","corrplot","ggplot2","FactoMineR","CircStats","missMDA", "pastecs","devtools","factoextra","gridExtra","survey","labelled","purrr","DT")
lapply(library_list, require, character.only = TRUE)


###########################
# generate input grid table
time_var<-c("wddur","wedur","slpdur","slp_short_wd","slp_long_wd","we_minus_wd")
disturb_var<-c("whiirs","restless","pill","backtosleep","earlywake","freqwake","fallasleep","ess","essgt10","snore")
ah_var<-c("slpa36","slpa54","slpa36gt5","slpa36gt15","slpa54gt5","slpa54gt15","slpa39","slpa42","slpa91","slpa92","slpa97","per90","event_length_sec","hypoxicburden_harmonized")
hr_var<-c("slpa111","slpa112","slpa113","slpa114")
sleep_var<-c(time_var,disturb_var,ah_var,hr_var)
sleep_trait_input_list<-c("slpdur","slpa54","per90","slpa36","slpa114")
circular_sleep_trait_input_list<-c("wdwake_rad","wewake_rad","wdbed_rad","webed_rad","wdmid_rad","wemid_rad")
model_list<-c("model_1","model_2")
model_covariates<-c("age gender center background bmi batch","age gender center background bmi batch alcohol_use cigarette_use gpaq_total_met ahei2010")
rank_normalization_input_list<-c(TRUE)
rescale_yes_input_list<-c(FALSE)
imp_data_input_list<-c("metacov")
statifier_list<-c("gender")

# create an expanded grid of sleep_var, predior_seq,covariate_seq,strata_seq, and dataset_seq
# linear sleep traits
atlas_input_grid<-expand.grid(imp_dataset=imp_data_input_list,
                              sleep_trait=sleep_var,
                              model=model_list,                              
                              rank_norm=rank_normalization_input_list,
                              rescale=rescale_yes_input_list,
                              stratifier=statifier_list)
atlas_input_grid$covariates<-ifelse(atlas_input_grid$model=="model_1",model_covariates[1],model_covariates[2])
write.csv(atlas_input_grid,file="data/atlast_input_grid.csv",row.names = T)
# circular sleep traits
circular_atlas_input_grid<-expand.grid(imp_dataset=imp_data_input_list,
                              sleep_trait=circular_sleep_trait_input_list,
                              model=model_list,
                              rank_norm=rank_normalization_input_list,
                              rescale=rescale_yes_input_list,
                              stratifier=statifier_list)
circular_atlas_input_grid$covariates<-ifelse(circular_atlas_input_grid$model=="model_1",model_covariates[1],model_covariates[2])
write.csv(circular_atlas_input_grid,file="data/atlast_cir_input_grid.csv",row.names = T)
