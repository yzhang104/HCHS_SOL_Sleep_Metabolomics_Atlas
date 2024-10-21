args <- commandArgs(trailingOnly=TRUE)
comb_number <- as.numeric(args[1])
# Source functions
source("sol_metabolomics_functions.R")
source("acat.R")
library("survey")
library("tidyverse")
library("dplyr")
options(survey.lonely.psu="remove")   
options(scipen = 999, stringsAsFactors = F)

cir_input_grid<-read.csv(file="data/atlast_cir_input_grid.csv",stringsAsFactors = F,header = T)
# read in the nonlinearity results table
nl_tbl<-read.csv("data/nonlinearity_model.csv")

dataset_code<-cir_input_grid[comb_number,"imp_dataset"]
sleep_trait<-cir_input_grid[comb_number,"sleep_trait"]
rank_normalization_yes<-cir_input_grid[comb_number,"rank_norm"]
rescale_yes<-cir_input_grid[comb_number,"rescale"]
model<-cir_input_grid[comb_number,"model"]
covars<-if (!is.null(cir_input_grid[comb_number,"covariates"])) {unlist(strsplit(cir_input_grid[comb_number,"covariates"],split = " "))} else {na}
stratifier<-cir_input_grid[comb_number,"stratifier"]

# read in annotation file
hchs_annot<-readxl::read_xlsx("data/hchs_batch2/batch2_data.xlsx",sheet="metabolites.info")%>%
  dplyr::mutate(chemid=paste0("chem_",CHEM_ID),
                CHEM_ID=as.character(CHEM_ID))
colnames(hchs_annot)<-tolower(colnames(hchs_annot))
xeno_metab_list<-unname(unlist(hchs_annot[which(hchs_annot$super_pathway=="Xenobiotics"),"chemid"]))
nonxeno_metab_list<-unname(unlist(hchs_annot[which(hchs_annot$super_pathway!="Xenobiotics"),"chemid"]))

# read in multiple imputed datasets
b1_imp_path_0="data/Imputed_metab_datasets/b1_zero.csv" # the original file (without imputed values)
b1_imp_path_1=paste0("data/Imputed_metab_datasets/b1_",dataset_code,"_1.csv")
b1_imp_path_2=paste0("data/Imputed_metab_datasets/b1_",dataset_code,"_2.csv")
b1_imp_path_3=paste0("data/Imputed_metab_datasets/b1_",dataset_code,"_3.csv")
b1_imp_path_4=paste0("data/Imputed_metab_datasets/b1_",dataset_code,"_4.csv")
b1_imp_path_5=paste0("data/Imputed_metab_datasets/b1_",dataset_code,"_5.csv")
b2_imp_path_0="data/Imputed_metab_datasets/b2_zero.csv" # the original file (without imputed values)
b2_imp_path_1=paste0("data/Imputed_metab_datasets/b2_",dataset_code,"_1.csv")
b2_imp_path_2=paste0("data/Imputed_metab_datasets/b2_",dataset_code,"_2.csv")
b2_imp_path_3=paste0("data/Imputed_metab_datasets/b2_",dataset_code,"_3.csv")
b2_imp_path_4=paste0("data/Imputed_metab_datasets/b2_",dataset_code,"_4.csv")
b2_imp_path_5=paste0("data/Imputed_metab_datasets/b2_",dataset_code,"_5.csv")

b1_metacov_0<-read.csv(paste0(b1_imp_path_0),header = TRUE, sep = ",")
b1_metacov_1<-read.csv(paste0(b1_imp_path_1),header = TRUE, sep = ",")
b1_metacov_2<-read.csv(paste0(b1_imp_path_2),header = TRUE, sep = ",")
b1_metacov_3<-read.csv(paste0(b1_imp_path_3),header = TRUE, sep = ",")
b1_metacov_4<-read.csv(paste0(b1_imp_path_4),header = TRUE, sep = ",")
b1_metacov_5<-read.csv(paste0(b1_imp_path_5),header = TRUE, sep = ",")
b2_metacov_0<-read.csv(paste0(b2_imp_path_0),header = TRUE, sep = ",")
b2_metacov_1<-read.csv(paste0(b2_imp_path_1),header = TRUE, sep = ",")
b2_metacov_2<-read.csv(paste0(b2_imp_path_2),header = TRUE, sep = ",")
b2_metacov_3<-read.csv(paste0(b2_imp_path_3),header = TRUE, sep = ",")
b2_metacov_4<-read.csv(paste0(b2_imp_path_4),header = TRUE, sep = ",")
b2_metacov_5<-read.csv(paste0(b2_imp_path_5),header = TRUE, sep = ",")

# min value imputed datasets
b1_metacov_0[b1_metacov_0 == 0] <- NA
b2_metacov_0[b2_metacov_0 == 0] <- NA

# assess missingness in the raw datasets and exclude obs with missing more than 25% metabolites
b1_ID_list <- b1_metacov_0[rowSums(is.na(b1_metacov_0)) < (ncol(b1_metacov_0)-2)*0.25, "SOL_ID"]
b2_ID_list <- b2_metacov_0[rowSums(is.na(b2_metacov_0)) < (ncol(b2_metacov_0)-2)*0.25, "SOL_ID"]

# drop non chemical name (except ID) columns
edit_list<-list(b1_metacov_1=b1_metacov_1,
                b1_metacov_2=b1_metacov_2,
                b1_metacov_3=b1_metacov_3,
                b1_metacov_4=b1_metacov_4,
                b1_metacov_5=b1_metacov_5,
                b2_metacov_1=b2_metacov_1,
                b2_metacov_2=b2_metacov_2,
                b2_metacov_3=b2_metacov_3,
                b2_metacov_4=b2_metacov_4,
                b2_metacov_5=b2_metacov_5)
edit_list <- lapply(edit_list, function(x) x %>% 
                      dplyr::select(-one_of(intersect(xeno_metab_list,names(.))))%>% 
                      dplyr::select(SOL_ID,tidyselect::any_of(intersect(hchs_annot$chemid,names(.))))%>%
                      dplyr::filter(SOL_ID%in%c(b1_ID_list,b2_ID_list)))

# write the modified data frames back to the global environment with the same names as the original data frames
list2env(edit_list, envir = globalenv())

b1_imp_list<-list(b1_1=b1_metacov_1,b1_2=b1_metacov_2,b1_3=b1_metacov_3,b1_4=b1_metacov_5,b1_5=b1_metacov_5)
b2_imp_list<-list(b2_1=b2_metacov_1,b2_2=b2_metacov_2,b2_3=b2_metacov_3,b2_4=b2_metacov_5,b2_5=b2_metacov_5)

# impute the xenobiotic metabolites with half value
b1_half<-impute_func(data=b1_metacov_0[,colnames(b1_metacov_0)%in%c("SOL_ID",xeno_metab_list)],method="half")%>% 
  dplyr::filter(SOL_ID%in%c(b1_ID_list,b2_ID_list))
b2_half<-impute_func(data=b2_metacov_0[,colnames(b2_metacov_0)%in%c("SOL_ID",xeno_metab_list)],method="half")%>% 
  dplyr::filter(SOL_ID%in%c(b1_ID_list,b2_ID_list))

# preprocessing the imputed metabolomics data
# 1. Compute the SDm of the imputed  metabolites (only for metabolites exclude all covariates and ID)
b1_imp_std_list<-list()
for (m in seq_along(b1_imp_list)){
  b1_imp_std_list[[m]]<-cal_matrix_sd(dataset=b1_imp_list[[m]],col_excluded=c(colnames(b1_imp_list[[m]])[!colnames(b1_imp_list[[m]])%in%colnames(b1_metacov_0)],"SOL_ID"))
  names(b1_imp_std_list)[[m]]<-paste0(names(b1_imp_list)[[m]],"_std")
}
b2_imp_std_list<-list()
for (m in seq_along(b2_imp_list)){
  b2_imp_std_list[[m]]<-cal_matrix_sd(dataset=b2_imp_list[[m]],col_excluded=c(colnames(b2_imp_list[[m]])[!colnames(b2_imp_list[[m]])%in%colnames(b2_metacov_0)],"SOL_ID"))
  names(b1_imp_std_list)[[m]]<-paste0(names(b1_imp_list)[[m]],"_std")
}

b1_half_std<-cal_matrix_sd(dataset=b1_half,col_excluded="SOL_ID")
names(b1_half_std)<-paste0(names(b1_half),"_std")

b2_half_std<-cal_matrix_sd(dataset=b2_half,col_excluded="SOL_ID")
names(b2_half_std)<-paste0(names(b2_half),"_std")

b1b2_half_list<-list(b1_half=b1_half,b2_half=b2_half)
b1b2_half_std_list<-list(b1_half_std=b1_half_std,b2_half_std=b2_half_std)

# 2. Rank normalization
# 2.1 Rank-normalize continuous metabolites
if (rank_normalization_yes==T){
  # apply rank normalization to a data frame from a list with selected columns
  # using the following function can't update the dataframes
  if (rescale_yes==TRUE){
    b1_imp_norm_list<-purrr::map2(b1_imp_list,b1_imp_std_list, ~{
      # select all columns except for ID
      cols_to_normalize <- names(.x)[names(.x) %in% colnames(b1_metacov_0[,-1])]
      # print(cols_to_normalize)
      # apply rank normalization to each column
      .x[cols_to_normalize] <- lapply(.x[cols_to_normalize], rank_normalisation)
      .x[cols_to_normalize] <-as.matrix(.x[cols_to_normalize])%*%.y
      # return the updated data frame
      return(.x)
    })
    b2_imp_norm_list<-purrr::map2(b2_imp_list,b2_imp_std_list, ~{
      # select all columns except for ID
      cols_to_normalize <- names(.x)[names(.x) %in% colnames(b2_metacov_0[,-1])]
      # print(cols_to_normalize)
      # apply rank normalization to each column
      .x[cols_to_normalize] <- lapply(.x[cols_to_normalize], rank_normalisation)
      .x[cols_to_normalize] <-as.matrix(.x[cols_to_normalize])%*%.y
      # return the updated data frame
      return(.x)
    })
    b1b2_half_norm_list<-purrr::map2(b1b2_half_list,b1b2_half_std_list, ~{
      # select all columns except for ID
      cols_to_normalize <- names(.x)[names(.x) %in% colnames(b1_metacov_0[,-1])]
      # print(cols_to_normalize)
      # apply rank normalization to each column
      .x[cols_to_normalize] <- lapply(.x[cols_to_normalize], rank_normalisation)
      .x[cols_to_normalize] <-as.matrix(.x[cols_to_normalize])%*%.y
      # return the updated data frame
      return(.x)
    })
    
    # 5.2 rescale with the original SDm (optional step)
    # metab_imp_cont_rescaled<-metab_imp_cont %*% metab_raw_std
    # colnames(metab_imp_cont_rescaled)<-colnames(metab_imp_cont)
    # metab_imp_cont_rescaled<-as.data.frame(metab_imp_cont_rescaled)
    # metab_imp_cont_rescaled$ID<-metab_imp_cont_ID$ID
    # metab_imp_cont<-metab_imp_cont_rescaled
  } else {
    # use the map()function from purrr package
    b1_imp_norm_list<-purrr::map2(b1_imp_list,b1_imp_std_list, ~{
      # select all columns except for ID
      cols_to_normalize <- names(.x)[names(.x) %in% colnames(b1_metacov_0[,-1])]
      # print(cols_to_normalize)
      # apply rank normalization to each column
      .x[cols_to_normalize] <- lapply(.x[cols_to_normalize], rank_normalisation)
      # return the updated data frame
      return(.x)
    })
    b2_imp_norm_list<-purrr::map2(b2_imp_list,b2_imp_std_list, ~{
      # select all columns except for ID
      cols_to_normalize <- names(.x)[names(.x) %in% colnames(b2_metacov_0[,-1])]
      # print(cols_to_normalize)
      # apply rank normalization to each column
      .x[cols_to_normalize] <- lapply(.x[cols_to_normalize], rank_normalisation)
      # return the updated data frame
      return(.x)
    })
    b1b2_half_norm_list<-purrr::map2(b1b2_half_list,b1b2_half_std_list, ~{
      # select all columns except for ID
      cols_to_normalize <- names(.x)[names(.x) %in% colnames(b1_metacov_0[,-1])]
      # print(cols_to_normalize)
      # apply rank normalization to each column
      .x[cols_to_normalize] <- lapply(.x[cols_to_normalize], rank_normalisation)
      # return the updated data frame
      return(.x)
    })
  }
} else {
  b1_imp_norm_list<-b1_imp_list
  b2_imp_norm_list<-b2_imp_list
  b1b2_half_norm_list<-b1b2_half_list
}

###########################
# read in phenotype data
pheno<-read.csv("data/metsleep_covariates_20220126.csv",header = T,sep=",") # add menopause variables
# recode empty cell value and "N" to missing
recode_cont_var<-c("LABA91","INSULIN_FAST","LABA69","LABA68","LABA66","SLEA2A_2401","SLEA1A_2401","SLEA2C_2401","SLEA1C_2401")
pheno[,recode_cont_var]<-sapply(pheno[,recode_cont_var], function(x){ifelse(x=="",NA,ifelse(x%in%c("N","L","Z","H"),NA,x))})
# recode missing values in categorical variables
recode_var<-c("OCEA13","SLEA4","SLEA5","SLEA6","SLEA7","SLEA8","SLEA11")
pheno[,recode_var]<-sapply(pheno[,recode_var], function(x){ifelse(is.na(x),NA,ifelse(x%in%c(1:7),x,NA))})
# Define the weight used in the analysis
pheno$WEIGHT<-pheno$WEIGHT_FINAL_NORM_OVERALL
pheno<-pheno%>% 
  unique()%>%
  dplyr::mutate(
    ALCOHOL_USE=factor(ALCOHOL_USE,levels=c(1,2,3),labels=c("never","former","current")),
    CIGARETTE_USE=factor(CIGARETTE_USE,levels=c(1,2,3),labels=c("never","former","current")),
    GENDER=factor(GENDER,levels = c("F","M"),labels = c("Female","Male")),
    OSA=ifelse(is.na(SLPA54)==0,ifelse(SLPA54>15,1,0),NA),
    OSA_status=factor(OSA,levels=c(0,1),labels = c("Non-OSA","OSA")),
    Per90=transformation(SLPA97),
    Per90_status=factor(Per90,levels=c(0,1),labels=c("Non-hypoxia","Hypoxia")),
    bkgrd1_clean=ifelse(BKGRD1_C7%in%c("","Q"),NA,as.character(BKGRD1_C7)),
    background=factor(bkgrd1_clean,levels=c(0:6),labels=c("Dominican","Central_American","Cuban","Mexican","Puerto_Rican", "South American","Multi")),
    # 0: Dominican
    # 1: Central American
    # 2: Cuban
    # 3: Mexican
    # 4: Puerto Rican
    # 5: South American
    # 6: More than one/Other heritage
    # Q: Unknown
    DYSLIPIDEMIA=factor(DYSLIPIDEMIA,levels=c(0,1),labels=c("No","Yes")),
    LABA91=ifelse(!is.na(LABA91),ifelse(LABA91%in%c("","N","L","Z","H"),NA,as.numeric(as.character(LABA91))),NA), # fix to as.numeric rounding up numbers
    INCIDENT_DM_V1V2=ifelse(is.na(DIABETES3),NA, # if baseline is missing, then incident DM is missing
                            ifelse(is.na(DIABETES3_INDICATOR_V2),# if follow up is missing, then incident DM is missing
                                   NA,
                                   ifelse(DIABETES3%in%c(1,2)&DIABETES3_INDICATOR_V2==1,1,0))), # if baseline is 1 or 2, while follow up indicator is 1, then incident =1, otherwise incident =0
    DIABETES3_INDICATOR=ifelse(!is.na(DIABETES3),ifelse(DIABETES3==3,1,0),NA),
    DIABETES3_INDICATOR=factor(DIABETES3_INDICATOR,levels=c(0,1),labels=c("No","Yes")),
    DIABETES3_INDICATOR_V2=factor(DIABETES3_INDICATOR_V2,levels=c(0,1),labels=c("No","Yes")),
    HTN=dplyr::case_when(
      HYPERTENSION2==1|SBPA6>90|SBPA5>140            ~"Yes",
      is.na(HYPERTENSION2)|is.na(SBPA6)|is.na(SBPA5) ~ NA_character_,
      TRUE                                          ~ "NO"
    ), # use hypertension2 instead
    OBSES=ifelse(!is.na(BMI),ifelse(BMI>=30,1,0),NA),
    OBSES_FCT=factor(OBSES,levels=c(0,1),labels=c("No","Yes")),
    HYPERTENSION=factor(HYPERTENSION2,levels=c(0,1),labels=c("No","Yes")), # overwrite hypertension with factor hypertension2
    HIGH_TOTAL_CHOL=factor(HIGH_TOTAL_CHOL,levels=c(0,1),labels=c("No","Yes")),
    MED_ANTIHYPERT=factor(MED_ANTIHYPERT,levels=c(0,1),labels=c("No","Yes")),
    MED_LLD=factor(MED_LLD,levels=c(0,1),labels=c("No","Yes")),
    INSULIN_FAST=ifelse(!is.na(INSULIN_FAST),ifelse(INSULIN_FAST%in%c("N","L",""),NA,as.numeric(INSULIN_FAST)),NA),
    LABA69=ifelse(!is.na(LABA69),ifelse(LABA69%in%c("N","L","","Z","H"),NA,as.numeric(as.character(LABA69))),NA), # fix to as.numeric rounding up numbers
    LABA68=ifelse(!is.na(LABA68),ifelse(LABA68%in%c("N","L","","Z","H"),NA,as.numeric(as.character(LABA68))),NA), # fix to as.numeric rounding up numbers
    LABA66=ifelse(!is.na(LABA66),ifelse(LABA66%in%c("N","L","","Z","H"),NA,as.numeric(as.character(LABA66))),NA), # fix to as.numeric rounding up numbers
    NO_MED=ifelse(!is.na(MUEA2),ifelse(MUEA2=="1",1,ifelse(MUEA2%in%c("S","2"),0,NA)),NA),
    MED_CURRENT_HRT=factor(MED_CURRENT_HRT,levels=c(0,1),labels=c("No","Yes")),
    MED_ANTIANGINAL=factor(MED_ANTIANGINAL,levels=c(0,1),labels=c("No","Yes")),
    MED_ANTIASTHMATICS=factor(MED_ANTIASTHMATICS,levels=c(0,1),labels=c("No","Yes")),
    MED_ANTIARRHYTHMICS=factor(MED_ANTIARRHYTHMICS,levels=c(0,1),labels=c("No","Yes")),
    MED_ANTICOAG=factor(MED_ANTICOAG,levels=c(0,1),labels=c("No","Yes")),
    MED_ANTIDIAB=factor(MED_ANTIDIAB,levels=c(0,1),labels=c("No","Yes")),
    MED_ASPIRIN=factor(MED_ASPIRIN,levels=c(0,1),labels=c("No","Yes")),
    MED_BB_OPHTHALMIC=factor(MED_BB_OPHTHALMIC,levels=c(0,1),labels=c("No","Yes")),
    MED_CARDIACGLYCOSIDES=factor(MED_CARDIACGLYCOSIDES,levels=c(0,1),labels=c("No","Yes")),
    MED_CHEMO=factor(MED_CHEMO,levels=c(0,1),labels=c("No","Yes")),
    MED_FIBARES_NICOACID=factor(MED_FIBARES_NICOACID,levels=c(0,1),labels=c("No","Yes")),
    MED_NSAID=factor(MED_NSAID,levels=c(0,1),labels=c("No","Yes")),
    MED_OI_STEROID=factor(MED_OI_STEROID,levels=c(0,1),labels=c("No","Yes")),
    MED_ANTIANXI=factor(MED_ANTIANXI,levels=c(0,1),labels=c("No","Yes")),
    MED_ANTIDEPRESS=factor(MED_ANTIDEPRESS,levels=c(0,1),labels=c("No","Yes")),
    HTN_NEW=ifelse(!is.na(HYPERTENSION2),
                   ifelse(!is.na(HYPERTENSION2_V2),
                          ifelse(HYPERTENSION2==0&HYPERTENSION2_V2==1,1,0),
                          0),
                   NA)
  )
#######
# 12.2.2020
# clean min_o2, avg_o2: if >100 set to missing, if <=0 set to missing
pheno[which(pheno$SLPA91>100),"SLPA91"]<-100
pheno[which(pheno$SLPA92<=0),"SLPA92"]<-NA

# Define the weight used in the analysis
pheno$WEIGHT<-pheno$WEIGHT_FINAL_NORM_OVERALL

#key to convert hours to radians
conv <- 2*pi/24 ## hrs -> radians

# derive sleep traits
pheno<-pheno%>%
  tidyr::separate(SLEA2A_2401,c("wdwake_hr","wdwake_mn"),":",remove = FALSE)%>%
  tidyr::separate(SLEA1A_2401,c("wdbed_hr","wdbed_mn"),":",remove = FALSE)%>%
  tidyr::separate(SLEA2C_2401,c("wewake_hr","wewake_mn"),":",remove = FALSE)%>%
  tidyr::separate(SLEA1C_2401,c("webed_hr","webed_mn"),":",remove = FALSE)%>%
  dplyr::mutate(wdwake_time=as.numeric(wdwake_hr)+as.numeric(wdwake_mn)/60, # convert character timet to numeric time (in hr)
                wdbed_time=as.numeric(wdbed_hr)+as.numeric(wdbed_mn)/60,    # convert character timet to numeric time (in hr)
                wewake_time=as.numeric(wewake_hr)+as.numeric(wewake_mn)/60, # convert character timet to numeric time (in hr)
                webed_time=as.numeric(webed_hr)+as.numeric(webed_mn)/60,    # convert character timet to numeric time (in hr)
                wddur=ifelse(is.na(wdbed_time)|is.na(wdwake_time),NA,ifelse(wdbed_time>wdwake_time,wdwake_time+24-wdbed_time,wdwake_time-wdbed_time)), # if wake time>bed time then add 24
                wedur=ifelse(is.na(webed_time)|is.na(wewake_time),NA,ifelse(webed_time>wewake_time,wewake_time+24-webed_time,wewake_time-webed_time)), # if wake time>bed time then add 24
                wdmid=ifelse(is.na(wdwake_time)|is.na(wddur),NA,ifelse(wdbed_time+0.5*wddur>24,wdbed_time+0.5*wddur-24,wdbed_time+0.5*wddur)),
                wemid=ifelse(is.na(wewake_time)|is.na(wedur),NA,ifelse(webed_time+0.5*wedur>24,webed_time+0.5*wedur-24,webed_time+0.5*wedur)),
                we_minus_wd=ifelse(is.na(wdbed_time)|is.na(wdwake_time)|is.na(webed_time)|is.na(wewake_time),NA,ifelse(wemid-wdmid<(-12),wemid-wdmid+24,ifelse(wemid-wdmid>12,wemid-wdmid-24,wemid-wdmid))),
                wdwake_hrs=(as.numeric(lubridate::hm(SLEA2A_2401))/(60*60)),
                wdbed_hrs=(as.numeric(lubridate::hm(SLEA1A_2401))/(60*60)),
                wewake_hrs=(as.numeric(lubridate::hm(SLEA2C_2401))/(60*60)),
                webed_hrs=(as.numeric(lubridate::hm(SLEA1C_2401))/(60*60)),
                wdwake_rad=circular::circular(conv*as.numeric(wdwake_hr)),
                wdbed_rad=circular::circular(conv*as.numeric(wdbed_hr)),
                wewake_rad=circular::circular(conv*as.numeric(wewake_hr)),
                webed_rad=circular::circular(conv*as.numeric(webed_hr)),
                wdmid_rad=circular::circular(conv*as.numeric(wdmid)),
                wemid_rad=circular::circular(conv*as.numeric(wemid)),
                we_minus_wd_rad=circular::circular(conv*as.numeric(we_minus_wd)),
                wdwake_time_shift=ifelse(!is.na(wdwake_time),ifelse(wdwake_time<12, (wdwake_time+24), wdwake_time),NA),
                wewake_time_shift=ifelse(!is.na(wewake_time),ifelse(wewake_time<12, (wewake_time+24), wewake_time),NA),
                wdbed_time_shift=ifelse(!is.na(wdbed_time),ifelse(wdbed_time<12, (wdbed_time+24), wdbed_time),NA),
                webed_time_shift=ifelse(!is.na(webed_time),ifelse(webed_time<12, (webed_time+24), webed_time),NA),
                wdmid_shift=ifelse(!is.na(wdmid),ifelse(wdmid<12, (wdmid+24), wdmid),NA),
                wemid_time_shift=ifelse(!is.na(wemid),ifelse(wemid<12, (wemid+24), wemid),NA),
                we_minus_wd_shift=ifelse(!is.na(we_minus_wd),ifelse(we_minus_wd<12, (we_minus_wd+24), we_minus_wd),NA),
                osa_ess=ifelse(is.na(SLPA54)|is.na(ESS),NA,ifelse(SLPA54>15&ESS>10,1,0)),
                slp_short_wd=ifelse(!is.na(wddur),ifelse(wddur<=5,1,0),NA),
                slp_long_wd=ifelse(!is.na(wddur),ifelse(wddur==9,1,0),NA),
                snore=ifelse(SLEA13==4,1,ifelse(SLEA13%in%c(1,2,3),0,NA)), #always or almost always (6-7 nights a week)=1, less than 6-7 nights a weeek=0,dont know and others =NA
                pill=ifelse(!is.na(SLEA8),ifelse(SLEA8%in%c("4","5",4,5),1,0),NA), # 3 or more times a week=1
                restless=ifelse(!is.na(SLEA11),ifelse(SLEA11%in%c("3","4",3,4),1,0),NA), # restless or very restless=1
                backtosleep=ifelse(!is.na(SLEA7),ifelse(SLEA7%in%c("4","5",4,5),1,0),NA), # 3 or more times a week=1
                earlywake=ifelse(!is.na(SLEA6),ifelse(SLEA6%in%c("4","5",4,5),1,0),NA), # 3 or more times a week=1
                freqwake=ifelse(!is.na(SLEA5),ifelse(SLEA5%in%c("4","5",4,5),1,0),NA), # 3 or more times a week=1
                fallasleep=ifelse(!is.na(SLEA4),ifelse(SLEA4%in%c("4","5",4,5),1,0),NA), # 3 or more times a week=1
                CKD=ifelse(!is.na(CKD),ifelse(CKD==5,1,0),NA),
                CKD_fct=factor(CKD,levels=c(0,1),labels=c("No","Yes")),
                ASTHMA_CURR=factor(ASTHMA_CURR,levels=c(0,1),labels=c("No","Yes")),
                COPD_EVER=factor(COPD_EVER,levels=c(0,1),labels=c("No","Yes")),
                short_sleep_wd=factor(slp_short_wd,levels=c(0,1),labels=c("No","Yes")),
                long_sleep_wd=factor(slp_long_wd,levels=c(0,1),labels=c("No","Yes"))
  )

# change variable names to uppercase
colnames(pheno)<-tolower(colnames(pheno))

# combine batch 1 and 2
b1b2_imp_norm_list<-purrr::map2(b1_imp_norm_list,b2_imp_norm_list, ~{
  .x["batch"]<-"batch1"
  .y["batch"]<-"batch2"
  combine_df<-rbind(.x,.y)
  return(combine_df)
})
b1b2_half_norm_list[[1]]$batch<-"batch1"
b1b2_half_norm_list[[2]]$batch<-"batch2"

# merge metabolomic data with phenotype data
b1b2_imp_pheno_list<-list()
for (m in seq_along(b1b2_imp_norm_list)){
  b1b2_imp_pheno_list[[m]]<-merge(b1b2_imp_norm_list[[m]],pheno,by.y="id",by.x="SOL_ID",all.y=T)
}
b1b2_half_pheno<-rbind(b1b2_half_norm_list[[1]],b1b2_half_norm_list[[2]])%>%
  merge(.,pheno,by.y="id",by.x="SOL_ID",all.y=T)

b1b2_imp_mdl<-list()
for (m in seq_along(b1b2_imp_pheno_list)){
  # multiple imputation datasets
  b1b2_imp_mdl[[m]]<-atlas_cir_svyreg_loop_nonlinear(data=b1b2_imp_pheno_list[[m]],covar=covars,end=ncol(b1_imp_norm_list[[m]]),metab_is_cont=T,metab_is_complete=T,trait=sleep_trait,trait_as_predictor=T,model_input=nl_tbl)
}
b1b2_imp_mdl_append<-dplyr::bind_rows(b1b2_imp_mdl, .id = "imp_dataset")
# unnest the list column to get the combined sin estimates in separate columns
b1b2_imp_mdl_sin_combine <- b1b2_imp_mdl_append %>%
  dplyr::group_by(metabolite) %>%
  dplyr::summarise(rubins_result = list(rubins_rule(ests = sin_beta, ses = sin_se, round_digit = 2))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(metabolite = factor(metabolite, levels = unique(metabolite))) %>%
  tidyr::unnest_wider(rubins_result) %>%
  dplyr::mutate(sin_beta=est,
                sin_se=SE,
                sin_p_val=pval,
                n=unique(b1b2_imp_mdl_append$n))
# unnest the list column to get the combined cos estimates in separate columns
b1b2_imp_mdl_cos_combine <- b1b2_imp_mdl_append %>%
  dplyr::group_by(metabolite) %>%
  dplyr::summarise(rubins_result = list(rubins_rule(ests = cos_beta, ses = cos_se, round_digit = 2))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(metabolite = factor(metabolite, levels = unique(metabolite))) %>%
  tidyr::unnest_wider(rubins_result) %>%
  dplyr::mutate(cos_beta=est,
                cos_se=SE,
                cos_p_val=pval,
                n=unique(b1b2_imp_mdl_append$n))
b1b2_imp_mdl_p_combine <- b1b2_imp_mdl_append %>%
  dplyr::group_by(metabolite) %>%
  dplyr::summarise(p_acat = ACAT(p_val)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(metabolite = factor(metabolite, levels = unique(metabolite))) 

b1b2_imp_mdl_combine<-merge(b1b2_imp_mdl_sin_combine[,c("metabolite","sin_beta","sin_se","sin_p_val","n")],
                            b1b2_imp_mdl_cos_combine[,c("metabolite","cos_beta","cos_se","cos_p_val")],by="metabolite")%>%
  merge(.,b1b2_imp_mdl_p_combine,by="metabolite")

# xenobiotic metabolites
b1b2_half_mdl<-atlas_cir_svyreg_loop_nonlinear(data=b1b2_half_pheno,covar=covars,end=ncol(b1b2_half_list[[1]]),metab_is_cont=T,metab_is_complete=T,trait=sleep_trait,trait_as_predictor=T,model_input=nl_tbl)%>%
  dplyr::mutate(imp_dataset="half",
                p_acat=p_val)
# combine xenobiotic and nonxenobiotic results
b1b2_imp_mdl_total<-rbind(b1b2_imp_mdl_combine,b1b2_half_mdl[,c("metabolite","n","sin_beta","sin_se","sin_p_val","cos_beta","cos_se","cos_p_val","p_acat")])%>%
  dplyr::mutate(p_val_fdr = p.adjust(p_acat, method = "BH"),
                trait=sleep_trait,
                model=model,
                covar=cir_input_grid[comb_number,"covariates"],
                strata="both")

# stratified analysis
# non-xenobiotic metabolites
b1b2_imp_strat_mdl<-list()
for (m in seq_along(b1b2_imp_pheno_list)){
  # multiple imputation datasets
  b1b2_imp_strat_mdl[[m]]<-atlas_cir_strat_svyreg_loop_nonlinear(data=b1b2_imp_pheno_list[[m]],covar=covars,end=ncol(b1_imp_norm_list[[m]]),metab_is_cont=T,metab_is_complete=T,
                                                       trait_for_model="original",trait=sleep_trait,trait_as_predictor=T,stratifier=stratifier,model_input=nl_tbl)
}
b1b2_imp_strat_mdl_append<-dplyr::bind_rows(b1b2_imp_strat_mdl, .id = "imp_dataset")

# unnest the list column to get the combined sin estimates in separate columns
# female
b1b2_imp_strat_female_mdl_sin_combine <- b1b2_imp_strat_mdl_append[which(b1b2_imp_strat_mdl_append$strata=="Female"),]  %>%
  dplyr::group_by(metabolite) %>%
  dplyr::summarise(rubins_result = list(rubins_rule(ests = sin_beta, ses = sin_se, round_digit = 2))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(metabolite = factor(metabolite, levels = unique(metabolite))) %>%
  tidyr::unnest_wider(rubins_result) %>%
  dplyr::mutate(sin_beta=est,
                sin_se=SE,
                sin_p_val=pval,
                strata="Female",
                n=unique(b1b2_imp_strat_mdl_append[which(b1b2_imp_strat_mdl_append$strata=="Female"),"n"]))
# unnest the list column to get the combined cos estimates in separate columns
b1b2_imp_strat_female_mdl_cos_combine <- b1b2_imp_strat_mdl_append[which(b1b2_imp_strat_mdl_append$strata=="Female"),] %>%
  dplyr::group_by(metabolite) %>%
  dplyr::summarise(rubins_result = list(rubins_rule(ests = cos_beta, ses = cos_se, round_digit = 2))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(metabolite = factor(metabolite, levels = unique(metabolite))) %>%
  tidyr::unnest_wider(rubins_result) %>%
  dplyr::mutate(cos_beta=est,
                cos_se=SE,
                cos_p_val=pval,
                strata="Female",
                n=unique(b1b2_imp_strat_mdl_append[which(b1b2_imp_strat_mdl_append$strata=="Female"),"n"]))
b1b2_imp_strat_female_mdl_p_combine <- b1b2_imp_strat_mdl_append[which(b1b2_imp_strat_mdl_append$strata=="Female"),]  %>%
  dplyr::group_by(metabolite) %>%
  dplyr::summarise(p_acat = ACAT(p_val)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(metabolite = factor(metabolite, levels = unique(metabolite))) 

b1b2_imp_strat_female_mdl_combine<-merge(b1b2_imp_strat_female_mdl_sin_combine[,c("metabolite","sin_beta","sin_se","sin_p_val","strata","n")],
                            b1b2_imp_strat_female_mdl_cos_combine[,c("metabolite","cos_beta","cos_se","cos_p_val")],by="metabolite")%>%
  merge(.,b1b2_imp_strat_female_mdl_p_combine,by="metabolite")
# Male
b1b2_imp_strat_male_mdl_sin_combine <- b1b2_imp_strat_mdl_append[which(b1b2_imp_strat_mdl_append$strata=="Male"),]  %>%
  dplyr::group_by(metabolite) %>%
  dplyr::summarise(rubins_result = list(rubins_rule(ests = sin_beta, ses = sin_se, round_digit = 2))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(metabolite = factor(metabolite, levels = unique(metabolite))) %>%
  tidyr::unnest_wider(rubins_result) %>%
  dplyr::mutate(sin_beta=est,
                sin_se=SE,
                sin_p_val=pval,
                strata="male",
                n=unique(b1b2_imp_strat_mdl_append[which(b1b2_imp_strat_mdl_append$strata=="Male"),"n"]))
# unnest the list column to get the combined cos estimates in separate columns
b1b2_imp_strat_male_mdl_cos_combine <- b1b2_imp_strat_mdl_append[which(b1b2_imp_strat_mdl_append$strata=="Male"),] %>%
  dplyr::group_by(metabolite) %>%
  dplyr::summarise(rubins_result = list(rubins_rule(ests = cos_beta, ses = cos_se, round_digit = 2))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(metabolite = factor(metabolite, levels = unique(metabolite))) %>%
  tidyr::unnest_wider(rubins_result) %>%
  dplyr::mutate(cos_beta=est,
                cos_se=SE,
                cos_p_val=pval,
                strata="male",
                n=unique(b1b2_imp_strat_mdl_append[which(b1b2_imp_strat_mdl_append$strata=="Male"),"n"]))
b1b2_imp_strat_male_mdl_p_combine <- b1b2_imp_strat_mdl_append[which(b1b2_imp_strat_mdl_append$strata=="Male"),]  %>%
  dplyr::group_by(metabolite) %>%
  dplyr::summarise(p_acat = ACAT(p_val)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(metabolite = factor(metabolite, levels = unique(metabolite))) 

b1b2_imp_strat_male_mdl_combine<-merge(b1b2_imp_strat_male_mdl_sin_combine[,c("metabolite","sin_beta","sin_se","sin_p_val","strata","n")],
                                       b1b2_imp_strat_male_mdl_cos_combine[,c("metabolite","cos_beta","cos_se","cos_p_val")],by="metabolite")%>%
  merge(.,b1b2_imp_strat_male_mdl_p_combine,by="metabolite")

# xenobiotic metabolites
b1b2_half_strat_mdl<-atlas_cir_strat_svyreg_loop_nonlinear(data=b1b2_half_pheno,covar=covars,end=ncol(b1b2_half_list[[1]]),metab_is_cont=T,metab_is_complete=T,
                                                 trait_for_model="original",trait=sleep_trait,trait_as_predictor=T,stratifier=stratifier,model_input=nl_tbl)%>%
  dplyr::mutate(imp_dataset="half",
                p_acat=p_val)

b1b2_imp_strat_female_mdl_total<-rbind(b1b2_imp_strat_female_mdl_combine[,c("metabolite","sin_beta","sin_se","sin_p_val","cos_beta","cos_se","cos_p_val","p_acat","strata","n")],b1b2_half_strat_mdl[which(b1b2_half_strat_mdl$strata=="Female"),c("metabolite","sin_beta","sin_se","sin_p_val","cos_beta","cos_se","cos_p_val","p_acat","strata","n")])%>%
  dplyr::mutate(p_val_fdr = p.adjust(p_acat, method = "BH"))
b1b2_imp_strat_male_mdl_total<-rbind(b1b2_imp_strat_male_mdl_combine[,c("metabolite","sin_beta","sin_se","sin_p_val","cos_beta","cos_se","cos_p_val","p_acat","strata","n")],b1b2_half_strat_mdl[which(b1b2_half_strat_mdl$strata=="Male"),c("metabolite","sin_beta","sin_se","sin_p_val","cos_beta","cos_se","cos_p_val","p_acat","strata","n")])%>%
  dplyr::mutate(p_val_fdr = p.adjust(p_acat, method = "BH"))

b1b2_imp_strat_mdl_total<-rbind(b1b2_imp_strat_female_mdl_total,b1b2_imp_strat_male_mdl_total)%>%
  dplyr::mutate(
    trait=sleep_trait,
    model=model,
    covar=cir_input_grid[comb_number,"covariates"])
b1b2_imp_mdl_total<-rbind(b1b2_imp_mdl_total,b1b2_imp_strat_mdl_total)
write.csv(b1b2_imp_mdl_total,file =paste0("output/circular/circular_main_",comb_number,".csv"),row.names = F)
