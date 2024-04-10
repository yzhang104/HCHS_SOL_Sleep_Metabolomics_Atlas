# function to check installed package, load packages
instant_pkgs<- function(pkgs) { 
  pkgs_miss <- pkgs[which(!pkgs %in% installed.packages()[, 1])]
  if (length(pkgs_miss) > 0) {
    install.packages(pkgs_miss,repos = "http://cran.us.r-project.org", method = "curl")
  }
  
  if (length(pkgs_miss) == 0) {
    message("\n ...Packages were already installed!\n")
  }
  
  # install packages not already loaded:
  pkgs_miss <- pkgs[which(!pkgs %in% installed.packages()[, 1])]
  if (length(pkgs_miss) > 0) {
    install.packages(pkgs_miss,repos = "http://cran.us.r-project.org", method = "curl")
  }
  
  # load packages not already loaded:
  attached <- search()
  attached_pkgs <- attached[grepl("package", attached)]
  need_to_attach <- pkgs[which(!pkgs %in% gsub("package:", "", attached_pkgs))]
  
  if (length(need_to_attach) > 0) {
    for (i in 1:length(need_to_attach)) require(need_to_attach[i], character.only = TRUE)
  }
  
  if (length(need_to_attach) == 0) {
    message("\n ...Packages were already loaded!\n")
  }
}

# Function to dichotomize trait continuous variables
dichotmize_func<-function(x,cutoff_value,case) {
  if (!is.null(cutoff_value)) {
    if (case=="<") {
      case_index<-which(x<cutoff_value)
    } else if (case==">") {
      case_index<-which(x>cutoff_value)
    } else if (case%in%c("<=","=<")) {
      case_index<-which(x<=cutoff_value)
    } else if (case%in%c(">=","=>")) {
      case_index<-which(x>=cutoff_value)
    } else if (case=="=") {
      case_index<-which(x=cutoff_value)
    } else stop('Only inequality symbols are allowed for trait_dich_case parameter')
    
    val <- rep(0, length(x))
    na_index<-which(is.na(x))
    val[na_index] <- NA
    val[case_index]<-1
    return(val)
  }
  else val<-x
  return(val)
}

# Function to transform oxygen saturation below 90%
transformation <- function(x){
  val <- rep(0, length(x))
  val[which(x >=5)] <- 1
  val[which(is.na(x))] <- NA
  return(val)
}

# Function to summarize metabolites distribution
distribution_func<-function(seed,data,notes) {
  set.seed(seed)
  random_col<-sample(ncol(data),50,replace=F)
  # boxplot of 50 random metabolites 
  plot1<-ggplot(stack(data[,random_col]),aes(x=ind,y=values))+theme(axis.text=element_text(size=7))+
    geom_boxplot()+ coord_flip()+
    scale_y_continuous(labels = function(x) format(x, scientific = TRUE))+
    labs(title=paste("Boxplot for metabolites concentrations\n",notes))
  # Density plot of all means
  plot2<-ggplot(data=data.frame(mean=apply(data, 2, mean, na.rm=TRUE)),aes(x=mean))+geom_density(aes(x=mean, y=..scaled..))+
    labs(title=paste("Density plot of average metabolites concentrations\n",notes))
  output_list<-list(boxplot=plot1,densityplot=plot2)
  # class(output_list)="three_plots"
  # multiplot(plot1,plot2,plot3,cols=1)
  return(output_list)
}

# Function to impute missing values
impute_func<-function(data,method){
  # mean: replace with mean
  # median: replace with median
  # min: replace with min
  # half: replace with half of the min value
  # knn: k-nearest neighbour - HCHS/SOL doesn't have any complete rows (all obs have some missing values), knn can't be implemented
  # rf: random forest
if (method=="mean"){
  imp_data<-data.frame(sapply(data,function(x) ifelse(is.na(x),mean(x, na.rm = TRUE),x)),check.names = F,stringsAsFactors = F)
} else if (method=="median") {
    imp_data<-data.frame(sapply(data,function(x) ifelse(is.na(x),median(x, na.rm = TRUE),x)),check.names = F,stringsAsFactors = F)
  }  else if (method=="min") {
    imp_data<-data.frame(sapply(data,function(x) ifelse(is.na(x),min(x, na.rm = TRUE),x)),check.names = F,stringsAsFactors = F)
  } else if (method=="rf") {
    imp_data<-missForest::missForest(data, maxiter = 10, ntree = 100)
  } else if (method=="half") {
    imp_data<-data.frame(sapply(data,function(x) ifelse(is.na(x),0.5*min(x, na.rm = TRUE),x)),check.names = F,stringsAsFactors = F)
  } else stop('Method not supported')
  # else if (method=="knn") {
  #   imp_data<-VIM::kNN(data, variable = colnames(data), k = 5, numFun = median)
  # }
  return(imp_data)
}

# Function to rank normalize var
rank_normalisation<-function(x){
  qnorm((rank(x,na.last="keep")-0.5)/length(x))
}


# function to export the coefficient and map to chemical names
export_lasso_coeff<-function(lasso_function_output,min_col_number,key_col,drop_covariates,metab_reference_table,export_path){
  if (length(lasso_function_output$coef_lasso)>3){
    lasso_coeff<-data.frame(lasso_function_output$coef_lasso)
    lasso_coeff$feature_id<-rownames(lasso_coeff)
    colnames(lasso_coeff)<-c("coeff",key_col)
    lasso_coeff<-lasso_coeff[which(!lasso_coeff[,key_col]%in%c("(Intercept)",drop_covariates)),] # drop coefficients for covariates
    lasso_coeff<-merge(lasso_coeff,metab_reference_table,by=key_col)
    write.csv(lasso_coeff,file=paste0(export_path),row.names = F)
    return(lasso_coeff)
  }
}

# function to calculate a summary score based on metabolite concentration and coefficient table
summary_score<-function(data, coeff){
  final<-data.frame()
  for (i in 1:nrow(coeff)){
    # print(i)
    component<-coeff[i,1]
    coefficient<-coeff[i,2]
    temp<-data.frame(metab=as.numeric(data[,colnames(data)%in%component]),coef=rep(coefficient,nrow(data)))
    temp$product<-temp$metab*temp$coef
    
    if (i==1){
      final<-data.frame(temp[,3])
    } else {
      final<-data.frame(cbind(final,temp[,3]))
    }
    colnames(final)[i]<-component
    # print(i)
  }
  if (ncol(final)==1) {
    score<-final
    colnames(score)[1]<-"summary"
    } else {score<-data.frame(summary=rowSums(final,na.rm = T))}
  return(score)
}

# Function to loop regression model and export model summary using survey weight
svyreg_loop<-function(data,covar,end,metab_is_cont,metab_is_complete,trait_original,trait_binary,trait_for_model,trait_as_predictor){
  dat<-data[,1:end]
  dat<-dat[,!colnames(dat)%in%c("LAB_ID","ID","SOL_ID","id")]
  out_nvar=ncol(dat)
  out_beta=rep(NA,out_nvar)
  out_se=rep(NA,out_nvar)
  out_pvalue=rep(NA,out_nvar)
  out_nobs=rep(NA,out_nvar)
  if (trait_for_model=="binary") {
    trait<-"TraitBinary"
  } else if (trait_for_model=="original"){
    trait<-trait_original
  } else { trait<-"Trait" }

    # if the original trait variable is binary or the user pick the binary version to use in the model
  trait_is_cont=ifelse(length(unique(data[!is.na(data[,trait]),trait]))==2,FALSE,
                       ifelse(trait_for_model=="binary",FALSE,TRUE))

  
  for(i in 1:(ncol(dat))) {
    # met_df<-cbind(data[,i+2],data[,c("ID",covar)])
    if (!"ID"%in%colnames(data)){
       data<-data%>%
        dplyr::rename("ID"="SOL_ID")
    }
    met_df<-cbind(data[,colnames(dat)[i]],data[,c("ID",covar)])
    met_df_id<-met_df[complete.cases(met_df),"ID"]
    require(survey)

    if (all(grep("[A-Z]", names(data[,-1])) == 0)) { # if column names are all in lower case
      survey_design=svydesign(id=~psu_id, strata=~strat, weights=~weight, data=data)
    } else {
      survey_design=svydesign(id=~PSU_ID, strata=~STRAT, weights=~WEIGHT, data=data)
    }
    
    if (trait_as_predictor==T){
      if (metab_is_cont==T){
        model<- svyglm(as.formula(paste0(colnames(dat)[i],"~",trait,"+",paste(covar,collapse= "+"))),design=subset(survey_design,ID%in%met_df_id))
      } else {
        model<- svyglm(as.formula(paste0(colnames(dat)[i],"~",trait,"+",paste(covar,collapse = "+"))),design=subset(survey_design,ID%in%met_df_id),family=quasipoisson(link='log'))
      }
    } else {
      if (trait_is_cont==T){
        model<- svyglm(as.formula(paste0(trait,"~",colnames(dat)[i],"+",paste(covar,collapse= "+"))),design=subset(survey_design,ID%in%met_df_id))
      } else {
        model<- svyglm(as.formula(paste0(trait,"~",colnames(dat)[i],"+",paste(covar,collapse = "+"))),design=subset(survey_design,ID%in%met_df_id),family=quasipoisson(link='log'))
      }
    }

    Vcov <- vcov(model, useScale = FALSE)
    beta<- coef(model)
    se<- sqrt(diag(vcov(model, useScale = FALSE)))
    zval<- beta / se
    pval<- 2 * pnorm(abs(zval), lower.tail = FALSE)
    out_beta[i]=as.numeric(beta[2])
    out_se[i] = as.numeric(se[2])
    out_pvalue[i] = as.numeric(pval[2])
    out_nobs[i]=as.numeric(nobs(model))
  }
  regress_output<-data.frame(metabolite=colnames(dat),
                             beta=out_beta,
                             se=out_se,
                             n=out_nobs,
                             p_val=out_pvalue,
                             p_val_Bonf=p.adjust(out_pvalue,method="bonferroni"),
                             p_val_fdr=p.adjust(out_pvalue,method="BH")
  )%>%
    dplyr::mutate(
      pval_Bonf_neglog=-log10(p_val_Bonf),
      pval_fdr_neglog=-log10(p_val_fdr),
      sig_Bonf=ifelse(p_val_Bonf<=0.05,"Bonf<0.05","Not Sig"),
      sig_fdr=ifelse(p_val_fdr<=0.05,"FDR<0.05","Not Sig"),
      is_continuous=metab_is_cont
    )
  identified_data<-regress_output%>%
    dplyr::filter(!stringr::str_detect(metabolite,"^X"))%>%
    dplyr::mutate(p_val_fdr_named=p.adjust(p_val,method="BH"),
           sig_fdr_named=ifelse(p_val_fdr_named<=0.05,"FDR<0.05","Not Sig"))
  unidentified_data<-regress_output%>%
    dplyr::filter(stringr::str_detect(metabolite,"^X"))%>%
    dplyr::mutate(p_val_fdr_named=NA,
           sig_fdr_named=NA)
  regress_output<-rbind(identified_data,unidentified_data)
  if (metab_is_complete==T){
    regress_output$is_continuous<-"Complete-cases"
  } else {
    regress_output$is_continuous<-factor(regress_output$is_continuous,levels=c(TRUE,FALSE),labels=c("Imputed","Dichotomized"))
    }
  
  return(regress_output)
}

# Function to loop the interaction regression model and export model summary using survey weight
svyreg_loop_interaction<-function(data,covar,end,metab_is_cont,metab_is_complete,trait_original,trait_binary,trait_for_model,trait_as_predictor,interaction){
  dat<-data[,1:end]
  dat<-dat[,!colnames(dat)%in%c("LAB_ID","ID")]
  if (is.vector(dat)){
    dat<-data.frame(dat)
    colnames(dat)[1]<-colnames(data)[2]
  }
  out_nvar=ncol(dat) 
  out_beta=rep(NA,out_nvar)
  out_se=rep(NA,out_nvar)
  out_pvalue=rep(NA,out_nvar)
  out_nobs=rep(NA,out_nvar)
  int_beta=rep(NA,out_nvar)
  int_se=rep(NA,out_nvar)
  int_pvalue=rep(NA,out_nvar)
  cross_beta=rep(NA,out_nvar)
  cross_se=rep(NA,out_nvar)
  cross_pvalue=rep(NA,out_nvar)
  
  # drop the interaction term from the covariate list
  covar2<-covar[!covar%in%interaction]
  
  if (trait_for_model=="binary") {
    trait<-"TraitBinary"
  }  else if (trait_for_model=="original"){
    trait<-trait_original
  } else { trait<-"Trait" }
  
  # if the original trait variable is binary or the user pick the binary version to use in the model
  trait_is_cont=ifelse(length(unique(data[!is.na(data[,trait]),trait]))==2,FALSE,
                       ifelse(trait_for_model=="binary",FALSE,TRUE))
  
  if (all(grep("[A-Z]", names(data[,-1])) == 0)) { # if column names are all in lower case
    survey_design=svydesign(id=~psu_id, strata=~strat, weights=~weight, data=data)
  } else {
    survey_design=svydesign(id=~PSU_ID, strata=~STRAT, weights=~WEIGHT, data=data)
  }

  for(i in 1:(ncol(dat))) {
    met_df<-cbind(data[,colnames(dat)[i]],data[,c("ID",covar)])
    met_df_id<-met_df[complete.cases(met_df),"ID"]
    if (trait_as_predictor==T){
      if (metab_is_cont==T){
        model<- svyglm(as.formula(paste0(colnames(dat)[i],"~",trait,"+",paste(covar2,collapse= "+"),"+",trait,"*",interaction)),design=subset(survey_design,ID%in%met_df_id))
      } else {
        model<- svyglm(as.formula(paste0(colnames(dat)[i],"~",trait,"+",paste(covar2,collapse = "+"),"+",trait,"*",interaction)),design=subset(survey_design,ID%in%met_df_id),family=quasipoisson(link='log'))
      }
    } else {
      if (trait_is_cont==T){
        model<- svyglm(as.formula(paste0(trait,"~",colnames(dat)[i],"+",paste(covar2,collapse= "+"),"+",colnames(dat)[i],"*",interaction)),design=subset(survey_design,ID%in%met_df_id))
      } else {
        model<- svyglm(as.formula(paste0(trait,"~",colnames(dat)[i],"+",paste(covar2,collapse = "+"),"+",colnames(dat)[i],"*",interaction)),design=subset(survey_design,ID%in%met_df_id),family=quasipoisson(link='log'))
      }
    }
    
    Vcov <- vcov(model, useScale = FALSE)
    beta<- coef(model)
    se<- sqrt(diag(vcov(model, useScale = FALSE)))
    zval<- beta / se
    pval<- 2 * pnorm(abs(zval), lower.tail = FALSE)
    out_beta[i]=as.numeric(beta[2])
    out_se[i] = round(as.numeric(se[2]),digits = 3)
    out_pvalue[i] = as.numeric(pval[2])
    out_nobs[i]=as.numeric(nobs(model))
    int_beta[i]=as.numeric(beta[length(beta)-1])
    int_se[i] = round(as.numeric(se[length(se)-1]),digits = 3)
    int_pvalue[i] = round(as.numeric(pval[length(pval)-1]),digits=9)
    cross_beta[i]=as.numeric(beta[length(beta)])
    cross_se[i] = round(as.numeric(se[length(se)]),digits = 3)
    cross_pvalue[i] = as.numeric(pval[length(pval)])
  }
  regress_output<-data.frame(metabolite=colnames(dat),
                             trait_beta=out_beta,
                             trait_se=out_se,
                             n=out_nobs,
                             trait_p_val=out_pvalue,
                             trait_p_val_Bonf=p.adjust(out_pvalue,method="bonferroni"),
                             trait_p_val_fdr=p.adjust(out_pvalue,method="BH"),
                             int_main_beta=int_beta,
                             int_main_se=int_se,
                             int_main_p_val=int_pvalue,
                             # int_main_p_val_Bonf=p.adjust(int_pvalue,method="bonferroni"),
                             # int_main_p_val_fdr=p.adjust(int_pvalue,method="BH"),
                             int_cross_beta=cross_beta,
                             int_cross_se=cross_se,
                             int_cross_p_val=cross_pvalue,
                             # int_cross_p_val_Bonf=p.adjust(cross_pvalue,method="bonferroni"),
                             # int_cross_p_val_fdr=p.adjust(cross_pvalue,method="BH"),
                             interaction_term=names(beta)[length(beta)-1],
                             cross_term=names(beta)[length(beta)]
                             
  )%>%
    dplyr::mutate(
      pval_Bonf_neglog=-log10(trait_p_val_Bonf),
      pval_fdr_neglog=-log10(trait_p_val_fdr),
      sig_Bonf=ifelse(trait_p_val_Bonf<=0.05,"Bonf<0.05","Not Sig"),
      sig_fdr=ifelse(trait_p_val_fdr<=0.05,"FDR<0.05","Not Sig"),
      is_continuous=metab_is_cont
    )
  identified_data<-regress_output%>%
    dplyr::filter(!stringr::str_detect(metabolite,"^X"))%>%
    dplyr::mutate(p_val_fdr_named=p.adjust(trait_p_val,method="BH"),
                  sig_fdr_named=ifelse(p_val_fdr_named<=0.05,"FDR<0.05","Not Sig"))
  unidentified_data<-regress_output%>%
    dplyr::filter(stringr::str_detect(metabolite,"^X"))%>%
    dplyr::mutate(p_val_fdr_named=NA,
                  sig_fdr_named=NA)
  regress_output<-rbind(identified_data,unidentified_data)
  if (metab_is_complete==T){
    regress_output$is_continuous<-"Complete-cases"
  } else {
    regress_output$is_continuous<-factor(regress_output$is_continuous,levels=c(TRUE,FALSE),labels=c("Imputed","Dichotomized"))
  }
  
  return(regress_output)
}


# Function to compare results from different models
compare_func<-function(compare_list,data_set,p_named_only){
  output_table<-data.frame()
  for (i in seq_along(compare_list)){
    dataset<-as.data.frame(compare_list[i])
    # output[[i]]<-DescTools::Freq(dataset[,col])[1,]%>%
    #   mutate(model=paste0("Model_",i))
    if (p_named_only==T) {
      sig_data <- dataset %>%
        # filter(p_val_fdr < 0.05) %>%
        dplyr::group_by(sig_fdr_named) %>%
        dplyr::summarise(freq = n())%>%
        ungroup()%>%
        tidyr::complete(sig_fdr_named,fill=list(freq=0))%>% #fill in the rows where the combination is missing (no significant results)
        dplyr::filter(sig_fdr_named=="FDR<0.05")%>%
        dplyr::select(freq)
    } else {
      sig_data <- dataset %>%
      # filter(p_val_fdr < 0.05) %>%
      dplyr::group_by(sig_fdr) %>%
      dplyr::summarise(freq = n())%>%
      ungroup()%>%
      tidyr::complete(sig_fdr,fill=list(freq=0))%>% #fill in the rows where the combination is missing (no significant results)
      dplyr::filter(sig_fdr=="FDR<0.05")%>%
      dplyr::select(freq)
    }

    if (nrow(sig_data)==0) {
      sig_data<-data.frame(freq=c(0))
    } else {sig_data<-sig_data}
    total_data <- dataset %>%
      dplyr::summarise(total = n())
    summary_data <- cbind(sig_data, total_data) %>%
      mutate(
        perc = round(freq / total, digits = 3),
        model = paste0("Model_", i),
        strata = "Both"
      ) %>%
      dplyr::select(model, strata, freq, perc)
    output_table <- rbind(output_table, summary_data)
    
  }
  # names(output)<-name
  # output_table<-data.table::rbindlist(output,use.names = T,fill = T)%>%
  output_table<-output_table%>%
    mutate(perc=round(perc,digits = 3),
           strata="All")%>%
    dplyr::select(c(model,strata,freq,perc))%>%
    dplyr::rename_at(vars(-c(model,strata)), ~ paste0(.,'_',data_set))
  return(output_table)
}

# Function to compare results from different models (using bonferroni correction)
compare_bonf_sig_func<-function(compare_list,data_set,p_named_only){
  output_table<-data.frame()
  for (i in seq_along(compare_list)){
    dataset<-as.data.frame(compare_list[i])
    # output[[i]]<-DescTools::Freq(dataset[,col])[1,]%>%
    #   mutate(model=paste0("Model_",i))
    if (p_named_only==T) {
      sig_data <- dataset %>%
        # filter(p_val_fdr < 0.05) %>%
        dplyr::group_by(sig_Bonf) %>%
        dplyr::summarise(freq = n())%>%
        ungroup()%>%
        tidyr::complete(sig_Bonf,fill=list(freq=0))%>% #fill in the rows where the combination is missing (no significant results)
        dplyr::filter(sig_Bonf=="Bonf<0.05")%>%
        dplyr::select(freq)
    } else {
      sig_data <- dataset %>%
        # filter(p_val_fdr < 0.05) %>%
        dplyr::group_by(sig_Bonf) %>%
        dplyr::summarise(freq = n())%>%
        ungroup()%>%
        tidyr::complete(sig_Bonf,fill=list(freq=0))%>% #fill in the rows where the combination is missing (no significant results)
        dplyr::filter(sig_Bonf=="Bonf<0.05")%>%
        dplyr::select(freq)
    }
    
    if (nrow(sig_data)==0) {
      sig_data<-data.frame(freq=c(0))
    } else {sig_data<-sig_data}
    total_data <- dataset %>%
      dplyr::summarise(total = n())
    summary_data <- cbind(sig_data, total_data) %>%
      mutate(
        perc = round(freq / total, digits = 3),
        model = paste0("Model_", i),
        strata = "Both"
      ) %>%
      dplyr::select(model, strata, freq, perc)
    output_table <- rbind(output_table, summary_data)
    
  }
  # names(output)<-name
  # output_table<-data.table::rbindlist(output,use.names = T,fill = T)%>%
  output_table<-output_table%>%
    mutate(perc=round(perc,digits = 3),
           strata="All")%>%
    dplyr::select(c(model,strata,freq,perc))%>%
    dplyr::rename_at(vars(-c(model,strata)), ~ paste0(.,'_',data_set))
  return(output_table)
}

# Functions to compare results from different models (stratification analysis)
compare_strat_func<-function(compare_list,p_named_only,data_set){
  output_table<-data.frame()
  for (i in seq_along(compare_list)){
    dataset<-as.data.frame(compare_list[i])
    if (p_named_only==T){
      sig_data<-dataset%>%
        # filter(p_val_fdr_named<0.05)%>%
        group_by(strata,sig_fdr_named)%>%
        dplyr::summarise(freq=n())%>%
        ungroup()%>%
        tidyr::complete(strata,sig_fdr_named,fill=list(freq=0))%>% #fill in the rows where the combination is missing (no significant results)
        dplyr::filter(sig_fdr_named=="FDR<0.05")%>%
        dplyr::select(strata,freq)
    } else {
      sig_data<-dataset%>%
        # filter(p_val_fdr<0.05)%>%
        group_by(strata,sig_fdr)%>%
        dplyr::summarise(freq=n())%>%
        ungroup()%>%
        tidyr::complete(strata,sig_fdr,fill=list(freq=0))%>% #fill in the rows where the combination is missing (no significant results)
        dplyr::filter(sig_fdr=="FDR<0.05")%>%
        dplyr::select(strata,freq)
    }
    total_data<-dataset%>%
      group_by(strata)%>%
      dplyr::summarise(total=n())
    if (nrow(sig_data)==0) {
     sig_data<-data.frame(strata=total_data[,1],freq=c(0,0)) 
    } else {sig_data<-sig_data}
    summary_data<-merge(sig_data,total_data,by="strata")%>%
      dplyr::mutate(perc=round(freq/total,digits = 3),
             model=paste0("Model_",i))%>%
      dplyr::select(model,strata,freq,perc)
    output_table<-rbind(output_table,summary_data)

  }
  output_table<-output_table%>%
    dplyr::rename_at(vars(-c(model,strata)), ~ paste0(.,'_',data_set))
  return(output_table)
}

# Functions to compare results from different models (stratification analysis) (bonferroni correction)
compare_bonf_sig_strat_func<-function(compare_list,p_named_only,data_set){
  output_table<-data.frame()
  for (i in seq_along(compare_list)){
    dataset<-as.data.frame(compare_list[i])
    if (p_named_only==T){
      sig_data<-dataset%>%
        # filter(p_val_fdr_named<0.05)%>%
        group_by(strata,sig_Bonf)%>%
        dplyr::summarise(freq=n())%>%
        ungroup()%>%
        tidyr::complete(strata,sig_Bonf,fill=list(freq=0))%>% #fill in the rows where the combination is missing (no significant results)
        dplyr::filter(sig_Bonf=="Bonf<0.05")%>%
        dplyr::select(strata,freq)
    } else {
      sig_data<-dataset%>%
        # filter(p_val_fdr<0.05)%>%
        group_by(strata,sig_Bonf)%>%
        dplyr::summarise(freq=n())%>%
        ungroup()%>%
        tidyr::complete(strata,sig_Bonf,fill=list(freq=0))%>% #fill in the rows where the combination is missing (no significant results)
        dplyr::filter(sig_Bonf=="Bonf<0.05")%>%
        dplyr::select(strata,freq)
    }
    total_data<-dataset%>%
      group_by(strata)%>%
      dplyr::summarise(total=n())
    if (nrow(sig_data)==0) {
      sig_data<-data.frame(strata=total_data[,1],freq=c(0,0)) 
    } else {sig_data<-sig_data}
    summary_data<-merge(sig_data,total_data,by="strata")%>%
      dplyr::mutate(perc=round(freq/total,digits = 3),
                    model=paste0("Model_",i))%>%
      dplyr::select(model,strata,freq,perc)
    output_table<-rbind(output_table,summary_data)
    
  }
  output_table<-output_table%>%
    dplyr::rename_at(vars(-c(model,strata)), ~ paste0(.,'_',data_set))
  return(output_table)
}

# Function to generate pvalue vs beta graphs
pval_beta_plot<-function(data,model_text,footnote_text,is_imputed,p_named_only,trait_label_original,trait_label_binary,trait_for_model){
  if (trait_for_model=="binary") {
    title_text<-trait_label_binary
  } else { title_text<-trait_label_original}
  footnote_text<-paste(footnote_text,collapse = ",")
if (p_named_only==T){
  data$pval_fdr_neglog<--log10(data$p_val_fdr_named)
  pval_top15<-quantile(data$pval_fdr_neglog,probs = 0.985,na.rm = T)
  text_label_pval<-ifelse(pval_top15>-log10(0.05),pval_top15,-log10(0.05))
  data<-data[which(!is.na(data$p_val_fdr_named)),]
  
  if (is_imputed==T) {
  plot<-ggplot(data=data,aes(x=beta,y=pval_fdr_neglog))+
    geom_point(aes(color=sig_fdr_named,shape = is_continuous))+
    scale_color_manual(values=c("red","grey"),name="Significant")+
    scale_shape_manual(values = c(6, 20),name="Imputed or Dichotomized")+
    geom_hline(aes(yintercept=-log10(0.05)))+
    labs(x="Beta",y="-log10(p value)",title = paste(model_text,":\nScatter plot of beta coefficient and FDR adjusted P value\n(-log10) of metabolites concentrations\nand",title_text),caption=paste("metabolites with pvalue<",formatC(exp(-text_label_pval), format = "e", digits = 2),"annotated in the graph\nmodel adjusted for",tolower(footnote_text)))+
    geom_text_repel(data=subset(data, pval_fdr_neglog>text_label_pval),aes(label=metabolite),size=3)
} else {
  plot<-ggplot(data=data,aes(x=beta,y=pval_fdr_neglog))+
    geom_point(aes(color=sig_fdr_named))+
    scale_color_manual(values=c("red","grey"),name="Significant")+
    geom_hline(aes(yintercept=-log10(0.05)))+
    labs(x="Beta",y="-log10(p value)",title = paste(model_text,":\nScatter plot of beta coefficient and FDR adjusted P value\n(-log10) of metabolites concentrations\nand",title_text),caption=paste("metabolites with pvalue<",formatC(exp(-text_label_pval), format = "e", digits = 2),"annotated in the graph\nmodel adjusted for",tolower(footnote_text)))+
    geom_text_repel(data=subset(data, pval_fdr_neglog>text_label_pval),aes(label=metabolite),size=3)
} 
  } else {
    pval_top15<-quantile(data$pval_fdr_neglog,probs = 0.985,na.rm = T)
    text_label_pval<-ifelse(pval_top15>-log10(0.05),pval_top15,-log10(0.05))
  if (is_imputed==T) {
    plot<-ggplot(data=data,aes(x=beta,y=pval_fdr_neglog))+
      geom_point(aes(color=sig_fdr,shape = is_continuous))+
      scale_color_manual(values=c("red","grey"),name="Significant")+
      scale_shape_manual(values = c(6, 20),name="Imputed or Dichotomized")+
      geom_hline(aes(yintercept=-log10(0.05)))+
      labs(x="Beta",y="-log10(p value)",title = paste(model_text,":\nScatter plot of beta coefficient and FDR adjusted P value\n(-log10) of metabolites concentrations\nand",title_text),caption=paste("metabolites with pvalue<",formatC(exp(-text_label_pval), format = "e", digits = 2),"annotated in the graph\nmodel adjusted for",tolower(footnote_text)))+
      geom_text_repel(data=subset(data, pval_fdr_neglog>text_label_pval),aes(label=metabolite),size=3)
  } else {
    plot<-ggplot(data=data,aes(x=beta,y=pval_fdr_neglog))+
      geom_point(aes(color=sig_fdr))+
      scale_color_manual(values=c("red","grey"),name="Significant")+
      geom_hline(aes(yintercept=-log10(0.05)))+
      labs(x="Beta",y="-log10(p value)",title = paste(model_text,":\nScatter plot of beta coefficient and FDR adjusted P value\n(-log10) of metabolites concentrations\nand",title_text),caption=paste("metabolites with pvalue<",formatC(exp(-text_label_pval), format = "e", digits = 2),"annotated in the graph\nmodel adjusted for",tolower(footnote_text)))+
      geom_text_repel(data=subset(data, pval_fdr_neglog>text_label_pval),aes(label=metabolite),size=3)
  }
}
  return(plot)
}

# function to compare results
compare_annot<-function(imp_output,comp_output,annot,p_named_only) {
  if (p_named_only==T){
    imp_sig_results<-imp_output[which(imp_output$p_val_fdr_named<=0.05),]%>%
      mutate(p_val_fdr=p_val_fdr_named)
    comp_sig_results<-comp_output[which(comp_output$p_val_fdr_named<=0.05),]%>%
      mutate(p_val_fdr=p_val_fdr_named)
  } else {
    imp_sig_results<-imp_output[which(imp_output$p_val_fdr<=0.05),]
    comp_sig_results<-comp_output[which(comp_output$p_val_fdr<=0.05),]
  }
  imp_candidate<-annot[which(annot$METABOLITE_STD%in%imp_sig_results$metabolite),]
  comp_candidate<-annot[which(annot$METABOLITE_STD%in%comp_sig_results$metabolite),]
  
  if (nrow(imp_candidate)==0&nrow(comp_candidate)==0){
    print("No metabolite candidate identified")
  } else if (nrow(imp_candidate)==0&nrow(comp_candidate)!=0) {
    comp_tbl<-merge(comp_sig_results,comp_candidate,by.x="metabolite",by.y="METABOLITE_STD")
    comp_unique<-comp_tbl%>%
      dplyr::mutate(unique="Complete-case only")%>%
      merge(.,imp_output[,c("metabolite","beta","se","p_val_fdr","p_val")],by="metabolite",suffixes=c("_comp","_imp"))
    union_table<-comp_unique%>%
      mutate_at(vars(beta_imp, beta_comp,se_imp, se_comp), funs(round(., 3)))%>%
      mutate_at(vars(p_val_imp,p_val_fdr_imp,p_val_comp,p_val_fdr_comp), funs(formatC(., format = "e", digits = 2)))
    return(as.data.frame(union_table[,c("metabolite","BIOCHEMICAL","SUPER_PATHWAY","SUB_PATHWAY","PLATFORM","MASS","HMDB","unique","beta_imp","se_imp","p_val_imp","p_val_fdr_imp","beta_comp","se_comp","p_val_comp","p_val_fdr_comp")]))
  } else if (nrow(imp_candidate)!=0&nrow(comp_candidate)==0) {
    imp_tbl<-merge(imp_sig_results,imp_candidate,by.x="metabolite",by.y="METABOLITE_STD")
    imp_unique<-imp_tbl%>%
      dplyr::mutate(unique="Imputed only")%>%
      merge(.,comp_output[,c("metabolite","beta","se","p_val_fdr","p_val")],by="metabolite",suffixes=c("_imp","_comp"))
    union_table<-imp_unique%>%
      mutate_at(vars(beta_imp, beta_comp,se_imp, se_comp), funs(round(., 3)))%>%
      mutate_at(vars(p_val_imp,p_val_fdr_imp,p_val_comp,p_val_fdr_comp), funs(formatC(., format = "e", digits = 2)))
    return(as.data.frame(union_table[,c("metabolite","BIOCHEMICAL","SUPER_PATHWAY","SUB_PATHWAY","PLATFORM","MASS","HMDB","unique","beta_imp","se_imp","p_val_imp","p_val_fdr_imp","beta_comp","se_comp","p_val_comp","p_val_fdr_comp")]))
  } else {
    imp_tbl<-merge(imp_sig_results,imp_candidate,by.x="metabolite",by.y="METABOLITE_STD")
    comp_tbl<-merge(comp_sig_results,comp_candidate,by.x="metabolite",by.y="METABOLITE_STD")
    
    imp_unique<-imp_tbl%>%
      dplyr::filter(!metabolite%in%comp_sig_results$metabolite)%>%
      dplyr::mutate(unique="Imputed only")%>%
      merge(.,comp_output[,c("metabolite","beta","se","p_val_fdr","p_val")],by="metabolite",suffixes=c("_imp","_comp"))
    
    comp_unique<-comp_tbl%>%
      dplyr::filter(!metabolite%in%imp_sig_results$metabolite)%>%
      dplyr::mutate(unique="Complete-case only")%>%
      merge(.,imp_output[,c("metabolite","beta","se","p_val_fdr","p_val")],by="metabolite",suffixes=c("_comp","_imp"))
    
    imp_shared<-imp_tbl%>%
      dplyr::filter(metabolite%in%comp_sig_results$metabolite)%>%
      mutate(unique="Shared")
    
    comp_shared<-comp_tbl%>%
      dplyr::filter(metabolite%in%imp_tbl$metabolite)%>%
      dplyr::mutate(unique="Shared")
    
    both_shared<-merge(imp_shared[,c("metabolite","BIOCHEMICAL","SUPER_PATHWAY","SUB_PATHWAY","PLATFORM","MASS","PUBCHEM","KEGG","HMDB","beta","se","p_val_fdr","p_val")],comp_shared[,c("metabolite","beta","se","p_val_fdr","p_val")],by="metabolite",suffixes = c("_imp","_comp"))%>%
      dplyr::mutate(unique="Shared")
    
    union_table<-rbind.fill(imp_unique,comp_unique,both_shared)
    # write.csv(union_table,file = paste0("output/metabolite_candidates_",params$trait_input,".csv"))
    
    union_table<-union_table%>%
      mutate_at(vars(beta_imp, beta_comp,se_imp, se_comp), funs(round(., 3)))%>%
      mutate_at(vars(p_val_imp,p_val_fdr_imp,p_val_comp,p_val_fdr_comp), funs(formatC(., format = "e", digits = 2)))
    return(as.data.frame(union_table[,c("metabolite","BIOCHEMICAL","SUPER_PATHWAY","SUB_PATHWAY","PLATFORM","MASS","HMDB","unique","beta_imp","se_imp","p_val_imp","p_val_fdr_imp","beta_comp","se_comp","p_val_comp","p_val_fdr_comp")]))
}
}

# combine all results into one table
combine_output<-function(data_list){
  output_list<-list()
  # all_cand<-lapply(data_list,FUN=function(data){data<-data[which(data$p_val_Bonf<0.05),]})
  for (i in seq_along(data_list)){
    # output_list[i]<-as.data.frame(data_list[i])%>%dplyr::filter(.,p_val_Bonf<0.05)
    output_list[i]<-as.data.frame(data_list[[i]])[which(data_list[[i]]$p_val_Bonf<0.05),]
  }
  all_cand<-purrr::reduce(union,data_list)
  return(all_cand)
}  

# univariate glm model with sampling weights
by_svyglm_uni<-function(data,outcome,exposure){
  svy_design=survey::svydesign(id=~PSU_ID, strata=~STRAT, weights=~WEIGHT, data=data)
  outcome_is_binary<-ifelse(length(unique(data[!is.na(data[,outcome]),outcome]))==2,TRUE,FALSE)
  if (outcome_is_binary==T) {
  svymodel<-survey::svyglm(as.formula(paste(outcome,"~",exposure)),design=svy_design,family=quasipoisson(link='log'))  
  } else {
  svymodel<-survey::svyglm(as.formula(paste(outcome,"~",exposure)),design=svy_design)  
  }
  return(svymodel)
}
# multivariate glm model with sampling weights
by_svyglm_multi<-function(data,outcome,exposure,covar,include_ID){
  data<-data[!is.na(data$WEIGHT),]
  svy_design=survey::svydesign(id=~PSU_ID, strata=~STRAT, weights=~WEIGHT, data=data)
  outcome_is_binary<-ifelse(length(unique(data[!is.na(data[,outcome]),outcome]))==2,TRUE,FALSE)
  include_df<-data[which(data$ID%in%include_ID),]
  include_df<-include_df[complete.cases(include_df[,c(outcome,exposure,covar)]),]
  
  if (outcome_is_binary==T) {
    svymodel<-survey::svyglm(as.formula(paste(outcome,"~",exposure,"+",paste(covar,collapse = "+"))),design=subset(svy_design,ID%in%include_df$ID),family=quasipoisson(link='log'))
  } else {
    svymodel<-survey::svyglm(as.formula(paste(outcome,"~",exposure,"+",paste(covar,collapse = "+"))),design=subset(svy_design,ID%in%include_df$ID))
  }
  return(svymodel)
}


# Stratified plot
pval_beta_strat_plot<-function(data,model_text,footnote_text,is_imputed,p_named_only,trait_label_original,trait_label_binary,trait_for_model,stratifier){
  if (trait_for_model=="binary") {
    title_text<-trait_label_binary
  } else { title_text<-trait_label_original
  } 
  footnote_text<-footnote_text[!footnote_text%in%stratifier]
  footnote_text<-paste(footnote_text,collapse = ",")
  if (p_named_only==T) {
      data$pval_fdr_neglog<--log10(data$p_val_fdr_named)
      pval_top15<-quantile(data$pval_fdr_neglog,probs = 0.985,na.rm = T)
      data$sig_fdr<-data$sig_fdr_named
    } else {pval_top15<-quantile(data$pval_fdr_neglog,probs = 0.985,na.rm = T)}
    text_label_pval<-ifelse(pval_top15>-log10(0.05),pval_top15,-log10(0.05))
    if (is_imputed==T) {
      plot<-ggplot(data=data,aes(x=beta,y=pval_fdr_neglog))+
        geom_point(aes(color=sig_fdr,shape = is_continuous))+
        scale_color_manual(values=c("red","grey"),name="Significant")+
        scale_shape_manual(values = c(6, 20),name="Imputed or Dichotomized")+
        geom_hline(aes(yintercept=-log10(0.05)))+
        labs(x="Beta",y="-log10(p value)",title = paste(model_text,":\nScatter plot of beta coefficient and FDR adjusted P value\n(-log10) of metabolites concentrations\nand",title_text),caption=paste("metabolites with pvalue<",formatC(exp(-text_label_pval), format = "e", digits = 2),"annotated in the graph\nmodel adjusted for",tolower(footnote_text)))+
        geom_text_repel(data=subset(data, pval_fdr_neglog>text_label_pval),aes(label=metabolite),size=3)
    plot2<-plot+facet_wrap(~strata)
    } else {
      plot<-ggplot(data=data,aes(x=beta,y=pval_fdr_neglog))+
        geom_point(aes(color=sig_fdr))+
        scale_color_manual(values=c("red","grey"),name="Significant")+
        geom_hline(aes(yintercept=-log10(0.05)))+
        labs(x="Beta",y="-log10(p value)",title = paste(model_text,":\nScatter plot of beta coefficient and FDR adjusted P value\n(-log10) of metabolites concentrations\nand",title_text),caption=paste("metabolites with pvalue<",formatC(exp(-text_label_pval), format = "e", digits = 2),"annotated in the graph\nmodel adjusted for",tolower(footnote_text)))+
        geom_text_repel(data=subset(data, pval_fdr_neglog>text_label_pval),aes(label=metabolite),size=3)
    plot2<-plot+facet_wrap(~strata)
    }
  return(plot2)
}


# Function to loop regression model and export model summary using survey weight
strat_svyreg_loop<-function(data,covar,end,metab_is_cont,metab_is_complete,trait_original,trait_binary,trait_for_model,trait_as_predictor,stratifier){
  if (trait_for_model=="binary") {
    trait<-"TraitBinary"
  }  else if (trait_for_model=="original"){
    trait<-trait_original
  } else { trait<-"Trait" }
  
  dat<-data[,1:end]
  dat<-dat[,!colnames(dat)%in%c("SOL_ID","ID","LAB_ID")]
  
  covar<-covar[!covar%in%stratifier]
  group_by_value<-unique(data[,stratifier])
  out_nvar<-ncol(dat)
  out_beta<-rep(NA,out_nvar)
  out_se<-rep(NA,out_nvar)
  out_pvalue<-rep(NA,out_nvar)
  out_nobs<-rep(NA,out_nvar)
  group<-rep(NA,out_nvar)
  regress_output<-data.frame()
  # if the original trait variable is binary or the user pick the binary version to use in the model
  trait_is_cont<-ifelse(length(unique(data[!is.na(data[,trait]),trait]))==2,FALSE,
                       ifelse(trait_for_model=="binary",FALSE,TRUE))
  
  
for (j in seq_along(group_by_value)){
  newdata<-data[which(data[,stratifier]==as.character(group_by_value[j])),]
  # survey_design=svydesign(id=~PSU_ID, strata=~STRAT, weights=~WEIGHT, data=data)
  if (!"ID"%in%colnames(data)){
    data<-data%>%
      dplyr::rename("ID"="SOL_ID")
    newdata<-newdata%>%
      dplyr::rename("ID"="SOL_ID")
  }
  dat<-newdata[,1:end]
  dat<-dat[,!colnames(dat)%in%c("ID","LAB_ID")]

  if (all(grep("[A-Z]", names(data[,-1])) == 0)) { # if column names are all in lower case
    survey_design=svydesign(id=~psu_id, strata=~strat, weights=~weight, data=data)
  } else {
    survey_design=svydesign(id=~PSU_ID, strata=~STRAT, weights=~WEIGHT, data=data)
  }
  
  for(i in 1:(ncol(dat))) {
    # print(i)
    met_df<-cbind(newdata[,colnames(newdata)%in%colnames(dat)[i]],newdata[,c("ID",covar)])
    met_df_id<-met_df[complete.cases(met_df),"ID"]
    if (trait_as_predictor==T){
      if (metab_is_cont==T){
        model<- svyglm(as.formula(paste0(colnames(dat)[i],"~",trait,"+",paste(covar,collapse= "+"))),design=subset(survey_design,ID%in%met_df_id))
      } else {
        model<- svyglm(as.formula(paste0(colnames(dat)[i],"~",trait,"+",paste(covar,collapse = "+"))),design=subset(survey_design,ID%in%met_df_id),family=quasipoisson(link='log'))
      }
    }else {
      if (trait_is_cont==T){
        model<- svyglm(as.formula(paste0(trait,"~",colnames(dat)[i],"+",paste(covar,collapse= "+"))),design=subset(survey_design,ID%in%met_df_id))
      }else{
        model<- svyglm(as.formula(paste0(trait,"~",colnames(dat)[i],"+",paste(covar,collapse = "+"))),design=subset(survey_design,ID%in%met_df_id),family=quasipoisson(link='log'))
      }
    }
    
    Vcov <- vcov(model, useScale = FALSE)
    beta<- coef(model)
    se<- sqrt(diag(vcov(model, useScale = FALSE)))
    zval<- beta / se
    pval<- 2 * pnorm(abs(zval), lower.tail = FALSE)
    out_beta[i]=as.numeric(beta[2])
    out_se[i] = as.numeric(se[2])
    out_pvalue[i] = as.numeric(pval[2])
    out_nobs[i]=as.numeric(nobs(model))
    group[i]=as.character(group_by_value[j])
    # print(i)
  }
  strata_output<-data.frame(metabolite=colnames(dat),
                             beta=out_beta,
                             se=out_se,
                             n=out_nobs,
                             p_val=out_pvalue,
                             p_val_Bonf=p.adjust(out_pvalue,method="bonferroni"),
                             p_val_fdr=p.adjust(out_pvalue,method="BH"),
                             strata=group)%>%
    dplyr::mutate(
      pval_Bonf_neglog=-log10(p_val_Bonf),
      pval_fdr_neglog=-log10(p_val_fdr),
      sig_Bonf=ifelse(p_val_Bonf<=0.05,"Bonf<0.05","Not Sig"),
      sig_fdr=ifelse(p_val_fdr<=0.05,"FDR<0.05","Not Sig"),
      is_continuous=metab_is_cont
    )
  regress_output<-rbind(regress_output,strata_output)
}
  identified_data<-regress_output%>%
    dplyr::filter(!stringr::str_detect(metabolite,"^X"))%>%
    dplyr::group_by(strata)%>%
    dplyr::mutate(p_val_fdr_named=p.adjust(p_val,method="BH"),
           sig_fdr_named=ifelse(p_val_fdr_named<=0.05,"FDR<0.05","Not Sig"))
  unidentified_data<-regress_output%>%
    dplyr::filter(stringr::str_detect(metabolite,"^X"))%>%
    dplyr::mutate(p_val_fdr_named=NA,
           sig_fdr_named=NA)
  regress_output<-plyr::rbind.fill(identified_data,unidentified_data)
  if (metab_is_complete==T){
    regress_output$is_continuous<-"Complete-cases"
  } else {
    regress_output$is_continuous<-factor(regress_output$is_continuous,levels=c(TRUE,FALSE),labels=c("Imputed","Dichotomized"))
  }
  
  return(regress_output)
}

# Simple non-templated version
# Function to loop regression model and export model summary using survey weight
strat_svyreg_loop_nontemplate<-function(data,covar,start,end,metab_is_cont,metab_is_complete,trait,trait_as_predictor,stratifier){
  covar<-covar[!covar%in%stratifier]
  group_by_value<-unique(data[,stratifier])
  out_nvar<-end-start+1
  out_beta<-rep(NA,out_nvar)
  out_se<-rep(NA,out_nvar)
  out_pvalue<-rep(NA,out_nvar)
  group<-rep(NA,out_nvar)
  regress_output<-data.frame()
  # if the original trait variable is binary or the user pick the binary version to use in the model
  trait_is_cont<-ifelse(length(unique(data[!is.na(data[,trait]),trait]))==2,FALSE,TRUE)
  
  for (j in seq_along(group_by_value)){
    newdata<-data[which(data[,stratifier]==as.character(group_by_value[j])),]
    dat<-newdata[,start:end]
    survey_design=svydesign(id=~PSU_ID, strata=~STRAT, weights=~WEIGHT, data=newdata)
    for(i in 1:(ncol(dat))) {
      if (trait_as_predictor==T){
        if (metab_is_cont==T){
          model<- svyglm(as.formula(paste0(colnames(dat)[i],"~",trait,"+",paste(covar,collapse= "+"))),design=survey_design)
        } else {
          model<- svyglm(as.formula(paste0(colnames(dat)[i],"~",trait,"+",paste(covar,collapse = "+"))),design=survey_design,family=quasipoisson(link='log'))
        }
      }else {
        if (trait_is_cont==T){
          model<- svyglm(as.formula(paste0(trait,"~",colnames(dat)[i],"+",paste(covar,collapse= "+"))),design=survey_design)
        }else{
          model<- svyglm(as.formula(paste0(trait,"~",colnames(dat)[i],"+",paste(covar,collapse = "+"))),design=survey_design,family=quasipoisson(link='log'))
        }
      }
      
      Vcov <- vcov(model, useScale = FALSE)
      beta<- coef(model)
      se<- sqrt(diag(vcov(model, useScale = FALSE)))
      zval<- beta / se
      pval<- 2 * pnorm(abs(zval), lower.tail = FALSE)
      out_beta[i]=as.numeric(beta[2])
      out_se[i] = as.numeric(se[2])
      out_pvalue[i] = as.numeric(pval[2])
      group[i]=as.character(group_by_value[j])
    }
    strata_output<-data.frame(metabolite=colnames(dat),
                              beta=out_beta,
                              se=out_se,
                              p_val=out_pvalue,
                              p_val_Bonf=p.adjust(out_pvalue,method="bonferroni"),
                              p_val_fdr=p.adjust(out_pvalue,method="BH"),
                              strata=group)%>%
      mutate(
        pval_Bonf_neglog=-log10(p_val_Bonf),
        pval_fdr_neglog=-log10(p_val_fdr),
        sig_Bonf=ifelse(p_val_Bonf<=0.05,"Bonf<0.05","Not Sig"),
        sig_fdr=ifelse(p_val_fdr<=0.05,"FDR<0.05","Not Sig"),
        is_continuous=metab_is_cont
      )
    regress_output<-rbind(regress_output,strata_output)
  }
  identified_data<-regress_output%>%
    dplyr::filter(!stringr::str_detect(metabolite,"^X"))%>%
    dplyr::group_by(strata)%>%
    dplyr::mutate(p_val_fdr_named=p.adjust(p_val,method="BH"),
                  sig_fdr_named=ifelse(p_val_fdr_named<=0.05,"FDR<0.05","Not Sig"))
  unidentified_data<-regress_output%>%
    dplyr::filter(stringr::str_detect(metabolite,"^X"))%>%
    dplyr::mutate(p_val_fdr_named=NA,
                  sig_fdr_named=NA)
  regress_output<-plyr::rbind.fill(identified_data,unidentified_data)
  if (metab_is_complete==T){
    regress_output$is_continuous<-"Complete-cases"
  } else {
    regress_output$is_continuous<-factor(regress_output$is_continuous,levels=c(TRUE,FALSE),labels=c("Imputed","Dichotomized"))
  }
  
  return(regress_output)
}
#################################
svyreg_loop_uni<-function(data,predictor,end){
  dat<-data[,1:end]
  dat<-dat[,!colnames(dat)%in%c("LAB_ID","ID","SOL_ID","id")]
    # dat<-data[,3:end]
  # out_nvar=end-3+1
  out_nvar=ncol(dat)
  out_beta=rep(NA,out_nvar)
  out_se=rep(NA,out_nvar)
  out_pvalue=rep(NA,out_nvar)
  
  survey_design=svydesign(id=~PSU_ID, strata=~STRAT, weights=~WEIGHT, data=data)
  
  for(i in 1:(ncol(dat))) {
    model<- svyglm(as.formula(paste0(colnames(dat)[i],"~",predictor)),design=survey_design)
    Vcov <- vcov(model, useScale = FALSE)
    beta<- coef(model)
    se<- sqrt(diag(vcov(model, useScale = FALSE)))
    zval<- beta / se
    pval<- 2 * pnorm(abs(zval), lower.tail = FALSE)
    out_beta[i]=as.numeric(beta[2])
    out_se[i] = round(as.numeric(se[2]),digits = 3)
    out_pvalue[i] = as.numeric(pval[2])
  }
  regress_output<-data.frame(metabolite=colnames(dat),
                             beta=out_beta,
                             se=out_se,
                             p_val=out_pvalue,
                             p_val_Bonf=p.adjust(out_pvalue,method="bonferroni"),
                             p_val_fdr=p.adjust(out_pvalue,method="BH")
  )%>%
    dplyr::mutate(
      pval_Bonf_neglog=-log10(p_val_Bonf),
      pval_fdr_neglog=-log10(p_val_fdr),
      sig_Bonf=ifelse(p_val_Bonf<=0.05,"Bonf<0.05","Not Sig"),
      sig_fdr=ifelse(p_val_fdr<=0.05,"FDR<0.05","Not Sig")
    )
  identified_data<-regress_output%>%
    dplyr::filter(!stringr::str_detect(metabolite,"^X"))%>%
    dplyr::mutate(p_val_fdr_named=p.adjust(p_val,method="BH"),
           sig_fdr_named=ifelse(p_val_fdr_named<=0.05,"FDR<0.05","Not Sig"))
  unidentified_data<-regress_output%>%
    dplyr::filter(stringr::str_detect(metabolite,"^X"))%>%
    dplyr::mutate(p_val_fdr_named=NA,
           sig_fdr_named=NA)
  regress_output<-rbind(identified_data,unidentified_data)
  return(regress_output)
}

##########################################################################
# functions for atlas

# Function to loop regression model and export model summary using survey weight
atlas_svyreg_loop<-function(data,covar,end,metab_is_cont,metab_is_complete,trait,trait_as_predictor,trait_for_model){
  dat<-data[,1:end]
  dat<-dat[,!colnames(dat)%in%c("LAB_ID","ID","SOL_ID","id")]
  # dat<-data[,3:end]
  # out_nvar=end-3+1
  out_nvar=ncol(dat)
  out_beta=rep(NA,out_nvar)
  out_se=rep(NA,out_nvar)
  out_pvalue=rep(NA,out_nvar)
  out_n=rep(NA,out_nvar)
  # if the original trait variable is binary or the user pick the binary version to use in the model
  trait_is_cont=ifelse(length(unique(data[!is.na(data[,trait]),trait]))==2,FALSE,TRUE)
  
  survey_design=svydesign(id=~PSU_ID, strata=~STRAT, weights=~WEIGHT, data=data)
  
  for(i in 1:(ncol(dat))) {
    met_df<-cbind(data[,i+2],data[,c("ID",covar)])
    met_df_id<-met_df[complete.cases(met_df),"ID"]
    # print(i)
    if (trait_as_predictor==T){
      if (metab_is_cont==T){
        model<- svyglm(as.formula(paste0(colnames(dat)[i],"~",trait,"+",paste(covar,collapse= "+"))),design=subset(survey_design,ID%in%met_df_id))
      } else {
        model<- svyglm(as.formula(paste0(colnames(dat)[i],"~",trait,"+",paste(covar,collapse = "+"))),design=subset(survey_design,ID%in%met_df_id),family=quasipoisson(link='log'))
      }
    } else {
      if (trait_is_cont==T){
        model<- svyglm(as.formula(paste0(trait,"~",colnames(dat)[i],"+",paste(covar,collapse= "+"))),design=subset(survey_design,ID%in%met_df_id))
      } else {
        model<- svyglm(as.formula(paste0("as.numeric(",trait,")~",colnames(dat)[i],"+",paste(covar,collapse = "+"))),design=subset(survey_design,ID%in%met_df_id),family=quasipoisson(link='log'))
      }

    }
    
    Vcov <- vcov(model, useScale = FALSE)
    beta<- coef(model)
    se<- sqrt(diag(vcov(model, useScale = FALSE)))
    zval<- beta / se
    pval<- 2 * pnorm(abs(zval), lower.tail = FALSE)
    out_beta[i]=as.numeric(beta[2])
    out_se[i] = as.numeric(se[2])
    out_pvalue[i] = as.numeric(pval[2])
    out_n[i]=as.numeric(nobs(model))
    # print(i)
  }
  
  regress_output<-data.frame(metabolite=colnames(dat),
                             beta=out_beta,
                             se=out_se,
                             p_val=out_pvalue,
                             p_val_Bonf=p.adjust(out_pvalue,method="bonferroni"),
                             p_val_fdr=p.adjust(out_pvalue,method="BH"),
                             n=out_n
  )%>%
    mutate(
      pval_Bonf_neglog=-log10(p_val_Bonf),
      pval_fdr_neglog=-log10(p_val_fdr),
      sig_Bonf=ifelse(p_val_Bonf<=0.05,"Bonf<0.05","Not Sig"),
      sig_fdr=ifelse(p_val_fdr<=0.05,"FDR<0.05","Not Sig"),
      is_continuous=metab_is_cont
    )
  # identified_data<-regress_output%>%
  #   filter(!str_detect(metabolite,"^X"))%>%
  #   mutate(p_val_fdr_named=p.adjust(p_val,method="BH"),
  #          sig_fdr_named=ifelse(p_val_fdr_named<=0.05,"FDR<0.05","Not Sig"))
  # unidentified_data<-regress_output%>%
  #   filter(str_detect(metabolite,"^X"))%>%
  #   mutate(p_val_fdr_named=NA,
  #          sig_fdr_named=NA)
  # regress_output<-rbind(identified_data,unidentified_data)
  if (metab_is_complete==T){
    regress_output$is_continuous<-"Complete-cases"
  } else {
    regress_output$is_continuous<-factor(regress_output$is_continuous,levels=c(TRUE,FALSE),labels=c("Imputed","Dichotomized"))
  }
  
  return(regress_output)
}

# Function to loop regression sex-interaction model and export model summary using survey weight
atlas_svyreg_loop_interaction<-function(data,covar,end,metab_is_cont,metab_is_complete,trait,trait_as_predictor,interaction_term){
  dat<-data[,1:end]
  dat<-dat[,!colnames(dat)%in%c("LAB_ID","ID","SOL_ID","id")]
  # dat<-data[,3:end]
  # out_nvar=end-3+1
  out_nvar=ncol(dat)
  out_beta=rep(NA,out_nvar)
  int_beta=rep(NA,out_nvar)
  sex_beta=rep(NA,out_nvar)
  out_se=rep(NA,out_nvar)
  int_se=rep(NA,out_nvar)
  sex_se=rep(NA,out_nvar)
  out_pvalue=rep(NA,out_nvar)
  int_pvalue=rep(NA,out_nvar)
  sex_pvalue=rep(NA,out_nvar)
  out_n=rep(NA,out_nvar)
  # if the original trait variable is binary or the user pick the binary version to use in the model
  trait_is_cont=ifelse(length(unique(data[!is.na(data[,trait]),trait]))==2,FALSE,TRUE)
 
   if (!"ID"%in%colnames(data)){
    data<-data%>%
      dplyr::rename("ID"="SOL_ID")

  }
  
  if (all(grep("[A-Z]", names(data[,-1])) == 0)) { # if column names are all in lower case
    survey_design=svydesign(id=~psu_id, strata=~strat, weights=~weight, data=data)
  } else {
    survey_design=svydesign(id=~PSU_ID, strata=~STRAT, weights=~WEIGHT, data=data)
  }
  # survey_design=svydesign(id=~PSU_ID, strata=~STRAT, weights=~WEIGHT, data=data)
  
  for(i in 1:(ncol(dat))) {
    met_df<-cbind(data[,i+2],data[,c("ID",covar)])
    met_df_id<-met_df[complete.cases(met_df),"ID"]
    # print(i)
    if (trait_as_predictor==T){
      if (metab_is_cont==T){
        model<- svyglm(as.formula(paste0(colnames(dat)[i],"~",trait,"+",paste(c(paste0(trait,"*",interaction_term),covar),collapse= "+"))),design=subset(survey_design,ID%in%met_df_id))
      } else {
        model<- svyglm(as.formula(paste0(colnames(dat)[i],"~",trait,"+",paste(c(paste0(trait,"*",interaction_term),covar),collapse= "+"))),design=subset(survey_design,ID%in%met_df_id),family=quasipoisson(link='log'))
      }
    } else {
      if (trait_is_cont==T){
        model<- svyglm(as.formula(paste0(trait,"~",colnames(dat)[i],"+",paste(c(paste0(trait,"*",interaction_term),covar),collapse= "+"))),design=subset(survey_design,ID%in%met_df_id))
      } else {
        model<- svyglm(as.formula(paste0("as.numeric(",trait,")~",colnames(dat)[i],"+",paste(c(paste0(trait,"*",interaction_term),covar),collapse= "+"))),design=subset(survey_design,ID%in%met_df_id),family=quasipoisson(link='log'))
      }
      
    }
    
    Vcov <- vcov(model, useScale = FALSE)
    beta<- coef(model)
    se<- sqrt(diag(vcov(model, useScale = FALSE)))
    zval<- beta / se
    pval<- 2 * pnorm(abs(zval), lower.tail = FALSE)
    out_beta[i]=as.numeric(beta[2])
    int_beta[i]=as.numeric(beta[length(beta)])
    sex_beta[i]=as.numeric(beta[3])
    out_se[i] = as.numeric(se[2])
    int_se[i] = as.numeric(se[length(se)])
    sex_se[i] = as.numeric(se[3])
    out_pvalue[i] = as.numeric(pval[2])
    int_pvalue[i] = as.numeric(pval[length(pval)])
    sex_pvalue[i]= as.numeric(pval[3])
    out_n[i]=as.numeric(nobs(model))
    # print(i)
  }
  
  regress_output<-data.frame(metabolite=colnames(dat),
                             beta=out_beta,
                             int_beta=int_beta,
                             sex_beta=sex_beta,
                             se=out_se,
                             int_se=int_se,
                             sex_se=sex_se,
                             p_val=out_pvalue,
                             int_p_val=int_pvalue,
                             sex_p_val=sex_pvalue,
                             p_val_Bonf=p.adjust(out_pvalue,method="bonferroni"),
                             p_val_fdr=p.adjust(out_pvalue,method="BH"),
                             int_p_val_fdr=p.adjust(int_pvalue,method="BH"),
                             n=out_n
  )%>%
    mutate(
      pval_Bonf_neglog=-log10(p_val_Bonf),
      pval_fdr_neglog=-log10(p_val_fdr),
      sig_Bonf=ifelse(p_val_Bonf<=0.05,"Bonf<0.05","Not Sig"),
      sig_fdr=ifelse(p_val_fdr<=0.05,"FDR<0.05","Not Sig"),
      is_continuous=metab_is_cont
    )
  # identified_data<-regress_output%>%
  #   filter(!str_detect(metabolite,"^X"))%>%
  #   mutate(p_val_fdr_named=p.adjust(p_val,method="BH"),
  #          sig_fdr_named=ifelse(p_val_fdr_named<=0.05,"FDR<0.05","Not Sig"))
  # unidentified_data<-regress_output%>%
  #   filter(str_detect(metabolite,"^X"))%>%
  #   mutate(p_val_fdr_named=NA,
  #          sig_fdr_named=NA)
  # regress_output<-rbind(identified_data,unidentified_data)
  if (metab_is_complete==T){
    regress_output$is_continuous<-"Complete-cases"
  } else {
    regress_output$is_continuous<-factor(regress_output$is_continuous,levels=c(TRUE,FALSE),labels=c("Imputed","Dichotomized"))
  }
  
  return(regress_output)
}

# Function to loop regression model and export model summary using survey weight
atlas_strat_svyreg_loop<-function(data,covar,end,metab_is_cont,metab_is_complete,trait,trait_for_model,trait_as_predictor,stratifier){
  if (stratifier=="Gender") { # if stratifier is meno_stat then also drop gender from the covaraite list
    covar<-covar[!covar%in%stratifier]
  } else {
    covar<-covar[!covar%in%c(stratifier,"GENDER")]
  }
  
  group_by_value<-unique(data[,stratifier])
  group_by_value<-group_by_value[!is.na(group_by_value)]
  out_nvar<-length(colnames(data[,1:end])[!colnames(data[,1:end])%in%c("LAB_ID","ID","SOL_ID","id")])
  # out_nvar<-end-3+1
  out_beta<-rep(NA,out_nvar)
  out_se<-rep(NA,out_nvar)
  out_pvalue<-rep(NA,out_nvar)
  out_n<-rep(NA,out_nvar)
  group<-rep(NA,out_nvar)
  regress_output<-data.frame()
  # if the original trait variable is binary or the user pick the binary version to use in the model
  trait_is_cont<-ifelse(length(unique(data[!is.na(data[,trait]),trait]))==2,FALSE,TRUE)
  
  
  for (j in seq_along(group_by_value)){
    newdata<-data[which(data[,stratifier]==as.character(group_by_value[j])),]
    # dat<-data[,3:end]
    # out_nvar=end-3+1
    dat<-newdata[,1:end]
    dat<-dat[,!colnames(dat)%in%c("LAB_ID","ID","SOL_ID","id")]
    out_nvar=ncol(dat)
    survey_design=svydesign(id=~PSU_ID, strata=~STRAT, weights=~WEIGHT, data=data)
    for(i in 1:(ncol(dat))) {
      met_df<-cbind(newdata[,i+2],newdata[,c("ID",covar)])
      met_df_id<-met_df[complete.cases(met_df),"ID"]
      if (trait_as_predictor==T){
        if (metab_is_cont==T){
          model<- svyglm(as.formula(paste0(colnames(dat)[i],"~",trait,"+",paste(covar,collapse= "+"))),design=subset(survey_design,ID%in%met_df_id))
        } else {
          model<- svyglm(as.formula(paste0(colnames(dat)[i],"~",trait,"+",paste(covar,collapse = "+"))),design=subset(survey_design,ID%in%met_df_id),family=quasipoisson(link='log'))
        }
      }else {
        if (trait_is_cont==T){
          model<- svyglm(as.formula(paste0(trait,"~",colnames(dat)[i],"+",paste(covar,collapse= "+"))),design=subset(survey_design,ID%in%met_df_id))
        }else{
          model<- svyglm(as.formula(paste0("as.numeric(",trait,")~",colnames(dat)[i],"+",paste(covar,collapse = "+"))),design=subset(survey_design,ID%in%met_df_id),family=quasipoisson(link='log'))
        }
      }
      
      Vcov <- vcov(model, useScale = FALSE)
      beta<- coef(model)
      se<- sqrt(diag(vcov(model, useScale = FALSE)))
      zval<- beta / se
      pval<- 2 * pnorm(abs(zval), lower.tail = FALSE)
      out_beta[i]=as.numeric(beta[2])
      out_se[i] = as.numeric(se[2])
      out_pvalue[i] = as.numeric(pval[2])
      out_n[i]=as.numeric(nobs(model))
      group[i]=as.character(group_by_value[j])
    }
    strata_output<-data.frame(metabolite=colnames(dat),
                              beta=out_beta,
                              se=out_se,
                              p_val=out_pvalue,
                              p_val_Bonf=p.adjust(out_pvalue,method="bonferroni"),
                              p_val_fdr=p.adjust(out_pvalue,method="BH"),
                              strata=group,
                              n=out_n)%>%
      mutate(
        pval_Bonf_neglog=-log10(p_val_Bonf),
        pval_fdr_neglog=-log10(p_val_fdr),
        sig_Bonf=ifelse(p_val_Bonf<=0.05,"Bonf<0.05","Not Sig"),
        sig_fdr=ifelse(p_val_fdr<=0.05,"FDR<0.05","Not Sig"),
        is_continuous=metab_is_cont
      )
    regress_output<-rbind(regress_output,strata_output)
  }
  if (metab_is_complete==T){
    regress_output$is_continuous<-"Complete-cases"
  } else {
    regress_output$is_continuous<-factor(regress_output$is_continuous,levels=c(TRUE,FALSE),labels=c("Imputed","Dichotomized"))
  }
  
  return(regress_output)
}

#######################################################################################################
# for circular variables
# # Function to loop regression model and export model summary using survey weight and circular variables
# atlas_cir_svyreg_loop<-function(data,covar,end,metab_is_cont,metab_is_complete,trait,trait_as_predictor){
#   dat<-data[,1:end]
#   dat<-dat[,!colnames(dat)%in%c("LAB_ID","ID","SOL_ID","id")]
#   # dat<-data[,3:end]
#   # out_nvar=end-3+1
#   out_nvar=ncol(dat)
#   out_beta<-rep(NA,out_nvar)
#   out_se<-rep(NA,out_nvar)
#   out_pvalue<-rep(NA,out_nvar)
#   out_nobs<-rep(NA,out_nvar)
#   out_transformation<-rep(NA,out_nvar)
#   # if the original trait variable is binary or the user pick the binary version to use in the model
#   # trait_is_cont=ifelse(length(unique(data[!is.na(data[,trait]),trait]))==2,FALSE,TRUE)
#   if (all(grep("[A-Z]", names(data[,-1])) == 0)) { # if column names are all in lower case
#     survey_design=svydesign(id=~psu_id, strata=~strat, weights=~weight, data=data)
#   } else {
#     survey_design=svydesign(id=~PSU_ID, strata=~STRAT, weights=~WEIGHT, data=data)
#   }
#   # survey_design=svydesign(id=~PSU_ID, strata=~STRAT, weights=~WEIGHT, data=data)
#   
#   
#   
#   for(i in 1:(ncol(dat))) {
#     if (trait_as_predictor==T){
#       if (metab_is_cont==T){
#         model<- svyglm(as.formula(paste0(colnames(dat)[i],"~sin(",trait,")+cos(",trait,")+",paste(unlist(covar),collapse= "+"))),design=survey_design)
#       } else {
#         model<- svyglm(as.formula(paste0(colnames(dat)[i],"~sin(",trait,")+cos(",trait,")+",paste(unlist(covar),collapse = "+"))),design=survey_design,family=quasipoisson(link='log'))
#       }
#       # extract the coefficients and their covariance matrix
#       Vcov_matrix <- vcov(model, useScale = FALSE)
#       coef_matrix<- coef(model)
#       # if the correlation betwee sin() and cos() is negative, use atan() transformation instead of sin() and cos()
#       if (Vcov_matrix[2,3]<0){
#         if (metab_is_cont==T){
#           model<- svyglm(as.formula(paste0(colnames(dat)[i],"~atan(",trait,")+",paste(unlist(covar),collapse= "+"))),design=survey_design)
#         } else {
#           model<- svyglm(as.formula(paste0(colnames(dat)[i],"~atan(",trait,")+",paste(unlist(covar),collapse = "+"))),design=survey_design,family=quasipoisson(link='log'))
#         }
#         Vcov <- vcov(model, useScale = FALSE)
#         beta<- coef(model)
#         se<- sqrt(diag(vcov(model, useScale = FALSE)))
#         zval<- beta / se
#         pval<- 2 * pnorm(abs(zval), lower.tail = FALSE)
#         out_beta[i]=as.numeric(beta[2])
#         out_se[i] = round(as.numeric(se[2]),digits = 3)
#         out_pvalue[i] = as.numeric(pval[2])
#         out_nobs[i]=as.numeric(nobs(model))
#         out_transformation[i]="atan"
#       } else {
#         # Compute the estimated coefficient and its standard error for the circular predictor (coef_matrix[2] is sin(x), coef_matrix[3] is cos(x))
#         coef_circ <- sqrt(coef_matrix[2]^2 + coef_matrix[3]^2)
#         se_circ <- sqrt((coef_matrix[2]/coef_circ)^2*Vcov_matrix[2, 2]+(coef_matrix[3]/coef_circ)^2*Vcov_matrix[3,3]+2*coef_matrix[2]*coef_matrix[3]/(coef_circ^3)*Vcov_matrix[2, 3])
#         
#         # Compute the t-statistic and p-value for the circular predictor
#         t_stat_circ <- coef_circ / se_circ
#         p_value_circ <- 2 * pt(abs(t_stat_circ), df = model$df.residual, lower.tail = FALSE)
#         
#         out_beta[i]=as.numeric(coef_circ)
#         out_se[i] = as.numeric(se_circ)
#         out_pvalue[i] = as.numeric(p_value_circ)
#         out_nobs[i]=as.numeric(nobs(model)) 
#         out_transformation[i]="sin_cos"
#       }
# 
#     } else {
#       if (trait_is_cont==T){
#         model<- svyglm(as.formula(paste0(trait,"~",colnames(dat)[i],"+",paste(covar,collapse= "+"))),design=subset(survey_design,ID%in%met_df_id))
#       } else {
#         model<- svyglm(as.formula(paste0(trait,"~",colnames(dat)[i],"+",paste(covar,collapse = "+"))),design=subset(survey_design,ID%in%met_df_id),family=quasipoisson(link='log'))
#       }
#        Vcov <- vcov(model, useScale = FALSE)
#     beta<- coef(model)
#     se<- sqrt(diag(vcov(model, useScale = FALSE)))
#     zval<- beta / se
#     pval<- 2 * pnorm(abs(zval), lower.tail = FALSE)
#     out_beta[i]=as.numeric(beta[2])
#     out_se[i] = round(as.numeric(se[2]),digits = 3)
#     out_pvalue[i] = as.numeric(pval[2])
#     out_nobs[i]=as.numeric(nobs(model))   
#     out_transformation[i]=NA
#     }
#   }
#   regress_output<-data.frame(metabolite=colnames(dat),
#                              beta=out_beta,
#                              se=out_se,
#                              p_val=out_pvalue,
#                              n=out_nobs,
#                              p_val_Bonf=p.adjust(out_pvalue,method="bonferroni"),
#                              p_val_fdr=p.adjust(out_pvalue,method="BH"),
#                              transformation=out_transformation
#   )%>%
#     mutate(
#       pval_Bonf_neglog=-log10(p_val_Bonf),
#       pval_fdr_neglog=-log10(p_val_fdr),
#       sig_Bonf=ifelse(p_val_Bonf<=0.05,"Bonf<0.05","Not Sig"),
#       sig_fdr=ifelse(p_val_fdr<=0.05,"FDR<0.05","Not Sig"),
#       sig_Bonf=ifelse(p_val_Bonf<=0.05,"Bonf<0.05","Not Sig"),
#       sig_fdr=ifelse(p_val_fdr<=0.05,"FDR<0.05","Not Sig"),
#       is_continuous=metab_is_cont
#     )
#   
#   if (metab_is_complete==T){
#     regress_output$is_continuous<-"Complete-cases"
#   } else {
#     regress_output$is_continuous<-factor(regress_output$is_continuous,levels=c(TRUE,FALSE),labels=c("Imputed","Dichotomized"))
#   }
#   
#   return(regress_output)
# }


# # # Function to loop regression model and export model summary using survey weight and circular variables
#  archive_atlas_cir_svyreg_loop<-function(data,covar,end,metab_is_cont,metab_is_complete,trait,trait_as_predictor){
#   dat<-data[,1:end]
#   dat<-dat[,!colnames(dat)%in%c("LAB_ID","ID","SOL_ID","id")]
#   # dat<-data[,3:end]
#   # out_nvar=end-3+1
#   out_nvar=ncol(dat)
#   out_sin_beta<-rep(NA,out_nvar)
#   out_sin_se<-rep(NA,out_nvar)
#   out_sin_pvalue<-rep(NA,out_nvar)
#   out_cos_beta<-rep(NA,out_nvar)
#   out_cos_se<-rep(NA,out_nvar)
#   out_cos_pvalue<-rep(NA,out_nvar)
#   out_pvalue<-rep(NA,out_nvar)
#   out_1_beta<-rep(NA,out_nvar)
#   out_1_sd<-rep(NA,out_nvar)
#   out_1_lb<-rep(NA,out_nvar)
#   out_1_ub<-rep(NA,out_nvar)
#   out_2_beta<-rep(NA,out_nvar)
#   out_2_sd<-rep(NA,out_nvar)
#   out_2_lb<-rep(NA,out_nvar)
#   out_2_ub<-rep(NA,out_nvar)
# 
#   # if the original trait variable is binary or the user pick the binary version to use in the model
#   # trait_is_cont=ifelse(length(unique(data[!is.na(data[,trait]),trait]))==2,FALSE,TRUE)
#   if (all(grep("[A-Z]", names(data[,-1])) == 0)) { # if column names are all in lower case
#     survey_design=svydesign(id=~psu_id, strata=~strat, weights=~weight, data=data)
#   } else {
#     survey_design=svydesign(id=~PSU_ID, strata=~STRAT, weights=~WEIGHT, data=data)
#   }
#   # survey_design=svydesign(id=~PSU_ID, strata=~STRAT, weights=~WEIGHT, data=data)
# 
#   
#   
#   for(i in 1:(ncol(dat))) {
#     if (trait_as_predictor==T){
#       if (metab_is_cont==T){
#         model<- svyglm(as.formula(paste0(colnames(dat)[i],"~sin(",trait,")+cos(",trait,")+",paste(unlist(covar),collapse= "+"))),design=survey_design)
#       } else {
#         model<- svyglm(as.formula(paste0(colnames(dat)[i],"~sin(",trait,")+cos(",trait,")+",paste(unlist(covar),collapse = "+"))),design=survey_design,family=quasipoisson(link='log'))
#       }
#       # extract the coefficients and their covariance matrix
#       Vcov_matrix <- vcov(model, useScale = FALSE)
#       coef_matrix<- coef(model)
#       # Compute the estimated coefficient and its standard error for the circular predictor
#       coef_circ <- sqrt(coef_matrix[2]^2 + coef_matrix[3]^2)
#       se_circ <- sqrt((coef_matrix[2]/coef_circ)^2*Vcov_matrix[2, 2]+(coef_matrix[3]/coef_circ)^2*Vcov_matrix[3,3]+2*coef_matrix[2]*coef_matrix[3]/(coef_circ^3)*Vcov_matrix[2, 3])
#       
#       # Compute the t-statistic and p-value for the circular predictor
#       t_stat_circ <- coef_circ / se_circ
#       p_value_circ <- 2 * pt(abs(t_stat_circ), df = model$df.residual, lower.tail = FALSE)
#       
#       # # if the correlation betwee sin() and cos() is negative, which indicates the circular predictor doesn't have uniform distribution
#       # # Use the von Mises-Fisher regression model instead
#       # if (Vcov_matrix[2,3]<0){
#       #   require(CircStats)
#       #   # Define the von Mises-Fisher regression model log-likelihood function
#       #   vmf_ll <- function(beta, y, x) {
#       #     # Compute the predicted values for the circular predictor x
#       #     mu <- exp(cos(beta[1]) * x + sin(beta[1]) * sqrt(1 - x^2) + beta[2] * age)
#       #     # Compute the log-likelihood of the von Mises-Fisher distribution
#       #     ll <- sum(log(dvmf(y, mu, kappa = beta[3], log = TRUE)))
#       #     return(-ll)
#       #   }
#       #   
#       #   # Fit the von Mises-Fisher regression model using maximum likelihood estimation
#       #   if (metab_is_cont==T){
#       #     model<- svymle(as.formula(paste0(colnames(dat)[i],"~sin(",trait,")+cos(",trait,")+",paste(unlist(covar),collapse= "+"))),design=survey_designll = vmf_ll, start = c(0, 0, 1))
#       #   } else {
#       #     model<- svymle(as.formula(paste0(colnames(dat)[i],"~sin(",trait,")+cos(",trait,")+",paste(unlist(covar),collapse = "+"))),design=survey_design,family=quasipoisson(link='log'),ll = vmf_ll, start = c(0, 0, 1))
#       #   }
#       # }
#       # compute the t-statistics and p-values for the circular predictor
#       # t_stat <- coef_matrix[2] / sqrt(Vcov_matrix[2,2])
#       # p_value <- 2 * pt(abs(t_stat), df = model$df.residual, lower.tail = FALSE)
#       
#       # if the correlation betwee sin() and cos() is negative, use atan() transformation instead of sin() and cos()
#       if (Vcov_matrix[2,3]<0){
#         if (metab_is_cont==T){
#           model<- svyglm(as.formula(paste0(colnames(dat)[i],"~atan(",trait,")+",paste(unlist(covar),collapse= "+"))),design=survey_design)
#         } else {
#           model<- svyglm(as.formula(paste0(colnames(dat)[i],"~atan(",trait,")+",paste(unlist(covar),collapse = "+"))),design=survey_design,family=quasipoisson(link='log'))
#         }
#       }
#      
#       # hm_sin<-paste(names(coef(model))[2],"=0")
#       # hm_cos<-paste(names(coef(model))[3],"=0")
#       # hm<-car::linearHypothesis(model, c(hm_sin,hm_cos),test=c("F"),white.adjust = "hc1")
# 
#       se<- sqrt(diag(vcov(model, useScale = FALSE)))
#       zval<- beta / se
#       pval<- 2 * pnorm(abs(zval), lower.tail = FALSE)
#       out_sin_beta[i]=as.numeric(beta[2])
#       out_sin_se[i] = as.numeric(se[2])
#       out_sin_pvalue[i] = as.numeric(pval[2])
#       out_cos_beta[i]=as.numeric(beta[3])
#       out_cos_se[i] = as.numeric(se[3])
#       out_cos_pvalue[i] = as.numeric(pval[3])
#       out_pvalue[i]=as.numeric(hm[2,4])
#       regress_output<-data.frame(metabolite=colnames(dat),
#                                  sin_beta=out_sin_beta,
#                                  sin_se=out_sin_se,
#                                  sin_p_val=out_sin_pvalue,
#                                  sin_p_val_Bonf=p.adjust(out_sin_pvalue,method="bonferroni"),
#                                  sin_p_val_fdr=p.adjust(out_sin_pvalue,method="BH"),
#                                  cos_beta=out_cos_beta,
#                                  cos_se=out_cos_se,
#                                  cos_p_val=out_cos_pvalue,
#                                  cos_p_val_Bonf=p.adjust(out_cos_pvalue,method="bonferroni"),
#                                  cos_p_val_fdr=p.adjust(out_cos_pvalue,method="BH"),
#                                  p_val=out_pvalue,
#                                  p_val_Bonf=p.adjust(out_pvalue,method="bonferroni"),
#                                  p_val_fdr=p.adjust(out_pvalue,method="BH")
#       )%>%
#         mutate(
#           sin_pval_Bonf_neglog=-log10(sin_p_val_Bonf),
#           sin_pval_fdr_neglog=-log10(sin_p_val_fdr),
#           sin_sig_Bonf=ifelse(sin_p_val_Bonf<=0.05,"Bonf<0.05","Not Sig"),
#           sin_sig_fdr=ifelse(sin_p_val_fdr<=0.05,"FDR<0.05","Not Sig"),
#           cos_pval_Bonf_neglog=-log10(cos_p_val_Bonf),
#           cos_pval_fdr_neglog=-log10(cos_p_val_fdr),
#           cos_sig_Bonf=ifelse(cos_p_val_Bonf<=0.05,"Bonf<0.05","Not Sig"),
#           cos_sig_fdr=ifelse(cos_p_val_fdr<=0.05,"FDR<0.05","Not Sig"),
#           pval_Bonf_neglog=-log10(p_val_Bonf),
#           pval_fdr_neglog=-log10(p_val_fdr),
#           sig_Bonf=ifelse(p_val_Bonf<=0.05,"Bonf<0.05","Not Sig"),
#           sig_fdr=ifelse(p_val_fdr<=0.05,"FDR<0.05","Not Sig"),
#           is_continuous=metab_is_cont
#         )
#     } else {
#         model<- bpnreg::bpnr(pred.I=as.formula(paste0(trait,"~",colnames(dat)[i],"+",paste(unlist(covar),collapse= "+"))),its = 5000, burn=500, n.lag=3, seed=101, data=data)
#         out_1_beta[i]=model$lin.coef.I[2,1]
#         out_1_sd[i]=model$lin.coef.I[2,3]
#         out_1_lb[i]=model$lin.coef.I[2,4]
#         out_1_ub[i]=model$lin.coef.I[2,5]
#         out_2_beta[i]=model$lin.coef.II[2,1]
#         out_2_sd[i]=model$lin.coef.II[2,3]
#         out_2_lb[i]=model$lin.coef.II[2,4]
#         out_2_ub[i]=model$lin.coef.II[2,5]
#         regress_output<-data.frame(metabolite=colnames(dat),
#                                    i_beta=out_1_beta,
#                                    i_sd=out_1_sd,
#                                    i_lb=out_1_lb,
#                                    i_ub=out_1_ub,
#                                    ii_beta=out_2_beta,
#                                    ii_sd=out_2_sd,
#                                    ii_lb=out_2_lb,
#                                    ii_ub=out_2_ub
#         )%>%
#           mutate(
#             bayeci_1=ifelse(i_ub<0&i_lb<0,-1,ifelse(i_ub>0&i_lb>0,1,0)),
#             bayeci_2=ifelse(ii_ub<0&ii_lb<0,-1,ifelse(ii_ub>0&ii_lb>0,1,0)),
#             is_continuous=metab_is_cont
#           )
#     }
# 
#   }
#   if (metab_is_complete==T){
#     regress_output$is_continuous<-"Complete-cases"
#   } else {
#     regress_output$is_continuous<-factor(regress_output$is_continuous,levels=c(TRUE,FALSE),labels=c("Imputed","Dichotomized"))
#   }
# 
#   return(regress_output)
# }


# Function to loop regression model and export model summary using survey weight and circular variables by modeling the circular response variables as sin(y) and cos(y)
# for the multi-response regression
atlas_cir_svyreg_loop<-function(data,covar,end,metab_is_cont,metab_is_complete,trait,trait_as_predictor){
  dat<-data[,1:end]
  dat<-dat[,!colnames(dat)%in%c("LAB_ID","ID","SOL_ID","id")]
  # dat<-data[,3:end]
  # out_nvar=end-3+1
  out_nvar=ncol(dat)
  out_sin_beta<-rep(NA,out_nvar)
  out_sin_se<-rep(NA,out_nvar)
  out_sin_pvalue<-rep(NA,out_nvar)
  out_cos_beta<-rep(NA,out_nvar)
  out_cos_se<-rep(NA,out_nvar)
  out_cos_pvalue<-rep(NA,out_nvar)
  out_pvalue<-rep(NA,out_nvar)
  out_intercept<-rep(NA,out_nvar)
  out_intercept_se<-rep(NA,out_nvar)
  out_nobs<-rep(NA,out_nvar)
  # if the original trait variable is binary or the user pick the binary version to use in the model
  # trait_is_cont=ifelse(length(unique(data[!is.na(data[,trait]),trait]))==2,FALSE,TRUE)
    if (all(grep("[A-Z]", names(data[,-1])) == 0)) { # if column names are all in lower case
      survey_design=svydesign(id=~psu_id, strata=~strat, weights=~weight, data=data)
    } else {
      survey_design=svydesign(id=~PSU_ID, strata=~STRAT, weights=~WEIGHT, data=data)
    }
  # survey_design=svydesign(id=~PSU_ID, strata=~STRAT, weights=~WEIGHT, data=data)

  for(i in 1:(ncol(dat))) {
    if (trait_as_predictor==T){
      if (metab_is_cont==T){
        require(survey)
        require(jtools)
        model<- svyglm(as.formula(paste0(colnames(dat)[i],"~sin(",trait,")+cos(",trait,")+",paste(unlist(covar),collapse= "+"))),design=survey_design)
        # corr<-svycor(as.formula(paste0("~sin(",trait,")+cos(",trait,")")), design = survey_design,na.rm = T)
      } else {
        require(survey)
        require(jtools)
        model<- svyglm(as.formula(paste0(colnames(dat)[i],"~sin(",trait,")+cos(",trait,")+",paste(unlist(covar),collapse = "+"))),design=survey_design,family=quasipoisson(link='log'))
        # corr<-svycor(as.formula(paste0("~sin(",trait,")+cos(",trait,")")), design = survey_design,na.rm = T)
      }
      hm_sin<-paste(names(coef(model))[2],"=0")
      hm_cos<-paste(names(coef(model))[3],"=0")
      hm<-car::linearHypothesis(model, c(hm_sin,hm_cos),test=c("F"),white.adjust = "hc1")
      Vcov <- vcov(model, useScale = FALSE)
      beta<- coef(model)
      se<- sqrt(diag(vcov(model, useScale = FALSE)))
      zval<- beta / se
      pval<- 2 * pnorm(abs(zval), lower.tail = FALSE)
      out_sin_beta[i]=as.numeric(beta[2])
      out_sin_se[i] = as.numeric(se[2])
      out_sin_pvalue[i] = as.numeric(pval[2])
      out_cos_beta[i]=as.numeric(beta[3])
      out_cos_se[i] = as.numeric(se[3])
      out_cos_pvalue[i] = as.numeric(pval[3])
      out_pvalue[i]=as.numeric(hm[2,4])
      out_intercept[i]=as.numeric(coef(model)[1])
      out_intercept_se[i]=as.numeric(se[1])
      out_nobs[i]=as.numeric(nobs(model))

    } else {
      # model<- svyglm(as.formula(paste0("cbind(",trait,".sin,",trait,".cos)","~",colnames(dat)[i],"+",paste(unlist(covar),collapse= "+"))),design=survey_design)
      # model<-svyglm(cbind(wdwake_time_cir.sin,wdwake_time_cir.cos)~uridine+AGE+GENDER+CENTER+background,design=survey_design)
      model<-estimatr::lm_robust(formula=as.formula(paste0("cbind(",trait,".sin,",trait,".cos)","~",colnames(dat)[i],"+",paste(unlist(covar),collapse= "+"))),data=data, weights = WEIGHT, clusters=STRAT)
      test<-manova(as.formula(paste0("cbind(",trait,".sin,",trait,".cos)","~",colnames(dat)[i],"+",paste(unlist(covar),collapse= "+"))),data=data, weights = WEIGHT)
      pval<-summary(test,test="Pillai")$stats[1,6]

      out_sin_beta[i]=as.numeric(model$coefficients[2,1])
      out_sin_se[i] = as.numeric(model$std.error[2,1])
      out_sin_pvalue[i] = as.numeric(model$p.value[2,1])
      out_cos_beta[i]=as.numeric(model$coefficients[2,2])
      out_cos_se[i] = as.numeric(model$std.error[2,2])
      out_cos_pvalue[i] = as.numeric(model$p.value[2,2])
      out_pvalue[i]=as.numeric(pval)
      out_intercept[i]=as.numeric(coef(model)[1])
      out_intercept_se[i]=as.numeric(se[1])
      out_nobs[i]=as.numeric(nobs(model))
      # out_1_beta[i]=model$lin.coef.I[2,1]
      # out_1_sd[i]=model$lin.coef.I[2,3]
      # out_1_lb[i]=model$lin.coef.I[2,4]
      # out_1_ub[i]=model$lin.coef.I[2,5]
      # out_2_beta[i]=model$lin.coef.II[2,1]
      # out_2_sd[i]=model$lin.coef.II[2,3]
      # out_2_lb[i]=model$lin.coef.II[2,4]
      # out_2_ub[i]=model$lin.coef.II[2,5]
      # regress_output<-data.frame(metabolite=colnames(dat),
      #                            i_beta=out_1_beta,
      #                            i_sd=out_1_sd,
      #                            i_lb=out_1_lb,
      #                            i_ub=out_1_ub,
      #                            ii_beta=out_2_beta,
      #                            ii_sd=out_2_sd,
      #                            ii_lb=out_2_lb,
      #                            ii_ub=out_2_ub
      # )%>%
      #   mutate(
      #     bayeci_1=ifelse(i_ub<0&i_lb<0,-1,ifelse(i_ub>0&i_lb>0,1,0)),
      #     bayeci_2=ifelse(ii_ub<0&ii_lb<0,-1,ifelse(ii_ub>0&ii_lb>0,1,0)),
      #     is_continuous=metab_is_cont
      #   )
    }
    regress_output<-data.frame(metabolite=colnames(dat),
                               n=out_nobs,
                               sin_beta=out_sin_beta,
                               sin_se=out_sin_se,
                               sin_p_val=out_sin_pvalue,
                               sin_p_val_Bonf=p.adjust(out_sin_pvalue,method="bonferroni"),
                               sin_p_val_fdr=p.adjust(out_sin_pvalue,method="BH"),
                               cos_beta=out_cos_beta,
                               cos_se=out_cos_se,
                               cos_p_val=out_cos_pvalue,
                               cos_p_val_Bonf=p.adjust(out_cos_pvalue,method="bonferroni"),
                               cos_p_val_fdr=p.adjust(out_cos_pvalue,method="BH"),
                               p_val=out_pvalue,
                               p_val_Bonf=p.adjust(out_pvalue,method="bonferroni"),
                               p_val_fdr=p.adjust(out_pvalue,method="BH"),
                               mesor=out_intercept,
                               mesor_se=out_intercept_se
    )%>%
      mutate(
        sin_pval_Bonf_neglog=-log10(sin_p_val_Bonf),
        sin_pval_fdr_neglog=-log10(sin_p_val_fdr),
        sin_sig_Bonf=ifelse(sin_p_val_Bonf<=0.05,"Bonf<0.05","Not Sig"),
        sin_sig_fdr=ifelse(sin_p_val_fdr<=0.05,"FDR<0.05","Not Sig"),
        cos_pval_Bonf_neglog=-log10(cos_p_val_Bonf),
        cos_pval_fdr_neglog=-log10(cos_p_val_fdr),
        cos_sig_Bonf=ifelse(cos_p_val_Bonf<=0.05,"Bonf<0.05","Not Sig"),
        cos_sig_fdr=ifelse(cos_p_val_fdr<=0.05,"FDR<0.05","Not Sig"),
        pval_Bonf_neglog=-log10(p_val_Bonf),
        pval_fdr_neglog=-log10(p_val_fdr),
        sig_Bonf=ifelse(p_val_Bonf<=0.05,"Bonf<0.05","Not Sig"),
        sig_fdr=ifelse(p_val_fdr<=0.05,"FDR<0.05","Not Sig"),
        is_continuous=metab_is_cont
      )

  }
  if (metab_is_complete==T){
    regress_output$is_continuous<-"Complete-cases"
  } else {
    regress_output$is_continuous<-factor(regress_output$is_continuous,levels=c(TRUE,FALSE),labels=c("Imputed","Dichotomized"))
  }

  return(regress_output)
}


# Function to loop regression model and export model summary using survey weight
atlas_cir_strat_svyreg_loop<-function(data,covar,end,metab_is_cont,metab_is_complete,trait,trait_for_model,trait_as_predictor,stratifier){
  covar<-unlist(covar)
  covar<-covar[!covar%in%stratifier]
  group_by_value<-unique(data[,stratifier])
  out_nvar<-length(colnames(data[,1:end])[!colnames(data[,1:end])%in%c("LAB_ID","ID","SOL_ID","id")])
  # out_nvar<-end-3+1
  out_sin_beta<-rep(NA,out_nvar)
  out_sin_se<-rep(NA,out_nvar)
  out_sin_pvalue<-rep(NA,out_nvar)
  out_cos_beta<-rep(NA,out_nvar)
  out_cos_se<-rep(NA,out_nvar)
  out_cos_pvalue<-rep(NA,out_nvar)
  out_pvalue<-rep(NA,out_nvar)
  out_intercept<-rep(NA,out_nvar)
  out_intercept_se<-rep(NA,out_nvar)
  out_nobs<-rep(NA,out_nvar)
  group<-rep(NA,out_nvar)
  regress_output<-data.frame()
  # if the original trait variable is binary or the user pick the binary version to use in the model
  # trait_is_cont<-ifelse(length(unique(data[!is.na(data[,trait]),trait]))==2,FALSE,TRUE)
  
  for (j in seq_along(group_by_value)){
    newdata<-data[which(data[,stratifier]==as.character(group_by_value[j])),]
    dat<-newdata[,1:end]
    dat<-dat[,!colnames(dat)%in%c("LAB_ID","ID","SOL_ID","id")]
    if (all(grep("[A-Z]", names(data[,-1])) == 0)) { # if column names are all in lower case
      survey_design=svydesign(id=~psu_id, strata=~strat, weights=~weight, data=newdata)
    } else {
      survey_design=svydesign(id=~PSU_ID, strata=~STRAT, weights=~WEIGHT, data=newdata)
    }
    for(i in 1:(ncol(dat))) {
      if (trait_as_predictor==T){
        if (metab_is_cont==T){
          require(survey)
          require(jtools)
          model<- svyglm(as.formula(paste0(colnames(dat)[i],"~sin(",trait,")+cos(",trait,")+",paste(covar,collapse= "+"))),design=survey_design)
        } else {
          require(survey)
          require(jtools)
          model<- svyglm(as.formula(paste0(colnames(dat)[i],"~sin(",trait,")+cos(",trait,")+",paste(covar,collapse = "+"))),design=survey_design,family=quasipoisson(link='log'))
        }
        hm_sin<-paste(names(coef(model))[2],"=0")
        hm_cos<-paste(names(coef(model))[3],"=0")
        hm<-car::linearHypothesis(model, c(hm_sin,hm_cos),test=c("F"),white.adjust = "hc1")
        Vcov <- vcov(model, useScale = FALSE)
        beta<- coef(model)
        se<- sqrt(diag(vcov(model, useScale = FALSE)))
        zval<- beta / se
        pval<- 2 * pnorm(abs(zval), lower.tail = FALSE)
        out_sin_beta[i]=as.numeric(beta[2])
        out_sin_se[i] = as.numeric(se[2])
        out_sin_pvalue[i] = as.numeric(pval[2])
        out_cos_beta[i]=as.numeric(beta[3])
        out_cos_se[i] = as.numeric(se[3])
        out_cos_pvalue[i] = as.numeric(pval[3])
        out_pvalue[i]=as.numeric(hm[2,4])
        out_intercept[i]=as.numeric(coef(model)[1])
        out_intercept_se[i]=as.numeric(se[1])
        out_nobs[i]=as.numeric(nobs(model))
        group[i]=as.character(group_by_value[j])

      }else {
        model<-estimatr::lm_robust(formula=as.formula(paste0("cbind(",trait,".sin,",trait,".cos)","~",colnames(dat)[i],"+",paste(unlist(covar),collapse= "+"))),data=newdata, weights = WEIGHT, clusters=STRAT)
        test<-manova(as.formula(paste0("cbind(",trait,".sin,",trait,".cos)","~",colnames(dat)[i],"+",paste(unlist(covar),collapse= "+"))),data=data, weights = WEIGHT)
        pval<-summary(test,test="Pillai")$stats[1,6]
        group[i]=as.character(group_by_value[j])
        out_sin_beta[i]=as.numeric(model$coefficients[2,1])
        out_sin_se[i] = as.numeric(model$std.error[2,1])
        out_sin_pvalue[i] = as.numeric(model$p.value[2,1])
        out_cos_beta[i]=as.numeric(model$coefficients[2,2])
        out_cos_se[i] = as.numeric(model$std.error[2,2])
        out_cos_pvalue[i] = as.numeric(model$p.value[2,2])
        out_pvalue[i]=as.numeric(pval)
        out_intercept[i]=as.numeric(coef(model)[1])
        out_intercept_se[i]=as.numeric(se[1])
        out_nobs[i]=as.numeric(nobs(model))
        # model<- bpnreg::bpnr(pred.I=as.formula(paste0(trait,"~",colnames(dat)[i],"+",paste(unlist(covar),collapse= "+"))),its = 5000, burn=500, n.lag=3, seed=101, data=newdata)
        # out_1_sd[i]=model$lin.coef.I[2,3]
        # out_1_lb[i]=model$lin.coef.I[2,4]
        # out_1_ub[i]=model$lin.coef.I[2,5]
        # out_2_beta[i]=model$lin.coef.II[2,1]
        # out_2_sd[i]=model$lin.coef.II[2,3]
        # out_2_lb[i]=model$lin.coef.II[2,4]
        # out_2_ub[i]=model$lin.coef.II[2,5]
        # strata_output<-data.frame(metabolite=colnames(dat),
        #                            i_beta=out_1_beta,
        #                            i_sd=out_1_sd,
        #                            i_lb=out_1_lb,
        #                            i_ub=out_1_ub,
        #                            ii_beta=out_2_beta,
        #                            ii_sd=out_2_sd,
        #                            ii_lb=out_2_lb,
        #                            ii_ub=out_2_ub,
        #                            strata=group
        # )%>%
        #   mutate(
        #     bayeci_1=ifelse(i_ub<0&i_lb<0,-1,ifelse(i_ub>0&i_lb>0,1,0)),
        #     bayeci_2=ifelse(ii_ub<0&ii_lb<0,-1,ifelse(ii_ub>0&ii_lb>0,1,0)),
        #     is_continuous=metab_is_cont
        #   )
      }
      strata_output<-data.frame(metabolite=colnames(dat),
                                n=out_nobs,
                                sin_beta=out_sin_beta,
                                sin_se=out_sin_se,
                                sin_p_val=out_sin_pvalue,
                                sin_p_val_Bonf=p.adjust(out_sin_pvalue,method="bonferroni"),
                                sin_p_val_fdr=p.adjust(out_sin_pvalue,method="BH"),
                                cos_beta=out_cos_beta,
                                cos_se=out_cos_se,
                                cos_p_val=out_cos_pvalue,
                                cos_p_val_Bonf=p.adjust(out_cos_pvalue,method="bonferroni"),
                                cos_p_val_fdr=p.adjust(out_cos_pvalue,method="BH"),
                                p_val=out_pvalue,
                                p_val_Bonf=p.adjust(out_pvalue,method="bonferroni"),
                                p_val_fdr=p.adjust(out_pvalue,method="BH"),
                                mesor=out_intercept,
                                mesor_se=out_intercept_se,
                                strata=group)%>%
        mutate(
          sin_pval_Bonf_neglog=-log10(sin_p_val_Bonf),
          sin_pval_fdr_neglog=-log10(sin_p_val_fdr),
          sin_sig_Bonf=ifelse(sin_p_val_Bonf<=0.05,"Bonf<0.05","Not Sig"),
          sin_sig_fdr=ifelse(sin_p_val_fdr<=0.05,"FDR<0.05","Not Sig"),
          cos_pval_Bonf_neglog=-log10(cos_p_val_Bonf),
          cos_pval_fdr_neglog=-log10(cos_p_val_fdr),
          cos_sig_Bonf=ifelse(cos_p_val_Bonf<=0.05,"Bonf<0.05","Not Sig"),
          cos_sig_fdr=ifelse(cos_p_val_fdr<=0.05,"FDR<0.05","Not Sig"),
          pval_Bonf_neglog=-log10(p_val_Bonf),
          pval_fdr_neglog=-log10(p_val_fdr),
          sig_Bonf=ifelse(p_val_Bonf<=0.05,"Bonf<0.05","Not Sig"),
          sig_fdr=ifelse(p_val_fdr<=0.05,"FDR<0.05","Not Sig"),
          is_continuous=metab_is_cont
        )
    }
    regress_output<-rbind(regress_output,strata_output)
  }
  # identified_data<-regress_output%>%
  #   filter(!str_detect(metabolite,"^X"))%>%
  #   group_by(strata)%>%
  #   mutate(p_val_fdr_named=p.adjust(p_val,method="BH"),
  #          sig_fdr_named=ifelse(p_val_fdr_named<=0.05,"FDR<0.05","Not Sig"))
  # unidentified_data<-regress_output%>%
  #   filter(str_detect(metabolite,"^X"))%>%
  #   mutate(p_val_fdr_named=NA,
  #          sig_fdr_named=NA)
  # regress_output<-plyr::rbind.fill(identified_data,unidentified_data)
  if (metab_is_complete==T){
    regress_output$is_continuous<-"Complete-cases"
  } else {
    regress_output$is_continuous<-factor(regress_output$is_continuous,levels=c(TRUE,FALSE),labels=c("Imputed","Dichotomized"))
  }
  
  return(regress_output)
}






# from htc502/ewrefxn documentation 
# function to convert long format data frame to matrix for pheatmap
tbl2hmap= function(tbl.long, rowVar, colVar, valueVar,
                   colAnnVars = NULL, rowAnnVars = NULL) {
  Mat0 = dplyr::select(tbl.long, one_of(c(rowVar, colVar,valueVar)))%>%
    tidyr::spread_(key = colVar, value=valueVar)
  Mat = select(Mat0, -one_of(rowVar)) %>% as.matrix()
  rownames(Mat) = unlist(select(Mat0,one_of(rowVar)))
  
  colAnn = rowAnn= NULL
  
  if(!is.null(colAnnVars)) {
    colAnn0 = dplyr::select(tbl.long,one_of( c(colVar,colAnnVars))) %>%
      unique()
    colAnn = select(colAnn0, - one_of(colVar)) %>% as.data.frame()
    rownames(colAnn) = unlist(select( colAnn0,one_of(colVar)))
    # row.names(colAnn) <-  make.unique(unlist(select( colAnn0,one_of(colVar))))
  }
  
  if(!is.null(rowAnnVars)) {
    rowAnn0 = dplyr::select(tbl.long,one_of(c( rowVar,rowAnnVars))) %>%
      unique()
    rowAnn = select(rowAnn0, -one_of( rowVar)) %>% as.data.frame()
    rownames(rowAnn) = unlist(select(rowAnn0,one_of(rowVar)))
    # row.names(rowAnn) <-  make.unique(unlist(select( rowAnn0,one_of(rowVar))))
  }
  
  output2<-list(mat = Mat, rowAnn =rowAnn, colAnn = colAnn)
}


# function to convert long format data frame to matrix for pheatmap
tbl2hmap_sort= function(tbl.long, rowVar, colVar, valueVar,
                   colAnnVars = NULL, rowAnnVars = NULL) {
  Mat0 = dplyr::select(tbl.long, one_of(c(rowVar, colVar,valueVar)))%>%
    tidyr::spread_(key = colVar, value=valueVar)
  Mat = select(Mat0, -one_of(rowVar)) %>% as.matrix()
  rownames(Mat) = unlist(select(Mat0,one_of(rowVar)))
  Mat_sort<-cbind(rownames(Mat),Mat)
  
  colAnn = rowAnn= NULL
  
  if(!is.null(colAnnVars)) {
    colAnn0 = dplyr::select(tbl.long,one_of( c(colVar,colAnnVars))) %>%
      unique()
    colAnn = select(colAnn0, - one_of(colVar)) %>% as.data.frame()
    # rownames(colAnn) = unlist(select( colAnn0,one_of(colVar)))
    row.names(colAnn) <-  make.unique(unlist(select( colAnn0,one_of(colVar))))
    col_sort<-tbl.long[,c(colVar,colAnnVars)]
    # col_sort<-col_sort[,order(colAnnVars,colVar)]
    Mat_sort<-Mat_sort[match(col_sort[,1],Mat_sort[,1]),] 
    rownames(Mat_sort)<-Mat_sort[,1]
  }
  
  if(!is.null(rowAnnVars)) {
    rowAnn0 = dplyr::select(tbl.long,one_of(c( rowVar,rowAnnVars))) %>%
      unique()
    rowAnn = select(rowAnn0, -one_of( rowVar)) %>% as.data.frame()
    # rownames(rowAnn) = unlist(select(rowAnn0,one_of(rowVar)))
    row.names(rowAnn) <-  make.unique(unlist(select( rowAnn0,one_of(rowVar))))
    row_sort<-tbl.long[,c(rowVar,rowAnnVars)]
    # row_sort<-row_sort[,order(rowAnnVars,rowVar)]
    Mat_sort<-Mat_sort[match(row_sort[,1],Mat_sort[,1]),] 
    rownames(Mat_sort)<-Mat_sort[,1]
  }
  # Mat_sorted<-Mat_sort[,-1]
  rn<-Mat_sort[,1]
  Mat_sorted<-apply(Mat_sort[,-1], 2, as.numeric)
  rownames(Mat_sorted)<-rn
  # Mat_sorted<-as.matrix(Mat_sort,rownames=TRUE)
  class(Mat_sorted) <- "numeric"
  Mat_sorted2<-sapply(Mat_sorted[,-1], as.numeric)
  output<-list(mat = Mat, rowAnn =rowAnn, colAnn = colAnn)
}


# calculate the time difference between the bed/get-up time and midnight in hr (positive value means after midnight)
hr_ps_midnight<-function(hr,mn){
  hr<-as.numeric(hr)
  mn<-as.numeric(mn)
  diff<-ifelse(!is.na(hr),ifelse(!is.na(mn),
                                 ifelse(hr>12,(24-hr)*60-mn,
                                        (0-hr)*60-mn),
                                 ifelse(hr>12,(24-hr)*60,
                                        (0-hr)*60)),NA)
  return(diff)
}
# functions to calculate mean squared prediction error and root mean squared error to validate different imputation method
validation_function<-function(dataset,trait,covar, trait_as_predictor,imp_method) {
  data_fit<-dataset[which(dataset[,"ID"]%in%random_id),]
  data_test<-dataset[which(!dataset[,"ID"]%in%random_id),]
  metab_var<-colnames(dataset)[!colnames(dataset)%in%c(colnames(pheno_complete),"LAB_ID")]
  trait_is_cont<-ifelse(length(unique(dataset[!is.na(dataset[,trait]),trait]))==2,FALSE,TRUE)
  survey_design=svydesign(id=~PSU_ID, strata=~STRAT, weights=~WEIGHT, data=data_fit)
  
  out_mspe=rep(NA,length(metab_var))
  out_mae=rep(NA,length(metab_var))
  out_rmse=rep(NA,length(metab_var))
  
  for(i in seq_along(metab_var)) {
    if (trait_as_predictor==T){
      model<- svyglm(as.formula(paste0(metab_var[i],"~",trait,"+",paste(covar,collapse= "+"))),design=survey_design)
      pred_test<-predict(model,newdata=data_test,type="response")
      mspe<-mean((data_test[,metab_var[i]]-pred_test[1])^2)
      mae<-mean(abs(data_test[,metab_var[i]]-pred_test[1]))
      rmse<-sqrt(mean((data_test[,metab_var[i]]-pred_test[1])^2))
      
    } else {
      if (trait_is_cont==T){
        model<- svyglm(as.formula(paste0(trait,"~",metab_var,"+",paste(covar,collapse= "+"))),design=survey_design)
        pred_test<-predict(model,newdata=data_test,type="response")
        mspe<-mean((data_test[,metab_var[i]]-pred_test[1])^2)
        mae<-mean(abs(data_test[,metab_var[i]]-pred_test[1]))
        rmse<-sqrt(mean((data_test[,metab_var[i]]-pred_test[1])^2))
      } else {
        model<- svyglm(as.formula(paste0(trait,"~",metab_var,"+",paste(covar,collapse = "+"))),design=survey_design,family=quasipoisson(link='log'))
        pred_test<-predict(model,newdata=data_test,type="response")
        mspe<-mean((data_test[,metab_var[i]]-pred_test[1])^2)
        mae<-mean(abs(data_test[,metab_var[i]]-pred_test[1]))
        rmse<-sqrt(mean((data_test[,metab_var[i]]-pred_test[1])^2))
      }
    }
    # Vcov <- vcov(model, useScale = FALSE)
    # beta<- coef(model)
    # se<- sqrt(diag(vcov(model, useScale = FALSE)))
    # zval<- beta / se
    # pval<- 2 * pnorm(abs(zval), lower.tail = FALSE)
    # out_beta[i]=as.numeric(beta[2])
    # out_se[i] = round(as.numeric(se[2]),digits = 3)
    # out_pvalue[i] = as.numeric(pval[2])
    out_mspe[i]=mspe
    out_mae[i]=mae
    out_rmse[i]=rmse
  }
  # regress_output<-data.frame(metabolite=metab_var,
  #                            beta=out_beta,
  #                            se=out_se,
  #                            p_val=out_pvalue,
  #                            p_val_Bonf=round(p.adjust(out_pvalue,method="bonferroni"),digits = 3),
  #                            p_val_fdr=p.adjust(out_pvalue,method="BH")
  # )%>%
  #   mutate(
  #     pval_Bonf_neglog=round(-log10(p_val_Bonf),digits = 3),
  #     pval_fdr_neglog=-log10(p_val_fdr),
  #     sig_Bonf=ifelse(p_val_Bonf<=0.05,"Bonf<0.05","Not Sig"),
  #     sig_fdr=ifelse(p_val_fdr<=0.05,"FDR<0.05","Not Sig"),
  #     is_continuous=metab_is_cont
  #   )
  # identified_data<-regress_output%>%
  #   filter(!str_detect(metabolite,"^X"))%>%
  #   mutate(p_val_fdr_named=p.adjust(p_val,method="BH"),
  #          sig_fdr_named=ifelse(p_val_fdr_named<=0.05,"FDR<0.05","Not Sig"))
  # unidentified_data<-regress_output%>%
  #   filter(str_detect(metabolite,"^X"))%>%
  #   mutate(p_val_fdr_named=NA,
  #          sig_fdr_named=NA)
  # regress_output<-rbind(identified_data,unidentified_data)
  accuracy_output<-data.frame(metabolite=metab_var,
                              mspe=out_mspe,
                              mae=out_mae,
                              rmse=out_rmse,
                              imp=imp_method)
  return(accuracy_output)
}

# function to read in all csv files in a folder and rbind
load_data <- function(path) { 
  files <- dir(path, pattern = '\\.csv', full.names = TRUE)
  tables <- lapply(files, read.csv)
  do.call(dplyr::bind_rows, tables)
}

# function to read in all csv files in a folder, rbind and use the file name as id
load_data_name<-function(path) { 
  files.list <- list.files(path=path,pattern='*.csv')
  files <- dir(path, pattern = '\\.csv', full.names = TRUE)
  df.list <- setNames(lapply(files,read.csv ), files.list)
  df <- dplyr::bind_rows(df.list, .id = "id")
}

# function to get the percentage of metabolite meeting the exclusion criteria
conditional_tally<-function(data,col,equal_sign,value){
  metab_count<-data%>%
    group_by(metabolite)%>%
    summarise(n=n())
  if (equal_sign=="<") {
    filter_data<-data[which(data[,col]<value),]
  } else if (equal_sign%in%c("=<","<=")) {
    filter_data<-data[which(data[,col]<=value),]
  } else if (equal_sign%in%c(">=","=>")) {
    filter_data<-data[which(data[,col]>=value),]
  } else if (equal_sign==">") {
    filter_data<-data[which(data[,col]>value),]
  } else if (equal_sign=="=") {
    filter_data<-data[which(data[,col]==value),]
  } else if (equal_sign=="!=") {
    filter_data<-data[which(data[,col]!=value),]
  } else stop("operator not accepted!")
  
  metab_neg<-filter_data%>%
    group_by(metabolite)%>%
    summarise(neg_n=n())
  metab_exclude<-merge(metab_neg,metab_count,by="metabolite")%>%
    mutate(perc_p=neg_n/n)
  return(metab_exclude)
}

# remove outliers
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 3 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

# identify the ID of outliers
identify_outliers <- function(data,id,var, na.rm = TRUE, ...) {
  outlier_id_sum<-NA
  for (i in seq_along(var)){
    # print(i)
    x<-var[i]
    qnt <- quantile(data[,x], probs=c(.25, .75), na.rm = na.rm)
    H <- 3 * IQR(data[,x], na.rm = na.rm)
    y<-data.frame(ID=data[,id],X=data[,x])
    miss_id<-y[is.na(y[,2]),1]
    y[which(y$X < (qnt[1] - H)),2] <- NA
    y[which(y$X > (qnt[2] + H)),2] <- NA
    out_id<-y[is.na(y[,2]&!y[,1]%in%miss_id),1]
    if (length(out_id)==0) {
      outlier_id_sum<-outlier_id_sum
    } else {
      outlier_id<-data.frame(ID=out_id,var=rep(x,length(out_id)))
      colnames(outlier_id)[2]<-paste0(x)
      if(is.na(outlier_id_sum)[1]){
        outlier_id_sum<-outlier_id
      } else {
        outlier_id_sum<-merge(outlier_id_sum,outlier_id,by="ID",all=TRUE)
      }
    }
    # print(i)
  }
  
  return(outlier_id_sum)
}

# convert factor to numeric without loosing the level value
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}

# add new column for sin and cos
add.cal<-function(data){
  sindata<-sin(data*pi/12)
  cosdata<-cos(data*pi/12)
  sinName<-paste0(names(data),"sin")
  cosName<-paste0(names(data),"cos")
  newdata<-cbind(sindata,cosdata)
  colnames(newdata)<-c(sinName,cosName)
  return(newdata)
}


# convert a long format data set into individual tables in which metabolites are row, and each column is the value (either p or beta) for each trait
long_to_wide<-function(data,var_included,value_ind,var_factor,trait_order){
  # var_value: variable of which the value will be used to populate the whole wide table and be used to generate heatmap
  # var_included: string of variable names that are going to be included in the wide table (usually metaboltie, var_value, super_pathway, sub_pathway)
  # var_ind: which varible in the var_included string that's var_value. If following the example above, should be 2.
  # var_factor: index of variable to be formatted to factor
  wide_data<-data.frame()
  by_var<-var_included[-value_ind]
  for (i in seq_along(unique(data$trait))) {
    temp_data<-data[which(data$trait==unique(data$trait)[i]),var_included]
    colnames(temp_data)[value_ind]<-paste0(unique(data$trait)[i])
    if (i==1) {
      wide_data<-temp_data
    } else {
      wide_data<-merge(wide_data,temp_data, by=as.vector(by_var))
    }
  }
  wide_data[,var_factor]<-lapply(wide_data[,var_factor], as.factor)
  wide_data<-wide_data[,trait_order]%>%
    tibble::column_to_rownames(.,var="metabolite")
  return(wide_data)
}
  
# functions to prepare data set for heatmap 
heatmap_comparison<-function(data, met_list,trait_list,trait_group_list, trait_index, main_row_order,met_font,trait_font,plot_height){
  # data: long format data (regress_results)
  # met_list: vector of metabolites that should be included in the heatmap
  # trait_group_list: a factor vector assining trait group to trait (e.g. trait_group_list<-factor(c(rep("AHI",2),rep("Disturbance",10),rep("Heart-rate",4),rep("Timing",6))))
  # trait_index: column index of the matrix to indicate which columns should be included in the heatmap (e.g. 3:24)
  # main_row_order: row order extracted from the main heatmap
  
  # derive beta_p=unadjusted p value times the +/- sign of the beta coefficient
  data$beta_p<-ifelse(data$beta>=0,-log10(data$p_val),-log10(data$p_val)*(-1))
  # include metabolites from the list
  data<-data[which(data[,"metabolite"]%in%met_list),]
  # sort by multiple columns
  data.table::setorder(data,trait,super_pathway,sub_pathway,metabolite)
  
  # convert a long format data set into individual tables in which metabolites are row, and each column is the value (either p or beta) for each trait
  both_beta<-long_to_wide(data=data[which(data[,"strata"]=="Both"),],  var_included=c("metabolite","beta","super_pathway","sub_pathway"),value_ind = 2,var_factor=c("super_pathway","sub_pathway"),trait_order=trait_list)
  both_p<-long_to_wide(data=data[which(data[,"strata"]=="Both"),], var_included=c("metabolite","beta_p","super_pathway","sub_pathway"),value_ind = 2,var_factor=c("super_pathway","sub_pathway"),trait_order=trait_list)
  female_beta<-long_to_wide(data=data[which(data[,"strata"]=="Female"),], var_included=c("metabolite","beta","super_pathway","sub_pathway"),value_ind = 2,var_factor=c("super_pathway","sub_pathway"),trait_order=trait_list)
  female_p<-long_to_wide(data=data[which(data[,"strata"]=="Female"),], var_included=c("metabolite","beta_p","super_pathway","sub_pathway"),value_ind = 2,var_factor=c("super_pathway","sub_pathway"),trait_order=trait_list)
  male_beta<-long_to_wide(data=data[which(data[,"strata"]=="Male"),], var_included=c("metabolite","beta","super_pathway","sub_pathway"),value_ind = 2,var_factor=c("super_pathway","sub_pathway"),trait_order=trait_list)
  male_p<-long_to_wide(data=data[which(data[,"strata"]=="Male"),], var_included=c("metabolite","beta_p","super_pathway","sub_pathway"),value_ind = 2,var_factor=c("super_pathway","sub_pathway"),trait_order=trait_list)

  # metabolite annotation color
  sp_color <- colorRampPalette(RColorBrewer::brewer.pal(12,"Paired"))(length(levels(both_p$super_pathway)))
  names(sp_color) <- levels(both_p$super_pathway)
  row_ha <- HeatmapAnnotation(df = data.frame(super_pathway=both_p$super_pathway), which="row", col = list(super_pathway = sp_color),show_annotation_name = FALSE)
  # trait annotation color
  # trait_group_list<-factor(c(rep("AHI",2),rep("Disturbance",10),rep("Heart-rate",4),rep("Timing",6)))
  trait_color <- colorRampPalette(RColorBrewer::brewer.pal(8,"Pastel2"))(length(levels(trait_group_list)))
  names(trait_color) <- levels(trait_group_list)
  col_ha <- HeatmapAnnotation(df = data.frame(trait_group=trait_group_list), which="column", col = list(trait_group = trait_color))
  
  both_p_ht<-Heatmap((as.matrix(both_p[,trait_index])), name = "Both:-log10(p)", row_title= NULL, column_title = "Both:Unadjusted p",
                     row_split=both_p$super_pathway,
                     row_gap = unit(1, "mm"),
                     column_split = trait_group_list,
                     top_annotation = col_ha,
                     left_annotation = row_ha,
                     border = TRUE,
                     cluster_rows = F,
                     cluster_columns = F,
                     show_column_dend = FALSE,
                     column_names_rot = 45,
                     column_names_max_height = max_text_width(
                       rownames(both_p[,trait_index]),
                       gp = gpar(fontsize = 0.5)),
                     row_names_side = "left",
                     row_dend_side="right",
                     column_names_gp = gpar(fontsize=trait_font),
                     row_names_max_width = unit(6, "cm"),
                     row_names_gp = gpar(fontsize = met_font),
                     row_order=unlist(main_row_order),
                     heatmap_height = unit(plot_height,"cm")
  )
  
  both_beta_ht<-Heatmap((as.matrix(both_beta[,trait_index])), name = "Both:beta", row_title= NULL, column_title = "Both:beta coefficient",
                         row_split=both_beta$super_pathway,
                         row_gap = unit(1, "mm"),
                         column_split = trait_group_list,
                         top_annotation = col_ha,
                         left_annotation = row_ha,
                         border = TRUE,
                         cluster_rows = F,
                         cluster_columns = F,
                         show_column_dend = FALSE,
                         column_names_rot = 45,
                         column_names_max_height = max_text_width(
                           rownames(both_beta[,trait_index]),
                           gp = gpar(fontsize = 0.5)),
                         row_names_side = "left",
                         row_dend_side="right",
                         column_names_gp = gpar(fontsize=trait_font),
                         row_names_max_width = unit(6, "cm"),
                         row_names_gp = gpar(fontsize = met_font),
                         row_order=unlist(main_row_order),
                        heatmap_height = unit(plot_height,"cm")
  )
  male_p_ht<-Heatmap((as.matrix(male_p[,trait_index])), name = "Male:-log10(p)", row_title= NULL, column_title = "Male:unadjusted p",
                     row_split=male_p$super_pathway,
                     row_gap = unit(1, "mm"),
                     column_split = trait_group_list,
                     top_annotation = col_ha,
                     left_annotation = row_ha,
                     border = TRUE,
                     cluster_rows = F,
                     cluster_columns = F,
                     show_column_dend = FALSE,
                     column_names_rot = 45,
                     column_names_max_height = max_text_width(
                       rownames(male_p[,trait_index]),
                       gp = gpar(fontsize = 0.5)),
                     row_names_side = "left",
                     row_dend_side="right",
                     column_names_gp = gpar(fontsize=trait_font),
                     row_names_max_width = unit(6, "cm"),
                     row_names_gp = gpar(fontsize = met_font),
                     row_order=unlist(main_row_order),
                     heatmap_height = unit(plot_height,"cm")
  )
  male_beta_ht<-Heatmap((as.matrix(male_beta[,trait_index])), name = "Male:beta", row_title= NULL, column_title = "Male:beta coefficient",
                     row_split=male_beta$super_pathway,
                     row_gap = unit(1, "mm"),
                     column_split = trait_group_list,
                     top_annotation = col_ha,
                     left_annotation = row_ha,
                     border = TRUE,
                     cluster_rows = F,
                     cluster_columns = F,
                     show_column_dend = FALSE,
                     column_names_rot = 45,
                     column_names_max_height = max_text_width(
                       rownames(male_beta[,trait_index]),
                       gp = gpar(fontsize = 0.5)),
                     row_names_side = "left",
                     row_dend_side="right",
                     column_names_gp = gpar(fontsize=trait_font),
                     row_names_max_width = unit(6, "cm"),
                     row_names_gp = gpar(fontsize = met_font),
                     row_order=unlist(main_row_order),
                     heatmap_height = unit(plot_height,"cm")
  )
  female_p_ht<-Heatmap((as.matrix(female_p[,trait_index])), name = "Female:-log10(p)", row_title= NULL, column_title = "Female:unadjusted p",
                       row_split=female_p$super_pathway,
                       row_gap = unit(1, "mm"),
                       column_split = trait_group_list,
                       top_annotation = col_ha,
                       left_annotation = row_ha,
                       border = TRUE,
                       cluster_rows = F,
                       cluster_columns = F,
                       show_column_dend = FALSE,
                       column_names_rot = 45,
                       column_names_max_height = max_text_width(
                         rownames(female_p[,trait_index]),
                         gp = gpar(fontsize = 0.5)),
                       row_names_side = "left",
                       row_dend_side="right",
                       column_names_gp = gpar(fontsize=trait_font),
                       row_names_max_width = unit(6, "cm"),
                       row_names_gp = gpar(fontsize = met_font),
                       row_order=unlist(main_row_order),
                       heatmap_height = unit(plot_height,"cm")
  )
  female_beta_ht<-Heatmap((as.matrix(female_beta[,trait_index])), name = "Female:beta", row_title= NULL, column_title = "Female: beta coefficient",
                          row_split=female_beta$super_pathway,
                          row_gap = unit(1, "mm"),
                          column_split = trait_group_list,
                          top_annotation = col_ha,
                          left_annotation = row_ha,
                          border = TRUE,
                          cluster_rows = F,
                          cluster_columns = F,
                          show_column_dend = FALSE,
                          column_names_rot = 45,
                          column_names_max_height = max_text_width(
                            rownames(female_beta[,trait_index]),
                            gp = gpar(fontsize = 0.5)),
                          row_names_side = "left",
                          row_dend_side="right",
                          column_names_gp = gpar(fontsize=trait_font),
                          row_names_max_width = unit(6, "cm"),
                          row_names_gp = gpar(fontsize = met_font),
                          row_order=unlist(main_row_order),
                          heatmap_height = unit(plot_height,"cm")
  )
  ht_list=both_p_ht+both_beta_ht+female_p_ht+female_beta_ht+male_p_ht+male_beta_ht
  return(ht_list)
}

# heatmap comparison functions that only takes one strata (either female or male) across three models
heatmap_comparison_pfilter<-function(data, met_list,trait_list,trait_group_list, trait_index, main_row_order,met_font,trait_font,plot_height){
  # data: long format data (regress_results)
  # met_list: vector of metabolites that should be included in the heatmap
  # trait_group_list: a factor vector assining trait group to trait (e.g. trait_group_list<-factor(c(rep("AHI",2),rep("Disturbance",10),rep("Heart-rate",4),rep("Timing",6))))
  # trait_index: column index of the matrix to indicate which columns should be included in the heatmap (e.g. 3:24)
  # main_row_order: row order extracted from the main heatmap
  
  # derive beta_p=unadjusted p value times the +/- sign of the beta coefficient
  # data$beta_p<-ifelse(data$beta>=0,-log10(data$p_val),-log10(data$p_val)*(-1))
  # include metabolites from the list
  data<-data[which(data[,"metabolite"]%in%met_list),]
  # sort by multiple columns
  data.table::setorder(data,trait,super_pathway,sub_pathway,metabolite)
  
  # convert a long format data set into individual tables in which metabolites are row, and each column is the value (either p or beta) for each trait
  md1_beta<-long_to_wide(data=data[which(data[,"model"]=="Model_1"),],  var_included=c("metabolite","beta","super_pathway","sub_pathway","pfilter_1"),value_ind = 2,var_factor=c("super_pathway","sub_pathway"),trait_order=trait_list)
  md1_p<-long_to_wide(data=data[which(data[,"model"]=="Model_1"),], var_included=c("metabolite","beta_p","super_pathway","sub_pathway","pfilter_1"),value_ind = 2,var_factor=c("super_pathway","sub_pathway"),trait_order=trait_list)
  md2_beta<-long_to_wide(data=data[which(data[,"model"]=="Model_2"),], var_included=c("metabolite","beta","super_pathway","sub_pathway","pfilter_1"),value_ind = 2,var_factor=c("super_pathway","sub_pathway"),trait_order=trait_list)
  md2_p<-long_to_wide(data=data[which(data[,"model"]=="Model_2"),], var_included=c("metabolite","beta_p","super_pathway","sub_pathway","pfilter_1"),value_ind = 2,var_factor=c("super_pathway","sub_pathway"),trait_order=trait_list)
  md3_beta<-long_to_wide(data=data[which(data[,"model"]=="Model_3"),], var_included=c("metabolite","beta","super_pathway","sub_pathway","pfilter_1"),value_ind = 2,var_factor=c("super_pathway","sub_pathway"),trait_order=trait_list)
  md3_p<-long_to_wide(data=data[which(data[,"model"]=="Model_3"),], var_included=c("metabolite","beta_p","super_pathway","sub_pathway","pfilter_1"),value_ind = 2,var_factor=c("super_pathway","sub_pathway"),trait_order=trait_list)
  
  # metabolite annotation color
  sp_color <- colorRampPalette(RColorBrewer::brewer.pal(12,"Paired"))(length(levels(md1_p$super_pathway)))
  names(sp_color) <- levels(md1_p$super_pathway)
  row_ha <- HeatmapAnnotation(df = data.frame(super_pathway=md1_p$super_pathway), which="row", col = list(super_pathway = sp_color),show_annotation_name = FALSE)
  # trait annotation color
  # trait_group_list<-factor(c(rep("AHI",2),rep("Disturbance",10),rep("Heart-rate",4),rep("Timing",6)))
  trait_color <- colorRampPalette(RColorBrewer::brewer.pal(8,"Pastel2"))(length(levels(trait_group_list)))
  names(trait_color) <- levels(trait_group_list)
  col_ha <- HeatmapAnnotation(df = data.frame(trait_group=trait_group_list), which="column", col = list(trait_group = trait_color))
  
  md1_p_ht<-Heatmap((as.matrix(md1_p[,trait_index])), name = "Mdl_1:-log10(p)", row_title= NULL, column_title = "Model 1:Unadjusted p",
                     row_split=md1_p$super_pathway,
                     row_gap = unit(1, "mm"),
                     column_split = trait_group_list,
                     top_annotation = col_ha,
                     left_annotation = row_ha,
                     border = TRUE,
                     cluster_rows = F,
                     cluster_columns = F,
                     show_column_dend = FALSE,
                     column_names_rot = 45,
                     column_names_max_height = max_text_width(
                       rownames(md1_p[,trait_index]),
                       gp = gpar(fontsize = 0.5)),
                     row_names_side = "left",
                     row_dend_side="right",
                     column_names_gp = gpar(fontsize=trait_font),
                     row_names_max_width = unit(6, "cm"),
                     row_names_gp = gpar(fontsize = met_font),
                     row_order=unlist(main_row_order),
                     heatmap_height = unit(plot_height,"cm")
  )
  
  md1_beta_ht<-Heatmap((as.matrix(md1_beta[,trait_index])), name = "Mdl_1:beta", row_title= NULL, column_title = "Model 1:beta coefficient",
                        row_split=md1_beta$super_pathway,
                        row_gap = unit(1, "mm"),
                        column_split = trait_group_list,
                        top_annotation = col_ha,
                        left_annotation = row_ha,
                        border = TRUE,
                        cluster_rows = F,
                        cluster_columns = F,
                        show_column_dend = FALSE,
                        column_names_rot = 45,
                        column_names_max_height = max_text_width(
                          rownames(md1_beta[,trait_index]),
                          gp = gpar(fontsize = 0.5)),
                        row_names_side = "left",
                        row_dend_side="right",
                        column_names_gp = gpar(fontsize=trait_font),
                        row_names_max_width = unit(6, "cm"),
                        row_names_gp = gpar(fontsize = met_font),
                        row_order=unlist(main_row_order),
                        heatmap_height = unit(plot_height,"cm")
  )
  md2_p_ht<-Heatmap((as.matrix(md2_p[,trait_index])), name = "Mdl_2:-log10(p)", row_title= NULL, column_title = "Model 2:unadjusted p",
                     row_split=md2_p$super_pathway,
                     row_gap = unit(1, "mm"),
                     column_split = trait_group_list,
                     top_annotation = col_ha,
                     left_annotation = row_ha,
                     border = TRUE,
                     cluster_rows = F,
                     cluster_columns = F,
                     show_column_dend = FALSE,
                     column_names_rot = 45,
                     column_names_max_height = max_text_width(
                       rownames(md2_p[,trait_index]),
                       gp = gpar(fontsize = 0.5)),
                     row_names_side = "left",
                     row_dend_side="right",
                     column_names_gp = gpar(fontsize=trait_font),
                     row_names_max_width = unit(6, "cm"),
                     row_names_gp = gpar(fontsize = met_font),
                     row_order=unlist(main_row_order),
                     heatmap_height = unit(plot_height,"cm")
  )
  md2_beta_ht<-Heatmap((as.matrix(md2_beta[,trait_index])), name = "Mdl_2:beta", row_title= NULL, column_title = "Model 2:beta coefficient",
                        row_split=md2_beta$super_pathway,
                        row_gap = unit(1, "mm"),
                        column_split = trait_group_list,
                        top_annotation = col_ha,
                        left_annotation = row_ha,
                        border = TRUE,
                        cluster_rows = F,
                        cluster_columns = F,
                        show_column_dend = FALSE,
                        column_names_rot = 45,
                        column_names_max_height = max_text_width(
                          rownames(md2_beta[,trait_index]),
                          gp = gpar(fontsize = 0.5)),
                        row_names_side = "left",
                        row_dend_side="right",
                        column_names_gp = gpar(fontsize=trait_font),
                        row_names_max_width = unit(6, "cm"),
                        row_names_gp = gpar(fontsize = met_font),
                        row_order=unlist(main_row_order),
                        heatmap_height = unit(plot_height,"cm")
  )
  md3_p_ht<-Heatmap((as.matrix(md3_p[,trait_index])), name = "Mdl_3:-log10(p)", row_title= NULL, column_title = "Model 3:unadjusted p",
                       row_split=md3_p$super_pathway,
                       row_gap = unit(1, "mm"),
                       column_split = trait_group_list,
                       top_annotation = col_ha,
                       left_annotation = row_ha,
                       border = TRUE,
                       cluster_rows = F,
                       cluster_columns = F,
                       show_column_dend = FALSE,
                       column_names_rot = 45,
                       column_names_max_height = max_text_width(
                         rownames(md3_p[,trait_index]),
                         gp = gpar(fontsize = 0.5)),
                       row_names_side = "left",
                       row_dend_side="right",
                       column_names_gp = gpar(fontsize=trait_font),
                       row_names_max_width = unit(6, "cm"),
                       row_names_gp = gpar(fontsize = met_font),
                       row_order=unlist(main_row_order),
                       heatmap_height = unit(plot_height,"cm")
  )
  md3_beta_ht<-Heatmap((as.matrix(md3_beta[,trait_index])), name = "Mdl_3:beta", row_title= NULL, column_title = "Model 3: beta coefficient",
                          row_split=md3_beta$super_pathway,
                          row_gap = unit(1, "mm"),
                          column_split = trait_group_list,
                          top_annotation = col_ha,
                          left_annotation = row_ha,
                          border = TRUE,
                          cluster_rows = F,
                          cluster_columns = F,
                          show_column_dend = FALSE,
                          column_names_rot = 45,
                          column_names_max_height = max_text_width(
                            rownames(md3_beta[,trait_index]),
                            gp = gpar(fontsize = 0.5)),
                          row_names_side = "left",
                          row_dend_side="right",
                          column_names_gp = gpar(fontsize=trait_font),
                          row_names_max_width = unit(6, "cm"),
                          row_names_gp = gpar(fontsize = met_font),
                          row_order=unlist(main_row_order),
                          heatmap_height = unit(plot_height,"cm")
  )
  ht_list=md1_p_ht+md1_beta_ht+md2_p_ht+md2_beta_ht+md3_p_ht+md3_beta_ht
  return(ht_list)
}

# functions to prepare data set for heatmap 
heatmap_beta_comparison<-function(data, met_list,trait_list,trait_group_list, trait_index, main_row_order,met_font,trait_font,plot_height){
  # data: long format data (regress_results)
  # met_list: vector of metabolites that should be included in the heatmap
  # trait_group_list: a factor vector assining trait group to trait (e.g. trait_group_list<-factor(c(rep("AHI",2),rep("Disturbance",10),rep("Heart-rate",4),rep("Timing",6))))
  # trait_index: column index of the matrix to indicate which columns should be included in the heatmap (e.g. 3:24)
  # main_row_order: row order extracted from the main heatmap
  # include metabolites from the list
  data<-data[which(data$metabolite%in%met_list),]
  # convert a long format data set into individual tables in which metabolites are row, and each column is the value (either p or beta) for each trait
  both_beta<-long_to_wide(data=data[which(data[,"strata"]=="Both"),],  var_included=c("metabolite","outbound","super_pathway","sub_pathway"),value_ind = 2,var_factor=c("super_pathway","sub_pathway"),trait_order=trait_list)
  female_beta<-long_to_wide(data=data[which(data[,"strata"]=="Female"),], var_included=c("metabolite","outbound","super_pathway","sub_pathway"),value_ind = 2,var_factor=c("super_pathway","sub_pathway"),trait_order=trait_list)
  male_beta<-long_to_wide(data=data[which(data[,"strata"]=="Male"),], var_included=c("metabolite","outbound","super_pathway","sub_pathway"),value_ind = 2,var_factor=c("super_pathway","sub_pathway"),trait_order=trait_list)
  
  # metabolite annotation color
  sp_color <- colorRampPalette(RColorBrewer::brewer.pal(12,"Paired"))(length(levels(both_beta$super_pathway)))
  names(sp_color) <- levels(both_beta$super_pathway)
  row_ha <- HeatmapAnnotation(df = data.frame(super_pathway=both_beta$super_pathway), which="row", col = list(super_pathway = sp_color),show_annotation_name = FALSE)
  # trait annotation color
  # trait_group_list<-factor(c(rep("AHI",2),rep("Disturbance",10),rep("Heart-rate",4),rep("Timing",6)))
  trait_color <- colorRampPalette(RColorBrewer::brewer.pal(8,"Pastel2"))(length(levels(trait_group_list)))
  names(trait_color) <- levels(trait_group_list)
  col_ha <- HeatmapAnnotation(df = data.frame(trait_group=trait_group_list), which="column", col = list(trait_group = trait_color))
  
  both_beta_ht<-Heatmap((as.matrix(both_beta[,trait_index])), name = "Both:change of beta", row_title= NULL, column_title = "Both:beta outside of 95% CI",
                        row_split=both_beta$super_pathway,
                        row_gap = unit(1, "mm"),
                        column_split = trait_group_list,
                        top_annotation = col_ha,
                        left_annotation = row_ha,
                        border = TRUE,
                        cluster_rows = F,
                        cluster_columns = F,
                        show_column_dend = FALSE,
                        column_names_rot = 45,
                        column_names_max_height = max_text_width(
                          rownames(both_beta[,trait_index]),
                          gp = gpar(fontsize = 0.5)),
                        row_names_side = "left",
                        row_dend_side="right",
                        column_names_gp = gpar(fontsize=trait_font),
                        row_names_max_width = unit(6, "cm"),
                        row_names_gp = gpar(fontsize = met_font),
                        row_order=unlist(main_row_order),
                        heatmap_height = unit(plot_height,"cm")
  )
  
  male_beta_ht<-Heatmap((as.matrix(male_beta[,trait_index])), name = "Male:change of beta", row_title= NULL, column_title = "Male:beta outside of 95% CI",
                        row_split=male_beta$super_pathway,
                        row_gap = unit(1, "mm"),
                        column_split = trait_group_list,
                        top_annotation = col_ha,
                        left_annotation = row_ha,
                        border = TRUE,
                        cluster_rows = F,
                        cluster_columns = F,
                        show_column_dend = FALSE,
                        column_names_rot = 45,
                        column_names_max_height = max_text_width(
                          rownames(male_beta[,trait_index]),
                          gp = gpar(fontsize = 0.5)),
                        row_names_side = "left",
                        row_dend_side="right",
                        column_names_gp = gpar(fontsize=trait_font),
                        row_names_max_width = unit(6, "cm"),
                        row_names_gp = gpar(fontsize = met_font),
                        row_order=unlist(main_row_order),
                        heatmap_height = unit(plot_height,"cm")
  )
  
  female_beta_ht<-Heatmap((as.matrix(female_beta[,trait_index])), name = "Female:change of beta", row_title= NULL, column_title = "Female:beta outside of 95% CI",
                          row_split=female_beta$super_pathway,
                          row_gap = unit(1, "mm"),
                          column_split = trait_group_list,
                          top_annotation = col_ha,
                          left_annotation = row_ha,
                          border = TRUE,
                          cluster_rows = F,
                          cluster_columns = F,
                          show_column_dend = FALSE,
                          column_names_rot = 45,
                          column_names_max_height = max_text_width(
                            rownames(female_beta[,trait_index]),
                            gp = gpar(fontsize = 0.5)),
                          row_names_side = "left",
                          row_dend_side="right",
                          column_names_gp = gpar(fontsize=trait_font),
                          row_names_max_width = unit(6, "cm"),
                          row_names_gp = gpar(fontsize = met_font),
                          row_order=unlist(main_row_order),
                          heatmap_height = unit(plot_height,"cm")
  )
  ht_list=both_beta_ht+female_beta_ht+male_beta_ht
  return(ht_list)
}

# Function to plot age, gender, bmi, physical activity, dietary score, HTN status and DM status of top 100 vs the rest for the top two hit metabolite for each sleep trait
demo_bymetab<-function(metdata,phenodata,metabolite){
  # metdata<-metdata[order(metdata[,metabolite],na.last = TRUE, decreasing = T),]
  # top_met<-metdata[,colnames(metdata)%in%c("ID",metabolite)]
  pheno_bymetab<-phenodata[,c("ID","AGE","GENDER","BMI","AHEI2010","GPAQ_TOTAL_MET","HYPERTENSION","DIABETES2_INDICATOR")]%>%
    merge(.,metdata[,colnames(metdata)%in%c("ID",metabolite)],by="ID")
  pheno_bymetab<-pheno_bymetab[order(metdata[,metabolite],na.last = TRUE, decreasing = T),]
  pheno_bymetab$top100<-"Others"
  pheno_bymetab[1:100,"top100"]<-"Top100"
  
  # age_plot<-ggplot(pheno_bymetab[which(!is.na(pheno_bymetab$AGE)),c("AGE","top100")], aes(x = AGE, fill = top100)) +
  #   geom_density(position="identity", alpha=0.5) +
  #   scale_x_continuous(name = "Age") +
  #   scale_y_continuous(name = "Density") +
  #   ggtitle("Age")
  bmi_plot<-ggplot(pheno_bymetab[which(!is.na(pheno_bymetab$BMI)),c("BMI","top100")], aes(x = BMI, fill = top100)) +
    geom_density(position="identity", alpha=0.5) +
    scale_x_continuous(name = "BMI") +
    scale_y_continuous(name = "Density") +
    ggtitle("BMI")
  phy_plot<-ggplot(pheno_bymetab[which(!is.na(pheno_bymetab$GPAQ_TOTAL_MET)),c("GPAQ_TOTAL_MET","top100")], aes(x = GPAQ_TOTAL_MET, fill = top100)) +
    geom_density(position="identity", alpha=0.5) +
    scale_x_continuous(name = "Physical Activity in Met-hr") +
    scale_y_continuous(name = "Density") +
    ggtitle("Physical activity (Met-hr)")
  ahei_plot<-ggplot(pheno_bymetab[which(!is.na(pheno_bymetab$AHEI2010)),c("AHEI2010","top100")], aes(x = AHEI2010, fill = top100)) +
    geom_density(position="identity", alpha=0.5) +
    scale_x_continuous(name = "Dietary Score (AHEI2010)") +
    scale_y_continuous(name = "Density") +
    ggtitle("Dietary Score (AHEI2010)")
  
  
  # working violin plot with boxplot overlay
  age_plot<-ggplot(pheno_bymetab[which(!is.na(pheno_bymetab$AGE)),c("AGE","top100")], aes(x = top100,y=AGE, fill=top100)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.2,hjust=0.5)+
    theme(legend.position = "none")+
    scale_y_continuous(name = "Age") +
    ggtitle("Age")
  
  bmi_plot<-ggplot(pheno_bymetab[which(!is.na(pheno_bymetab$BMI)),c("BMI","top100")], aes(x = top100,y=BMI, fill=top100)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.2,hjust=0.5)+
    theme(legend.position = "none")+
    scale_y_continuous(name = "BMI") +
    ggtitle("BMI")
  # 
  # phy_plot<-ggplot(pheno_bymetab[which(!is.na(pheno_bymetab$GPAQ_TOTAL_MET)),c("GPAQ_TOTAL_MET","top100")], aes(x = top100,y=GPAQ_TOTAL_MET, fill=top100)) +
  #   geom_violin(trim = FALSE) +
  #   geom_boxplot(width = 0.2,hjust=0.5)+
  #   theme(legend.position = "none")+
  #   scale_y_continuous(name = "Physical Activity in Met-hr") +
  #   ggtitle("Physical Activity in Met-hr")
  # 
  ahei_plot<-ggplot(pheno_bymetab[which(!is.na(pheno_bymetab$AHEI2010)),c("AHEI2010","top100")], aes(x = top100,y=AHEI2010, fill=top100)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.2,hjust=0.5)+
    theme(legend.position = "none")+
    scale_y_continuous(name = "AHEI2010") +
    ggtitle("AHEI2010")
  
  gender_temp<-pheno_bymetab[which(!is.na(pheno_bymetab$GENDER)),c("GENDER","top100")]
  gender_temp<-data.table::setDT(gender_temp)[,list(count = .N), by = .(top100,GENDER)][,list(GENDER = GENDER, count = count,
                                                                                              percent_fmt = paste0(formatC(count*100/sum(count), digits = 2), "%"),
                                                                                              percent_num = count/sum(count)
  ), by = top100]
  gender_plot<-ggplot(data=gender_temp, aes(x=top100, y= percent_num, fill=GENDER)) +   
    geom_bar(stat = "identity", width=0.7) +
    geom_text(aes(label = percent_fmt),position = position_stack(vjust = 0.5))+
    ggtitle("Gender")
  
  htn_temp<-pheno_bymetab[which(!is.na(pheno_bymetab$HYPERTENSION)),c("HYPERTENSION","top100")]
  htn_temp<-data.table::setDT(htn_temp)[,list(count = .N), by = .(top100,HYPERTENSION)][,list(HYPERTENSION = HYPERTENSION, count = count,
                                                                                              percent_fmt = paste0(formatC(count*100/sum(count), digits = 2), "%"),
                                                                                              percent_num = count/sum(count)
  ), by = top100]
  htn_plot<-ggplot(data=htn_temp, aes(x=top100, y= percent_num, fill=HYPERTENSION)) +   
    geom_bar( stat = "identity", width=0.7) +
    geom_text(aes(label = percent_fmt),position = position_stack(vjust = 0.5))+
    ggtitle("Hypertension")+ labs(fill = "HTN")
  
  dm_temp<-pheno_bymetab[which(!is.na(pheno_bymetab$DIABETES2_INDICATOR)),c("DIABETES2_INDICATOR","top100")]
  dm_temp<-data.table::setDT(dm_temp)[,list(count = .N), by = .(top100,DIABETES2_INDICATOR)][,list(DIABETES2_INDICATOR = DIABETES2_INDICATOR, count = count,
                                                                                                   percent_fmt = paste0(formatC(count*100/sum(count), digits = 2), "%"),
                                                                                                   percent_num = count/sum(count)
  ), by = top100]
  dm_plot<-ggplot(data=dm_temp, aes(x=top100, y= percent_num, fill=DIABETES2_INDICATOR)) +   
    geom_bar(stat = "identity", width=0.7) +
    geom_text(aes(label = percent_fmt),position = position_stack(vjust = 0.5))+
    ggtitle("Diabetes")+ labs(fill = "DM")
  
  
  grid.arrange(
    age_plot,bmi_plot,phy_plot,ahei_plot,gender_plot,htn_plot,dm_plot,
    nrow = 2,
    top = paste0("Comparison between top 100 obs\nwith the hightest concentration of ", metabolite," and the rest")
  )
}

#########################################
# functions to regress against covariates and take the residuals as the new value
residual_val<-function(met_data,pheno_data,dep_list,cov){
  data.residuals<-met_data[,c("ID","LAB_ID")]
  data<-merge(met_data,pheno_data,by="ID",all.y=T)%>%
    droplevels()
  # survey_design=svydesign(id=~PSU_ID, strata=~STRAT, weights=~WEIGHT, data=data)
  survey_design=svydesign(id=~PSU_ID, strata=~STRAT, weights=~WEIGHT_FINAL_NORM_OVERALL, data=data) # when regress sleep trait and covariates, should use baseline weights
  
  for (k in seq_along(dep_list)) {
    # if the trait variable is binary
    # print(paste("k=",k))
    var<-dep_list[k]
    trait_is_cont=ifelse(length(unique(data[!is.na(data[,var]),var]))==2,FALSE,TRUE)
    output<-data.frame()
    comp_id<-data[complete.cases(data[,c(dep_list[k],cov)]),"ID"]
    if (trait_is_cont==T){
      model<- svyglm(as.formula(paste(dep_list[k],"~",paste(cov,collapse= "+"))),design=subset(survey_design,ID%in%comp_id))
    } else {
      model<- svyglm(as.formula(paste("as.numeric(",dep_list[k],")~",paste(cov,collapse= "+"))),design=subset(survey_design,ID%in%comp_id),family=quasipoisson(link='log'))
    }
    output<-cbind(comp_id,model$residuals)
    colnames(output)<-c("ID",paste(dep_list[k]))
    data.residuals<-merge(x=data.residuals, y=output, by="ID", all = TRUE)
    # print(paste("k=",k))
  }
  return(data.residuals)
  # data.residuals
}

############################################
# function to loop over multiple outcomes and exposures using sampling weight
# functions using sampling weight - multivariate regression model (binary outcome)
fun_svyglm_multi<-function(dat,out_start,out_end,exp_start,exp_end,covar,is_factor,interaction){
  # define the column number for the function
  # outcome
  out_nvar=out_end-out_start+1
  
  out_variable=rep(NA,out_nvar)
  out_beta=rep(NA,out_nvar)
  out_se=rep(NA,out_nvar)
  out_pvalue=rep(NA,out_nvar)
  out_nobs=rep(NA,out_nvar)
  # exposure
  exp_nvar=exp_end-exp_start+1
  
  exp_mdl=rep(NA,exp_nvar)
  exp_variable=rep(NA,exp_nvar)
  exp_beta=rep(NA,exp_nvar)
  exp_se=rep(NA,exp_nvar)
  exp_pvalue=rep(NA,exp_nvar)
  exp_nobs=rep(NA,out_nvar)
  number=1
  
  # define the survey design
  survey_design=survey::svydesign(id=~PSU_ID, strata=~STRAT, weights=~WEIGHT,data=dat)
  for (i in out_start:out_end){
    # print(i)
    outcome=colnames(dat)[i]
    for (j in exp_start:exp_end){
      # print(j)
      exposure=colnames(dat)[j]
      # model<-glm(as.formula(paste(outcome,"~",exposure)), family=binomial(link = "logit"), data=dat)
      # model<- svyglm(as.formula(paste(outcome,"~",exposure,"+",paste(covar,collapse = "+"))),
      #                design=survey_design,family=quasipoisson(link='log'))
      if (is.null(interaction)){
        if (is_factor==T){
          model<- survey::svyglm(as.formula(paste(outcome,"~",exposure,"+",paste(covar,collapse = "+"))),
                                 design=survey_design,family=quasipoisson(link='log'))
        }else{
          model<- survey::svyglm(as.formula(paste(outcome,"~as.numeric(",exposure,")","+",paste(covar,collapse = "+"))),
                                 design=survey_design,family=quasipoisson(link='log'))
        }
      }else {
        if (is_factor==T){
          model<- survey::svyglm(as.formula(paste(outcome,"~",exposure,"+",paste(covar,collapse = "+"),"+",exposure,"*",interaction)),
                                 design=survey_design,family=quasipoisson(link='log'))
        }else{
          model<- survey::svyglm(as.formula(paste(outcome,"~as.numeric(",exposure,")","+",paste(covar,collapse = "+"),"+",exposure,"*",interaction)),
                                 design=survey_design,family=quasipoisson(link='log'))
        }
      }
      
      Vcov <- vcov(model, useScale = FALSE)
      beta <- coef(model)
      se <- sqrt(diag(Vcov))
      zval <- beta / se
      pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
      noobs<-nobs(model)
      
      for (k in 2:length(beta)){
        out_beta[number]=as.numeric(beta[k])
        out_se[number] = as.numeric(se[k])
        out_pvalue[number] = as.numeric(pval[k])
        out_variable[number] = outcome
        out_nobs[number]=noobs[1]
        number = number + 1
        exp_beta[number] = as.numeric(beta[k])
        exp_se[number] = as.numeric(se[k])
        exp_pvalue[number] = as.numeric(pval[k])
        exp_variable[number]=names(beta)[k]
        exp_nobs[number]=noobs[1]
        exp_mdl[number]=exposure
        number = number + 1
      }
      # print(j)
    }
    # print(i)
  }
  outcome = data.frame(out_variable, out_beta, out_se, out_pvalue,out_nobs)
  exposure = data.frame(exp_variable, exp_beta, exp_se, exp_pvalue,exp_mdl,exp_nobs)
  outcome = outcome %>% 
    dplyr::rename(
      variable = out_variable,
      beta = out_beta,
      se = out_se,
      pvalue = out_pvalue,
      obs = out_nobs
    )
  exposure = exposure %>% 
    dplyr::rename(
      variable = exp_variable,
      beta = exp_beta,
      se = exp_se,
      pvalue = exp_pvalue,
      obs = exp_nobs,
      exposure_model=exp_mdl
    )
  all=merge(outcome,exposure,by=c("beta","se","pvalue","obs"))
  # all = rbind(outcome, exposure)
  all = na.omit(all)
  data = all %>%
    dplyr::rename(outcome=variable.x,
           exposure=variable.y)%>%
    mutate (
      beta = round(beta, 5),
      se = round(se, 5),
      pvalue = round(pvalue, 5)
    ) %>%
    mutate(or=exp(beta),
           lower95=exp(beta-1.96*se),
           upper95=exp(beta+1.96*se),
           flag_pvalue_high=ifelse(pvalue<0.05,1,0),
           flag_pvalue_low=ifelse(pvalue<0.2,1,0))%>%
    dplyr::select(outcome,exposure_model,exposure, beta, se, or,lower95,upper95,pvalue,flag_pvalue_high,flag_pvalue_low,obs)%>%
    arrange(outcome,exposure)
  return(data)
}

# functions using sampling weight - multivariate regression model (continous outcome)
fun_svylm_multi<-function(dat,out_start,out_end,exp_start,exp_end,covar,is_factor,interaction){
  # define the column number for the function
  # outcome
  out_nvar=out_end-out_start+1
  
  out_variable=rep(NA,out_nvar)
  out_beta=rep(NA,out_nvar)
  out_se=rep(NA,out_nvar)
  out_pvalue=rep(NA,out_nvar)
  out_nobs=rep(NA,out_nvar)
  # exposure
  exp_nvar=exp_end-exp_start+1
  
  exp_mdl=rep(NA,exp_nvar)
  exp_variable=rep(NA,exp_nvar)
  exp_beta=rep(NA,exp_nvar)
  exp_se=rep(NA,exp_nvar)
  exp_pvalue=rep(NA,exp_nvar)
  exp_nobs=rep(NA,out_nvar)
  number=1
  
  # define the survey design
  survey_design=survey::svydesign(id=~PSU_ID, strata=~STRAT, weights=~WEIGHT,data=dat)
  for (i in out_start:out_end){
    # print(i)
    outcome=colnames(dat)[i]
    for (j in exp_start:exp_end){
      # print(j)
      exposure=colnames(dat)[j]
      if (!is.null(covar)){
        if (is.null(interaction)){
          if (is_factor==T){
            model<- survey::svyglm(as.formula(paste(outcome,"~",exposure,"+",paste(covar,collapse = "+"))),
                                   design=survey_design)
          }else{
            model<- survey::svyglm(as.formula(paste(outcome,"~as.numeric(",exposure,")","+",paste(covar,collapse = "+"))),
                                   design=survey_design)
          }
        }else {
          if (is_factor==T){
            model<- survey::svyglm(as.formula(paste(outcome,"~",exposure,"+",paste(covar,collapse = "+"),"+",exposure,"*",interaction)),
                                   design=survey_design)
          }else{
            model<- survey::svyglm(as.formula(paste(outcome,"~as.numeric(",exposure,")","+",paste(covar,collapse = "+"),"+",exposure,"*",interaction)),
                                   design=survey_design)
          }
        }
      } else {
        if (is_factor==T){
          model<- survey::svyglm(as.formula(paste(outcome,"~",exposure)),
                                 design=survey_design)
        } else{
          model<- survey::svyglm(as.formula(paste(outcome,"~as.numeric(",exposure,")")),
                                 design=survey_design)
        }
        
        
      }
      
      Vcov <- vcov(model, useScale = FALSE)
      beta <- coef(model)
      se <- sqrt(diag(Vcov))
      zval <- beta / se
      pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
      noobs<-nobs(model)
      
      for (k in 2:length(beta)){
        out_beta[number]=as.numeric(beta[k])
        out_se[number] = as.numeric(se[k])
        out_pvalue[number] = as.numeric(pval[k])
        out_variable[number] = outcome
        out_nobs[number]=noobs[1]
        number = number + 1
        exp_beta[number] = as.numeric(beta[k])
        exp_se[number] = as.numeric(se[k])
        exp_pvalue[number] = as.numeric(pval[k])
        exp_variable[number]=names(beta)[k]
        exp_nobs[number]=noobs[1]
        exp_mdl[number]=exposure
        number = number + 1
      }
      # print(j)
    }
    # print(i)
  }
  outcome = data.frame(out_variable, out_beta, out_se, out_pvalue,out_nobs)
  exposure = data.frame(exp_variable, exp_beta, exp_se, exp_pvalue,exp_mdl,exp_nobs)
  outcome = outcome %>% 
    dplyr::rename(
      variable = out_variable,
      beta = out_beta,
      se = out_se,
      pvalue = out_pvalue,
      obs = out_nobs
    )
  exposure = exposure %>% 
    dplyr::rename(
      variable = exp_variable,
      beta = exp_beta,
      se = exp_se,
      pvalue = exp_pvalue,
      obs = exp_nobs,
      exposure_model=exp_mdl
    )
  all=merge(outcome,exposure,by=c("beta","se","pvalue","obs"))
  # all = rbind(outcome, exposure)
  all = na.omit(all)
  data = all %>%
    dplyr::rename(outcome=variable.x,
                  exposure=variable.y)%>%
    dplyr::mutate (
      beta = round(beta, 5),
      se = round(se, 5),
      pvalue = round(pvalue, 5)
    ) %>%
    dplyr::mutate(
      # or=exp(beta),
      #      lower95=exp(beta-1.96*se),
      #      upper95=exp(beta+1.96*se),
           flag_pvalue_high=ifelse(pvalue<0.05,1,0),
           flag_pvalue_low=ifelse(pvalue<0.2,1,0))%>%
    dplyr::select(outcome,exposure_model,exposure, beta, se, pvalue,flag_pvalue_high,flag_pvalue_low,obs)%>%
    arrange(outcome,exposure)
  return(data)
}

########################################
# elastic net regression
enet_function<-function(x,y,binary_outcome,var_not_penalized,standardized){
  penalize_vector<-data.frame(factor=rep(1,ncol(x)),predictor=colnames(x))
  penalize_vector$penalize<-ifelse(penalize_vector$predictor%in%var_not_penalized,0,1)
  if (binary_outcome==TRUE) {
    glmmod<-glmnet(x,y,alpha=0.5,family="binomial",penalty.factor = penalize_vector$penalize,standardize=standardized)
    lambdas_to_try <- 10^seq(-3, 5, length.out = 100)
    
    a <- seq(0.1, 0.9, 0.5)
    search <- foreach(i = a, .combine = rbind) %dopar% {
      set.seed(123)
      glmmod_cv <- cv.glmnet(x,y,type.measure="class", family = "binomial", nfold = 10, paralle = TRUE, alpha = i,penalty.factor = penalize_vector$penalize,standardize=standardized)
      data.frame(cvm = glmmod_cv$cvm[glmmod_cv$lambda == glmmod_cv$lambda.min], lambda.1se = glmmod_cv$lambda.1se,lambda.min = glmmod_cv$lambda.min, alpha = i)
    }
    cv3 <- search[search$cvm == min(search$cvm), ]
    glmmod_cv <- cv.glmnet(x,y,type.measure="class", family = "binomial", nfold = 10, paralle = TRUE, alpha = cv3$alpha,penalty.factor = penalize_vector$penalize,standardize=standardized)
    # md3 <- glmnet(x, y, family = "binomial", lambda = cv3$lambda.min, alpha = cv3$alpha,penalty.factor = penalize_vector$penalize,standardize=standardized)
    # Plot cross-validation results
    plot(glmmod_cv)
    plot_cv<-recordPlot()
    best_lambda<-glmmod_cv$lambda.min # find the lambda that minimizes the CV error
    mean_lambda<-glmmod_cv$cvm #cvm:mean cross-validated error
    se1_lambda<-glmmod_cv$lambda.1se
    plot(glmmod, xvar="lambda")
    abline(v = log(glmmod_cv$lambda.min), col = "red", lty = "dashed")
    abline(v = log(glmmod_cv$lambda.1se), col = "blue", lty = "dashed")
    plot_coef_lambda<-recordPlot()
    plot_influ<-coef(glmmod_cv, s = "lambda.min") %>%
      broom::tidy() %>%
      dplyr::filter(row != "(Intercept)") %>%
      ggplot(aes(value, reorder(row, value), color = value > 0)) +
      geom_point(show.legend = FALSE) +
      ggtitle("Influential variables") +
      xlab("Coefficient") +
      ylab(NULL) 
    coef<-coef(glmmod_cv, s = "lambda.min")
    coef_lasso<-coef[which(coef[,1]!=0),1]
    # glmmod_final<-glmnet(x,y,alpha=1,family="binomial",lambda = best_lambda,penalty.factor = penalize_vector$penalize)
    # output<-list(glmmod_final=glmmod_final,plot_coef_lambda=plot_coef_lambda,plot_cv=plot_cv, 
    #              best_lambda=best_lambda,mean_lambda=mean_lambda,se1_lambda=se1_lambda,coef_lasso=coef_lasso,plot_influ=plot_influ)
    output<-list(plot_coef_lambda=plot_coef_lambda,plot_cv=plot_cv, 
                 best_lambda=best_lambda,mean_lambda=mean_lambda,se1_lambda=se1_lambda,coef_lasso=coef_lasso,plot_influ=plot_influ)
  }else {
    glmmod<-glmnet(x,y,alpha=1,penalty.factor = penalize_vector$penalize,standardize=standardized)
    lambdas_to_try <- 10^seq(-3, 5, length.out = 100)
    a <- seq(0.1, 0.9, 0.5)
    search <- foreach(i = a, .combine = rbind) %dopar% {
      set.seed(123)
      glmmod_cv <- cv.glmnet(x,y,type.measure="mse", nfold = 10, paralle = TRUE, alpha = i,penalty.factor = penalize_vector$penalize,standardize=standardized)
      data.frame(cvm = glmmod_cv$cvm[glmmod_cv$lambda == glmmod_cv$lambda.min], lambda.1se = glmmod_cv$lambda.1se,lambda.min = glmmod_cv$lambda.min, alpha = i)
    }
    cv3 <- search[search$cvm == min(search$cvm), ]
    glmmod_cv <- cv.glmnet(x,y,type.measure="mse", nfold = 10, paralle = TRUE, alpha = cv3$alpha,penalty.factor = penalize_vector$penalize,standardize=standardized)
    
    # set.seed(123)
    # glmmod_cv <- cv.glmnet(x, y, alpha = 1, lambda = lambdas_to_try,nfolds = 10,penalty.factor = penalize_vector$penalize,standardize=standardized)
    # Plot cross-validation results
    plot(glmmod_cv)
    plot_cv<-recordPlot()
    best_lambda<-glmmod_cv$lambda.min # find the lambda that minimizes the CV error
    mean_lambda<-glmmod_cv$cvm #cvm:mean cross-validated error
    se1_lambda<-glmmod_cv$lambda.1se
    plot(glmmod, xvar="lambda")
    abline(v = log(glmmod_cv$lambda.min), col = "red", lty = "dashed")
    abline(v = log(glmmod_cv$lambda.1se), col = "blue", lty = "dashed")
    plot_coef_lambda<-recordPlot()
    plot_influ<-coef(glmmod_cv, s = "lambda.min") %>%
      broom::tidy() %>%
      dplyr::filter(row != "(Intercept)") %>%
      ggplot(aes(value, reorder(row, value), color = value > 0)) +
      geom_point(show.legend = FALSE) +
      ggtitle("Influential variables") +
      xlab("Coefficient") +
      ylab(NULL)+ theme(axis.text.y = element_text(size = 7))
    coef<-coef(glmmod_cv, s = "lambda.min")
    coef_lasso<-coef[which(coef[,1]!=0),1]
    # glmmod_final<-glmnet(x,y,alpha=1,lambda = best_lambda,penalty.factor = penalize_vector$penalize)
    # output<-list(glmmod_final=glmmod_final,plot_coef_lambda=plot_coef_lambda,plot_cv=plot_cv, best_lambda=best_lambda,mean_lambda=mean_lambda,
    # se1_lambda=se1_lambda,coef_lasso=coef_lasso,plot_influ=plot_influ)
    output<-list(plot_coef_lambda=plot_coef_lambda,plot_cv=plot_cv, best_lambda=best_lambda,mean_lambda=mean_lambda,se1_lambda=se1_lambda,
                 coef_lasso=coef_lasso,plot_influ=plot_influ)
    
  }
  return(output)
}
# LASSO Regression
lasso_function<-function(x,y,binary_outcome,var_not_penalized,standardized){
  # Note alpha=1 for lasso only and can blend with ridge penalty down to
  # alpha=0 ridge only.
  penalize_vector<-data.frame(factor=rep(1,ncol(x)),predictor=colnames(x))
  penalize_vector$penalize<-ifelse(penalize_vector$predictor%in%var_not_penalized,0,1)

  if (binary_outcome==TRUE) {
    glmmod<-glmnet(x,y,alpha=1,family="binomial",penalty.factor = penalize_vector$penalize,standardize=standardized)
    lambdas_to_try <- 10^seq(-3, 5, length.out = 100)
    set.seed(123)
    glmmod_cv <- cv.glmnet(x, y, alpha = 1, type.measure="auc",lambda = lambdas_to_try, nfolds = 10,family="binomial",penalty.factor = penalize_vector$penalize,standardize=standardized)
    # Plot cross-validation results
    plot(glmmod_cv)
    plot_cv<-recordPlot()
    best_lambda<-glmmod_cv$lambda.min # find the lambda that minimizes the CV error
    mean_lambda<-glmmod_cv$cvm #cvm:mean cross-validated error
    se1_lambda<-glmmod_cv$lambda.1se
    plot(glmmod, xvar="lambda")
    abline(v = log(glmmod_cv$lambda.min), col = "red", lty = "dashed")
    abline(v = log(glmmod_cv$lambda.1se), col = "blue", lty = "dashed")
    plot_coef_lambda<-recordPlot()
    plot_influ<-data.frame(as.matrix(coef(glmmod_cv, s = "lambda.min"))) %>%
      dplyr::mutate(row=row.names(.))%>%
      dplyr::rename(value="s1")%>%
      # broom::tidy() %>%
      dplyr::filter(row != "(Intercept)"&value!=0) %>%
      ggplot(aes(value, reorder(row, value), color = value > 0)) +
      geom_point(show.legend = FALSE) +
      ggtitle("Influential variables") +
      xlab("Coefficient") +
      ylab(NULL) 
    coef<-coef(glmmod_cv, s = "lambda.min")
    coef_lasso_all<-coef
    coef_lasso<-coef[which(coef[,1]!=0),1]
    if (length(coef_lasso)==(length(var_not_penalized)+1)){
      set.seed(1234)
      glmmod_cv <- cv.glmnet(x, y, alpha = 1, type.measure="auc",lambda = lambdas_to_try, nfolds = 10,family="binomial",penalty.factor = penalize_vector$penalize,standardize=standardized)
      # Plot cross-validation results
      plot(glmmod_cv)
      plot_cv<-recordPlot()
      best_lambda<-glmmod_cv$lambda.min # find the lambda that minimizes the CV error
      mean_lambda<-glmmod_cv$cvm #cvm:mean cross-validated error
      se1_lambda<-glmmod_cv$lambda.1se
      plot(glmmod, xvar="lambda")
      abline(v = log(glmmod_cv$lambda.min), col = "red", lty = "dashed")
      abline(v = log(glmmod_cv$lambda.1se), col = "blue", lty = "dashed")
      plot_coef_lambda<-recordPlot()
      plot_influ<-data.frame(as.matrix(coef(glmmod_cv, s = "lambda.min"))) %>%
        dplyr::mutate(row=row.names(.))%>%
        dplyr::rename(value="s1")%>%
        # broom::tidy() %>%
        dplyr::filter(row != "(Intercept)"&value!=0) %>%
        ggplot(aes(value, reorder(row, value), color = value > 0)) +
        geom_point(show.legend = FALSE) +
        ggtitle("Influential variables") +
        xlab("Coefficient") +
        ylab(NULL) 
      coef<-coef(glmmod_cv, s = "lambda.min")
      coef_lasso_all<-coef
      coef_lasso<-coef[which(coef[,1]!=0),1]
    }
    # glmmod_final<-glmnet(x,y,alpha=1,family="binomial",lambda = best_lambda,penalty.factor = penalize_vector$penalize)
    # output<-list(glmmod_final=glmmod_final,plot_coef_lambda=plot_coef_lambda,plot_cv=plot_cv, 
    #              best_lambda=best_lambda,mean_lambda=mean_lambda,se1_lambda=se1_lambda,coef_lasso=coef_lasso,plot_influ=plot_influ)
    output<-list(plot_coef_lambda=plot_coef_lambda,plot_cv=plot_cv, 
                 best_lambda=best_lambda,mean_lambda=mean_lambda,se1_lambda=se1_lambda,coef_lasso=coef_lasso,coef_lasso_all=coef_lasso_all,plot_influ=plot_influ)
  } else {
    glmmod<-glmnet(x,y,alpha=1,penalty.factor = penalize_vector$penalize,standardize=standardized)
    lambdas_to_try <- 10^seq(-3, 5, length.out = 100)
    set.seed(123)
    glmmod_cv <- cv.glmnet(x, y, alpha = 1, lambda = lambdas_to_try,nfolds = 10,penalty.factor = penalize_vector$penalize,standardize=standardized)
    # Plot cross-validation results
    plot(glmmod_cv)
    plot_cv<-recordPlot()
    best_lambda<-glmmod_cv$lambda.min # find the lambda that minimizes the CV error
    mean_lambda<-glmmod_cv$cvm #cvm:mean cross-validated error
    se1_lambda<-glmmod_cv$lambda.1se
    plot(glmmod, xvar="lambda")
    abline(v = log(glmmod_cv$lambda.min), col = "red", lty = "dashed")
    abline(v = log(glmmod_cv$lambda.1se), col = "blue", lty = "dashed")
    plot_coef_lambda<-recordPlot()
    plot_influ<-data.frame(as.matrix(coef(glmmod_cv, s = "lambda.min"))) %>%
      dplyr::mutate(row=row.names(.))%>%
      dplyr::rename(value="s1")%>%
      # broom::tidy() %>%
      dplyr::filter(row != "(Intercept)"&value!=0) %>%
      ggplot(aes(value, reorder(row, value), color = value > 0)) +
      geom_point(show.legend = FALSE) +
      ggtitle("Influential variables") +
      xlab("Coefficient") +
      ylab(NULL)+ theme(axis.text.y = element_text(size = 7))
    coef<-coef(glmmod_cv, s = "lambda.min")
    coef_lasso<-coef[which(coef[,1]!=0),1]
    coef_lasso_all<-coef
    if (length(coef_lasso)==(length(var_not_penalized)+1)){
      # set.seed(100)
      set.seed(1026)
      glmmod_cv <- cv.glmnet(x, y, alpha = 1, lambda = lambdas_to_try,nfolds = 10,penalty.factor = penalize_vector$penalize,standardize=standardized)
      # Plot cross-validation results
      plot(glmmod_cv)
      plot_cv<-recordPlot()
      best_lambda<-glmmod_cv$lambda.min # find the lambda that minimizes the CV error
      mean_lambda<-glmmod_cv$cvm #cvm:mean cross-validated error
      se1_lambda<-glmmod_cv$lambda.1se
      plot(glmmod, xvar="lambda")
      abline(v = log(glmmod_cv$lambda.min), col = "red", lty = "dashed")
      abline(v = log(glmmod_cv$lambda.1se), col = "blue", lty = "dashed")
      plot_coef_lambda<-recordPlot()
      plot_influ<-data.frame(as.matrix(coef(glmmod_cv, s = "lambda.min"))) %>%
        dplyr::mutate(row=row.names(.))%>%
        dplyr::rename(value="s1")%>%
        # broom::tidy() %>%
        dplyr::filter(row != "(Intercept)"&value!=0) %>%
        ggplot(aes(value, reorder(row, value), color = value > 0)) +
        geom_point(show.legend = FALSE) +
        ggtitle("Influential variables") +
        xlab("Coefficient") +
        ylab(NULL) 
      coef<-coef(glmmod_cv, s = "lambda.min")
      coef_lasso_all<-coef
      coef_lasso<-coef[which(coef[,1]!=0),1]
      if (length(coef_lasso)==(length(var_not_penalized)+1)){
        set.seed(10000)
        glmmod_cv <- cv.glmnet(x, y, alpha = 1, lambda = lambdas_to_try,nfolds = 10,penalty.factor = penalize_vector$penalize,standardize=standardized)
        # Plot cross-validation results
        plot(glmmod_cv)
        plot_cv<-recordPlot()
        best_lambda<-glmmod_cv$lambda.min # find the lambda that minimizes the CV error
        mean_lambda<-glmmod_cv$cvm #cvm:mean cross-validated error
        se1_lambda<-glmmod_cv$lambda.1se
        plot(glmmod, xvar="lambda")
        abline(v = log(glmmod_cv$lambda.min), col = "red", lty = "dashed")
        abline(v = log(glmmod_cv$lambda.1se), col = "blue", lty = "dashed")
        plot_coef_lambda<-recordPlot()
        plot_influ<-data.frame(as.matrix(coef(glmmod_cv, s = "lambda.min"))) %>%
          dplyr::mutate(row=row.names(.))%>%
          dplyr::rename(value="s1")%>%
          # broom::tidy() %>%
          dplyr::filter(row != "(Intercept)"&value!=0) %>%
          ggplot(aes(value, reorder(row, value), color = value > 0)) +
          geom_point(show.legend = FALSE) +
          ggtitle("Influential variables") +
          xlab("Coefficient") +
          ylab(NULL) 
        coef<-coef(glmmod_cv, s = "lambda.min")
        coef_lasso_all<-coef
        coef_lasso<-coef[which(coef[,1]!=0),1]
        if (length(coef_lasso)==(length(var_not_penalized)+1)){
          set.seed(100)
          glmmod_cv <- cv.glmnet(x, y, alpha = 1, lambda = lambdas_to_try,nfolds = 10,penalty.factor = penalize_vector$penalize,standardize=standardized)
          # Plot cross-validation results
          plot(glmmod_cv)
          plot_cv<-recordPlot()
          best_lambda<-glmmod_cv$lambda.min # find the lambda that minimizes the CV error
          mean_lambda<-glmmod_cv$cvm #cvm:mean cross-validated error
          se1_lambda<-glmmod_cv$lambda.1se
          plot(glmmod, xvar="lambda")
          abline(v = log(glmmod_cv$lambda.min), col = "red", lty = "dashed")
          abline(v = log(glmmod_cv$lambda.1se), col = "blue", lty = "dashed")
          plot_coef_lambda<-recordPlot()
          plot_influ<-data.frame(as.matrix(coef(glmmod_cv, s = "lambda.min"))) %>%
            dplyr::mutate(row=row.names(.))%>%
            dplyr::rename(value="s1")%>%
            # broom::tidy() %>%
            dplyr::filter(row != "(Intercept)"&value!=0) %>%
            ggplot(aes(value, reorder(row, value), color = value > 0)) +
            geom_point(show.legend = FALSE) +
            ggtitle("Influential variables") +
            xlab("Coefficient") +
            ylab(NULL) 
          coef<-coef(glmmod_cv, s = "lambda.min")
          coef_lasso_all<-coef
          coef_lasso<-coef[which(coef[,1]!=0),1]
        }
      }
    }
    # glmmod_final<-glmnet(x,y,alpha=1,lambda = best_lambda,penalty.factor = penalize_vector$penalize)
    # output<-list(glmmod_final=glmmod_final,plot_coef_lambda=plot_coef_lambda,plot_cv=plot_cv, best_lambda=best_lambda,mean_lambda=mean_lambda,
    # se1_lambda=se1_lambda,coef_lasso=coef_lasso,plot_influ=plot_influ)
    output<-list(plot_coef_lambda=plot_coef_lambda,plot_cv=plot_cv, best_lambda=best_lambda,mean_lambda=mean_lambda,se1_lambda=se1_lambda,
                 coef_lasso=coef_lasso,coef_lasso_all=coef_lasso_all,plot_influ=plot_influ)
    
  }
  return(output)
}

# LASSO regression performance evaluation: continous outcome
# Compute R^2 from true and predicted values
eval_results <- function(true, predicted, df) {
  SSE <- sum((predicted - true)^2)
  SST <- sum((true - mean(true))^2)
  R_square <- 1 - SSE / SST
  RMSE = sqrt(SSE/nrow(df))
  
  # Model performance metrics
  data.frame(
    RMSE = RMSE,
    Rsquare = R_square
  )
}

# Permutation function to assess false positive rate
permute_resids<-function(data,covar,dep_list){
  data.residuals<-data[,c("ID","LAB_ID")]
  survey_design=svydesign(id=~PSU_ID, strata=~STRAT, weights=~WEIGHT, data=data)
  for (i in seq_along(dep_list)) {
    # print(i)
    output<-data.frame()
    comp_id<-data[complete.cases(data[,dep_list[i]]),"ID"]
    model<- svyglm(as.formula(paste(dep_list[i],"~",paste(covar,collapse= "+"))),design=survey_design)
    resid<-model$residuals
    resid_permuted <- sample(resid)
    value_permuted<-resid_permuted + fitted.values(model)	
    output<-cbind(comp_id,value_permuted)
    colnames(output)<-c("ID",paste(dep_list[i]))
    data.residuals<-merge(x=data.residuals, y=output, by="ID", all = TRUE)
    # print(i)
  }
  return(data.residuals)
}



# modified loop regression models using sampling weight and bonferroni correction

  fun_svyglm_multi_correction<-function(dat,out_start,out_end,exp_start,exp_end,covar,is_factor,interaction,n_compare){
    # define the Bonferronic corrected p value cutoff
    corrected_p<-0.05/n_compare
    # define the column number for the function
    # outcome
    exp_end=ncol(dat)
    out_nvar=out_end-out_start+1
    
    out_variable=rep(NA,out_nvar)
    out_beta=rep(NA,out_nvar)
    out_se=rep(NA,out_nvar)
    out_pvalue=rep(NA,out_nvar)
    out_nobs=rep(NA,out_nvar)
    # exposure
    exp_nvar=exp_end-exp_start+1
    
    exp_mdl=rep(NA,exp_nvar)
    exp_variable=rep(NA,exp_nvar)
    exp_beta=rep(NA,exp_nvar)
    exp_se=rep(NA,exp_nvar)
    exp_pvalue=rep(NA,exp_nvar)
    exp_nobs=rep(NA,out_nvar)
    number=1
    
    # define the survey design
    survey_design=survey::svydesign(id=~PSU_ID, strata=~STRAT, weights=~WEIGHT,data=dat)
    for (i in out_start:out_end){
      # print(i)
      outcome=colnames(dat)[i]
      for (j in exp_start:exp_end){
        # print(j)
        exposure=colnames(dat)[j]
        if (is.null(interaction)){
          if (is_factor==T){
            model<- survey::svyglm(as.formula(paste(outcome,"~",exposure,"+",paste(covar,collapse = "+"))),
                                   design=survey_design,family=quasipoisson(link='log'))
          }else{
            model<- survey::svyglm(as.formula(paste(outcome,"~as.numeric(",exposure,")","+",paste(covar,collapse = "+"))),
                                   design=survey_design,family=quasipoisson(link='log'))
          }
        }else {
          if (is_factor==T){
            model<- survey::svyglm(as.formula(paste(outcome,"~",exposure,"+",paste(covar,collapse = "+"),"+",exposure,"*",interaction)),
                                   design=survey_design,family=quasipoisson(link='log'))
          }else{
            model<- survey::svyglm(as.formula(paste(outcome,"~as.numeric(",exposure,")","+",paste(covar,collapse = "+"),"+",exposure,"*",interaction)),
                                   design=survey_design,family=quasipoisson(link='log'))
          }
        }
        
        Vcov <- vcov(model, useScale = FALSE)
        beta <- coef(model)
        se <- sqrt(diag(Vcov))
        zval <- beta / se
        pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
        noobs<-nobs(model)
        
        for (k in 2:length(beta)){
          out_beta[number]=as.numeric(beta[k])
          out_se[number] = as.numeric(se[k])
          out_pvalue[number] = as.numeric(pval[k])
          out_variable[number] = outcome
          out_nobs[number]=noobs[1]
          number = number + 1
          exp_beta[number] = as.numeric(beta[k])
          exp_se[number] = as.numeric(se[k])
          exp_pvalue[number] = as.numeric(pval[k])
          exp_variable[number]=names(beta)[k]
          exp_nobs[number]=noobs[1]
          exp_mdl[number]=exposure
          number = number + 1
        }
        # print(j)
      }
      # print(i)
    }
    outcome = data.frame(out_variable, out_beta, out_se, out_pvalue,out_nobs)
    exposure = data.frame(exp_variable, exp_beta, exp_se, exp_pvalue,exp_mdl,exp_nobs)
    outcome = outcome %>% 
      dplyr::rename(
        variable = out_variable,
        beta = out_beta,
        se = out_se,
        pvalue = out_pvalue,
        obs = out_nobs
      )
    exposure = exposure %>% 
      dplyr::rename(
        variable = exp_variable,
        beta = exp_beta,
        se = exp_se,
        pvalue = exp_pvalue,
        obs = exp_nobs,
        exposure_model=exp_mdl
      )
    all=merge(outcome,exposure,by=c("beta","se","pvalue","obs"))
    # all = rbind(outcome, exposure)
    all = na.omit(all)
    data = all %>%
      dplyr::rename(outcome=variable.x,
                    exposure=variable.y)%>%
      mutate (
        beta = round(beta, 5),
        se = round(se, 5),
        pvalue = round(pvalue, 5)
      ) %>%
      mutate(or=exp(beta),
             lower95=exp(beta-1.96*se),
             upper95=exp(beta+1.96*se),
             flag_pvalue_high=ifelse(pvalue<corrected_p,1,0),
             flag_pvalue_low=ifelse(pvalue<0.05,1,0))%>%
      dplyr::select(outcome,exposure_model,exposure, beta, se, or,lower95,upper95,pvalue,flag_pvalue_high,flag_pvalue_low,obs)%>%
      arrange(outcome,exposure)
    return(data)
  }
  
  ################################################
  # The p-filter: multilayer FDR control for grouped hypotheses 
  pfilter = function(P,alphas,group){
    # P in [0,1]^n = vector of p-values
    # alphas in [0,1]^M = vector of target FDR levels
    # groups is a n-by-M matrix; 
    #	groups[i,m] = which group does P[i] belong to,
    #		for the m-th grouping
    
    # change groups to a matrix of col_1: number of rows, col_2: convert trait and metabolties to numeric value
    if (is.null(dim(group))) {
      groups<-as.matrix(data.frame(n_row=c(1:length(group)),n_group=as.numeric(as.factor(group))))
      alphas<-alphas[1:2]
    } else {
      groups<-as.matrix(data.frame(n_row=c(1:nrow(group)),n_group1=as.numeric(as.factor(group[,1])),n_group2=as.numeric(as.factor(group[,2]))))
    }
    
    n = length(P)
    M = length(alphas)
    G = apply(groups,2,max) # G[m] = # groups, for grouping m
    
    
    Simes = list()
    for(m in 1:M){
      Simes[[m]]=rep(0,G[m])
      for(g in 1:G[m]){
        group = which(groups[,m]==g)
        Simes[[m]][g]=min(sort(P[group])*length(group)/(1:length(group)))
      }
    }
    
    
    # initialize
    thresh = alphas
    Sh = 1:n
    for(m in 1:M){
      pass_Simes_m = which(is.element(groups[,m],which(Simes[[m]]<=thresh[m])))
      Sh = intersect(Sh,pass_Simes_m)
    }
    done = FALSE
    
    
    while(!done){
      thresh_old = thresh
      for(m in 1:M){
        # which groups, for the m-th grouping, 
        #	have any potential discoveries?
        Shm = sort(unique(groups[Sh,m]))
        
        # run BH on Simes[[m]], constraining to Shm
        Pvals_m = rep(1.01,G[m]); # >1 for groups not in Dm
        Pvals_m[Shm] = Simes[[m]][Shm]
        khatm = max(0,which(sort(Pvals_m)<=(1:G[m])/G[m]*alphas[m]))
        thresh[m] = khatm/G[m]*alphas[m]
        Sh = intersect(Sh,
                       which(is.element(groups[,m],which(Pvals_m<=thresh[m]))))
      }
      if(all(thresh_old==thresh)){done = TRUE}
    }
    
    Sh_temp = Sh;
    Sh = rep(0,n); Sh[Sh_temp] = 1
    Sh
    
  }
  
  # function to loop pfilter over multiple models and strata
  pfilter_loop<-function(data, model_input, strata_input,alpha_threshold){
    output_data<-data.frame()
    for (i in seq_along(model_input)){
      strata_output<-data.frame()
      for (j in seq_along(strata_input)){
        temp_data<-data%>%
          dplyr::filter(model==model_input[i],strata==strata_input[j])
        temp_data$pfilter_1<-pfilter(P=temp_data$p_val,alphas=alpha_threshold,group=temp_data$trait)
        temp_data$pfilter_2<-pfilter(P=temp_data$p_val,alphas=alpha_threshold,group=temp_data[,c("trait","metabolite")])
        strata_output<-rbind(strata_output,temp_data)
      }
      output_data<-rbind(output_data,strata_output)
    }
    return(output_data)
  }
  
  #################################################
  # modified loop regression models using sampling weight and pfilter correction
  
  fun_svyglm_multi_pfilter<-function(dat,out_start,out_end,exp_start,exp_end,covar,is_factor,interaction,trait_list,alpha_threshold){
    
    # define the column number for the function
    # outcome
    exp_end=ncol(dat)
    out_nvar=out_end-out_start+1
    
    out_variable=rep(NA,out_nvar)
    out_beta=rep(NA,out_nvar)
    out_se=rep(NA,out_nvar)
    out_pvalue=rep(NA,out_nvar)
    out_nobs=rep(NA,out_nvar)
    # exposure
    exp_nvar=exp_end-exp_start+1
    
    exp_mdl=rep(NA,exp_nvar)
    exp_variable=rep(NA,exp_nvar)
    exp_beta=rep(NA,exp_nvar)
    exp_se=rep(NA,exp_nvar)
    exp_pvalue=rep(NA,exp_nvar)
    exp_nobs=rep(NA,out_nvar)
    number=1
    
    # define the survey design
    survey_design=survey::svydesign(id=~PSU_ID, strata=~STRAT, weights=~WEIGHT,data=dat)
    for (i in out_start:out_end){
      # print(i)
      outcome=colnames(dat)[i]
      for (j in exp_start:exp_end){
        # print(j)
        exposure=colnames(dat)[j]
        if (is.null(interaction)){
          if (is_factor==T){
            model<- survey::svyglm(as.formula(paste(outcome,"~",exposure,"+",paste(covar,collapse = "+"))),
                                   design=survey_design,family=quasipoisson(link='log'))
          }else{
            model<- survey::svyglm(as.formula(paste(outcome,"~as.numeric(",exposure,")","+",paste(covar,collapse = "+"))),
                                   design=survey_design,family=quasipoisson(link='log'))
          }
        }else {
          if (is_factor==T){
            model<- survey::svyglm(as.formula(paste(outcome,"~",exposure,"+",paste(covar,collapse = "+"),"+",exposure,"*",interaction)),
                                   design=survey_design,family=quasipoisson(link='log'))
          }else{
            model<- survey::svyglm(as.formula(paste(outcome,"~as.numeric(",exposure,")","+",paste(covar,collapse = "+"),"+",exposure,"*",interaction)),
                                   design=survey_design,family=quasipoisson(link='log'))
          }
        }
        
        Vcov <- vcov(model, useScale = FALSE)
        beta <- coef(model)
        se <- sqrt(diag(Vcov))
        zval <- beta / se
        pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
        noobs<-nobs(model)
        
        for (k in 2:length(beta)){
          out_beta[number]=as.numeric(beta[k])
          out_se[number] = as.numeric(se[k])
          out_pvalue[number] = as.numeric(pval[k])
          out_variable[number] = outcome
          out_nobs[number]=noobs[1]
          number = number + 1
          exp_beta[number] = as.numeric(beta[k])
          exp_se[number] = as.numeric(se[k])
          exp_pvalue[number] = as.numeric(pval[k])
          exp_variable[number]=names(beta)[k]
          exp_nobs[number]=noobs[1]
          exp_mdl[number]=exposure
          number = number + 1
        }
        # print(j)
      }
      # print(i)
    }
    outcome = data.frame(out_variable, out_beta, out_se, out_pvalue,out_nobs)
    exposure = data.frame(exp_variable, exp_beta, exp_se, exp_pvalue,exp_mdl,exp_nobs)
    outcome = outcome %>% 
      dplyr::rename(
        variable = out_variable,
        beta = out_beta,
        se = out_se,
        pvalue = out_pvalue,
        obs = out_nobs
      )
    exposure = exposure %>% 
      dplyr::rename(
        variable = exp_variable,
        beta = exp_beta,
        se = exp_se,
        pvalue = exp_pvalue,
        obs = exp_nobs,
        exposure_model=exp_mdl
      )
    all=merge(outcome,exposure,by=c("beta","se","pvalue","obs"))
    # all = rbind(outcome, exposure)
    all = na.omit(all)
    data = all %>%
      dplyr::rename(outcome=variable.x,
                    exposure=variable.y)%>%
      mutate (
        beta = round(beta, 5),
        se = round(se, 5)
        # pvalue = round(pvalue, 5)
      ) %>%
      mutate(or=exp(beta),
             lower95=exp(beta-1.96*se),
             upper95=exp(beta+1.96*se)
             # flag_pvalue_high=ifelse(pvalue<corrected_p,1,0),
             # flag_pvalue_low=ifelse(pvalue<0.05,1,0)
             )
    # divide the output into covariates output and trait output
    data_cov<-data[which(!data$exposure%in%trait_list),]%>%
      dplyr::mutate(pfilter_1=NA,
             pfilter_2=NA)
    data_pfilter<-data[which(data$exposure%in%trait_list),]
    # apply pfilter on the trait output
    data_pfilter$pfilter_1<-pfilter(P=data_pfilter$pvalue,alphas=alpha_threshold,group=droplevels(data_pfilter$exposure))
    data_pfilter$pfilter_2<-pfilter(P=data_pfilter$pvalue,alphas=alpha_threshold,group=droplevels(data_pfilter[,c("exposure","outcome")]))
    data<-rbind(data_cov,data_pfilter)%>%
      dplyr::select(outcome,exposure_model,exposure, beta, se, or,lower95,upper95,pvalue,pfilter_1,pfilter_2,obs)%>%
      arrange(outcome,exposure)
    return(data)
  }
  
  
  #######################################
  # Functions for construct trait specific MRS (different metabolites for each trait)
  # construct a list of vectors, name of which is the sleep trait
  list_MetabByTrait<-function(dat,colname_by,colname_vector,filter_by,filter_val){
    data<-dat[which(dat[,filter_by]==filter_val),]
    list_vector<-list()
    n_trait<-length(unique(data[,colname_by]))
    for (i in 1:n_trait){
      trait_by<-unique(data[,colname_by])[i]
      list_vector[[i]]<-unique(data[which(data[,colname_by]==trait_by),colname_vector])
    }
    names(list_vector)<-unique(data[,colname_by])
    return(list_vector)
  }
  
  # # Loop MRS analysis (use stratifier)
  # loop_mrs<-function(metab_data,phenotype_data,sleep_trait,metab_list,covar,outcome_var,covariates_cog,admin_var,exclude_var,alpha_threshold_strict,id_bin,id_cont,stratifier){
    # mrs_results<-data.frame()
    # mrs_scores<-data.frame()
    # for (i in seq_along(sleep_trait)){
    #   # print(paste("i=",i))
    #   
    #   # add ECGA7 to the covariate list if sleep trait is hr_dim1 or hr_dim2
    #   if (sleep_trait[i]%in%c("hr_pc1","hr_pc2","hr_dim1","hr_dim2","avg_hr","std_hr","min_hr","max_hr")){
    #     covar_res<-c(covar,"ECGA7")
    #   } else covar_res<-covar
    #   # Regress the rank-normalized metabolites and sleep trait of interest over covariates, obtain residuals.
    #   res_SDm<-residual_val(met_data=metab_data,pheno_data=phenotype_data,dep_list=c(sleep_trait[i],metab_list[[i]]),cov=covar_res)
    #   # Obtain SDr of residuals (including traits)
    #   res_sdr_inv_SDm<-diag(1/sqrt(diag(cov(res_SDm[,!colnames(res_SDm)%in%c("ID","LAB_ID")],use="pairwise.complete.obs"))))
    #   # obtain the ID for the complete obs
    #   res_sdr_inv_SDm_completeID<-res_SDm[complete.cases(res_SDm[,!colnames(res_SDm)%in%c("ID","LAB_ID")]),"ID"]
    #   # replace NAs with 0, multiply with the SDm and then revert back to NAs
    #   res_0_SDm<- res_SDm[,!colnames(res_SDm)%in%c("ID","LAB_ID")]
    #   res_0_SDm[is.na(res_SDm[,!colnames(res_SDm)%in%c("ID","LAB_ID")])] <- 0
    #   # Divide metabolite residuals by SDr for CCA (because with CCA we want variables to have the same variance). 
    #   res_rescaled_SDm <- as.matrix(res_0_SDm) %*% res_sdr_inv_SDm
    #   res_rescaled_SDm[is.na(res_SDm[,!colnames(res_SDm)%in%c("ID","LAB_ID")])] <- NA
    #   colnames(res_rescaled_SDm)<-colnames(res_SDm[,!colnames(res_SDm)%in%c("ID","LAB_ID")])
    #   res_rescaled_SDm<-cbind(res_SDm[,c("ID","LAB_ID")],res_rescaled_SDm)
    #   
    #   # res_rescaled_SDm_trait<-res_rescaled_SDm[which(!is.na(res_rescaled_SDm[,sleep_trait[i]])),] # drop obs with missing sleep trait for rcc
    #   res_rescaled_SDm_trait<-res_rescaled_SDm[complete.cases(res_rescaled_SDm),] # drop any missing obs
    #   
    #   # Perform CCA. Obtain coefficients: every trait (single column) and all the metabolites
      # metab_matrix_SDm<-res_rescaled_SDm[,metab_list[[i]]]
      # trait_matrix_SDm<-res_rescaled_SDm[,sleep_trait[i]]
      # # non-missing value index for the sleep trait
      # trait_matrix_SDm_rowindex<-which(!is.na(trait_matrix_SDm))
      # # non-missing value index for the metabolites
      # metab_matrix_SDm_rowindex<-rownames(metab_matrix_SDm[which(rowSums(is.na(metab_matrix_SDm))==0),])
      # # row index for non-missing value for both matrix
      # combined_rowindex<-as.numeric(metab_matrix_SDm_rowindex[metab_matrix_SDm_rowindex%in%trait_matrix_SDm_rowindex])
      # 
      # # subset sdr for only the metabolites
      # if (length(metab_list[[i]])==1){
      #   res_metab_sdr_SDm<-sd(res_SDm[,colnames(res_SDm)%in%metab_list[[i]]], na.rm=T)
      # } else {
      #   res_metab_sdr_SDm<-sqrt(diag(cov(res_SDm[,colnames(res_SDm)%in%metab_list[[i]]],use="pairwise.complete.obs")))
      # }
      # # Perform CCA and multiply coefficients by SDr to go back to pre-CCA scale.(If only one metabolites then no need to multiply SD)
      #   # regul.par <- CCA::estim.regul(as.matrix(metab_matrix_SDm), as.matrix(trait_matrix_SDm), plt = T) # l=50
      # # cc_output_SDm<-CCA::rcc(as.matrix(metab_matrix_SDm),as.matrix(trait_matrix_SDm),regul.par$lambda1, regul.par$lambda2)
      # # cc_output_SDm<-mixOmics::rcc(as.matrix(metab_matrix_SDm),as.matrix(trait_matrix_SDm), ncomp = 3, method = 'shrinkage')
      # 
      # # basic CCA that takes y as a vector, missing value included
      # # cc_output_SDm<-CCA::cc(as.matrix(metab_matrix_SDm[trait_matrix_SDm_rowindex,]),as.matrix(trait_matrix_SDm[trait_matrix_SDm_rowindex])) # basic CCA to get beta_s
      # # basic CCA that takes y as a vector, missing value from either matrix excluded
      # cc_output_SDm<-CCA::cc(as.matrix(metab_matrix_SDm[combined_rowindex,]),as.matrix(trait_matrix_SDm[combined_rowindex])) # basic CCA to get beta_s
      # 
      # 
      # # cc_output_SDm<-CCA::rcc(as.matrix(metab_matrix_SDm),as.matrix(trait_matrix_SDm),0.1,0.2) # basic CCA to get beta_s
      # metab_coeff_SDm<-cc_output_SDm$xcoef/res_metab_sdr_SDm # SD of the residual metabolites
      # 
      # # comput the canonical correlation coefficient and other parameters
      # # metab_cc_corr<-comput(as.matrix(metab_matrix_SDm[combined_rowindex,]),as.matrix(trait_matrix_SDm[combined_rowindex]), cc_output_SDm)
      # 
      # # metab_coeff_SDm<-data.frame(coefficient=cc_output_SDm$xcoef/res_metab_sdr_SDm,trait=sleep_trait[i]) # SD of the residual metabolites
      # # metab_coeff_SDm$trait<-sleep_trait[i]
      # 
      # # Construct the score using rank-normalized metabolites and the coefficients from the previous step
      # x_SDm<-t(t(metab_data[,metab_list[[i]]])*as.vector(metab_coeff_SDm))
      # y_SDm<-rowSums(x_SDm,na.rm = TRUE)
      # # standardize the score
      # score_SDm<-as.vector(scale(y_SDm))
      # # score_SDm<-cbind(score_SDm,res_SDm$ID)
      # score_SDm<-cbind(score_SDm,metab_data$ID)
      # colnames(score_SDm)<-c(sleep_trait[i],"ID")
      # # output the score
      # if (i==1){
      #   mrs_scores<-score_SDm
      # } else {
      #   mrs_scores<-merge(mrs_scores,score_SDm,by="ID")
      # }
      # 
      # # merge with the phenotype data
      # phenotype_data$WEIGHT<-phenotype_data$WEIGHT_NORM_OVERALL_V2
      # pheno_mrs_SDm<-merge(phenotype_data[,c(outcome_var,covar_res,covariates_cog,admin_var,exclude_var,stratifier)],score_SDm,by="ID",all.x=T)
      # # pheno_mrs_SDm_ID<-pheno_mrs_SDm[which(!is.na(pheno_mrs_SDm$WEIGHT)&pheno_mrs_SDm$DIABETES3_INDICATOR!="Yes"&pheno_mrs_SDm$HYPERTENSION!="Yes"),]
      # pheno_mrs_SDm_ID<-pheno_mrs_SDm[which(!is.na(pheno_mrs_SDm$WEIGHT)),]
      # 
      # # run the association analysis
      #   out_nvar=length(outcome_var)
      #   out_beta=rep(NA,out_nvar)
      #   out_se=rep(NA,out_nvar)
      #   out_pvalue=rep(NA,out_nvar)
      #   out_noobs=rep(NA,out_nvar)     
      #   
        # out_strat1_nvar=length(outcome_var)
        # out_strat1_beta=rep(NA,out_nvar)
        # out_strat1_se=rep(NA,out_nvar)
        # out_strat1_pvalue=rep(NA,out_nvar)
        # out_strat1_noobs=rep(NA,out_nvar)
        # out_strat2_nvar=length(outcome_var)
        # out_strat2_beta=rep(NA,out_nvar)
        # out_strat2_se=rep(NA,out_nvar)
        # out_strat2_pvalue=rep(NA,out_nvar)
        # out_strat2_noobs=rep(NA,out_nvar)
        # out_strat3_nvar=length(outcome_var)
        # out_strat3_beta=rep(NA,out_nvar)
        # out_strat3_se=rep(NA,out_nvar)
        # out_strat3_pvalue=rep(NA,out_nvar)
        # out_strat3_noobs=rep(NA,out_nvar)
        # number=1
      # results=data.frame()
      # for (l in seq_along(outcome_var)){
      #   # add time elapse between baseline and follow up for MCI regression
      #   if (outcome_var[l]%in%c("mci")){
      #     cov_mrs<-c(cov_mrs,covariates_cog)
      #   } else {
      #     cov_mrs<-covar_res
      #   }
      #   # print(paste("l",l))
      #   if (!is.na(stratifier)){
      #     model_strat<-list()
      #     model<-by_svyglm_multi(data=pheno_mrs_SDm,outcome=outcome_var[l],exposure=sleep_trait[i],covar=cov_mrs,include_ID=pheno_mrs_SDm_ID$ID)
      #     strat_value<-unique(pheno_mrs_SDm[,stratifier])
      #     for (k in seq_along(strat_value[!is.na(strat_value)])){
      #       model_strat[[k]]<-by_svyglm_multi(data=pheno_mrs_SDm[which(pheno_mrs_SDm[,stratifier]==unique(pheno_mrs_SDm[,stratifier])[k]),],outcome=outcome_var[l],exposure=sleep_trait[i],covar=cov_mrs[!cov_mrs%in%c(stratifier,"GENDER")],include_ID=pheno_mrs_SDm_ID$ID)
      #     }
      #     # model_strat1<-by_svyglm_multi(data=subset(pheno_mrs_SDm,GENDER=="Female"),outcome=outcome_var[l],exposure=sleep_trait[i],covar=cov_mrs[!cov_mrs%in%c("GENDER")],include_ID=pheno_mrs_SDm_ID$ID)
      #     # model_male<-by_svyglm_multi(data=subset(pheno_mrs_SDm,GENDER=="Male"),outcome=outcome_var[l],exposure=sleep_trait[i],covar=cov_mrs[!cov_mrs%in%c("GENDER")],include_ID=pheno_mrs_SDm_ID$ID)
      #     Vcov <- vcov(model, useScale = FALSE)
      #     beta<- coef(model)
          # se<- sqrt(diag(vcov(model, useScale = FALSE)))
          # zval<- beta / se
          # pval<- 2 * pnorm(abs(zval), lower.tail = FALSE)
          # out_beta[number]=as.numeric(beta[2])
          # out_se[number] = as.numeric(se[2])
          # out_pvalue[number] = as.numeric(pval[2])
          # out_noobs[number]=nobs(model)
          # 
          # Vcov_strat1 <- vcov(model_strat[[1]], useScale = FALSE)
          # beta_strat1<- coef(model_strat[[1]])
          # se_strat1<- sqrt(diag(vcov(model_strat[[1]], useScale = FALSE)))
          # zval_strat1<- beta_strat1 / se_strat1
          # pval_strat1<- 2 * pnorm(abs(zval_strat1), lower.tail = FALSE)
          # out_strat1_beta[number]=as.numeric(beta_strat1[2])
          # out_strat1_se[number] = as.numeric(se_strat1[2])
          # out_strat1_pvalue[number] = as.numeric(pval_strat1[2])
          # out_strat1_noobs[number]=nobs(model_strat[[1]])
          # 
          # Vcov_strat2 <- vcov(model_strat[[2]], useScale = FALSE)
          # beta_strat2<- coef(model_strat[[2]])
          # se_strat2<- sqrt(diag(vcov(model_strat[[2]], useScale = FALSE)))
          # zval_strat2<- beta_strat2 / se_strat2
          # pval_strat2<- 2 * pnorm(abs(zval_strat2), lower.tail = FALSE)
          # out_strat2_beta[number]=as.numeric(beta_strat2[2])
          # out_strat2_se[number] = as.numeric(se_strat2[2])
          # out_strat2_pvalue[number] = as.numeric(pval_strat2[2])
          # out_strat2_noobs[number]=nobs(model_strat[[2]]) 
          # 
          # if (length(unique(pheno_mrs_SDm[,stratifier]))==3){
          #   Vcov_strat3 <- vcov(model_strat[[3]], useScale = FALSE)
          #   beta_strat3<- coef(model_strat[[3]])
          #   se_strat2<- sqrt(diag(vcov(model_strat[[3]], useScale = FALSE)))
          #   zval_strat3<- beta_strat3 / se_strat3
          #   pval_strat3<- 2 * pnorm(abs(zval_strat3), lower.tail = FALSE)
        #     out_strat3_beta[number]=as.numeric(beta_strat3[2])
        #     out_strat3_se[number] = as.numeric(se_strat3[2])
        #     out_strat3_pvalue[number] = as.numeric(pval_strat3[2])
        #     out_strat3_noobs[number]=nobs(model_strat[[3]]) 
        #   }
        #   number=number+1
        #   
        # } else {
        #   model<-by_svyglm_multi(data=pheno_mrs_SDm,outcome=outcome_var[l],exposure=sleep_trait[i],covar=cov_mrs,include_ID=pheno_mrs_SDm_ID$ID)
        # Vcov <- vcov(model, useScale = FALSE)
        # beta<- coef(model)
        # se<- sqrt(diag(vcov(model, useScale = FALSE)))
      #   zval<- beta / se
      #   pval<- 2 * pnorm(abs(zval), lower.tail = FALSE)
      #   out_beta[number]=as.numeric(beta[2])
      #   out_se[number] = as.numeric(se[2])
      #   out_pvalue[number] = as.numeric(pval[2])
      #   out_noobs[number]=nobs(model)
      #   number=number+1
      #   # print(paste("l=",l))
      #   }
      # }
      # if (!is.na(stratifier)){
      #   results_both<-data.frame(outcome=outcome_var,
      #                       beta=out_beta,
      #                       se=out_se,
      #                       or=exp(out_beta),
      #                       lower95=exp(out_beta-1.96*out_se),
      #                       upper95=exp(out_beta+1.96*out_se),
      #                       p_val=out_pvalue,
      #                       n=out_noobs,
      #                       sleep_traits=sleep_trait[i],
        #                     covar=paste(cov_mrs,collapse=","),
        #                     strata="Both"
        # )
        # results_strat1<-data.frame(outcome=outcome_var,
        #                     beta=out_strat1_beta,
        #                     se=out_strat1_se,
        #                     or=exp(out_strat1_beta),
        #                     lower95=exp(out_strat1_beta-1.96*out_strat1_se),
        #                     upper95=exp(out_strat1_beta+1.96*out_strat1_se),
        #                     p_val=out_strat1_pvalue,
        #                     n=out_strat1_noobs,
        #                     sleep_traits=sleep_trait[i],
        #                     covar=paste(cov_mrs,collapse=","),
        #                     strata=paste(strat_value[1])
        # )
        # results_strat2<-data.frame(outcome=outcome_var,
        #                            beta=out_strat2_beta,
        #                            se=out_strat2_se,
        #                            or=exp(out_strat2_beta),
        #                            lower95=exp(out_strat2_beta-1.96*out_strat2_se),
        #                            upper95=exp(out_strat2_beta+1.96*out_strat2_se),
        #                            p_val=out_strat2_pvalue,
        #                            n=out_strat2_noobs,
        #                            sleep_traits=sleep_trait[i],
        #                            covar=paste(cov_mrs,collapse=","),
        #                            strata=paste(strat_value[2])
        # )
        # results<-rbind(results_both,results_strat1,results_strat2)
        # if (length(unique(pheno_mrs_SDm[,stratifier]))==3){
        #   results_strat3<-data.frame(outcome=outcome_var,
        #                              beta=out_strat3_beta,
        #                              se=out_strat3_se,
        #                              or=exp(out_strat3_beta),
        #                              lower95=exp(out_strat3_beta-1.96*out_strat3_se),
        #                              upper95=exp(out_strat3_beta+1.96*out_strat3_se),
        #                              p_val=out_strat3_pvalue,
        #                              n=out_strat3_noobs,
        #                              sleep_traits=sleep_trait[i],
        #                              covar=paste(cov_mrs,collapse=","),
        #                              strata=paste(strat_value[3])
        #   )
    #       results<-rbind(results_both,results_strat1,results_strat2,results_strat3)
    #     }
    #     mrs_results<-rbind(mrs_results,results)
    #   } else {      results<-data.frame(outcome=outcome_var,
    #                       beta=out_beta,
    #                       se=out_se,
    #                       or=exp(out_beta),
    #                       lower95=exp(out_beta-1.96*out_se),
    #                       upper95=exp(out_beta+1.96*out_se),
    #                       p_val=out_pvalue,
    #                       n=out_noobs,
    #                       sleep_traits=sleep_trait[i],
    #                       covar=paste(cov_mrs,collapse=",")
    #   )
    #   mrs_results<-rbind(mrs_results,results)
    #   }
    # 
    #   # print(paste("i=",i))
    # }
    # 
    
  #   # apply pfilter on the MRS output
  #   # alpha_threshold_strict = c(0.05, 0.05,0.05)
  #   # we would like to control overall FDR at level 0.05,
  #   #	 	and first group-level(within trait) FDR at level 0.05
  #   # and second group-level(within metabolites) FDR at level 0.05
  #   
  #   # mrs_results$pfilter_1<-pfilter(P=mrs_results$p_val,alphas=alpha_threshold_strict,group=droplevels(mrs_results$outcome))
  #   # mrs_results$pfilter_2<-pfilter(P=mrs_results$p_val,alphas=alpha_threshold_strict,group=droplevels(mrs_results[,c("outcome","sleep_traits")]))
  #   # if ("GENDER"%in%cov_mrs){
  #   #   mrs_results$pfilter_1<-pfilter(P=mrs_results$p_val,alphas=alpha_threshold_strict,group=mrs_results$outcome)
  #   #   mrs_results$pfilter_2<-pfilter(P=mrs_results$p_val,alphas=alpha_threshold_strict,group=mrs_results[,c("outcome","sleep_traits")])
  #   # } else {
  #   mrs_results$pfilter_1<-pfilter(P=mrs_results$p_val,alphas=alpha_threshold_strict,group=mrs_results$outcome)
  #   mrs_results$pfilter_2<-pfilter(P=mrs_results$p_val,alphas=alpha_threshold_strict,group=mrs_results[,c("outcome","sleep_traits")])
  #   # }
  #   output=list(mrs_results=mrs_results,mrs_scores=mrs_scores)
  #   return(output)
  # }
  
  # Loop MRS analysis (only have gender stratification)
  loop_mrs<-function(metab_data,phenotype_data,sleep_trait,metab_list,covar,outcome_var,covariates_cog,admin_var,exclude_var,alpha_threshold_strict,id_bin,id_cont,stratifier){
    mrs_results<-data.frame()
    mrs_scores<-data.frame()
    for (i in seq_along(sleep_trait)){
      # print(paste("i=",i))
      
      # add ECGA7 to the covariate list if sleep trait is hr_dim1 or hr_dim2
      if (sleep_trait[i]%in%c("hr_pc1","hr_pc2","hr_dim1","hr_dim2","avg_hr","std_hr","min_hr","max_hr")){
        covar_res<-c(covar,"ECGA7")
      } else covar_res<-covar
      # Regress the rank-normalized metabolites and sleep trait of interest over covariates, obtain residuals.
      res_SDm<-residual_val(met_data=metab_data,pheno_data=phenotype_data,dep_list=c(sleep_trait[i],metab_list[[i]]),cov=covar_res)
      # Obtain SDr of residuals (including traits)
      res_sdr_inv_SDm<-diag(1/sqrt(diag(cov(res_SDm[,!colnames(res_SDm)%in%c("ID","LAB_ID")],use="pairwise.complete.obs"))))
      # obtain the ID for the complete obs
      res_sdr_inv_SDm_completeID<-res_SDm[complete.cases(res_SDm[,!colnames(res_SDm)%in%c("ID","LAB_ID")]),"ID"]
      # replace NAs with 0, multiply with the SDm and then revert back to NAs
      res_0_SDm<- res_SDm[,!colnames(res_SDm)%in%c("ID","LAB_ID")]
      res_0_SDm[is.na(res_SDm[,!colnames(res_SDm)%in%c("ID","LAB_ID")])] <- 0
      # Divide metabolite residuals by SDr for CCA (because with CCA we want variables to have the same variance). 
      res_rescaled_SDm <- as.matrix(res_0_SDm) %*% res_sdr_inv_SDm
      res_rescaled_SDm[is.na(res_SDm[,!colnames(res_SDm)%in%c("ID","LAB_ID")])] <- NA
      colnames(res_rescaled_SDm)<-colnames(res_SDm[,!colnames(res_SDm)%in%c("ID","LAB_ID")])
      res_rescaled_SDm<-cbind(res_SDm[,c("ID","LAB_ID")],res_rescaled_SDm)
      
      # res_rescaled_SDm_trait<-res_rescaled_SDm[which(!is.na(res_rescaled_SDm[,sleep_trait[i]])),] # drop obs with missing sleep trait for rcc
      res_rescaled_SDm_trait<-res_rescaled_SDm[complete.cases(res_rescaled_SDm),] # drop any missing obs
      
      # Perform CCA. Obtain coefficients: every trait (single column) and all the metabolites
      metab_matrix_SDm<-res_rescaled_SDm[,metab_list[[i]]]
      trait_matrix_SDm<-res_rescaled_SDm[,sleep_trait[i]]
      # non-missing value index for the sleep trait
      trait_matrix_SDm_rowindex<-which(!is.na(trait_matrix_SDm))
      # non-missing value index for the metabolites
      metab_matrix_SDm_rowindex<-rownames(metab_matrix_SDm[which(rowSums(is.na(metab_matrix_SDm))==0),])
      # row index for non-missing value for both matrix
      combined_rowindex<-as.numeric(metab_matrix_SDm_rowindex[metab_matrix_SDm_rowindex%in%trait_matrix_SDm_rowindex])
      
      # subset sdr for only the metabolites
      if (length(metab_list[[i]])==1){
        res_metab_sdr_SDm<-sd(res_SDm[,colnames(res_SDm)%in%metab_list[[i]]], na.rm=T)
      } else {
        res_metab_sdr_SDm<-sqrt(diag(cov(res_SDm[,colnames(res_SDm)%in%metab_list[[i]]],use="pairwise.complete.obs")))
      }
      # Perform CCA and multiply coefficients by SDr to go back to pre-CCA scale.(If only one metabolites then no need to multiply SD)
      # regul.par <- CCA::estim.regul(as.matrix(metab_matrix_SDm), as.matrix(trait_matrix_SDm), plt = T) # l=50
      # cc_output_SDm<-CCA::rcc(as.matrix(metab_matrix_SDm),as.matrix(trait_matrix_SDm),regul.par$lambda1, regul.par$lambda2)
      # cc_output_SDm<-mixOmics::rcc(as.matrix(metab_matrix_SDm),as.matrix(trait_matrix_SDm), ncomp = 3, method = 'shrinkage')
      
      # basic CCA that takes y as a vector, missing value included
      # cc_output_SDm<-CCA::cc(as.matrix(metab_matrix_SDm[trait_matrix_SDm_rowindex,]),as.matrix(trait_matrix_SDm[trait_matrix_SDm_rowindex])) # basic CCA to get beta_s
      # basic CCA that takes y as a vector, missing value from either matrix excluded
      cc_output_SDm<-CCA::cc(as.matrix(metab_matrix_SDm[combined_rowindex,]),as.matrix(trait_matrix_SDm[combined_rowindex])) # basic CCA to get beta_s
      
      
      # cc_output_SDm<-CCA::rcc(as.matrix(metab_matrix_SDm),as.matrix(trait_matrix_SDm),0.1,0.2) # basic CCA to get beta_s
      metab_coeff_SDm<-cc_output_SDm$xcoef/res_metab_sdr_SDm # SD of the residual metabolites
      
      # comput the canonical correlation coefficient and other parameters
      # metab_cc_corr<-comput(as.matrix(metab_matrix_SDm[combined_rowindex,]),as.matrix(trait_matrix_SDm[combined_rowindex]), cc_output_SDm)
      
      # metab_coeff_SDm<-data.frame(coefficient=cc_output_SDm$xcoef/res_metab_sdr_SDm,trait=sleep_trait[i]) # SD of the residual metabolites
      # metab_coeff_SDm$trait<-sleep_trait[i]
      
      # Construct the score using rank-normalized metabolites and the coefficients from the previous step
      x_SDm<-t(t(metab_data[,metab_list[[i]]])*as.vector(metab_coeff_SDm))
      y_SDm<-rowSums(x_SDm,na.rm = TRUE)
      # standardize the score
      score_SDm<-as.vector(scale(y_SDm))
      # score_SDm<-cbind(score_SDm,res_SDm$ID)
      score_SDm<-cbind(score_SDm,metab_data$ID)
      colnames(score_SDm)<-c(sleep_trait[i],"ID")
      # output the score
      if (i==1){
        mrs_scores<-score_SDm
      } else {
        mrs_scores<-merge(mrs_scores,score_SDm,by="ID")
      }
      
      # merge with the phenotype data
      phenotype_data$WEIGHT<-phenotype_data$WEIGHT_NORM_OVERALL_V2
      pheno_mrs_SDm<-merge(phenotype_data[,c(outcome_var,covar_res,covariates_cog,admin_var,exclude_var)],score_SDm,by="ID",all.x=T)
      # pheno_mrs_SDm_ID<-pheno_mrs_SDm[which(!is.na(pheno_mrs_SDm$WEIGHT)&pheno_mrs_SDm$DIABETES3_INDICATOR!="Yes"&pheno_mrs_SDm$HYPERTENSION!="Yes"),]
      pheno_mrs_SDm_ID<-pheno_mrs_SDm[which(!is.na(pheno_mrs_SDm$WEIGHT)),]
      
      # run the association analysis
      out_nvar=length(outcome_var)
      out_beta=rep(NA,out_nvar)
      out_se=rep(NA,out_nvar)
      out_pvalue=rep(NA,out_nvar)
      out_noobs=rep(NA,out_nvar)     
      out_female_nvar=length(outcome_var)
      out_female_beta=rep(NA,out_nvar)
      out_female_se=rep(NA,out_nvar)
      out_female_pvalue=rep(NA,out_nvar)
      out_female_noobs=rep(NA,out_nvar)
      out_male_nvar=length(outcome_var)
      out_male_beta=rep(NA,out_nvar)
      out_male_se=rep(NA,out_nvar)
      out_male_pvalue=rep(NA,out_nvar)
      out_male_noobs=rep(NA,out_nvar)
      number=1
      results=data.frame()
      for (l in seq_along(outcome_var)){
        # add time elapse between baseline and follow up for MCI regression
        if (outcome_var[l]%in%c("mci")){
          cov_mrs<-c(cov_mrs,covariates_cog)
        } else {
          cov_mrs<-covar_res
        }
        # print(paste("l",l))
        if ("GENDER"%in%cov_mrs){
          model<-by_svyglm_multi(data=pheno_mrs_SDm,outcome=outcome_var[l],exposure=sleep_trait[i],covar=cov_mrs,include_ID=pheno_mrs_SDm_ID$ID)
          model_female<-by_svyglm_multi(data=subset(pheno_mrs_SDm,GENDER=="Female"),outcome=outcome_var[l],exposure=sleep_trait[i],covar=cov_mrs[!cov_mrs%in%c("GENDER")],include_ID=pheno_mrs_SDm_ID$ID)
          model_male<-by_svyglm_multi(data=subset(pheno_mrs_SDm,GENDER=="Male"),outcome=outcome_var[l],exposure=sleep_trait[i],covar=cov_mrs[!cov_mrs%in%c("GENDER")],include_ID=pheno_mrs_SDm_ID$ID)
          Vcov <- vcov(model, useScale = FALSE)
          beta<- coef(model)
          se<- sqrt(diag(vcov(model, useScale = FALSE)))
          zval<- beta / se
          pval<- 2 * pnorm(abs(zval), lower.tail = FALSE)
          out_beta[number]=as.numeric(beta[2])
          out_se[number] = as.numeric(se[2])
          out_pvalue[number] = as.numeric(pval[2])
          out_noobs[number]=nobs(model)
          
          Vcov_female <- vcov(model_female, useScale = FALSE)
          beta_female<- coef(model_female)
          se_female<- sqrt(diag(vcov(model_female, useScale = FALSE)))
          zval_female<- beta_female / se_female
          pval_female<- 2 * pnorm(abs(zval_female), lower.tail = FALSE)
          out_female_beta[number]=as.numeric(beta_female[2])
          out_female_se[number] = as.numeric(se_female[2])
          out_female_pvalue[number] = as.numeric(pval_female[2])
          out_female_noobs[number]=nobs(model_female)
          Vcov_male <- vcov(model_male, useScale = FALSE)
          beta_male<- coef(model_male)
          se_male<- sqrt(diag(vcov(model_male, useScale = FALSE)))
          zval_male<- beta_male / se_male
          pval_male<- 2 * pnorm(abs(zval_male), lower.tail = FALSE)
          out_male_beta[number]=as.numeric(beta_male[2])
          out_male_se[number] = as.numeric(se_male[2])
          out_male_pvalue[number] = as.numeric(pval_male[2])
          out_male_noobs[number]=nobs(model_male)          
          number=number+1
        } else {
          model<-by_svyglm_multi(data=pheno_mrs_SDm,outcome=outcome_var[l],exposure=sleep_trait[i],covar=cov_mrs,include_ID=pheno_mrs_SDm_ID$ID)
          Vcov <- vcov(model, useScale = FALSE)
          beta<- coef(model)
          se<- sqrt(diag(vcov(model, useScale = FALSE)))
          zval<- beta / se
          pval<- 2 * pnorm(abs(zval), lower.tail = FALSE)
          out_beta[number]=as.numeric(beta[2])
          out_se[number] = as.numeric(se[2])
          out_pvalue[number] = as.numeric(pval[2])
          out_noobs[number]=nobs(model)
          number=number+1
          # print(paste("l=",l))
        }
      }
      if ("GENDER"%in%cov_mrs){
        results_both<-data.frame(outcome=outcome_var,
                                 beta=out_beta,
                                 se=out_se,
                                 or=exp(out_beta),
                                 lower95=exp(out_beta-1.96*out_se),
                                 upper95=exp(out_beta+1.96*out_se),
                                 p_val=out_pvalue,
                                 n=out_noobs,
                                 sleep_traits=sleep_trait[i],
                                 covar=paste(cov_mrs,collapse=","),
                                 strata="Both"
        )
        results_female<-data.frame(outcome=outcome_var,
                                   beta=out_female_beta,
                                   se=out_female_se,
                                   or=exp(out_female_beta),
                                   lower95=exp(out_female_beta-1.96*out_female_se),
                                   upper95=exp(out_female_beta+1.96*out_female_se),
                                   p_val=out_female_pvalue,
                                   n=out_female_noobs,
                                   sleep_traits=sleep_trait[i],
                                   covar=paste(cov_mrs,collapse=","),
                                   strata="Female"
        )
        results_male<-data.frame(outcome=outcome_var,
                                 beta=out_male_beta,
                                 se=out_male_se,
                                 or=exp(out_male_beta),
                                 lower95=exp(out_male_beta-1.96*out_male_se),
                                 upper95=exp(out_male_beta+1.96*out_male_se),
                                 p_val=out_male_pvalue,
                                 n=out_male_noobs,
                                 sleep_traits=sleep_trait[i],
                                 covar=paste(cov_mrs,collapse=","),
                                 strata="Male"
        )
        results<-rbind(results_both,results_female,results_male)
        mrs_results<-rbind(mrs_results,results)
      } else {      results<-data.frame(outcome=outcome_var,
                                        beta=out_beta,
                                        se=out_se,
                                        or=exp(out_beta),
                                        lower95=exp(out_beta-1.96*out_se),
                                        upper95=exp(out_beta+1.96*out_se),
                                        p_val=out_pvalue,
                                        n=out_noobs,
                                        sleep_traits=sleep_trait[i],
                                        covar=paste(cov_mrs,collapse=",")
      )
      mrs_results<-rbind(mrs_results,results)
      }
      
      # print(paste("i=",i))
    }
    
    
    # apply pfilter on the MRS output
    # alpha_threshold_strict = c(0.05, 0.05,0.05)
    # we would like to control overall FDR at level 0.05,
    #	 	and first group-level(within trait) FDR at level 0.05
    # and second group-level(within metabolites) FDR at level 0.05
    
    # mrs_results$pfilter_1<-pfilter(P=mrs_results$p_val,alphas=alpha_threshold_strict,group=droplevels(mrs_results$outcome))
    # mrs_results$pfilter_2<-pfilter(P=mrs_results$p_val,alphas=alpha_threshold_strict,group=droplevels(mrs_results[,c("outcome","sleep_traits")]))
    # if ("GENDER"%in%cov_mrs){
    #   mrs_results$pfilter_1<-pfilter(P=mrs_results$p_val,alphas=alpha_threshold_strict,group=mrs_results$outcome)
    #   mrs_results$pfilter_2<-pfilter(P=mrs_results$p_val,alphas=alpha_threshold_strict,group=mrs_results[,c("outcome","sleep_traits")])
    # } else {
    mrs_results$pfilter_1<-pfilter(P=mrs_results$p_val,alphas=alpha_threshold_strict,group=mrs_results$outcome)
    mrs_results$pfilter_2<-pfilter(P=mrs_results$p_val,alphas=alpha_threshold_strict,group=mrs_results[,c("outcome","sleep_traits")])
    # }
    output=list(mrs_results=mrs_results,mrs_scores=mrs_scores)
    return(output)
  }
  
  
  #####################
  # Split violin plot
  # GeomSplitViolin <- ggplot2::ggproto("GeomSplitViolin", GeomViolin,
  #                            draw_group = function(self, data, ..., draw_quantiles = NULL) {
  #                              data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
  #                              grp <- data[1, "group"]
  #                              newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
  #                              newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
  #                              newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
  # 
  #                              if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
  #                                stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
  #                                                                          1))
  #                                quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
  #                                aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
  #                                aesthetics$alpha <- rep(1, nrow(quantiles))
  #                                both <- cbind(quantiles, aesthetics)
  #                                quantile_grob <- GeomPath$draw_panel(both, ...)
  #                                ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
  #                              }
  #                              else {
  #                                ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
  #                              }
  #                            })
  
  geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                                draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                                show.legend = NA, inherit.aes = TRUE) {
    layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
          position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
          params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
  }
  
  # functions to construct split violin plots for each trait
  gen_split_violin<-function(mrs_data,pheno_data,trait,trait_label,group,group_label,gender,ind_title,ind_legend,legen_position_var){
    score_df<-mrs_data[,c("ID",trait)]
    combined_df<-merge(score_df,pheno_data[,c("ID",trait,group)],by="ID",suffixes=c(".mrs",".trait"),all.x = T)%>%
      tidyr::drop_na()%>%
      tidyr::gather(., "variable", "value", 2:3)
      # dplyr::filter(!is.na(background))
      combined_df<-combined_df[which(!is.na(combined_df[,group])),]
    
    combined_df$variable<-substring(combined_df$variable, regexpr("[.]", combined_df$variable)+1, nchar(combined_df$variable))
    
    colnames(combined_df)[2]<-"group"
    
    if (ind_title==T){
      if (ind_legend==T){
        ggplot(combined_df, aes(x=group,y= value, fill = variable)) + 
          geom_split_violin(trim = TRUE) + 
          geom_boxplot(width = 0.25, notch = FALSE, notchwidth = .4, outlier.shape = NA, coef=0) +
          labs(title=paste("Distribution of Trait and Trait-Specific MRS by",group_label,"\namong",gender),
               x =paste(group_label), y = "Density",fill = paste(trait_label)) +
          theme_classic() +
          theme(text = element_text(size = 12))
      }else {
        ggplot(combined_df, aes(x=group,y= value, fill = variable)) + 
      geom_split_violin(trim = TRUE) + 
      geom_boxplot(width = 0.25, notch = FALSE, notchwidth = .4, outlier.shape = NA, coef=0) +
      labs(title=paste("Distribution of Trait and Trait-Specific MRS by",group_label,"\namong",gender),
           x =paste(group_label), y = "Density",fill = paste(trait_label)) +
      theme_classic() + theme(legend.position = "none")+
      theme(text = element_text(size = 12))
      }
       
    } else {
      if (ind_legend==T){
        ggplot(combined_df, aes(x=group,y= value, fill = variable)) + 
          geom_split_violin(trim = TRUE) + 
          geom_boxplot(width = 0.25, notch = FALSE, notchwidth = .4, outlier.shape = NA, coef=0) +
          labs(title=paste(gender),
               x =paste(group_label), y = "Density",fill = paste(trait_label)) +
          theme_classic() + theme(legend.position=legen_position_var)+
          theme(text = element_text(size = 12)) 
      } else {
        ggplot(combined_df, aes(x=group,y= value, fill = variable)) + 
        geom_split_violin(trim = TRUE) + 
        geom_boxplot(width = 0.25, notch = FALSE, notchwidth = .4, outlier.shape = NA, coef=0) +
        labs(title=paste(gender),
             x =paste(group_label), y = "Density",fill = paste(trait_label)) +
        theme_classic() + theme(legend.position="none")+
        theme(text = element_text(size = 12)) 
      }
      
    }

  }
  
  
  # prediction and performance functions are from package ROCR
  compute_auc <- function(p, labels) {
    pred <- ROCR::prediction(predictions=p, labels=labels)
    auc <- ROCR::performance(pred, 'auc')
    auc <- unlist(slot(auc, 'y.values'))
    auc
  }
  
  
  merg_confmx<-function(data_list,col_name){
    # data<-lapply(data_list,function(x) x[!x$term%in%"mcnemar",col_name])
    data<-lapply(data_list,function(x) x[,col_name])
    
    data_frame<-dplyr::bind_cols(data)
    colnames(data_frame)<-names(data_list)
    data_frame<-cbind(data_list[[1]][,1],data_frame)
    return(data_frame)
  }
  
  
  # Multiple the predictor coefficient with the original STD
  multi_std_fun<-function(data,coeff){
    coeff$std<-NA
    coeff$std_coeff<-NA
    for (i in 1:nrow(coeff)){
      coeff[i,"std"]<-sd(data[,coeff[i,1]])
      coeff[i,"std_coeff"]<-sd(data[,coeff[i,1]])*coeff[i,2]
    }
    return(coeff)
  }
  
  
  ############################
  # modify tukey's power transformation
  # the original function (https://rdrr.io/cran/rcompanion/src/R/transformTukey.r) can't handle samples >5000 due to the limitation of Shapiro-Wilks test
  # the original function calculates S-W test even if Anderson-Darling is selected (statistic=2)
  # the modification simply takes S-W test statistics out of A-D test to avoid error
  transformTukey_mod = 
    function(x, start=-10, end=10, int=0.025,
             plotit=TRUE, verbose=FALSE, quiet=FALSE, statistic=1, 
             returnLambda=FALSE)
    {
      n=(end-start)/int
      lambda=as.numeric(rep(0.00, n))
      W=as.numeric(rep(0.00, n))
      Shapiro.p.value=as.numeric(rep(0.00, n))
      if(statistic==2){
        A=as.numeric(rep(1000.00, n))
        Anderson.p.value=as.numeric(rep(0.00, n))
      }
      for(i in (1:n)){
        lambda[i] = signif(start+(i-1)*int, digits=4)
        if (lambda[i]>0) {TRANS = x ^ lambda[i]}
        if (lambda[i]==0){TRANS = log(x)}
        if (lambda[i]<0) {TRANS = -1 * x ^ lambda[i]}
        W[i]=NA
        if(statistic==2){A[i]=NA}
        if (any(is.infinite(TRANS))==FALSE & any(is.nan(TRANS))==FALSE){
          if(statistic==1){
            W[i]=signif(shapiro.test(TRANS)$statistic, digits=4)
            Shapiro.p.value[i]=signif(shapiro.test(TRANS)$p.value, digits=4)
            }else if(statistic==2){
            A[i]=signif(ad.test(TRANS)$statistic, digits=4)
            Anderson.p.value[i]=signif(ad.test(TRANS)$p.value, digits=4)
          }
        } 
      }
      if(statistic==2){
        df=data.frame(lambda, A, Anderson.p.value)
      }
      if(statistic==1){
        df=data.frame(lambda, W, Shapiro.p.value)
      }
      if(verbose==TRUE){print(df)}
      if(plotit==TRUE){
        if(statistic==1){plot(lambda, W, col="black")}
        if(statistic==2){plot(lambda, A, col="blue")}
      }
      if(statistic==1){df2 = df[with(df, order(-W)),]}
      if(statistic==2){df2 = df[with(df, order(A)),]}
      if(quiet==FALSE){
        cat("\n")
        print(df2[1,])
        cat("\n")
        cat("if (lambda >  0){TRANS = x ^ lambda}","\n")
        cat("if (lambda == 0){TRANS = log(x)}","\n")
        cat("if (lambda <  0){TRANS = -1 * x ^ lambda}","\n")
        cat("\n")
      }
      lambda = df2[1,"lambda"]
      if (lambda>0) {TRANS = x ^ lambda}
      if (lambda==0){TRANS = log(x)}
      if (lambda<0) {TRANS = -1 * x ^ lambda}
      if(plotit==TRUE){
        plotNormalHistogram(TRANS, xlab="Transformed variable", 
                            linecol="red", col="lightgray")
      }
      if(plotit==TRUE){
        qqnorm(TRANS)
        qqline(TRANS, col="red")
      }
      if(returnLambda==FALSE){return(TRANS)}
      if(returnLambda==TRUE){
        names(lambda)="lambda"
        return(lambda)
      }
    }
  
  # function to merge two table 1 (weighted and nonweighted) and to write into a csv
  merge_table1<-function(tbl_nonweight,tbl_weight,out_path){
    # options(stringsAsFactors = FALSE)
    tbl_weight<-as.data.frame(tbl_weight, stringsAsFactors = FALSE)
    tbl_nonweight<-as.data.frame(tbl_nonweight, stringsAsFactors = FALSE)
    new_tbl<-tbl_weight # use the weighted table as base (including the p, test, missing)
    rownames(new_tbl)<-rownames(tbl_nonweight)
    if (ncol(tbl_weight)==2){
      end_col<-1
    } else {
      end_col<-ncol(tbl_weight)-2
    }
    for (i in 1:end_col){ # loop over each strata
      for (j in 1:nrow(tbl_weight)){ # loop over each row
        if (stringr::str_detect(rownames(tbl_weight)[j],"mean..SD")){ # if "mean  (SD)" is detected in the rowname
          new_tbl[j,i]<-new_tbl[j,i] # then don't modify the weighted table
        } else if (stringr::str_detect(new_tbl[j,i],"^\\s+$")) { # if the cell is empty (only contains white spaces)
          new_tbl[j,i]<-new_tbl[j,i]
        } else { 
          if (str_split_fixed(tbl_nonweight[j,i], '[(]', 2)[,2]==""){ # if no percentage is calculated (for the "n" row)
            new_tbl[j,i]<-tbl_nonweight[j,i] # then paste only the nonweighted count number
          } else { # otherwise stitch the count from the nonweighted and the percentage from the weighted together
            new_tbl[j,i]<-paste0(str_split_fixed(tbl_nonweight[j,i], '[(]', 2)[,1],"(",str_split_fixed(tbl_nonweight[j,i], '[(]', 2)[,2])
          }
          
        }
      }
    }
    write.csv(new_tbl, file = out_path)
  }
  
  
  
  # function to compare results from OSA and AHI
  compare_annot_trait<-function(trait1_output,trait2_output,trait1,trait2,annot) {
      trait1_sig_results<-trait1_output[which(trait1_output$p_val_fdr_named<=0.05),]
      trait2_sig_results<-trait2_output[which(trait2_output$p_val_fdr_named<=0.05),]

      trait1_candidate<-annot[which(annot$metabolite%in%trait1_sig_results$metabolite),]
      trait2_candidate<-annot[which(annot$metabolite%in%trait2_sig_results$metabolite),]
    
    if (nrow(trait1_candidate)==0&nrow(trait2_candidate)==0){
      print("No metabolite candidate identified")
    } else if (nrow(trait1_candidate)==0&nrow(trait2_candidate)!=0) {
      trait2_tbl<-merge(trait2_sig_results,trait2_candidate,by="metabolite")
      trait2_unique<-trait2_tbl%>%
        dplyr::mutate(unique=paste0(trait2))%>%
        merge(.,trait1_output[,c("metabolite","beta","se","p_val_fdr","p_val")],by="metabolite",suffixes=c(paste0("_",trait2),paste0("_",trait1)))
      union_table<-trait2_unique%>%
        mutate_at(vars(contains('beta_')),funs(round(., 3)))%>%
        mutate_at(vars(contains('se_')),funs(round(., 3)))%>%
        mutate_at(vars(contains('p_val_')),funs(formatC(., format = "e", digits = 2)))
      return(as.data.frame(union_table[,c(1,12,4,18:22,25:28,30,29,2,3,5,7)]))
    } else if (nrow(trait1_candidate)!=0&nrow(trait2_candidate)==0) {
      trait1_tbl<-merge(trait1_sig_results,trait1_candidate,by="metabolite")
      trait1_unique<-trait1_tbl%>%
        dplyr::mutate(unique=paste0(trait1))%>%
        merge(.,trait1_output[,c("metabolite","beta","se","p_val_fdr","p_val")],by="metabolite",suffixes=c(paste0("_",trait2),paste0("_",trait1)))
      union_table<-trait1_unique%>%
        mutate_at(vars(contains('beta_')),funs(round(., 3)))%>%
        mutate_at(vars(contains('se_')),funs(round(., 3)))%>%
        mutate_at(vars(contains('p_val_')),funs(formatC(., format = "e", digits = 2)))
      return(as.data.frame(union_table[,c(1,12,4,18:22,25:28,30,29,2,3,5,7)]))
      } else {
      trait1_tbl<-merge(trait1_sig_results,trait1_candidate,by="metabolite")
      trait2_tbl<-merge(trait2_sig_results,trait2_candidate,by="metabolite")
      
      trait1_unique<-trait1_tbl%>%
        dplyr::filter(!metabolite%in%trait2_sig_results$metabolite)%>%
        dplyr::mutate(unique=paste0(trait1))%>%
        merge(.,trait2_output[,c("metabolite","beta","se","p_val_fdr","p_val")],by="metabolite",suffixes=c(paste0("_",trait2),paste0("_",trait1)))
      
      trait2_unique<-trait2_tbl%>%
        dplyr::filter(!metabolite%in%trait1_sig_results$metabolite)%>%
        dplyr::mutate(unique=paste0(trait2))%>%
        merge(.,trait2_output[,c("metabolite","beta","se","p_val_fdr","p_val")],by="metabolite",suffixes=c(paste0("_",trait2),paste0("_",trait1)))
      
      trait1_shared<-trait1_tbl%>%
        dplyr::filter(metabolite%in%trait2_sig_results$metabolite)%>%
        mutate(unique="Shared")
      
      trait2_shared<-trait2_tbl%>%
        dplyr::filter(metabolite%in%trait1_tbl$metabolite)%>%
        dplyr::mutate(unique="Shared")
      
      both_shared<-merge(trait1_shared[,c("metabolite","biochemical","super_pathway","sub_pathway","platform","mass","pubchem","kegg","hmdb","beta","se","p_val_fdr","p_val")],trait2_shared[,c("metabolite","beta","se","p_val_fdr","p_val")],by="metabolite",suffixes = c(paste0("_",trait2),paste0("_",trait1)))%>%
        dplyr::mutate(unique="Shared")
      
      union_table<-rbind.fill(trait1_unique,trait2_unique,both_shared)
      # write.csv(union_table,file = paste0("output/metabolite_candidates_",params$trait_input,".csv"))
      
      union_table<-union_table%>%
        mutate_at(vars(contains('beta_')),funs(round(., 3)))%>%
        mutate_at(vars(contains('se_')),funs(round(., 3)))%>%
        mutate_at(vars(contains('p_val_')),funs(formatC(., format = "e", digits = 2)))
      return(as.data.frame(union_table[,c(1,12,4,18:22,25:28,30,29,2,3,5,7)]))
    }
  }
  
  # Function to loop regression model and export model summary 
  reg_loop<-function(data,covar,start,end,metab_is_cont,metab_is_complete,trait_original,trait_binary,trait_for_model,trait_as_predictor){
    if (trait_for_model=="binary") {
      trait<-"TraitBinary"
    } else if (trait_for_model=="original"){
      trait<-trait_original
    } else { trait<-"Trait" }
    dat<-data[,start:end]
    dat<-dat[,!colnames(dat)%in%c("LAB_ID","ID",trait)]
    out_nvar=ncol(dat)
    out_beta=rep(NA,out_nvar)
    out_se=rep(NA,out_nvar)
    out_pvalue=rep(NA,out_nvar)
    out_nobs=rep(NA,out_nvar)

    
    # if the original trait variable is binary or the user pick the binary version to use in the model
    trait_is_cont=ifelse(length(unique(data[!is.na(data[,trait]),trait]))==2,FALSE,
                         ifelse(trait_for_model=="binary",FALSE,TRUE))
    
    for(i in 1:(ncol(dat))) {
      met_df<-cbind(data[,colnames(dat)[i]],data[,c("ID",covar)])
      met_df_id<-met_df[complete.cases(met_df),"ID"]
      if (trait_as_predictor==T){
        if (metab_is_cont==T){
          model<- glm(as.formula(paste0(colnames(dat)[i],"~",trait,"+",paste(covar,collapse= "+"))),data=subset(data,ID%in%met_df_id))
        } else {
          model<- glm(as.formula(paste0(colnames(dat)[i],"~",trait,"+",paste(covar,collapse = "+"))),data=subset(data,ID%in%met_df_id),family=quasipoisson(link='log'))
        }
      } else {
        if (trait_is_cont==T){
          model<- glm(as.formula(paste0(trait,"~",colnames(dat)[i],"+",paste(covar,collapse= "+"))),data=subset(data,ID%in%met_df_id))
        } else {
          model<- glm(as.formula(paste0(trait,"~",colnames(dat)[i],"+",paste(covar,collapse = "+"))),data=subset(data,ID%in%met_df_id),family=quasipoisson(link='log'))
        }
      }
      
      Vcov <- vcov(model, useScale = FALSE)
      beta<- coef(model)
      se<- sqrt(diag(vcov(model, useScale = FALSE)))
      zval<- beta / se
      pval<- 2 * pnorm(abs(zval), lower.tail = FALSE)
      out_beta[i]=as.numeric(beta[2])
      out_se[i] = round(as.numeric(se[2]),digits = 3)
      out_pvalue[i] = as.numeric(pval[2])
      out_nobs[i]=as.numeric(nobs(model))
    }
    regress_output<-data.frame(metabolite=colnames(dat),
                               beta=out_beta,
                               se=out_se,
                               n=out_nobs,
                               p_val=out_pvalue,
                               p_val_Bonf=p.adjust(out_pvalue,method="bonferroni"),
                               p_val_fdr=p.adjust(out_pvalue,method="BH")
    )%>%
      mutate(
        pval_Bonf_neglog=-log10(p_val_Bonf),
        pval_fdr_neglog=-log10(p_val_fdr),
        sig_Bonf=ifelse(p_val_Bonf<=0.05,"Bonf<0.05","Not Sig"),
        sig_fdr=ifelse(p_val_fdr<=0.05,"FDR<0.05","Not Sig"),
        is_continuous=metab_is_cont
      )
    identified_data<-regress_output%>%
      dplyr::filter(!stringr::str_detect(metabolite,"^X"))%>%
      dplyr::mutate(p_val_fdr_named=p.adjust(p_val,method="BH"),
                    sig_fdr_named=ifelse(p_val_fdr_named<=0.05,"FDR<0.05","Not Sig"))
    unidentified_data<-regress_output%>%
      dplyr::filter(stringr::str_detect(metabolite,"^X"))%>%
      dplyr::mutate(p_val_fdr_named=NA,
                    sig_fdr_named=NA)
    regress_output<-rbind(identified_data,unidentified_data)
    if (metab_is_complete==T){
      regress_output$is_continuous<-"Complete-cases"
    } else {
      regress_output$is_continuous<-factor(regress_output$is_continuous,levels=c(TRUE,FALSE),labels=c("Imputed","Dichotomized"))
    }
    
    return(regress_output)
  }

  # Function to loop regression model and export model summary 
  mesa_reg_loop<-function(data,covar,start,end,metab_is_cont,metab_is_complete,trait_original,trait_binary,trait_for_model,trait_as_predictor){
    if (trait_for_model=="binary") {
      trait<-"TraitBinary"
    } else if (trait_for_model=="original"){
      trait<-trait_original
    } else { trait<-"Trait" }
    dat<-data[,start:end]
    dat<-dat[,!colnames(dat)%in%c("LAB_ID","ID",trait)]
    out_nvar=ncol(dat)
    out_beta=rep(NA,out_nvar)
    out_se=rep(NA,out_nvar)
    out_pvalue=rep(NA,out_nvar)
    out_nobs=rep(NA,out_nvar)
    
    
    # if the original trait variable is binary or the user pick the binary version to use in the model in MESA with incident outcomes
    trait_is_cont=ifelse(length(unique(data[!is.na(data[,trait]),trait]))==2,FALSE,
                         ifelse(trait_for_model=="binary",FALSE,TRUE))
    
    for(i in 1:(ncol(dat))) {
      met_df<-cbind(data[,colnames(dat)[i]],data[,c("ID",covar)])
      met_df_id<-met_df[complete.cases(met_df),"ID"]
      if (trait_as_predictor==T){
        if (metab_is_cont==T){
          model<- glm(as.formula(paste0(colnames(dat)[i],"~",trait,"+",paste(covar,collapse= "+"))),data=subset(data,ID%in%met_df_id))
        } else {
          if (colnames(dat)[i]=="inc_htn") {
            model<- glm(as.formula(paste0(colnames(dat)[i],"~",trait,"+",paste(covar,collapse = "+"))),data=subset(data,ID%in%met_df_id&prev_htn!=1),family=quasipoisson(link='log'))
          } else if (colnames(dat)[i]=="inc_cvd"){
            model<- glm(as.formula(paste0(colnames(dat)[i],"~",trait,"+",paste(covar,collapse = "+"))),data=subset(data,ID%in%met_df_id&prev_cvd!=1),family=quasipoisson(link='log'))
          } else if (colnames(dat)[i]=="inc_dm") {
            model<- glm(as.formula(paste0(colnames(dat)[i],"~",trait,"+",paste(covar,collapse = "+"))),data=subset(data,ID%in%met_df_id&prev_dm!=1),family=quasipoisson(link='log'))
          } else {
            model<- glm(as.formula(paste0(colnames(dat)[i],"~",trait,"+",paste(covar,collapse = "+"))),data=subset(data,ID%in%met_df_id),family=quasipoisson(link='log'))
          }
        }
      } else {
        if (trait_is_cont==T){
          model<- glm(as.formula(paste0(trait,"~",colnames(dat)[i],"+",paste(covar,collapse= "+"))),data=subset(data,ID%in%met_df_id))
        } else {
          if (colnames(dat)[i]=="inc_htn") {
            model<- glm(as.formula(paste0(trait,"~",colnames(dat)[i],"+",paste(covar,collapse = "+"))),data=subset(data,ID%in%met_df_id&prev_htn!=1),family=quasipoisson(link='log'))
          } else if (colnames(dat)[i]=="inc_cvd"){
            model<- glm(as.formula(paste0(trait,"~",colnames(dat)[i],"+",paste(covar,collapse = "+"))),data=subset(data,ID%in%met_df_id&prev_cvd!=1),family=quasipoisson(link='log'))
          } else if (colnames(dat)[i]=="inc_dm") {
            model<- glm(as.formula(paste0(trait,"~",colnames(dat)[i],"+",paste(covar,collapse = "+"))),data=subset(data,ID%in%met_df_id&prev_dm!=1),family=quasipoisson(link='log'))
          } else {
            model<- glm(as.formula(paste0(trait,"~",colnames(dat)[i],"+",paste(covar,collapse = "+"))),data=subset(data,ID%in%met_df_id),family=quasipoisson(link='log'))
          }
        }
      }
      
      Vcov <- vcov(model, useScale = FALSE)
      beta<- coef(model)
      se<- sqrt(diag(vcov(model, useScale = FALSE)))
      zval<- beta / se
      pval<- 2 * pnorm(abs(zval), lower.tail = FALSE)
      out_beta[i]=as.numeric(beta[2])
      out_se[i] = round(as.numeric(se[2]),digits = 3)
      out_pvalue[i] = as.numeric(pval[2])
      out_nobs[i]=as.numeric(nobs(model))
    }
    regress_output<-data.frame(metabolite=colnames(dat),
                               beta=out_beta,
                               se=out_se,
                               n=out_nobs,
                               p_val=out_pvalue,
                               p_val_Bonf=p.adjust(out_pvalue,method="bonferroni"),
                               p_val_fdr=p.adjust(out_pvalue,method="BH")
    )%>%
      mutate(
        pval_Bonf_neglog=-log10(p_val_Bonf),
        pval_fdr_neglog=-log10(p_val_fdr),
        sig_Bonf=ifelse(p_val_Bonf<=0.05,"Bonf<0.05","Not Sig"),
        sig_fdr=ifelse(p_val_fdr<=0.05,"FDR<0.05","Not Sig"),
        is_continuous=metab_is_cont
      )
    identified_data<-regress_output%>%
      dplyr::filter(!stringr::str_detect(metabolite,"^X"))%>%
      dplyr::mutate(p_val_fdr_named=p.adjust(p_val,method="BH"),
                    sig_fdr_named=ifelse(p_val_fdr_named<=0.05,"FDR<0.05","Not Sig"))
    unidentified_data<-regress_output%>%
      dplyr::filter(stringr::str_detect(metabolite,"^X"))%>%
      dplyr::mutate(p_val_fdr_named=NA,
                    sig_fdr_named=NA)
    regress_output<-rbind(identified_data,unidentified_data)
    if (metab_is_complete==T){
      regress_output$is_continuous<-"Complete-cases"
    } else {
      regress_output$is_continuous<-factor(regress_output$is_continuous,levels=c(TRUE,FALSE),labels=c("Imputed","Dichotomized"))
    }
    
    return(regress_output)
  }  
  
  # Function to loop regression model and export model summary
  strat_reg_loop<-function(data,covar,start,end,metab_is_cont,metab_is_complete,trait_original,trait_binary,trait_for_model,trait_as_predictor,stratifier){
    if (trait_for_model=="binary") {
      trait<-"TraitBinary"
    }  else if (trait_for_model=="original"){
      trait<-trait_original
    } else { trait<-"Trait" }
    
    dat<-data[,start:end]
    dat<-dat[,!colnames(dat)%in%c("ID","LAB_ID",trait)]
    
    covar<-covar[!covar%in%stratifier]
    group_by_value<-as.character(unique(data[,stratifier]))
    out_nvar<-ncol(dat)
    out_beta<-rep(NA,out_nvar)
    out_se<-rep(NA,out_nvar)
    out_pvalue<-rep(NA,out_nvar)
    out_nobs<-rep(NA,out_nvar)
    group<-rep(NA,out_nvar)
    regress_output<-data.frame()
    # if the original trait variable is binary or the user pick the binary version to use in the model
    trait_is_cont<-ifelse(length(unique(data[!is.na(data[,trait]),trait]))==2,FALSE,
                          ifelse(trait_for_model=="binary",FALSE,TRUE))
    
    
    for (j in seq_along(group_by_value)){
      newdata<-data[which(data[,stratifier]==as.character(group_by_value[j])),]
      dat<-newdata[,start:end]
      dat<-dat[,!colnames(dat)%in%c("ID","LAB_ID",trait)]
      for(i in 1:(ncol(dat))) {
        met_df<-cbind(dat[,colnames(dat)[i]],newdata[,c("ID",covar)])
        met_df_id<-met_df[complete.cases(met_df),"ID"]
        if (trait_as_predictor==T){
          if (metab_is_cont==T){
            model<- glm(as.formula(paste0(colnames(dat)[i],"~",trait,"+",paste(covar,collapse= "+"))),data=subset(data,ID%in%met_df_id))
          } else {
            if (colnames(dat)[i]=="inc_htn") {
              model<- glm(as.formula(paste0(colnames(dat)[i],"~",trait,"+",paste(covar,collapse = "+"))),data=subset(data,ID%in%met_df_id&prev_htn!=1),family=quasipoisson(link='log'))
            } else if (colnames(dat)[i]=="inc_cvd"){
              model<- glm(as.formula(paste0(colnames(dat)[i],"~",trait,"+",paste(covar,collapse = "+"))),data=subset(data,ID%in%met_df_id&prev_cvd!=1),family=quasipoisson(link='log'))
            } else if (colnames(dat)[i]=="inc_dm") {
              model<- glm(as.formula(paste0(colnames(dat)[i],"~",trait,"+",paste(covar,collapse = "+"))),data=subset(data,ID%in%met_df_id&prev_dm!=1),family=quasipoisson(link='log'))
            } else {
              model<- glm(as.formula(paste0(colnames(dat)[i],"~",trait,"+",paste(covar,collapse = "+"))),data=subset(data,ID%in%met_df_id),family=quasipoisson(link='log'))
            }
          }
        }else {
          if (trait_is_cont==T){
            model<- glm(as.formula(paste0(trait,"~",colnames(dat)[i],"+",paste(covar,collapse= "+"))),data=subset(data,ID%in%met_df_id))
          }else{
            if (colnames(dat)[i]=="inc_htn") {
              model<- glm(as.formula(paste0(trait,"~",colnames(dat)[i],"+",paste(covar,collapse = "+"))),data=subset(data,ID%in%met_df_id&prev_htn!=1),family=quasipoisson(link='log'))
            } else if (colnames(dat)[i]=="inc_cvd"){
              model<- glm(as.formula(paste0(trait,"~",colnames(dat)[i],"+",paste(covar,collapse = "+"))),data=subset(data,ID%in%met_df_id&prev_cvd!=1),family=quasipoisson(link='log'))
            } else if (colnames(dat)[i]=="inc_dm") {
              model<- glm(as.formula(paste0(trait,"~",colnames(dat)[i],"+",paste(covar,collapse = "+"))),data=subset(data,ID%in%met_df_id&prev_dm!=1),family=quasipoisson(link='log'))
            } else {
              model<- glm(as.formula(paste0(trait,"~",colnames(dat)[i],"+",paste(covar,collapse = "+"))),data=subset(data,ID%in%met_df_id),family=quasipoisson(link='log'))
            }
          }
        }
        
        Vcov <- vcov(model, useScale = FALSE)
        beta<- coef(model)
        se<- sqrt(diag(vcov(model, useScale = FALSE)))
        zval<- beta / se
        pval<- 2 * pnorm(abs(zval), lower.tail = FALSE)
        out_beta[i]=as.numeric(beta[2])
        out_se[i] = as.numeric(se[2])
        out_pvalue[i] = as.numeric(pval[2])
        out_nobs[i]=as.numeric(nobs(model))
        group[i]=group_by_value[j]
      }
      strata_output<-data.frame(metabolite=colnames(dat),
                                beta=out_beta,
                                se=out_se,
                                n=out_nobs,
                                p_val=out_pvalue,
                                p_val_Bonf=p.adjust(out_pvalue,method="bonferroni"),
                                p_val_fdr=p.adjust(out_pvalue,method="BH"),
                                strata=group)%>%
        mutate(
          pval_Bonf_neglog=-log10(p_val_Bonf),
          pval_fdr_neglog=-log10(p_val_fdr),
          sig_Bonf=ifelse(p_val_Bonf<=0.05,"Bonf<0.05","Not Sig"),
          sig_fdr=ifelse(p_val_fdr<=0.05,"FDR<0.05","Not Sig"),
          is_continuous=metab_is_cont
        )
      regress_output<-rbind(regress_output,strata_output)
    }
    identified_data<-regress_output%>%
      dplyr::filter(!stringr::str_detect(metabolite,"^X"))%>%
      group_by(strata)%>%
      dplyr::mutate(p_val_fdr_named=p.adjust(p_val,method="BH"),
                    sig_fdr_named=ifelse(p_val_fdr_named<=0.05,"FDR<0.05","Not Sig"))
    unidentified_data<-regress_output%>%
      dplyr::filter(stringr::str_detect(metabolite,"^X"))%>%
      dplyr::mutate(p_val_fdr_named=NA,
                    sig_fdr_named=NA)
    regress_output<-plyr::rbind.fill(identified_data,unidentified_data)
    if (metab_is_complete==T){
      regress_output$is_continuous<-"Complete-cases"
    } else {
      regress_output$is_continuous<-factor(regress_output$is_continuous,levels=c(TRUE,FALSE),labels=c("Imputed","Dichotomized"))
    }
    
    return(regress_output)
  }
  

  # Function to loop regression model and export model summary in MESA with incident outcomes
  mesa_strat_reg_loop<-function(data,covar,start,end,metab_is_cont,metab_is_complete,trait_original,trait_binary,trait_for_model,trait_as_predictor,stratifier){
    if (trait_for_model=="binary") {
      trait<-"TraitBinary"
    }  else if (trait_for_model=="original"){
      trait<-trait_original
    } else { trait<-"Trait" }
    
    dat<-data[,start:end]
    dat<-dat[,!colnames(dat)%in%c("ID","LAB_ID",trait)]
    
    covar<-covar[!covar%in%stratifier]
    group_by_value<-as.character(unique(data[,stratifier]))
    out_nvar<-ncol(dat)
    out_beta<-rep(NA,out_nvar)
    out_se<-rep(NA,out_nvar)
    out_pvalue<-rep(NA,out_nvar)
    out_nobs<-rep(NA,out_nvar)
    group<-rep(NA,out_nvar)
    regress_output<-data.frame()
    # if the original trait variable is binary or the user pick the binary version to use in the model
    trait_is_cont<-ifelse(length(unique(data[!is.na(data[,trait]),trait]))==2,FALSE,
                          ifelse(trait_for_model=="binary",FALSE,TRUE))
    
    
    for (j in seq_along(group_by_value)){
      newdata<-data[which(data[,stratifier]==as.character(group_by_value[j])),]
      dat<-newdata[,start:end]
      dat<-dat[,!colnames(dat)%in%c("ID","LAB_ID",trait)]
      for(i in 1:(ncol(dat))) {
        met_df<-cbind(dat[,colnames(dat)[i]],newdata[,c("ID",covar)])
        met_df_id<-met_df[complete.cases(met_df),"ID"]
        if (trait_as_predictor==T){
          if (metab_is_cont==T){
            model<- glm(as.formula(paste0(colnames(dat)[i],"~",trait,"+",paste(covar,collapse= "+"))),data=subset(data,ID%in%met_df_id))
          } else {
            model<- glm(as.formula(paste0(colnames(dat)[i],"~",trait,"+",paste(covar,collapse = "+"))),data=subset(data,ID%in%met_df_id),family=quasipoisson(link='log'))
          }
        }else {
          if (trait_is_cont==T){
            model<- glm(as.formula(paste0(trait,"~",colnames(dat)[i],"+",paste(covar,collapse= "+"))),data=subset(data,ID%in%met_df_id))
          }else{
            model<- glm(as.formula(paste0(trait,"~",colnames(dat)[i],"+",paste(covar,collapse = "+"))),data=subset(data,ID%in%met_df_id),family=quasipoisson(link='log'))
          }
        }
        
        Vcov <- vcov(model, useScale = FALSE)
        beta<- coef(model)
        se<- sqrt(diag(vcov(model, useScale = FALSE)))
        zval<- beta / se
        pval<- 2 * pnorm(abs(zval), lower.tail = FALSE)
        out_beta[i]=as.numeric(beta[2])
        out_se[i] = as.numeric(se[2])
        out_pvalue[i] = as.numeric(pval[2])
        out_nobs[i]=as.numeric(nobs(model))
        group[i]=group_by_value[j]
      }
      strata_output<-data.frame(metabolite=colnames(dat),
                                beta=out_beta,
                                se=out_se,
                                n=out_nobs,
                                p_val=out_pvalue,
                                p_val_Bonf=p.adjust(out_pvalue,method="bonferroni"),
                                p_val_fdr=p.adjust(out_pvalue,method="BH"),
                                strata=group)%>%
        mutate(
          pval_Bonf_neglog=-log10(p_val_Bonf),
          pval_fdr_neglog=-log10(p_val_fdr),
          sig_Bonf=ifelse(p_val_Bonf<=0.05,"Bonf<0.05","Not Sig"),
          sig_fdr=ifelse(p_val_fdr<=0.05,"FDR<0.05","Not Sig"),
          is_continuous=metab_is_cont
        )
      regress_output<-rbind(regress_output,strata_output)
    }
    identified_data<-regress_output%>%
      dplyr::filter(!stringr::str_detect(metabolite,"^X"))%>%
      group_by(strata)%>%
      dplyr::mutate(p_val_fdr_named=p.adjust(p_val,method="BH"),
                    sig_fdr_named=ifelse(p_val_fdr_named<=0.05,"FDR<0.05","Not Sig"))
    unidentified_data<-regress_output%>%
      dplyr::filter(stringr::str_detect(metabolite,"^X"))%>%
      dplyr::mutate(p_val_fdr_named=NA,
                    sig_fdr_named=NA)
    regress_output<-plyr::rbind.fill(identified_data,unidentified_data)
    if (metab_is_complete==T){
      regress_output$is_continuous<-"Complete-cases"
    } else {
      regress_output$is_continuous<-factor(regress_output$is_continuous,levels=c(TRUE,FALSE),labels=c("Imputed","Dichotomized"))
    }
    
    return(regress_output)
  }
  
  
  # Function to loop the interaction regression model and export model summary
  reg_loop_interaction<-function(data,covar,start,end,metab_is_cont,metab_is_complete,trait_original,trait_binary,trait_for_model,trait_as_predictor,interaction){
    if (trait_for_model=="binary") {
      trait<-"TraitBinary"
    }  else if (trait_for_model=="original"){
      trait<-trait_original
    } else { trait<-"Trait" }    
    dat<-data[,start:end]
    dat<-dat[,!colnames(dat)%in%c("LAB_ID","ID",trait)]
    out_nvar=ncol(dat)
    out_beta=rep(NA,out_nvar)
    out_se=rep(NA,out_nvar)
    out_pvalue=rep(NA,out_nvar)
    out_nobs=rep(NA,out_nvar)
    int_beta=rep(NA,out_nvar)
    int_se=rep(NA,out_nvar)
    int_pvalue=rep(NA,out_nvar)
    cross_beta=rep(NA,out_nvar)
    cross_se=rep(NA,out_nvar)
    cross_pvalue=rep(NA,out_nvar)
    
    # drop the interaction term from the covariate list
    covar2<-covar[!covar%in%interaction]

    
    # if the original trait variable is binary or the user pick the binary version to use in the model
    trait_is_cont=ifelse(length(unique(data[!is.na(data[,trait]),trait]))==2,FALSE,
                         ifelse(trait_for_model=="binary",FALSE,TRUE))
    
    for(i in 1:(ncol(dat))) {
      met_df<-cbind(data[,colnames(dat)[i]],data[,c("ID",covar)])
      met_df_id<-met_df[complete.cases(met_df),"ID"]
      if (trait_as_predictor==T){
        if (metab_is_cont==T){
          model<- glm(as.formula(paste0(colnames(dat)[i],"~",trait,"+",paste(covar2,collapse= "+"),"+",trait,"*",interaction)),data=subset(data,ID%in%met_df_id))
        } else {
          model<- glm(as.formula(paste0(colnames(dat)[i],"~",trait,"+",paste(covar2,collapse = "+"),"+",trait,"*",interaction)),data=subset(data,ID%in%met_df_id),family=quasipoisson(link='log'))
        }
      } else {
        if (trait_is_cont==T){
          model<- glm(as.formula(paste0(trait,"~",colnames(dat)[i],"+",paste(covar2,collapse= "+"),"+",colnames(dat)[i],"*",interaction)),data=subset(data,ID%in%met_df_id))
        } else {
          model<- glm(as.formula(paste0(trait,"~",colnames(dat)[i],"+",paste(covar2,collapse = "+"),"+",colnames(dat)[i],"*",interaction)),data=subset(data,ID%in%met_df_id),family=quasipoisson(link='log'))
        }
      }
      
      Vcov <- vcov(model, useScale = FALSE)
      beta<- coef(model)
      se<- sqrt(diag(vcov(model, useScale = FALSE)))
      zval<- beta / se
      pval<- 2 * pnorm(abs(zval), lower.tail = FALSE)
      out_beta[i]=as.numeric(beta[2])
      out_se[i] = round(as.numeric(se[2]),digits = 3)
      out_pvalue[i] = as.numeric(pval[2])
      out_nobs[i]=as.numeric(nobs(model))
      int_beta[i]=as.numeric(beta[length(beta)-1])
      int_se[i] = round(as.numeric(se[length(se)-1]),digits = 3)
      int_pvalue[i] = round(as.numeric(pval[length(pval)-1]),digits=9)
      cross_beta[i]=as.numeric(beta[length(beta)])
      cross_se[i] = round(as.numeric(se[length(se)]),digits = 3)
      cross_pvalue[i] = as.numeric(pval[length(pval)])
    }
    regress_output<-data.frame(metabolite=colnames(dat),
                               trait_beta=out_beta,
                               trait_se=out_se,
                               n=out_nobs,
                               trait_p_val=out_pvalue,
                               trait_p_val_Bonf=p.adjust(out_pvalue,method="bonferroni"),
                               trait_p_val_fdr=p.adjust(out_pvalue,method="BH"),
                               int_main_beta=int_beta,
                               int_main_se=int_se,
                               int_main_p_val=int_pvalue,
                               # int_main_p_val_Bonf=p.adjust(int_pvalue,method="bonferroni"),
                               # int_main_p_val_fdr=p.adjust(int_pvalue,method="BH"),
                               int_cross_beta=cross_beta,
                               int_cross_se=cross_se,
                               int_cross_p_val=cross_pvalue,
                               # int_cross_p_val_Bonf=p.adjust(cross_pvalue,method="bonferroni"),
                               # int_cross_p_val_fdr=p.adjust(cross_pvalue,method="BH"),
                               interaction_term=names(beta)[length(beta)-1],
                               cross_term=names(beta)[length(beta)]
                               
    )%>%
      dplyr::mutate(
        pval_Bonf_neglog=-log10(trait_p_val_Bonf),
        pval_fdr_neglog=-log10(trait_p_val_fdr),
        sig_Bonf=ifelse(trait_p_val_Bonf<=0.05,"Bonf<0.05","Not Sig"),
        sig_fdr=ifelse(trait_p_val_fdr<=0.05,"FDR<0.05","Not Sig"),
        is_continuous=metab_is_cont
      )
    identified_data<-regress_output%>%
      dplyr::filter(!stringr::str_detect(metabolite,"^X"))%>%
      dplyr::mutate(p_val_fdr_named=p.adjust(trait_p_val,method="BH"),
                    sig_fdr_named=ifelse(p_val_fdr_named<=0.05,"FDR<0.05","Not Sig"))
    unidentified_data<-regress_output%>%
      dplyr::filter(stringr::str_detect(metabolite,"^X"))%>%
      dplyr::mutate(p_val_fdr_named=NA,
                    sig_fdr_named=NA)
    regress_output<-rbind(identified_data,unidentified_data)
    if (metab_is_complete==T){
      regress_output$is_continuous<-"Complete-cases"
    } else {
      regress_output$is_continuous<-factor(regress_output$is_continuous,levels=c(TRUE,FALSE),labels=c("Imputed","Dichotomized"))
    }
    
    return(regress_output)
  }
  
  
  # Function to generate pvalue vs beta graphs (show  names for all significant metabolites)
  pval_beta_plot_label<-function(data,model_text,footnote_text,is_imputed,p_named_only,trait_label_original,trait_label_binary,trait_for_model){
    if (trait_for_model=="binary") {
      title_text<-trait_label_binary
    } else { title_text<-trait_label_original}
    footnote_text<-paste(footnote_text,collapse = ",")
    if (p_named_only==T){
      data$pval_fdr_neglog<--log10(data$p_val_fdr_named)
      pval_top15<-quantile(data$pval_fdr_neglog,probs = 0.985,na.rm = T)
      text_label_pval<-ifelse(pval_top15>-log10(0.05),pval_top15,-log10(0.05))
      data<-data[which(!is.na(data$p_val_fdr_named)),]
      
      if (is_imputed==T) {
        plot<-ggplot(data=data,aes(x=beta,y=pval_fdr_neglog))+
          geom_point(aes(color=sig_fdr_named,shape = is_continuous))+
          scale_color_manual(values=c("red","grey"),name="Significant")+
          scale_shape_manual(values = c(6, 20),name="Imputed or Dichotomized")+
          geom_hline(aes(yintercept=-log10(0.05)))+
          labs(x="Beta",y="-log10(p value)",title = paste(model_text,":\nScatter plot of beta coefficient and FDR adjusted P value\n(-log10) of metabolites concentrations\nand",title_text),caption=paste("metabolites with pvalue<",formatC(exp(-text_label_pval), format = "e", digits = 2),"annotated in the graph\nmodel adjusted for",tolower(footnote_text)))+
          geom_text_repel(data=subset(data, pval_fdr_neglog>-log10(0.05)),aes(label=metabolite),size=3)
      } else {
        plot<-ggplot(data=data,aes(x=beta,y=pval_fdr_neglog))+
          geom_point(aes(color=sig_fdr_named))+
          scale_color_manual(values=c("red","grey"),name="Significant")+
          geom_hline(aes(yintercept=-log10(0.05)))+
          labs(x="Beta",y="-log10(p value)",title = paste(model_text,":\nScatter plot of beta coefficient and FDR adjusted P value\n(-log10) of metabolites concentrations\nand",title_text),caption=paste("metabolites with pvalue<",formatC(exp(-text_label_pval), format = "e", digits = 2),"annotated in the graph\nmodel adjusted for",tolower(footnote_text)))+
          geom_text_repel(data=subset(data, pval_fdr_neglog>-log10(0.05)),aes(label=metabolite),size=3)
      } 
    } else {
      pval_top15<-quantile(data$pval_fdr_neglog,probs = 0.985,na.rm = T)
      text_label_pval<-ifelse(pval_top15>-log10(0.05),pval_top15,-log10(0.05))
      if (is_imputed==T) {
        plot<-ggplot(data=data,aes(x=beta,y=pval_fdr_neglog))+
          geom_point(aes(color=sig_fdr,shape = is_continuous))+
          scale_color_manual(values=c("red","grey"),name="Significant")+
          scale_shape_manual(values = c(6, 20),name="Imputed or Dichotomized")+
          geom_hline(aes(yintercept=-log10(0.05)))+
          labs(x="Beta",y="-log10(p value)",title = paste(model_text,":\nScatter plot of beta coefficient and FDR adjusted P value\n(-log10) of metabolites concentrations\nand",title_text),caption=paste("metabolites with pvalue<",formatC(exp(-text_label_pval), format = "e", digits = 2),"annotated in the graph\nmodel adjusted for",tolower(footnote_text)))+
          geom_text_repel(data=subset(data, pval_fdr_neglog>-log10(0.05)),aes(label=metabolite),size=3)
      } else {
        plot<-ggplot(data=data,aes(x=beta,y=pval_fdr_neglog))+
          geom_point(aes(color=sig_fdr))+
          scale_color_manual(values=c("red","grey"),name="Significant")+
          geom_hline(aes(yintercept=-log10(0.05)))+
          labs(x="Beta",y="-log10(p value)",title = paste(model_text,":\nScatter plot of beta coefficient and FDR adjusted P value\n(-log10) of metabolites concentrations\nand",title_text),caption=paste("metabolites with pvalue<",formatC(exp(-text_label_pval), format = "e", digits = 2),"annotated in the graph\nmodel adjusted for",tolower(footnote_text)))+
          geom_text_repel(data=subset(data, pval_fdr_neglog>-log10(0.05)),aes(label=metabolite),size=3)
      }
    }
    return(plot)
  }
  
  
  # Stratified plot (show names of all significant metaboiltes)
  pval_beta_strat_plot_label<-function(data,model_text,footnote_text,is_imputed,p_named_only,trait_label_original,trait_label_binary,trait_for_model,stratifier){
    if (trait_for_model=="binary") {
      title_text<-trait_label_binary
    } else { title_text<-trait_label_original
    } 
    footnote_text<-footnote_text[!footnote_text%in%stratifier]
    footnote_text<-paste(footnote_text,collapse = ",")
    if (p_named_only==T) {
      data$pval_fdr_neglog<--log10(data$p_val_fdr_named)
      pval_top15<-quantile(data$pval_fdr_neglog,probs = 0.985,na.rm = T)
      data$sig_fdr<-data$sig_fdr_named
    } else {pval_top15<-quantile(data$pval_fdr_neglog,probs = 0.985,na.rm = T)}
    # text_label_pval<-ifelse(pval_top15>-log10(0.05),pval_top15,-log10(0.05))
    text_label_pval<--log10(0.05)
    if (is_imputed==T) {
      plot<-ggplot(data=data,aes(x=beta,y=pval_fdr_neglog))+
        geom_point(aes(color=sig_fdr,shape = is_continuous))+
        scale_color_manual(values=c("red","grey"),name="Significant")+
        scale_shape_manual(values = c(6, 20),name="Imputed or Dichotomized")+
        geom_hline(aes(yintercept=-log10(0.05)))+
        labs(x="Beta",y="-log10(p value)",title = paste(model_text,":\nScatter plot of beta coefficient and FDR adjusted P value\n(-log10) of metabolites concentrations\nand",title_text),caption=paste("metabolites with pvalue<",formatC(exp(-text_label_pval), format = "e", digits = 2),"annotated in the graph\nmodel adjusted for",tolower(footnote_text)))+
        geom_text_repel(data=subset(data, pval_fdr_neglog>text_label_pval),aes(label=metabolite),size=3)
      plot2<-plot+facet_wrap(~strata)
    } else {
      plot<-ggplot(data=data,aes(x=beta,y=pval_fdr_neglog))+
        geom_point(aes(color=sig_fdr))+
        scale_color_manual(values=c("red","grey"),name="Significant")+
        geom_hline(aes(yintercept=-log10(0.05)))+
        labs(x="Beta",y="-log10(p value)",title = paste(model_text,":\nScatter plot of beta coefficient and FDR adjusted P value\n(-log10) of metabolites concentrations\nand",title_text),caption=paste("metabolites with pvalue<",formatC(exp(-text_label_pval), format = "e", digits = 2),"annotated in the graph\nmodel adjusted for",tolower(footnote_text)))+
        geom_text_repel(data=subset(data, pval_fdr_neglog>text_label_pval),aes(label=metabolite),size=3)
      plot2<-plot+facet_wrap(~strata)
    }
    return(plot2)
  }

  # extract model output of quartile regression of single exposure
  quartile_association_output<-function(model_list,dataset){
    qt2_beta<-rep(NA,length(model_list))
    qt2_se<-rep(NA,length(model_list))
    qt2_pvalue<-rep(NA,length(model_list))
    qt3_beta<-rep(NA,length(model_list))
    qt3_se<-rep(NA,length(model_list))
    qt3_pvalue<-rep(NA,length(model_list))
    qt4_beta<-rep(NA,length(model_list))
    qt4_se<-rep(NA,length(model_list))
    qt4_pvalue<-rep(NA,length(model_list))
    number<-1
    for (i in 1:length(model_list)){
      Vcov <- vcov(model_list[[i]], useScale = FALSE)
      beta<- coef(model_list[[i]])
      se<- sqrt(diag(vcov(model_list[[i]], useScale = FALSE)))
      zval<- beta / se
      pval<- 2 * pnorm(abs(zval), lower.tail = FALSE)
      qt2_beta[number]=as.numeric(beta[length(beta)-2])
      qt2_se[number] = as.numeric(se[length(beta)-2])
      qt2_pvalue[number] = as.numeric(pval[length(beta)-2])
      qt3_beta[number]=as.numeric(beta[length(beta)-1])
      qt3_se[number] = as.numeric(se[length(beta)-1])
      qt3_pvalue[number] = as.numeric(pval[length(beta)-1])
      qt4_beta[number]=as.numeric(beta[length(beta)])
      qt4_se[number] = as.numeric(se[length(beta)])
      qt4_pvalue[number] = as.numeric(pval[length(beta)])
      number<-number+1
    }
    qt2_output<-data.frame(beta=qt2_beta,se=qt2_se,pval=qt2_pvalue,Quartile=rep("2nd",length(model_list)))
    qt3_output<-data.frame(beta=qt3_beta,se=qt3_se,pval=qt3_pvalue,Quartile=rep("3rd",length(model_list)))        
    qt4_output<-data.frame(beta=qt4_beta,se=qt4_se,pval=qt4_pvalue,Quartile=rep("4th",length(model_list)))        
    qt_output<-rbind(qt2_output,qt3_output,qt4_output)%>%
      dplyr::mutate(OR=exp(beta),
                    lowerCI=exp(beta-1.96*se),
                    upperCI=exp(beta+1.96*se),
                    p_label=dplyr::case_when(
                      pval > 0.05 ~ "",
                      pval > 0.01 ~ "*",
                      pval > 0.001 ~ "**",
                      !is.na(pval) ~ "***",
                      TRUE ~ NA_character_
                    ),
                    model=rep(c("Model 1","Model 2","Model 3"),length(model_list)),
                    dataset=rep(paste(dataset),3*length(model_list)))
    
  }
  
  # extract model output of quartile regression of multiple exposures
  quartile_association_exposures_output<-function(model_list,dataset){
    qt2_beta<-rep(NA,length(model_list))
    qt2_se<-rep(NA,length(model_list))
    qt2_pvalue<-rep(NA,length(model_list))
    qt3_beta<-rep(NA,length(model_list))
    qt3_se<-rep(NA,length(model_list))
    qt3_pvalue<-rep(NA,length(model_list))
    qt4_beta<-rep(NA,length(model_list))
    qt4_se<-rep(NA,length(model_list))
    qt4_pvalue<-rep(NA,length(model_list))
    exposure<-rep(NA,length(model_list))
    out_nobs<-rep(NA,length(model_list))
    
    number<-1
    for (i in 1:length(model_list)){
      Vcov <- vcov(model_list[[i]], useScale = FALSE)
      beta<- coef(model_list[[i]])
      se<- sqrt(diag(vcov(model_list[[i]], useScale = FALSE)))
      zval<- beta / se
      pval<- 2 * pnorm(abs(zval), lower.tail = FALSE)
      qt2_beta[number]=as.numeric(beta[length(beta)-2])
      qt2_se[number] = as.numeric(se[length(beta)-2])
      qt2_pvalue[number] = as.numeric(pval[length(beta)-2])
      qt3_beta[number]=as.numeric(beta[length(beta)-1])
      qt3_se[number] = as.numeric(se[length(beta)-1])
      qt3_pvalue[number] = as.numeric(pval[length(beta)-1])
      qt4_beta[number]=as.numeric(beta[length(beta)])
      qt4_se[number] = as.numeric(se[length(beta)])
      qt4_pvalue[number] = as.numeric(pval[length(beta)])
      exposure[number]=substr(names(beta)[length(names(beta))],1,nchar(names(beta)[length(names(beta))])-1)
      out_nobs[number]=as.numeric(nobs(model_list[[i]]))
      
      number<-number+1
    }
    qt2_output<-data.frame(beta=qt2_beta,se=qt2_se,pval=qt2_pvalue,Quartile=rep("2nd",length(model_list)),Exposure=exposure,n=out_nobs)
    qt3_output<-data.frame(beta=qt3_beta,se=qt3_se,pval=qt3_pvalue,Quartile=rep("3rd",length(model_list)),Exposure=exposure,n=out_nobs)        
    qt4_output<-data.frame(beta=qt4_beta,se=qt4_se,pval=qt4_pvalue,Quartile=rep("4th",length(model_list)),Exposure=exposure,n=out_nobs)        
    qt_output<-rbind(qt2_output,qt3_output,qt4_output)%>%
      dplyr::mutate(OR=exp(beta),
                    lowerCI=exp(beta-1.96*se),
                    upperCI=exp(beta+1.96*se),
                    p_label=dplyr::case_when(
                      pval > 0.05 ~ "",
                      pval > 0.01 ~ "*",
                      pval > 0.001 ~ "**",
                      !is.na(pval) ~ "***",
                      TRUE ~ NA_character_
                    ),
                    model=rep(c("Model 1","Model 2","Model 3"),length(model_list)),
                    dataset=rep(paste(dataset),3*length(model_list)))
    
  }
  # extract model output of quartile regression of multiple exposures and models
  quartile_association_models_output<-function(model_list,dataset){
    qt2_beta<-rep(NA,length(model_list))
    qt2_se<-rep(NA,length(model_list))
    qt2_pvalue<-rep(NA,length(model_list))
    qt3_beta<-rep(NA,length(model_list))
    qt3_se<-rep(NA,length(model_list))
    qt3_pvalue<-rep(NA,length(model_list))
    qt4_beta<-rep(NA,length(model_list))
    qt4_se<-rep(NA,length(model_list))
    qt4_pvalue<-rep(NA,length(model_list))
    exposure<-rep(NA,length(model_list))
    covariates<-rep(NA,length(model_list))
    number<-1
    for (i in 1:length(model_list)){
      Vcov <- vcov(model_list[[i]], useScale = FALSE)
      beta<- coef(model_list[[i]])
      se<- sqrt(diag(vcov(model_list[[i]], useScale = FALSE)))
      zval<- beta / se
      pval<- 2 * pnorm(abs(zval), lower.tail = FALSE)
      qt2_beta[number]=as.numeric(beta[length(beta)-2])
      qt2_se[number] = as.numeric(se[length(beta)-2])
      qt2_pvalue[number] = as.numeric(pval[length(beta)-2])
      qt3_beta[number]=as.numeric(beta[length(beta)-1])
      qt3_se[number] = as.numeric(se[length(beta)-1])
      qt3_pvalue[number] = as.numeric(pval[length(beta)-1])
      qt4_beta[number]=as.numeric(beta[length(beta)])
      qt4_se[number] = as.numeric(se[length(beta)])
      qt4_pvalue[number] = as.numeric(pval[length(beta)])
      exposure[number]=substr(names(beta)[length(names(beta))],1,nchar(names(beta)[length(names(beta))])-1)
      predictors=paste(model_list[[i]]$formula[3])
      if (grepl("+",predictors)==T){
        covariates[number]=stringr::str_extract(predictors,".*(?=[+])")
        # pattern <- "^(.*?)\\+(.*?\\+[^+]+)$"
        # covariates[number]=stringr::str_extract(predictors,"^[^+]+(\\+[^+]+){1}$")
      } else {
        covariates[number]="None"
      }
      number<-number+1
    }
    qt2_output<-data.frame(beta=qt2_beta,se=qt2_se,pval=qt2_pvalue,Quartile=rep("2nd",length(model_list)),Exposure=exposure,Covariate=covariates)
    qt3_output<-data.frame(beta=qt3_beta,se=qt3_se,pval=qt3_pvalue,Quartile=rep("3rd",length(model_list)),Exposure=exposure,Covariate=covariates)        
    qt4_output<-data.frame(beta=qt4_beta,se=qt4_se,pval=qt4_pvalue,Quartile=rep("4th",length(model_list)),Exposure=exposure,Covariate=covariates)        
    qt_output<-rbind(qt2_output,qt3_output,qt4_output)%>%
      dplyr::mutate(OR=exp(beta),
                    lowerCI=exp(beta-1.96*se),
                    upperCI=exp(beta+1.96*se),
                    p_label=dplyr::case_when(
                      pval > 0.05 ~ "",
                      pval > 0.01 ~ "*",
                      pval > 0.001 ~ "**",
                      !is.na(pval) ~ "***",
                      TRUE ~ NA_character_
                    ),
                    # model=rep(c("Model 1","Model 2","Model 3"),length(model_list)),
                    dataset=paste(dataset))
    model_names<-paste0("model_",seq_along(unique(qt_output$Covariate)))
    qt_output<-qt_output%>%
      dplyr::mutate(model=rep(model_names,(nrow(qt_output)/length(model_names))))
  }  
  # extract model output of quartile regression of multiple exposures and models
  quartile_association_models_output_offset<-function(model_list,dataset){
    qt2_beta<-rep(NA,length(model_list))
    qt2_se<-rep(NA,length(model_list))
    qt2_pvalue<-rep(NA,length(model_list))
    qt3_beta<-rep(NA,length(model_list))
    qt3_se<-rep(NA,length(model_list))
    qt3_pvalue<-rep(NA,length(model_list))
    qt4_beta<-rep(NA,length(model_list))
    qt4_se<-rep(NA,length(model_list))
    qt4_pvalue<-rep(NA,length(model_list))
    exposure<-rep(NA,length(model_list))
    covariates<-rep(NA,length(model_list))
    out_nobs<-rep(NA,length(model_list))
    
    number<-1
    for (i in 1:length(model_list)){
      Vcov <- vcov(model_list[[i]], useScale = FALSE)
      beta<- coef(model_list[[i]])
      se<- sqrt(diag(vcov(model_list[[i]], useScale = FALSE)))
      zval<- beta / se
      pval<- 2 * pnorm(abs(zval), lower.tail = FALSE)
      qt2_beta[number]=as.numeric(beta[length(beta)-2])
      qt2_se[number] = as.numeric(se[length(beta)-2])
      qt2_pvalue[number] = as.numeric(pval[length(beta)-2])
      qt3_beta[number]=as.numeric(beta[length(beta)-1])
      qt3_se[number] = as.numeric(se[length(beta)-1])
      qt3_pvalue[number] = as.numeric(pval[length(beta)-1])
      qt4_beta[number]=as.numeric(beta[length(beta)])
      qt4_se[number] = as.numeric(se[length(beta)])
      qt4_pvalue[number] = as.numeric(pval[length(beta)])
      exposure[number]=substr(names(beta)[length(names(beta))],1,nchar(names(beta)[length(names(beta))])-1)
      predictors=paste(model_list[[i]]$formula[3])
      out_nobs[number]=as.numeric(nobs(model_list[[i]]))
      
        # covariates[number]=stringr::str_extract(predictors,".*(?=[+])")
        pattern <- "^(.*?)\\+(.*?\\+[^+]+)$"
        covariates[number]=sub(pattern, "\\1", predictors)
        if (covariates[number] == predictors) {
          covariates[number] <- "None"
        }
      number<-number+1
    }
    qt2_output<-data.frame(beta=qt2_beta,se=qt2_se,pval=qt2_pvalue,Quartile=rep("2nd",length(model_list)),Exposure=exposure,Covariate=covariates,n=out_nobs)
    qt3_output<-data.frame(beta=qt3_beta,se=qt3_se,pval=qt3_pvalue,Quartile=rep("3rd",length(model_list)),Exposure=exposure,Covariate=covariates,n=out_nobs)        
    qt4_output<-data.frame(beta=qt4_beta,se=qt4_se,pval=qt4_pvalue,Quartile=rep("4th",length(model_list)),Exposure=exposure,Covariate=covariates,n=out_nobs)        
    qt_output<-rbind(qt2_output,qt3_output,qt4_output)%>%
      dplyr::mutate(OR=exp(beta),
                    lowerCI=exp(beta-1.96*se),
                    upperCI=exp(beta+1.96*se),
                    p_label=dplyr::case_when(
                      pval > 0.05 ~ "",
                      pval > 0.01 ~ "*",
                      pval > 0.001 ~ "**",
                      !is.na(pval) ~ "***",
                      TRUE ~ NA_character_
                    ),
                    # model=rep(c("Model 1","Model 2","Model 3"),length(model_list)),
                    dataset=paste(dataset))
    model_names<-paste0("model_",seq_along(unique(qt_output$Covariate)))
    qt_output<-qt_output%>%
      dplyr::mutate(model=rep(model_names,(nrow(qt_output)/length(model_names))))
  }   
  
  ################################
  # Loop over association analysis gender stratified model using the gender-combined MRS for incident htn
  loop_association_analysis<-function(exposure_list,covariate_list,survey_design,outcome_var,dataset){
    output_mdl_list<-list()
    number=1
    model<-list()
    for (j in seq_along(exposure_list)){
      # print(paste("j=",j))
      for (k in seq_along(covariate_list)){
        # print(paste("k=",k))
        if (!is.na(covariate_list[k])){
          if (!is.na(exposure_list[j])){
            if (!is.null(survey_design)){
              if (length(unique(survey_design$variables[,outcome_var]))==2){
                model<-survey::svyglm(as.formula(paste0(outcome_var,"~",paste0(covariate_list[[k]],collapse="+"),"+",exposure_list[j])),design = survey_design, family=quasipoisson(link='log'))
              } else {
                model<-survey::svyglm(as.formula(paste0(outcome_var,"~",paste0(covariate_list[[k]],collapse="+"),"+",exposure_list[j])),design = survey_design)
              }
            } else {
              if (length(unique(dataset[,outcome_var]))==2) {
                model<-glm(as.formula(paste0(outcome_var,"~",paste0(covariate_list[[k]],collapse="+"),"+",exposure_list[j])),data = dataset, family=binomial(link='logit'))
              } else {
                model<-lm(as.formula(paste0(outcome_var,"~",paste0(covariate_list[[k]],collapse="+"),"+",exposure_list[j])),data = dataset)
              }
            }
          } else {
            if (!is.null(survey_design)){
              if (length(unique(survey_design$variables[,outcome_var]))==2){
                model<-survey::svyglm(as.formula(paste0(outcome_var,"~",paste0(covariate_list[[k]],collapse="+"))),design = survey_design, family=quasipoisson(link='log'))
              } else {
                model<-survey::svyglm(as.formula(paste0(outcome_var,"~",paste0(covariate_list[[k]],collapse="+"))),design = survey_design)
              }
            } else {
              if (length(unique(dataset[,outcome_var]))==2) {
                model<-glm(as.formula(paste0(outcome_var,"~",paste0(covariate_list[[k]],collapse="+"))),data = dataset, family=binomial(link='logit'))
              } else {
                model<-lm(as.formula(paste0(outcome_var,"~",paste0(covariate_list[[k]],collapse="+"))),data = dataset)
              }
            }
          }
        } else {
          if (!is.na(exposure_list[j])){
            if (!is.null(survey_design)){
              if (length(unique(survey_design$variables[,outcome_var]))==2) {
                model<-survey::svyglm(as.formula(paste0(outcome_var,"~",exposure_list[j])),design = survey_design, family=quasipoisson(link='log'))
              } else {
                model<-survey::svyglm(as.formula(paste0(outcome_var,"~",exposure_list[j])),design = survey_design)
              }
            } else {
              if (length(unique(dataset[,outcome_var]))==2){
                model<-glm(as.formula(paste0(outcome_var,"~",exposure_list[j])),data = dataset, family=binomial(link='logit'))
              } else (
                model<-lm(as.formula(paste0(outcome_var,"~",exposure_list[j])),data = dataset)
                
              )
            }
          } 
        }
        # print(paste("k=",k))
        output_mdl_list[[number]]<-model
        number=number+1
      }
      # print(paste("j=",j))
    }
    return(output_mdl_list)
  }
  
  # extract model summary in the format of a data frame after loop_association_analysis
  extract_mdl_output<-function(model_output_list,model_index_matrix,model_outcome,model_strata) {
    out_beta<-rep(NA,length(model_output_list))
    out_se<-rep(NA,length(model_output_list))
    out_pvalue<-rep(NA,length(model_output_list))
    out_nobs<-rep(NA,length(model_output_list))
    out_index<-rep(NA,length(model_output_list))
    out_trait<-rep(NA,length(model_output_list))
    out_model<-rep(NA,length(model_output_list))
    for (i in seq_along(model_output_list)){
      # print(i)
      Vcov <- vcov(model_output_list[[i]], useScale = FALSE)
      beta<- coef(model_output_list[[i]])
      se<- sqrt(diag(vcov(model_output_list[[i]], useScale = FALSE)))
      zval<- beta / se
      pval<- 2 * pnorm(abs(zval), lower.tail = FALSE)
      out_beta[i]=as.numeric(beta[length(beta)])
      out_se[i] = round(as.numeric(se[length(se)]),digits = 3)
      out_pvalue[i] = as.numeric(pval[length(pval)])
      out_nobs[i]=as.numeric(nobs(model_output_list[[i]]))
      out_index[i]=paste(i)
      out_trait[i]=colnames(model_index_matrix)[ceiling(i/nrow(model_index_matrix))]
      if (i%%nrow(model_index_matrix!=0)){
        out_model[i]=rownames(model_index_matrix)[i%%nrow(model_index_matrix)]
      } else {
        out_model[i]=rownames(model_index_matrix)[nrow(model_index_matrix)]
      }
      # print(i)
    }
    extracted_output<-data.frame(trait=out_trait,
                                 model=out_model,
                                 index=out_index,
                                 beta=out_beta,
                                 se=out_se,
                                 n=out_nobs,
                                 p_val=out_pvalue,
                                 outcome=model_outcome,
                                 strata=model_strata
    )
  }
  
  # extract model summary in the format of a data frame after loop_association_analysis
  extract_mdl_sig_output<-function(model_output_list,model_index_matrix,model_outcome,model_strata) {
    out_beta<-rep(NA,length(model_output_list))
    out_se<-rep(NA,length(model_output_list))
    out_pvalue<-rep(NA,length(model_output_list))
    out_nobs<-rep(NA,length(model_output_list))
    out_index<-rep(NA,length(model_output_list))
    out_trait<-rep(NA,length(model_output_list))
    out_model<-rep(NA,length(model_output_list))
    for (i in seq_along(model_output_list)){
      # print(i)
      Vcov <- vcov(model_output_list[[i]], useScale = FALSE)
      beta<- coef(model_output_list[[i]])
      se<- sqrt(diag(vcov(model_output_list[[i]], useScale = FALSE)))
      zval<- beta / se
      pval<- 2 * pnorm(abs(zval), lower.tail = FALSE)
      out_beta[i]=as.numeric(beta[length(beta)])
      out_se[i] = round(as.numeric(se[length(se)]),digits = 3)
      out_pvalue[i] = as.numeric(pval[length(pval)])
      out_nobs[i]=as.numeric(nobs(model_output_list[[i]]))
      out_index[i]=paste(i)
      out_trait[i]=colnames(model_index_matrix)[ceiling(i/nrow(model_index_matrix))]
      if (i%%nrow(model_index_matrix!=0)){
        out_model[i]=rownames(model_index_matrix)[i%%nrow(model_index_matrix)]
      } else {
        out_model[i]=rownames(model_index_matrix)[nrow(model_index_matrix)]
      }
      # print(i)
    }
    
    extracted_output<-data.frame(trait=out_trait,
                                 model=out_model,
                                 index=out_index,
                                 beta=out_beta,
                                 or=exp(out_beta),
                                 se=out_se,
                                 lower95=out_beta-1.96*out_se,
                                 upper95=out_beta+1.96*out_se,
                                 lowerCI=exp(out_beta-1.96*out_se),
                                 upperCI=exp(out_beta+1.96*out_se),
                                 n=out_nobs,
                                 p_val=out_pvalue,
                                 outcome=model_outcome,
                                 strata=model_strata
    )%>%
      dplyr::mutate(p_label=dplyr::case_when(
        p_val > 0.05 ~ "",
        p_val > 0.01 ~ "*",
        p_val > 0.001 ~ "**",
        !is.na(p_val) ~ "***",
        TRUE ~ NA_character_
      ))
  }
  
  # Loop over association analysis gender stratified model using the gender-combined MRS for incident htn
  loop_association_analysis_offset<-function(exposure_list,covariate_list,survey_design,outcome_var,dataset,offset_var){
    output_mdl_list<-list()
    number=1
    model<-list()
    for (j in seq_along(exposure_list)){
      # print(paste("j=",j))
      for (k in seq_along(covariate_list)){
        # print(paste("k=",k))
        if (!is.na(covariate_list[k])){
          if (!is.na(exposure_list[j])){
            if (!is.null(survey_design)){
              if (length(unique(survey_design$variables[,outcome_var]))==2){
                model<-survey::svyglm(as.formula(paste0(outcome_var,"~",paste0(covariate_list[[k]],collapse="+"),"+",exposure_list[j],"+offset(log(",offset_var,"))")),design = survey_design, family=quasipoisson(link='log'))
              } else {
                model<-survey::svyglm(as.formula(paste0(outcome_var,"~",paste0(covariate_list[[k]],collapse="+"),"+",exposure_list[j],"+offset(log(",offset_var,"))")),design = survey_design)
              }
            } else {
              if (length(unique(dataset[,outcome_var]))==2) {
                model<-glm(as.formula(paste0(outcome_var,"~",paste0(covariate_list[[k]],collapse="+"),"+",exposure_list[j])),data = dataset, family=binomial(link='logit'))
              } else {
                model<-lm(as.formula(paste0(outcome_var,"~",paste0(covariate_list[[k]],collapse="+"),"+",exposure_list[j])),data = dataset)
              }
            }
          } else {
            if (!is.null(survey_design)){
              if (length(unique(survey_design$variables[,outcome_var]))==2){
                model<-survey::svyglm(as.formula(paste0(outcome_var,"~",paste0(covariate_list[[k]],collapse="+"),"+offset(log(",offset_var,"))")),design = survey_design, family=quasipoisson(link='log'))
              } else {
                model<-survey::svyglm(as.formula(paste0(outcome_var,"~",paste0(covariate_list[[k]],collapse="+"),"+offset(log(",offset_var,"))")),design = survey_design)
              }
            } else {
              if (length(unique(dataset[,outcome_var]))==2) {
                model<-glm(as.formula(paste0(outcome_var,"~",paste0(covariate_list[[k]],collapse="+"))),data = dataset, family=binomial(link='logit'))
              } else {
                model<-lm(as.formula(paste0(outcome_var,"~",paste0(covariate_list[[k]],collapse="+"))),data = dataset)
              }
            }
          }
        } else {
          if (!is.na(exposure_list[j])){
            if (!is.null(survey_design)){
              if (length(unique(survey_design$variables[,outcome_var]))==2) {
                model<-survey::svyglm(as.formula(paste0(outcome_var,"~",exposure_list[j],"+offset(log(",offset_var,"))")),design = survey_design, family=quasipoisson(link='log'))
              } else {
                model<-survey::svyglm(as.formula(paste0(outcome_var,"~",exposure_list[j],"+offset(log(",offset_var,"))")),design = survey_design)
              }
            } else {
              if (length(unique(dataset[,outcome_var]))==2){
                model<-glm(as.formula(paste0(outcome_var,"~",exposure_list[j])),data = dataset, family=binomial(link='logit'))
              } else (
                model<-lm(as.formula(paste0(outcome_var,"~",exposure_list[j])),data = dataset)
                
              )
            }
          } 
        }
        # print(paste("k=",k))
        output_mdl_list[[number]]<-model
        number=number+1
      }
      # print(paste("j=",j))
    }
    return(output_mdl_list)
  }
  
  # extract model summary in the format of a data frame after loop_association_analysis
  extract_mdl_output<-function(model_output_list,model_index_matrix,model_outcome,model_strata) {
    out_beta<-rep(NA,length(model_output_list))
    out_se<-rep(NA,length(model_output_list))
    out_pvalue<-rep(NA,length(model_output_list))
    out_nobs<-rep(NA,length(model_output_list))
    out_index<-rep(NA,length(model_output_list))
    out_trait<-rep(NA,length(model_output_list))
    out_model<-rep(NA,length(model_output_list))
    for (i in seq_along(model_output_list)){
      # print(i)
      Vcov <- vcov(model_output_list[[i]], useScale = FALSE)
      beta<- coef(model_output_list[[i]])
      se<- sqrt(diag(vcov(model_output_list[[i]], useScale = FALSE)))
      zval<- beta / se
      pval<- 2 * pnorm(abs(zval), lower.tail = FALSE)
      out_beta[i]=as.numeric(beta[length(beta)])
      out_se[i] = round(as.numeric(se[length(se)]),digits = 3)
      out_pvalue[i] = as.numeric(pval[length(pval)])
      out_nobs[i]=as.numeric(nobs(model_output_list[[i]]))
      out_index[i]=paste(i)
      out_trait[i]=colnames(model_index_matrix)[ceiling(i/nrow(model_index_matrix))]
      if (i%%nrow(model_index_matrix!=0)){
        out_model[i]=rownames(model_index_matrix)[i%%nrow(model_index_matrix)]
      } else {
        out_model[i]=rownames(model_index_matrix)[nrow(model_index_matrix)]
      }
      # print(i)
    }
    extracted_output<-data.frame(trait=out_trait,
                                 model=out_model,
                                 index=out_index,
                                 beta=out_beta,
                                 se=out_se,
                                 n=out_nobs,
                                 p_val=out_pvalue,
                                 outcome=model_outcome,
                                 strata=model_strata
    )
  }
  
  # extract model summary in the format of a data frame after loop_association_analysis
  extract_mdl_sig_output<-function(model_output_list,model_index_matrix,model_outcome,model_strata) {
    out_beta<-rep(NA,length(model_output_list))
    out_se<-rep(NA,length(model_output_list))
    out_pvalue<-rep(NA,length(model_output_list))
    out_nobs<-rep(NA,length(model_output_list))
    out_index<-rep(NA,length(model_output_list))
    out_trait<-rep(NA,length(model_output_list))
    out_model<-rep(NA,length(model_output_list))
    for (i in seq_along(model_output_list)){
      # print(i)
      Vcov <- vcov(model_output_list[[i]], useScale = FALSE)
      beta<- coef(model_output_list[[i]])
      se<- sqrt(diag(vcov(model_output_list[[i]], useScale = FALSE)))
      zval<- beta / se
      pval<- 2 * pnorm(abs(zval), lower.tail = FALSE)
      out_beta[i]=as.numeric(beta[length(beta)])
      out_se[i] = round(as.numeric(se[length(se)]),digits = 3)
      out_pvalue[i] = as.numeric(pval[length(pval)])
      out_nobs[i]=as.numeric(nobs(model_output_list[[i]]))
      out_index[i]=paste(i)
      out_trait[i]=colnames(model_index_matrix)[ceiling(i/nrow(model_index_matrix))]
      if (i%%nrow(model_index_matrix!=0)){
        out_model[i]=rownames(model_index_matrix)[i%%nrow(model_index_matrix)]
      } else {
        out_model[i]=rownames(model_index_matrix)[nrow(model_index_matrix)]
      }
      # print(i)
    }
    
    extracted_output<-data.frame(trait=out_trait,
                                 model=out_model,
                                 index=out_index,
                                 beta=out_beta,
                                 or=exp(out_beta),
                                 se=out_se,
                                 lower95=out_beta-1.96*out_se,
                                 upper95=out_beta+1.96*out_se,
                                 lowerCI=exp(out_beta-1.96*out_se),
                                 upperCI=exp(out_beta+1.96*out_se),
                                 n=out_nobs,
                                 p_val=out_pvalue,
                                 outcome=model_outcome,
                                 strata=model_strata
    )%>%
      dplyr::mutate(p_label=dplyr::case_when(
        p_val > 0.05 ~ "",
        p_val > 0.01 ~ "*",
        p_val > 0.001 ~ "**",
        !is.na(p_val) ~ "***",
        TRUE ~ NA_character_
      ))
  }
  

  # can't use the pool function from the MI package
  # instead write a function that fit a model on the 5 datasets, and applies Rubin's rule.
  # for simplicity, assume insomnia is outcome, assume 5 imputations (from Tamar Sofer)
  
  rubins_rule <- function(ests, ses,round_digit = 2){
    vars <- ses^2
    m <- length(ses)
    
    pooled_est <- mean(ests)
    within_imp_var <- mean(vars)
    between_imp_var <- sum((vars - within_imp_var)^2)/(m -1)
    
    pooled_var <- within_imp_var + (1+1/m)*between_imp_var
    pooled_se <- sqrt(pooled_var)
    
    F_test_stat <- pooled_est^2/pooled_var
    F_df1 <- 1
    r <- (1 + 1/m)*between_imp_var/pooled_est
    F_df2 <- (m - 1)*(1 + 1/r)^2
    
    F_pval <- pf(F_test_stat, F_df1, F_df2, lower.tail = FALSE)
    
    CI <- paste0("(", round(pooled_est - 1.96*pooled_se, round_digit), ",", round(pooled_est + 1.96*pooled_se, round_digit), ")")
    
    return(data.frame(est = pooled_est, CI = CI,SE=pooled_se, pval = F_pval))
    
  }
# use rubin's rule to pool correlation coefficients
  rubins_rho <- function(rho_vec) {
    n_vec<-length(rho_vec)
    # Compute within-imputation variances
    s2_vec <- (1 - rho_vec^2) / (n_vec - 2)
    
    # Compute between-imputation variance
    B <- sum((n_vec - 1) / sum(n_vec - 1) * (rho_vec - mean(rho_vec))^2)
    
    # Compute combined estimate and standard error
    rho_bar <- mean(rho_vec)
    s2_w <- mean(s2_vec)
    s2_b <- B / sum(n_vec - 1)
    s2_comb <- s2_w + s2_b
    se_comb <- sqrt(s2_comb)
    
    # Compute t-statistic and p-value
    t_stat <- rho_bar / se_comb
    p_value <- 2 * (1 - pt(abs(t_stat), df = sum(n_vec) - 2))
    
    # Return results
    result <- data.frame(rho_bar = rho_bar,
                   se_comb = se_comb,
                   t_stat = t_stat,
                   p_value = p_value)
    return(result)
  }
  # use rubin's rule to pool intercept
  rubin_intercept <- function(b_vec, se_vec) {
    n_vec<-length(b_vec)
    # Compute within-imputation variances
    s2_vec <- se_vec^2
    
    # Compute between-imputation variance
    B <- sum((n_vec - 1) / sum(n_vec - 1) * (b_vec - mean(b_vec))^2)
    
    # Compute combined estimate and standard error
    b_bar <- mean(b_vec)
    s2_w <- mean(s2_vec)
    s2_b <- B / sum(n_vec - 1)
    s2_comb <- s2_w + s2_b
    se_comb <- sqrt(s2_comb)
    
    # Compute t-statistic and p-value
    t_stat <- b_bar / se_comb
    p_value <- 2 * (1 - pt(abs(t_stat), df = sum(n_vec) - 2))
    
    # Return results
    result <- list(mesor_bar = b_bar,
                   se_comb = se_comb,
                   t_stat = t_stat,
                   p_value = p_value)
    return(result)
  }
  
  # function to calculate matrix sd while excluding some columns  
  cal_matrix_sd<-function(dataset,col_excluded){
    dataset_std<-diag(sqrt(diag(cov(dataset[,!colnames(dataset)%in%col_excluded],use="pairwise.complete.obs"))))
  }
  # function to unnest a list of columns, similar to unnest_wider function from tidyr package
  unnest_wider_custom <- function(data, col_name, level = 1) {
    unnested <- lapply(data[[col_name]], function(x) unlist(x, recursive = FALSE, use.names = FALSE)[[level]])
    new_data <- cbind(data, do.call("rbind", unnested))
    new_data[[col_name]] <- NULL
    return(new_data)
  }
  
  
# function to convert long format to wide format table for complexheatmap  
# strata input is the column name of the variable for stratifying results (e.g., model, strata)
  
  convert_df_complexheatmap<-function(dataset,param_seq,param,strata,included_metab){
    outcome_list<-unique(dataset$trait)
   strata_list<-unique(dataset[,strata])
    combine_ind<-list()
    for (j in seq_along(strata_list)){
      param_list<-data.frame()
      for (i in seq_along(outcome_list)){
        all_metab<-dataset[which(dataset$trait==outcome_list[i]&dataset[,strata]==strata_list[j]),param_seq]
        param_list_temp<-all_metab[which(all_metab$metabolite%in%included_metab),c("metabolite",param)]
        colnames(param_list_temp)[2]<-outcome_list[i]
        if (i==1){
          param_list<-param_list_temp
        } else {
          param_list<-merge(param_list,param_list_temp,by="metabolite")
        }
      }
      param_list<-merge(param_list,dataset[which(dataset$trait==outcome_list[i]&dataset[,strata]==strata_list[j]),c("metabolite","super_pathway","sub_pathway")],by="metabolite",all.x=T)%>%
        data.table::setorder(.,super_pathway,sub_pathway)%>%
        dplyr::mutate(super_pathway=factor(super_pathway))%>%
        tibble::remove_rownames()%>%
        tibble::column_to_rownames(., var = "metabolite")%>%
      data.table::setorder(.,super_pathway,sub_pathway)
      combine_ind[[j]]<-param_list
      names(combine_ind)[j]<-strata_list[j]
      
    }
    return(combine_ind)
  }

  # generalize the convert function
  convert_df_complexheatmap_general<-function(dataset,param_seq,param,strata,included_x,heatmap_x,heatmap_y){
    outcome_list<-unlist(unique(dataset[,heatmap_y]))
    strata_list<-unlist(unique(dataset[,strata]))
    combine_ind<-list()
    for (j in seq_along(strata_list)){
      # print(paste("j=",j))
      param_list<-data.frame()
      for (i in seq_along(outcome_list)){
        # print(paste("i=",i))
        all_metab<-dataset[which(dataset[,heatmap_y]==outcome_list[i]&dataset[,strata]==strata_list[j]),param_seq]
        param_list_temp<-all_metab[which(all_metab[[heatmap_x]]%in%included_x),c(heatmap_x,param)]
        colnames(param_list_temp)[2]<-outcome_list[i]
        if (i==1){
          param_list<-param_list_temp
        } else {
          param_list<-merge(param_list,param_list_temp,by=paste(heatmap_x))
        }
        # print(paste("i=",i))
      }
      param_list <- replace(param_list, is.na(param_list), "unknown")
      param_list<-merge(param_list,dataset[which(dataset[,heatmap_y]==outcome_list[i]&dataset[,strata]==strata_list[j]),unique(c(heatmap_x,"super_pathway","sub_pathway"))],by=paste(heatmap_x),all.x=T)%>%
        data.table::setorder(.,super_pathway,sub_pathway)%>%
        dplyr::mutate(super_pathway=factor(super_pathway))%>%
        tibble::remove_rownames()%>%
        tibble::column_to_rownames(., var = paste(heatmap_x))
      if (heatmap_x=="sub_pathway"){
        param_list<-param_list%>%data.table::setorder(.,super_pathway)
      } else {
        param_list<-param_list%>%data.table::setorder(.,super_pathway,sub_pathway)
      }
      combine_ind[[j]]<-param_list
      names(combine_ind)[j]<-strata_list[j]
      # print("j=",j)
    }
    return(combine_ind)
  }
  
# create annotation for two cols in a dataset and order col_2 (subpathway), within each col_1(superpathway)  
  color_order_complexheatmap<-function(dataset,col_1,col_2){
    sp_color_col1 <- colorRampPalette(RColorBrewer::brewer.pal(12,"Paired"))(length(levels(dataset[,col_1])))
    names(sp_color_col1) <- levels(dataset[,col_1])
    # col2_order<-data.table::setorder(dataset,c(col_1,col_2))
    sortCols<-c(col_1,col_2)
    col2_order<-dataset[do.call("order", dataset[sortCols]), ]
    sp_color_col2 <- colorRampPalette(RColorBrewer::brewer.pal(12,"Paired"))(nrow(col2_order))
    names(sp_color_col2) <- col2_order[,col_2]
    # Reorder sub_pathway according to the order they appear in the table (first sort by super_pathway)
    dataset$col_2_ordered<-factor(dataset[,col_2],as.character(unique(dataset[,col_2])))
    require("ComplexHeatmap")
    # combine the super pathway and sub pathway annotation
    row_ha_total <- HeatmapAnnotation(df = data.frame(super_pathway=dataset[,col_1],
                                                      sub_pathway=dataset[,"col_2_ordered"]),
                                      which="row", col = list(
                                        super_pathway = sp_color_col1,
                                        sub_pathway=sp_color_col2
                                      ),
                                      show_annotation_name = FALSE)  
    return(row_ha_total)
  }
  
  # Insert empty cols to a data frame while carrying the same group name
  add_empty_rows<-function(df,empty_bar,col_name){
    to_add<-df[1:empty_bar,]
    to_add[]<-NA
    colnames(to_add)<-colnames(df)
    to_add[[col_name]]<-df[[col_name]][1]
    output_tbl<-rbind(df,to_add)
  }
  
  # Calculate Dice Coefficient for each pair of categories
  dice_coefficient <- function(cat1, cat2, data) {
    met_cat1 <- unique(data$metabolite[data$category == cat1])
    met_cat2 <- unique(data$metabolite[data$category == cat2])
    shared_metabolites <- length(intersect(met_cat1, met_cat2))
    total_metabolites <- length(unique(c(met_cat1, met_cat2)))
    dice <- 2 * shared_metabolites / total_metabolites
    return(dice)
  }