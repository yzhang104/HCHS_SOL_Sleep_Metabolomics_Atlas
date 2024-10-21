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
    out_mspe[i]=mspe
    out_mae[i]=mae
    out_rmse[i]=rmse
  }

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


  
  #####################
  # Split violin plot
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

# svyglm loop over multiple exposures and outcomes accmondating nonlinearity
  svyreg_loop_nonlinear<-function(data,covar,end,metab_is_cont,metab_is_complete,trait_original,trait_binary,trait_for_model,trait_as_predictor,model_input){
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
    # add nonlinear model to covaraites based on input table
    covar_nl<-covar
    for (j in seq_along(covar)){
      nl_col_number<-which(colnames(model_input) == paste(covar[j]))
      covar_nl[j]<-ifelse(covar[j]%in%colnames(model_input),
                       ifelse(model_input[which(model_input$metabolite==paste0(colnames(dat)[i])),nl_col_number]==0,covar[j],
                              paste0("ns(",covar[j],",df=",model_input[which(model_input$metabolite==paste0(colnames(dat)[i])),nl_col_number]+2,")")),covar[j])
    }
    # print(i)
    if (trait_as_predictor==T){
      if (metab_is_cont==T){
        model<- svyglm(as.formula(paste0(colnames(dat)[i],"~",trait,"+",paste(covar_nl,collapse= "+"))),design=subset(survey_design,ID%in%met_df_id))
      } else {
        model<- svyglm(as.formula(paste0(colnames(dat)[i],"~",trait,"+",paste(covar_nl,collapse = "+"))),design=subset(survey_design,ID%in%met_df_id),family=quasipoisson(link='log'))
      }
    } else {
      if (trait_is_cont==T){
        model<- svyglm(as.formula(paste0(trait,"~",colnames(dat)[i],"+",paste(covar_nl,collapse= "+"))),design=subset(survey_design,ID%in%met_df_id))
      } else {
        model<- svyglm(as.formula(paste0(trait,"~",colnames(dat)[i],"+",paste(covar_nl,collapse = "+"))),design=subset(survey_design,ID%in%met_df_id),family=quasipoisson(link='log'))
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
  
  
  # Function to loop regression model and export model summary using survey weight incorperating nonlinear model
  strat_svyreg_loop_nonlinear<-function(data,covar,end,metab_is_cont,metab_is_complete,trait_original,trait_binary,trait_for_model,trait_as_predictor,stratifier,model_input){
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
        # add nonlinear model to covaraites based on input table
        for (k in seq_along(covar)){
          covar_nl<-covar
          nl_col_number<-which(colnames(model_input) == paste(covar[k]))
          covar_nl[k]<-ifelse(covar[k]%in%colnames(model_input),
                           ifelse(model_input[which(model_input$metabolite==paste0(colnames(dat)[i])),nl_col_number]==0,covar[k],
                                  paste0("ns(",covar[k],",df=",model_input[which(model_input$metabolite==paste0(colnames(dat)[i])),nl_col_number]+2,")")),covar[k])
        }
        if (trait_as_predictor==T){
          if (metab_is_cont==T){
            model<- svyglm(as.formula(paste0(colnames(dat)[i],"~",trait,"+",paste(covar_nl,collapse= "+"))),design=subset(survey_design,ID%in%met_df_id))
          } else {
            model<- svyglm(as.formula(paste0(colnames(dat)[i],"~",trait,"+",paste(covar_nl,collapse = "+"))),design=subset(survey_design,ID%in%met_df_id),family=quasipoisson(link='log'))
          }
        }else {
          if (trait_is_cont==T){
            model<- svyglm(as.formula(paste0(trait,"~",colnames(dat)[i],"+",paste(covar_nl,collapse= "+"))),design=subset(survey_design,ID%in%met_df_id))
          }else{
            model<- svyglm(as.formula(paste0(trait,"~",colnames(dat)[i],"+",paste(covar_nl,collapse = "+"))),design=subset(survey_design,ID%in%met_df_id),family=quasipoisson(link='log'))
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
  
  # Function to loop the interaction regression model and export model summary using survey weight incorperating nonlinear model
  svyreg_loop_interaction_nonlinear<-function(data,covar,end,metab_is_cont,metab_is_complete,trait_original,trait_binary,trait_for_model,trait_as_predictor,interaction,model_input){
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
      for (j in seq_along(covar)){
        covar2_nl<-covar2
        nl_col_number<-which(colnames(model_input) == paste(covar2[j]))
        covar2_nl[j]<-ifelse(covar2[j]%in%colnames(model_input),
                         ifelse(model_input[which(model_input$metabolite==paste0(colnames(dat)[i])),nl_col_number]==0,covar2[j],
                                paste0("ns(",covar2[j],",df=",model_input[which(model_input$metabolite==paste0(colnames(dat)[i])),nl_col_number]+2,")")),covar2[j])
      }
      
      if (trait_as_predictor==T){
        if (metab_is_cont==T){
          model<- svyglm(as.formula(paste0(colnames(dat)[i],"~",trait,"+",paste(covar2_nl,collapse= "+"),"+",trait,"*",interaction)),design=subset(survey_design,ID%in%met_df_id))
        } else {
          model<- svyglm(as.formula(paste0(colnames(dat)[i],"~",trait,"+",paste(covar2_nl,collapse = "+"),"+",trait,"*",interaction)),design=subset(survey_design,ID%in%met_df_id),family=quasipoisson(link='log'))
        }
      } else {
        if (trait_is_cont==T){
          model<- svyglm(as.formula(paste0(trait,"~",colnames(dat)[i],"+",paste(covar2_nl,collapse= "+"),"+",colnames(dat)[i],"*",interaction)),design=subset(survey_design,ID%in%met_df_id))
        } else {
          model<- svyglm(as.formula(paste0(trait,"~",colnames(dat)[i],"+",paste(covar2_nl,collapse = "+"),"+",colnames(dat)[i],"*",interaction)),design=subset(survey_design,ID%in%met_df_id),family=quasipoisson(link='log'))
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
  
  
  # Function to loop regression model and export model summary using survey weight and circular variables by modeling the circular response variables as sin(y) and cos(y)
  # for the multi-response regression
  # incorperating nonlinear model
  atlas_cir_svyreg_loop_nonlinear<-function(data,covar,end,metab_is_cont,metab_is_complete,trait,trait_as_predictor,model_input){
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
      # add nonlinear model to covaraites based on input table
      covar_nl<-covar
      for (j in seq_along(covar)){
        nl_col_number<-which(colnames(model_input) == paste(covar[j]))
        covar_nl[j]<-ifelse(covar[j]%in%colnames(model_input),
                            ifelse(model_input[which(model_input$metabolite==paste0(colnames(dat)[i])),nl_col_number]==0,covar[j],
                                   paste0("ns(",covar[j],",df=",model_input[which(model_input$metabolite==paste0(colnames(dat)[i])),nl_col_number]+2,")")),covar[j])
      }
      
      if (trait_as_predictor==T){
        if (metab_is_cont==T){
          require(survey)
          require(jtools)
          model<- svyglm(as.formula(paste0(colnames(dat)[i],"~sin(",trait,")+cos(",trait,")+",paste(covar_nl,collapse= "+"))),design=survey_design)
          # corr<-svycor(as.formula(paste0("~sin(",trait,")+cos(",trait,")")), design = survey_design,na.rm = T)
        } else {
          require(survey)
          require(jtools)
          model<- svyglm(as.formula(paste0(colnames(dat)[i],"~sin(",trait,")+cos(",trait,")+",paste(covar_nl,collapse = "+"))),design=survey_design,family=quasipoisson(link='log'))
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
        model<-estimatr::lm_robust(formula=as.formula(paste0("cbind(",trait,".sin,",trait,".cos)","~",colnames(dat)[i],"+",paste(unlist(covar_nl),collapse= "+"))),data=data, weights = WEIGHT, clusters=STRAT)
        test<-manova(as.formula(paste0("cbind(",trait,".sin,",trait,".cos)","~",colnames(dat)[i],"+",paste(unlist(covar_nl),collapse= "+"))),data=data, weights = WEIGHT)
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
  
  
  # Function to loop regression model and export model summary using survey weight incoperating nonlinear model
  atlas_cir_strat_svyreg_loop_nonlinear<-function(data,covar,end,metab_is_cont,metab_is_complete,trait,trait_for_model,trait_as_predictor,stratifier,model_input){
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
          # add nonlinear model to covaraites based on input table
          covar_nl<-covar
          for (k in seq_along(covar)){
            nl_col_number<-which(colnames(model_input) == paste(covar[k]))
            covar_nl[k]<-ifelse(covar[k]%in%colnames(model_input),
                                ifelse(model_input[which(model_input$metabolite==paste0(colnames(dat)[i])),nl_col_number]==0,covar[k],
                                       paste0("ns(",covar[k],",df=",model_input[which(model_input$metabolite==paste0(colnames(dat)[i])),nl_col_number]+2,")")),covar[k])
          }
        if (trait_as_predictor==T){
          if (metab_is_cont==T){
            require(survey)
            require(jtools)
            model<- svyglm(as.formula(paste0(colnames(dat)[i],"~sin(",trait,")+cos(",trait,")+",paste(covar_nl,collapse= "+"))),design=survey_design)
          } else {
            require(survey)
            require(jtools)
            model<- svyglm(as.formula(paste0(colnames(dat)[i],"~sin(",trait,")+cos(",trait,")+",paste(covar_nl,collapse = "+"))),design=survey_design,family=quasipoisson(link='log'))
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
          model<-estimatr::lm_robust(formula=as.formula(paste0("cbind(",trait,".sin,",trait,".cos)","~",colnames(dat)[i],"+",paste(unlist(covar_nl),collapse= "+"))),data=newdata, weights = WEIGHT, clusters=STRAT)
          test<-manova(as.formula(paste0("cbind(",trait,".sin,",trait,".cos)","~",colnames(dat)[i],"+",paste(unlist(covar_nl),collapse= "+"))),data=data, weights = WEIGHT)
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
    if (metab_is_complete==T){
      regress_output$is_continuous<-"Complete-cases"
    } else {
      regress_output$is_continuous<-factor(regress_output$is_continuous,levels=c(TRUE,FALSE),labels=c("Imputed","Dichotomized"))
    }
    
    return(regress_output)
  }
  
  
  # Function to loop regression sex-interaction model and export model summary using survey weight with nonlinear model
  atlas_svyreg_loop_interaction_nonlinear<-function(data,covar,end,metab_is_cont,metab_is_complete,trait,trait_as_predictor,interaction_term,model_input){
    dat<-data[,1:end]
    dat<-dat[,!colnames(dat)%in%c("LAB_ID","ID","SOL_ID","id")]
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
      # add nonlinear model to covaraites based on input table
      covar_nl<-covar
      for (j in seq_along(covar)){
        nl_col_number<-which(colnames(model_input) == paste(covar[j]))
        covar_nl[j]<-ifelse(covar[j]%in%colnames(model_input),
                            ifelse(model_input[which(model_input$metabolite==paste0(colnames(dat)[i])),nl_col_number]==0,covar[j],
                                   paste0("ns(",covar[j],",df=",model_input[which(model_input$metabolite==paste0(colnames(dat)[i])),nl_col_number]+2,")")),covar[j])
      }
      # print(i)
      if (trait_as_predictor==T){
        if (metab_is_cont==T){
          model<- svyglm(as.formula(paste0(colnames(dat)[i],"~",trait,"+",paste(c(paste0(trait,"*",interaction_term),covar_nl),collapse= "+"))),design=subset(survey_design,ID%in%met_df_id))
        } else {
          model<- svyglm(as.formula(paste0(colnames(dat)[i],"~",trait,"+",paste(c(paste0(trait,"*",interaction_term),covar_nl),collapse= "+"))),design=subset(survey_design,ID%in%met_df_id),family=quasipoisson(link='log'))
        }
      } else {
        if (trait_is_cont==T){
          model<- svyglm(as.formula(paste0(trait,"~",colnames(dat)[i],"+",paste(c(paste0(trait,"*",interaction_term),covar_nl),collapse= "+"))),design=subset(survey_design,ID%in%met_df_id))
        } else {
          model<- svyglm(as.formula(paste0("as.numeric(",trait,")~",colnames(dat)[i],"+",paste(c(paste0(trait,"*",interaction_term),covar_nl),collapse= "+"))),design=subset(survey_design,ID%in%met_df_id),family=quasipoisson(link='log'))
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
    if (metab_is_complete==T){
      regress_output$is_continuous<-"Complete-cases"
    } else {
      regress_output$is_continuous<-factor(regress_output$is_continuous,levels=c(TRUE,FALSE),labels=c("Imputed","Dichotomized"))
    }
    
    return(regress_output)
  }
  