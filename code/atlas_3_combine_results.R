library("ComplexHeatmap")
library("circlize")
source("sol_metabolomics_functions.R")
library("tidyr")
library("dplyr")
library("igraph")
library("ggplot2")
library("ggraph")
library("tidygraph")
library("RColorBrewer")
library("bipartite")
# combine output files
directory_path<-getwd()
# get the file names of the run includes unknown metabolites
file_names <- dir("output/main")
cir_file_names<-dir("output/circular")
# read in annotation list
hchs_annot<-readxl::read_xlsx("data/hchs_batch2/batch2_data.xlsx",sheet="metabolites.info")%>%
  dplyr::mutate(metabolite=paste0("chem_",CHEM_ID),
                CHEM_ID=as.character(CHEM_ID),
                SUPER_PATHWAY=ifelse(!is.na(SUPER_PATHWAY),SUPER_PATHWAY,"Unnamed"))
colnames(hchs_annot)<-tolower(colnames(hchs_annot))
# readin metabolon correction table
lookup_table<-read.csv("data/metabolon_correction.csv") #https://www.nature.com/articles/s42255-024-01018-7/tables/1
# replace mislabelled metabolite chemical names
hchs_annot<-hchs_annot %>%
  mutate(
    chemical_name = case_when(
      chemical_name %in% lookup_table$unknown ~ 
        lookup_table$correct_id[match(chemical_name, lookup_table$unknown)],
      chemical_name %in% lookup_table$incorrect_id ~ 
        lookup_table$correct_id[match(chemical_name, lookup_table$incorrect_id)],
      TRUE ~ chemical_name  # Keep original value if no match
    )
  )
# get the input grid
input_grid<-read.csv("data/atlast_input_grid.csv",stringsAsFactors = F,header = T)
cir_input_grid<-read.csv(file="data/atlast_cir_input_grid.csv",stringsAsFactors = F,header = T)
sleep_trait_order<-c(unique(input_grid$sleep_trait),unique(cir_input_grid$sleep_trait))
# replace the trait name we_minus_wd with social_jetlag
sleep_trait_order[sleep_trait_order=="we_minus_wd"]<-"social_jetlag"

# read in each csv file and rbind
# including unknown metabolites
comb_df <- do.call(rbind,lapply(paste0("output/main/",file_names),read.csv))
# replace the trait name we_minus_wd with social_jetlag
comb_df[which(comb_df$trait=="we_minus_wd"),"trait"]<-"social_jetlag"
comb_df_annot<-merge(comb_df,hchs_annot[,c("metabolite","super_pathway","sub_pathway","hmdb")],by="metabolite",all.x=T)%>%
  dplyr::mutate(beta_p=ifelse(beta>=0,-log10(p_val),-log10(p_val)*(-1)),
                strata=tolower(strata))
# circular regression results
cir_comb_df<- do.call(rbind,lapply(paste0("output/circular/",cir_file_names),read.csv))%>%
    dplyr::mutate(amp=sqrt(sin_beta^2+cos_beta^2),
                se_amp=sqrt(sin_se^2+cos_se^2),
                acrophase=atan2(sin_beta,cos_beta))
cir_comb_df$strata<-tolower(cir_comb_df$strata)
cir_comb_df_annot<-merge(cir_comb_df,hchs_annot[,c("metabolite","super_pathway","sub_pathway","hmdb")],by="metabolite",all.x=T)%>%
  dplyr::mutate(beta=amp,
                se=se_amp,
                p_val=p_acat,
                beta_p=ifelse(beta>=0,-log10(p_val),-log10(p_val)*(-1)))
# combine linear and circular results
comb_df_annot<-plyr::rbind.fill(comb_df_annot,cir_comb_df_annot[,c("metabolite","beta","se","n","p_val","p_val_fdr" , "trait","model" ,"strata","covar","super_pathway","sub_pathway","hmdb" ,"beta_p")])
# ONLY include traits with significant results from model 1 (model 2 only trait excluded)
sig_trait_mdl1<-unique(comb_df_annot[comb_df_annot$p_val_fdr<0.05&comb_df_annot$model=="model_1","trait"])

# including unknown metabolites, including circular
ref_table<-data.frame(trait=c(unique(comb_df$trait),unique(cir_comb_df$trait))[match(sleep_trait_order,c(unique(comb_df$trait),unique(cir_comb_df$trait)))],
                      label=c("weekday sleep duration", "weekend sleep duration", "sleep duration", "weekday short sleep", "weekday long sleep",
                              "social jetlag","whiirs", "restless sleep", "sleep pill", "difficulty back to sleep", 
                              "early wake", "frequent wake", "difficulty fall asleep", "ess","ess>10", 
                              "snore", "rei0", "rei3", "rei0>5", "rei0>15", 
                              "rei3>5", "rei3>15", "total event count", "total event duration", "min o2", 
                              "avg o2", "perlt90", "per90", "avg event duration", "hypoxic burden", 
                              "min hr", "max hr", "avg hr", "std hr", "weekday wake time",
                              "weekend wake time", "weekday bed time", "weekend bed time", "weekday midpoint time", "weekend midpoint time"),
                      category=c(rep("duration",5),"timing",rep("insomnia",9),rep("sdb",15),rep("hr",4),rep("timing",6)),
                      dichotomized=c(FALSE,FALSE,FALSE,TRUE,TRUE,
                                     FALSE,FALSE,TRUE,TRUE,TRUE,
                                     TRUE,TRUE,TRUE,FALSE,TRUE,
                                     TRUE,FALSE,FALSE,TRUE,TRUE,
                                     TRUE,TRUE,FALSE,FALSE,FALSE,
                                     FALSE,FALSE,TRUE,FALSE,FALSE,
                                     FALSE,FALSE,FALSE,FALSE,FALSE,
                                     FALSE,FALSE,FALSE,FALSE,FALSE))

# subset model 1 with combined sex analysis
both_md1<-comb_df_annot[which(comb_df_annot$strata=="both"&comb_df_annot$model=="model_1"),]%>%
  dplyr::mutate(full_pathway=factor(paste(super_pathway,sub_pathway,metabolite)))%>%
  merge(.,ref_table,by="trait")

both_md1 <- both_md1[order(both_md1$category,both_md1$trait),]

######################
# Figure 2 - Figure 6
 # Figure 2
 # Plot the count of significant associations between subpathways and sleep trait categories (>=50% significant)
 subpathway_sig_list<-unique(as.vector(unlist(summary_subpathway_category_tbl[which(summary_subpathway_category_tbl$percent_significant>=25),"sub_pathway"])))
 sig_subpathway_df<-comb_df_annot[which(comb_df_annot$model=="model_1"&comb_df_annot$sub_pathway%in%subpathway_sig_list),]
 # Calculating the percentage of significant associations (both)
 sig_subpathway_df_summary <- sig_subpathway_df %>%
   dplyr::filter(strata=="both")%>%
   # dplyr::filter(strata=="male")%>%
   merge(.,ref_table,by="trait",all.x = T)%>%
   group_by(super_pathway, sub_pathway, category) %>%
   dplyr::summarise(
     num_total = n(),
     num_significant = sum(p_val_fdr < 0.05),
     percent_significant = (num_significant / num_total) * 100
   )%>%
   ungroup() %>%
   arrange(super_pathway, sub_pathway, category)
 
 # Creating a cumulative sum for stacking
 sig_subpathway_df_summary_2 <- sig_subpathway_df_summary %>%
   group_by(super_pathway, sub_pathway) %>%
   mutate(cumulative = cumsum(percent_significant)) %>%
   ungroup()
 
 # Transform data to long format
 sig_subpathway_df_long <- sig_subpathway_df_summary_2 %>%
   pivot_longer(cols = percent_significant, names_to = "observation", values_to = "value")

  # Set a number of 'empty bar' to add at the end of each super_pathway
 empty_bar <- 2
 sig_subpathway_df_long<-sig_subpathway_df_long%>%
   group_by(super_pathway) %>%
   group_modify(~add_empty_rows(.,empty_bar=empty_bar,col_name="super_pathway"))%>%
   ungroup() %>%
   arrange(super_pathway, sub_pathway, category)%>%
   mutate(id = ifelse(is.na(sub_pathway),paste0(super_pathway,"_ZZZZZ"),paste(super_pathway, sub_pathway, sep = "_")),
          id=as.numeric(factor(id))
   )

 # Prepare label and base data
 label_data <- sig_subpathway_df_long %>% 
   group_by(id,sub_pathway) %>% 
   dplyr::summarise(tot = sum(value, na.rm = TRUE)) %>%
   ungroup()
 number_of_bar <- nrow(label_data)
 label_data <- label_data %>%
   dplyr::mutate(angle = 90 - 360 * (id-0.5) / number_of_bar,
          hjust = ifelse(angle < -90, 1, 0),
          angle = ifelse(angle < -90, angle+180, angle)
          )
 label_data <- label_data %>%
   dplyr::mutate(wrapped_label = stringr::str_wrap(sub_pathway, width = 35)) # Adjust width as needed
 
 empty_bar_1=1
 base_data <- sig_subpathway_df_long %>% 
   group_by(super_pathway) %>% 
   dplyr::summarise(start = min(id), end = max(id)) %>% 
   dplyr::rowwise() %>% 
   dplyr::mutate(title = mean(c(start, end)))
 # Create the circular stacked bar plot
 circular_stacked_bar_plot <-
   ggplot(sig_subpathway_df_long, aes(x = as.factor(id), y = value, fill = category)) +
   geom_bar(stat = "identity", alpha = 0.5) +
   scale_fill_viridis(discrete = TRUE) +
   coord_polar(start = 0) +
   theme_minimal() +
   theme(legend.position = "right",
         axis.text = element_blank(),
         axis.title = element_blank(),
         panel.grid = element_blank(),
         # plot.margin = unit(rep(-1, 4), "cm"),
         plot.margin = margin(t = -10, r = 10, b = -5, l = -10, unit = "mm"),
         plot.background = element_rect(fill = "white", colour = "white")) +
   ylim(-80, max(label_data$tot+20, na.rm = TRUE)) +
   geom_text(data = label_data, aes(x = id, y = tot + 5, label = wrapped_label,  hjust = hjust, angle = angle), 
             color = "black", alpha = 0.5, size = 3, inherit.aes = FALSE, lineheight = 0.7) +
   geom_segment(data = base_data, aes(x = start-0.5, y = -5, xend = end-0.5, yend = -5), 
                colour = "black", alpha = 0.8, size = 0.6, inherit.aes = FALSE) +
   geom_text(data = base_data, aes(x = title, y = -18, label = super_pathway), 
             # hjust = c(0.5, 1,1,0,0.3, 0.5), colour = "black", alpha = 0.8, size = 3, # cutoff >=50%
             hjust = c(0.5, 1,1,0.5,0.3, 0.5,0.5), colour = "black", alpha = 0.8, size = 3, # cutoff >=25%
             fontface = "bold", inherit.aes = FALSE)
 ggsave(circular_stacked_bar_plot, file = "output/figure2.png", width = 12, height = 10, units = "in") # Adjust dimensions as needed

# Figure 3
# Heatmap for metabolite super pathway by sleep trait categories
  heatmap_category<-ggplot(summary_superpathway_category_tbl, aes(x = super_pathway, y = category, fill = percent_significant)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "blue") +
    labs(x = "Super Pathway", y = "Sleep Trait Category",
         fill = "Percentage Significant (%)") +
    theme(text = element_text(size = 13),legend.position="bottom")+
   coord_flip()
  png("output/figure3.png", width = 500, height = 500)
  print(heatmap_category)
  dev.off()

 # Figure 4
 # identify metabolites with significant associations with most sleep traits
 # Filter for significant associations
 metab_significant_associations <- comb_df_annot %>%
   filter(strata=="both"&model=="model_1")%>%
   filter(p_val_fdr < 0.05)%>%
   merge(.,ref_table,by="trait",all.x = T)%>%
   merge(.,hchs_annot[,c("metabolite","chemical_name")],all.x = T)%>%
   arrange(super_pathway, sub_pathway, category)
 
 # Count the number of unique outcomes and outcome categories for each predictor
 metab_significant_associations_counts <- metab_significant_associations %>%
   group_by(chemical_name, category) %>%
   dplyr::summarise(count = n(), .groups = 'drop') %>%
   merge(.,hchs_annot[,c("super_pathway","sub_pathway","chemical_name")],all.x = T)%>%
   arrange(super_pathway, sub_pathway,chemical_name, category)
 # identify the number of metabolites associated with more than one sleep trait
 length(unique(metab_significant_associations_counts[which(metab_significant_associations_counts$total_count>1),"chemical_name"]))
 num_metab_gt2_category<-metab_significant_associations_counts%>%
   group_by(chemical_name)%>%
   dplyr::summarise(category_count = n_distinct(category))
 dim(unique(num_metab_gt2_category[which(num_metab_gt2_category$category_count>1),"chemical_name"]))
  
top_metab_count<-metab_significant_associations_counts%>%
   group_by(chemical_name)%>%
   dplyr::summarise(total_count=sum(count))%>%
   arrange(desc(total_count))
 metab_levels <- top_metab_count$chemical_name
 
 metab_significant_associations_counts<-merge(metab_significant_associations_counts,top_metab_count,all.x=T)%>%
   arrange(desc(total_count)) %>%
   mutate(rank = min_rank(desc(total_count)))
 metab_significant_associations_counts$chemical_name <- factor(metab_significant_associations_counts$chemical_name, levels = metab_levels)
 
 metab_significant_associations_top_count<-metab_significant_associations_counts%>%
   filter(rank <= n() * 0.1)%>%
   arrange(desc(total_count))%>%
   dplyr::mutate(chemical_name_sort=as.factor(chemical_name))
  
 # Prepare a color mapping for domains and subdomains
super_pathway_colors <- setNames(brewer.pal(nlevels(as.factor(metab_significant_associations_top_count$super_pathway)), "Set3"), levels(as.factor(metab_significant_associations_top_count$super_pathwayn)))
 names(super_pathway_colors)<-levels(as.factor(metab_significant_associations_top_count$super_pathway))
 category_colors <- brewer.pal(nlevels(as.factor(metab_significant_associations_top_count$category)), "Set1")
 
 # Create an annotation data frame
 annotation_df <- metab_significant_associations_top_count %>%
   distinct(super_pathway, chemical_name) %>%
   mutate(super_pathway_color = super_pathway_colors[super_pathway]
          )

 metab_top_count_plot<-ggplot(metab_significant_associations_top_count, aes(x = chemical_name_sort, y = count)) +
   geom_bar(aes(fill = category), stat = "identity") +
   scale_fill_brewer(palette = "Dark2", name = "Sleep Phenotype Domain") +  # Using "Set3" palette for fill colors
   theme_minimal() +
   labs(title = "Top 10% Metabolites by Significant Associations",
        x = "Metabolites",
        y = "Count of Significant Associations") +
   theme(axis.text.x = element_text(angle = 0, hjust = 1),
         plot.background = element_rect(fill = "white", colour = "white")) + coord_flip()+
   # Use ggnewscale to add a new fill scale
   ggnewscale::new_scale_fill() +
   geom_tile(data = metab_significant_associations_top_count,
             aes(x = chemical_name_sort, y = -0.5, fill = super_pathway), # Adjust y to place below bars
             height = 0.2) +  # Adjust height as needed for the annotation bar
   scale_fill_brewer(palette = "Set3", name = "Super Pathway") +  
   scale_y_continuous(expand = c(0, 0), limits = c(-1, max(metab_significant_associations_top_count$total_count)))
 ggsave(metab_top_count_plot, file = "output/figure4.png", width = 12, height = 10, units = "in") 
 
# Figure 5
# Create Dice coefficient matrix for sleep traits category and metabolites
both_md1_sig<-comb_df_annot[which(comb_df_annot$strata=="both"&comb_df_annot$model=="model_1"&comb_df_annot$p_val_fdr<0.05),]%>%
  dplyr::mutate(full_pathway=factor(paste(super_pathway,sub_pathway,metabolite)))%>%
  merge(.,ref_table,by="trait")
both_md1_sig <- both_md1_sig[order(both_md1_sig$category,both_md1_sig$trait),]
# Create a matrix for the Dice coefficients (sleep traits by metabolites)
categories <- unique(both_md1_sig$category)
dice_matrix <- matrix(nrow = length(categories), ncol = length(categories),
                      dimnames = list(categories, categories))

for (i in 1:length(categories)) {
  for (j in 1:length(categories)) {
    dice_matrix[i, j] <- dice_coefficient(categories[i], categories[j], both_md1_sig)
  }
}

# Convert the matrix to a long format for plotting
dice_long <- reshape2::melt(dice_matrix)%>%
  dplyr::mutate(Var1_char=as.character(Var1),
                Var2_char=as.character(Var2))

dice_long_filtered <- dice_long[dice_long$Var1_char >= dice_long$Var2_char, ]
# Plotting the correlation matrix
dice_matrix_category_metab<-ggplot(dice_long_filtered, aes(Var1, Var2, fill = value)) +
  geom_tile() +geom_text(aes(label=round(value,digits=2)))+
  scale_fill_gradient2(low = "white", high = "blue", limit = c(0,0.5)) +
  theme_minimal() +
  theme(plot.background =element_rect(fill = "white"))+
  xlab("Domain") +
  ylab("Domain") +
  ggtitle("Dice Coefficient Matrix")
ggsave(dice_matrix_category_metab, file = "output/figure5.png", width = 6, height = 5, units = "in") 

# Network analysis
# Create the incidence matrix
incidence_matrix <- reshape2::acast(both_md1_sig, trait ~ metabolite, value.var = "trait", fun.aggregate = length, fill = 0)
#######################
# consolidate nodes into domains (metabolite->sub_pathway, sleep traits->sleep trait domains)
# Create the incidence matrix
incidence_matrix_consolidate <- reshape2::acast(both_md1_sig, category ~ sub_pathway, value.var = "category", fun.aggregate = length, fill = 0)
# replace value>0 to 1
incidence_matrix_consolidate_convert<-incidence_matrix_consolidate
incidence_matrix_consolidate_convert[incidence_matrix_consolidate_convert>1]=1
# This function generates a layout using a variant of Fruchterman and Reingold's force-directed placement algorithm.
n_consolidate <- ggnetwork(incidence_matrix_consolidate_convert,layout = "fruchtermanreingold", cell.jitter = 0.75,niter=10000)
# ggnetwork generated multiple rows for each sub_pathway

# single out the primary mode 
# add super_pathway col to the mode
n_consolidate$label<-ifelse(n_consolidate$vertex.name %in% both_md1_sig$category, 
                paste("Sleep:",n_consolidate$vertex.name),
                paste("Metab:",both_md1_sig$super_pathway[match(n_consolidate$vertex.name, both_md1_sig$sub_pathway)]))
n_consolidate$group<-ifelse(n_consolidate$vertex.name %in% both_md1_sig$category,"sleep","metabolite")
n_consolidate$definition<-ifelse(n_consolidate$vertex.name %in% both_md1_sig$category,both_md1_sig$label[match(n_consolidate$vertex.name, both_md1_sig$category)],n_consolidate$vertex.nam)
n_consolidate$category <- n_consolidate$vertex.names %in% unique(both_md1_sig$category)

# Convert the incidence matrix back to a long data frame format to get the counts
link_counts <- reshape2::melt(incidence_matrix_consolidate)
link_counts <- link_counts[link_counts$value > 0, ]  # Keep only rows with links
link_counts$value<-link_counts$value*10/(max(link_counts$value))
link_counts$Var1_sleep<-paste("Sleep:",link_counts$Var1)
# Compute the degree (number of links) for each node
category_degree <- rowSums(incidence_matrix_consolidate)
sub_pathway_degree <- colSums(incidence_matrix_consolidate)
# Add this degree information to the n object
n_consolidate$degree <- ifelse(n_consolidate$vertex.names %in% names(category_degree), 
                   category_degree[n_consolidate$vertex.names], 
                   sub_pathway_degree[n_consolidate$vertex.names])


# Generate color palettes
category_colors <- RColorBrewer::brewer.pal(5, "Dark2")
other_colors <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(9)

# Assign colors to nodes
n_consolidate$color <- ifelse(n_consolidate$category, 
                              category_colors[as.numeric(as.factor(subset(n_consolidate, category)$label))],
                              other_colors[as.numeric(as.factor(subset(n_consolidate, !category)$label))])

# Assign edge colors based on the start node color
link_counts$edgecolor <- n_consolidate$color[match(link_counts$Var1_sleep, n_consolidate$vertex.names)]

# Create a unique label dataset for non-category nodes
unique_labels_metabolite <- n_consolidate %>% 
  filter(group == "metabolite") %>%
  distinct(vertex.names, .keep_all = TRUE)

# Figure 6
Fruchterman_Reingold_network<-ggplot(n_consolidate, aes(x, y, xend = xend, yend = yend)) +
  # Edges from link_counts
geom_segment(data = link_counts,
             aes(x = n_consolidate$x[match(Var1, n_consolidate$vertex.names)], 
                 y = n_consolidate$y[match(Var1, n_consolidate$vertex.names)],
                 xend = n_consolidate$x[match(Var2, n_consolidate$vertex.names)], 
                 yend = n_consolidate$y[match(Var2, n_consolidate$vertex.names)],color = edgecolor),
                 size = link_counts$value, alpha = 0.5) +  
  geom_node_point(
    aes(x = x, y = y, color = label,shape = group),
    size = 4) +
  geom_point(data = n_consolidate, aes(x = x, y = y, color = label, shape = group, size = degree)) +
  geom_text_repel(data = subset(n_consolidate, category),aes(x = x, y = y, label = vertex.names),
                  box.padding = 0.5, point.padding = 0.5,
                  size = 5, fontface = "bold", color = "black",max.overlaps = Inf) +
  # Use unique_labels_metabolite for the non-category node labels
  geom_text_repel(data = unique_labels_metabolite,
                  aes(x = x, y = y, label = vertex.names),
                  box.padding = 0.1, point.padding = 0.1,
                  size = 3, fontface = "plain", color = "darkred", max.overlaps = Inf) +  
  theme_void()+
  theme(plot.background = element_rect(fill = "white"),
        panel.grid = element_blank())
ggsave(Fruchterman_Reingold_network, file = "output/figure6.png", width = 14, height = 8, units = "in") # Adjust dimensions as needed

#########################
# Table 2 - 3
# Table 2
summary_tbl<-comb_df_annot[which(comb_df_annot$model=="model_1"),] %>%
  dplyr::group_by(model, strata, trait) %>%
  dplyr::summarise(
    num_significant = sum(p_val_fdr < 0.05, na.rm = TRUE),
    percent_significant = sum(p_val_fdr < 0.05, na.rm = TRUE)* 100/768 
  )%>%
  merge(.,ref_table,by="trait",all.x = T)%>%
  dplyr::mutate(full_trait=paste(category,trait))
summary_tbl <- summary_tbl[order(summary_tbl$full_trait),]
summary_trait_tbl<-summary_tbl[which(summary_tbl$strata=="both"),]%>%dplyr::mutate(num_significant_div_metab=num_significant*100/768)
write.csv(summary_trait_tbl,file="output/table2.csv",row.names=F)

# Table 3
# export the top metabolites in a table
 write.csv(metab_significant_associations_top_count,file= "output/table3.csv",row.names = F)
 
#######################
# Supplementary Figure S2,S3, S5 - S8
# Supplemental Figure S2
# create manhattan style plot
manhattan_plot<-ggplot(both_md1) +
  geom_point(aes(full_pathway, -log10(p_val_fdr), colour = as.factor(super_pathway)), size = 1)+
  geom_hline(yintercept=-log10(0.05))+
  # facet_wrap(vars(label))+
  facet_wrap(vars(factor(label,levels=unique(label))))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  labs(y="-log10(FDR-p)",
       color="super pathway")
png("output/figure2.png", width = 1000, height = 800)
print(manhattan_plot)
dev.off()

# Supplemental Figure S3
# create a point plot for the summary table of the number of significant associations by trait and by strata
# summary table for the number of associations
summary_table<-data.frame(definition=NA,
                          median=NA,
                          mean=NA,
                          min=NA,
                          max=NA,
                          median_perc=NA,
                          mean_perc=NA,
                          min_perc=NA,
                          max_perc=NA,
                          percentage_note=NA)
# model 1 only
summary_tbl<-comb_df_annot[which(comb_df_annot$model=="model_1"),] %>%
  dplyr::group_by(model, strata, trait) %>%
  dplyr::summarise(
    num_significant = sum(p_val_fdr < 0.05, na.rm = TRUE),
    percent_significant = sum(p_val_fdr < 0.05, na.rm = TRUE)* 100/768 
  )%>%
  merge(.,ref_table,by="trait",all.x = T)%>%
  dplyr::mutate(full_trait=paste(category,trait))
summary_tbl <- summary_tbl[order(summary_tbl$full_trait),]

# scatter plot by sleep trait categories
scatterplot_category<-ggplot(summary_tbl[which(summary_tbl$model=="model_1"),])+
  geom_point(aes(label, num_significant , colour = as.factor(strata)), size = 2, position = position_dodge(0.5))+
  # geom_point(aes(label, num_significant , colour = as.factor(strata)),size = 2)+
  theme_bw()+ coord_flip()+
  facet_wrap(~category,scales="free_y")+
 labs(y="Number of significant associations",color='Strata')+
  theme(text = element_text(size = 15))
png("output/figure3.png", width = 1000, height = 800)
print(scatterplot_category)
dev.off()

# Supplementary Figure S5 
# Create correlation matrix of effect estimates of metabolites by sleep traits in both gender in model 1  
# Create matrix with effect estimates: columns are sleep-related phenotypes, rows are metabolites, 
# and we only use non-cyclical phenotypes. Compute correlation matrix between the phenotypes 
# (i.e. between the columns of this matrix). 
# Both sexes
linear_both_md1<-both_md1[which(both_md1$category!="Timing"),]%>%
  dplyr::mutate(beta2=beta^2)
beta_both_md1<-tidyr::spread(data=linear_both_md1[,c("label","metabolite","beta2")],key="label",value="beta2")
row.names(beta_both_md1)<-beta_both_md1[,1]
beta_both_md1<-beta_both_md1[,-1]
# Calculate the correlation matrix between phenotypes
beta_corr_both_md1 <- cor(beta_both_md1, use = "pairwise.complete.obs", method = "spearman")
corrplot(beta_corr_both_md1,type = "lower",order ="hclust")

pdf(file = "output/supplementary_figures5.pdf")
corrplot(beta_corr_both_md1,type = "lower",order ="hclust")
dev.off()

# Supplementary Figure S6
# Females only
linear_female_md1<-comb_df_annot[which(comb_df_annot$strata=="female"&comb_df_annot$model=="model_1"),]%>%
  dplyr::mutate(full_pathway=factor(paste(super_pathway,sub_pathway,metabolite)))%>%
  merge(.,ref_table,by="trait")%>%
  filter(category!="Timing")%>%
  dplyr::mutate(beta2=beta^2)
beta_female_md1<-tidyr::spread(data=linear_female_md1[,c("label","metabolite","beta2")],key="label",value="beta2")
row.names(beta_female_md1)<-beta_female_md1[,1]
beta_female_md1<-beta_female_md1[,-1]
# Calculate the correlation matrix between phenotypes
beta_corr_female_md1 <- cor(beta_female_md1, use = "pairwise.complete.obs", method = "spearman")
corrplot(beta_corr_female_md1,type = "lower",order ="hclust")

pdf(file = "output/supplementary_figures6.pdf")
corrplot::corrplot(beta_corr_female_md1,type = "lower",order ="hclust")
dev.off()

# Supplementary Figure S7
# Males only
linear_male_md1<-comb_df_annot[which(comb_df_annot$strata=="male"&comb_df_annot$model=="model_1"),]%>%
  dplyr::mutate(full_pathway=factor(paste(super_pathway,sub_pathway,metabolite)))%>%
  merge(.,ref_table,by="trait")%>%
  filter(category!="Timing")%>%
  dplyr::mutate(beta2=beta^2)
beta_male_md1<-tidyr::spread(data=linear_male_md1[,c("label","metabolite","beta2")],key="label",value="beta2")
row.names(beta_male_md1)<-beta_male_md1[,1]
beta_male_md1<-beta_male_md1[,-1]
# Calculate the correlation matrix between phenotypes
beta_corr_male_md1 <- cor(beta_male_md1, use = "pairwise.complete.obs", method = "spearman")
corrplot::corrplot(beta_corr_male_md1,type = "lower",order ="hclust")

pdf(file = "output/supplementary_figures7.pdf")
corrplot::corrplot(beta_corr_male_md1,type = "lower",order ="hclust")
dev.off()

 # supplementary Figure S8
 # circular stacked bar chart by sex
 # female
 sig_subpathway_df_summary_female <- sig_subpathway_df %>%
   dplyr::filter(strata=="female")%>%
   merge(.,ref_table,by="trait",all.x = T)%>%
   group_by(super_pathway, sub_pathway, category) %>%
   dplyr::summarise(
     num_total = n(),
     num_significant = sum(p_val_fdr < 0.05),
     percent_significant = (num_significant / num_total) * 100
   )%>%
   ungroup() %>%
   arrange(super_pathway, sub_pathway, category)
 
 # Creating a cumulative sum for stacking
 sig_subpathway_df_summary_female_2 <- sig_subpathway_df_summary_female %>%
   group_by(super_pathway, sub_pathway) %>%
   mutate(cumulative = cumsum(percent_significant)) %>%
   ungroup()
 
 # Transform data to long format
 sig_subpathway_df_female_long <- sig_subpathway_df_summary_female_2 %>%
   pivot_longer(cols = percent_significant, names_to = "observation", values_to = "value")
 
 # Set a number of 'empty bar' to add at the end of each super_pathway
 empty_bar <- 2
 sig_subpathway_df_female_long<-sig_subpathway_df_female_long%>%
   group_by(super_pathway) %>%
   group_modify(~add_empty_rows(.,empty_bar=empty_bar,col_name="super_pathway"))%>%
   ungroup() %>%
   arrange(super_pathway, sub_pathway, category)%>%
   mutate(id = ifelse(is.na(sub_pathway),paste0(super_pathway,"_ZZZZZ"),paste(super_pathway, sub_pathway, sep = "_")),
          id=as.numeric(factor(id))
   )
 
 # Prepare label and base data
 label_data_female <- sig_subpathway_df_female_long %>% 
   group_by(id,sub_pathway) %>% 
   dplyr::summarise(tot = sum(value, na.rm = TRUE)) %>%
   ungroup()
 number_of_bar <- nrow(label_data_female)
 label_data_female <- label_data_female %>%
   dplyr::mutate(angle = 90 - 360 * (id-0.5) / number_of_bar,
                 hjust = ifelse(angle < -90, 1, 0),
                 angle = ifelse(angle < -90, angle+180, angle)
   )
 label_data_female <- label_data_female %>%
   dplyr::mutate(wrapped_label = stringr::str_wrap(sub_pathway, width = 35)) # Adjust width as needed
 
 empty_bar_1=1
 base_data_female <- sig_subpathway_df_female_long %>% 
   group_by(super_pathway) %>% 
   dplyr::summarise(start = min(id), end = max(id)) %>% 
   dplyr::rowwise() %>% 
   dplyr::mutate(title = mean(c(start, end)))
 # Create the circular stacked bar plot
 circular_stacked_bar_plot_female <-
   ggplot(sig_subpathway_df_female_long, aes(x = as.factor(id), y = value, fill = category)) +
   geom_bar(stat = "identity", alpha = 0.5) +
   scale_fill_viridis(discrete = TRUE) +
   coord_polar(start = 0) +
   theme_minimal() +
   theme(legend.position = "right",
         axis.text = element_blank(),
         axis.title = element_blank(),
         panel.grid = element_blank(),
         # plot.margin = unit(rep(-1, 4), "cm"),
         plot.margin = margin(t = -10, r = 10, b = -5, l = -10, unit = "mm"),
         plot.background = element_rect(fill = "white", colour = "white")) +
   ylim(-80, max(label_data_female$tot+20, na.rm = TRUE)) +
   geom_text(data = label_data_female, aes(x = id, y = tot + 5, label = wrapped_label,  hjust = hjust, angle = angle), 
             color = "black", alpha = 0.5, size = 3, inherit.aes = FALSE, lineheight = 0.7) +
   geom_segment(data = base_data_female, aes(x = start-0.5, y = -5, xend = end-0.5, yend = -5), 
                colour = "black", alpha = 0.8, size = 0.6, inherit.aes = FALSE) +
   geom_text(data = base_data_female, aes(x = title, y = -18, label = super_pathway), 
             # hjust = c(0.5, 1,1,0,0.3, 0.5), colour = "black", alpha = 0.8, size = 3, # cutoff >=50%
             hjust = c(0.5, 1,1,0.5,0.3, 0.5,0.5), colour = "black", alpha = 0.8, size = 3, # cutoff >=25%
             fontface = "bold", inherit.aes = FALSE)+
   ggtitle("Female")
 # male
 sig_subpathway_df_summary_male <- sig_subpathway_df %>%
   dplyr::filter(strata=="male")%>%
   merge(.,ref_table,by="trait",all.x = T)%>%
   group_by(super_pathway, sub_pathway, category) %>%
   dplyr::summarise(
     num_total = n(),
     num_significant = sum(p_val_fdr < 0.05),
     percent_significant = (num_significant / num_total) * 100
   )%>%
   ungroup() %>%
   arrange(super_pathway, sub_pathway, category)
 
 # Creating a cumulative sum for stacking
 sig_subpathway_df_summary_male_2 <- sig_subpathway_df_summary_male %>%
   group_by(super_pathway, sub_pathway) %>%
   mutate(cumulative = cumsum(percent_significant)) %>%
   ungroup()
 
 # Transform data to long format
 sig_subpathway_df_male_long <- sig_subpathway_df_summary_male_2 %>%
   pivot_longer(cols = percent_significant, names_to = "observation", values_to = "value")
 
 # Set a number of 'empty bar' to add at the end of each super_pathway
 empty_bar <- 2
 sig_subpathway_df_male_long<-sig_subpathway_df_male_long%>%
   group_by(super_pathway) %>%
   group_modify(~add_empty_rows(.,empty_bar=empty_bar,col_name="super_pathway"))%>%
   ungroup() %>%
   arrange(super_pathway, sub_pathway, category)%>%
   mutate(id = ifelse(is.na(sub_pathway),paste0(super_pathway,"_ZZZZZ"),paste(super_pathway, sub_pathway, sep = "_")),
          id=as.numeric(factor(id))
   )
 
 # Prepare label and base data
 label_data_male <- sig_subpathway_df_male_long %>% 
   group_by(id,sub_pathway) %>% 
   dplyr::summarise(tot = sum(value, na.rm = TRUE)) %>%
   ungroup()
 number_of_bar <- nrow(label_data_male)
 label_data_male <- label_data_male %>%
   dplyr::mutate(angle = 90 - 360 * (id-0.5) / number_of_bar,
                 hjust = ifelse(angle < -90, 1, 0),
                 angle = ifelse(angle < -90, angle+180, angle)
   )
 label_data_male <- label_data_male %>%
   dplyr::mutate(wrapped_label = stringr::str_wrap(sub_pathway, width = 35)) # Adjust width as needed
 
 empty_bar_1=1
 base_data_male <- sig_subpathway_df_male_long %>% 
   group_by(super_pathway) %>% 
   dplyr::summarise(start = min(id), end = max(id)) %>% 
   dplyr::rowwise() %>% 
   dplyr::mutate(title = mean(c(start, end)))
 # Create the circular stacked bar plot
 circular_stacked_bar_plot_male <-
   ggplot(sig_subpathway_df_male_long, aes(x = as.factor(id), y = value, fill = category)) +
   geom_bar(stat = "identity", alpha = 0.5) +
   scale_fill_viridis(discrete = TRUE) +
   coord_polar(start = 0) +
   theme_minimal() +
   theme(legend.position = "right",
         axis.text = element_blank(),
         axis.title = element_blank(),
         panel.grid = element_blank(),
         # plot.margin = unit(rep(-1, 4), "cm"),
         plot.margin = margin(t = -10, r = 10, b = -5, l = -10, unit = "mm"),
         plot.background = element_rect(fill = "white", colour = "white")) +
   ylim(-80, max(label_data_male$tot+20, na.rm = TRUE)) +
   geom_text(data = label_data_male, aes(x = id, y = tot + 5, label = wrapped_label,  hjust = hjust, angle = angle), 
             color = "black", alpha = 0.5, size = 3, inherit.aes = FALSE, lineheight = 0.7) +
   geom_segment(data = base_data_male, aes(x = start-0.5, y = -5, xend = end-0.5, yend = -5), 
                colour = "black", alpha = 0.8, size = 0.6, inherit.aes = FALSE) +
   geom_text(data = base_data_male, aes(x = title, y = -18, label = super_pathway), 
             # hjust = c(0.5, 1,1,0,0.3, 0.5), colour = "black", alpha = 0.8, size = 3, # cutoff >=50%
             hjust = c(0.5, 1,1,0.5,0.3, 0.5,0.5), colour = "black", alpha = 0.8, size = 3, # cutoff >=25%
             fontface = "bold", inherit.aes = FALSE)+
   ggtitle("Male")
 ggarrange(circular_stacked_bar_plot_male,circular_stacked_bar_plot_female,common.legend=T,legend="bottom",labels=c("Male","Female"))
 ggsave(ggarrange(circular_stacked_bar_plot_male,circular_stacked_bar_plot_female,common.legend=T), file = "output/supplementary_figures8.png", width = 24, height = 10, units = "in") 
 
########################
# Supplementary Table S3, S4, S6 - S9
# Supplementary Table S3
summary_table[1,]=c("model_1,both",
                    median(summary_trait_tbl$num_significant),
                    mean(summary_trait_tbl$num_significant),
                    min(summary_trait_tbl$num_significant),
                    max(summary_trait_tbl$num_significant),
                    median(summary_trait_tbl$percent_significant),
                    mean(summary_trait_tbl$percent_significant),
                    min(summary_trait_tbl$percent_significant),
                    max(summary_trait_tbl$percent_significant),
                    "num_significant/768")
summary_table[2,]=c("model_1,both,dichotomized",
                    median(summary_trait_tbl[which(summary_trait_tbl$dichotomized==T),"num_significant"]),
                    mean(summary_trait_tbl[which(summary_trait_tbl$dichotomized==T),"num_significant"]),
                    min(summary_trait_tbl[which(summary_trait_tbl$dichotomized==T),"num_significant"]),
                    max(summary_trait_tbl[which(summary_trait_tbl$dichotomized==T),"num_significant"]),
                    median(summary_trait_tbl[which(summary_trait_tbl$dichotomized==T),"percent_significant"]),
                    mean(summary_trait_tbl[which(summary_trait_tbl$dichotomized==T),"percent_significant"]),
                    min(summary_trait_tbl[which(summary_trait_tbl$dichotomized==T),"percent_significant"]),
                    max(summary_trait_tbl[which(summary_trait_tbl$dichotomized==T),"percent_significant"]),
                    "num_significant/768")
summary_table[3,]=c("model_1,both,nondichotomized",
                    median(summary_trait_tbl[which(summary_trait_tbl$dichotomized==F),"num_significant"]),
                    mean(summary_trait_tbl[which(summary_trait_tbl$dichotomized==F),"num_significant"]),
                    min(summary_trait_tbl[which(summary_trait_tbl$dichotomized==F),"num_significant"]),
                    max(summary_trait_tbl[which(summary_trait_tbl$dichotomized==F),"num_significant"]),
                    median(summary_trait_tbl[which(summary_trait_tbl$dichotomized==F),"percent_significant"]),
                    mean(summary_trait_tbl[which(summary_trait_tbl$dichotomized==F),"percent_significant"]),
                    min(summary_trait_tbl[which(summary_trait_tbl$dichotomized==F),"percent_significant"]),
                    max(summary_trait_tbl[which(summary_trait_tbl$dichotomized==F),"percent_significant"]),
                    "num_significant/768")
summary_table[4,]=c("model_1,both,timing",
                    median(summary_trait_tbl[which(summary_trait_tbl$category=="timing"),"num_significant"]),
                    mean(summary_trait_tbl[which(summary_trait_tbl$category=="timing"),"num_significant"]),
                    min(summary_trait_tbl[which(summary_trait_tbl$category=="timing"),"num_significant"]),
                    max(summary_trait_tbl[which(summary_trait_tbl$category=="timing"),"num_significant"]),
                    median(summary_trait_tbl[which(summary_trait_tbl$category=="timing"),"percent_significant"]),
                    mean(summary_trait_tbl[which(summary_trait_tbl$category=="timing"),"percent_significant"]),
                    min(summary_trait_tbl[which(summary_trait_tbl$category=="timing"),"percent_significant"]),
                    max(summary_trait_tbl[which(summary_trait_tbl$category=="timing"),"percent_significant"]),
                    "num_significant/768")
summary_table[5,]=c("model_1,both,hr",
                    median(summary_trait_tbl[which(summary_trait_tbl$category=="hr"),"num_significant"]),
                    mean(summary_trait_tbl[which(summary_trait_tbl$category=="hr"),"num_significant"]),
                    min(summary_trait_tbl[which(summary_trait_tbl$category=="hr"),"num_significant"]),
                    max(summary_trait_tbl[which(summary_trait_tbl$category=="hr"),"num_significant"]),
                    median(summary_trait_tbl[which(summary_trait_tbl$category=="hr"),"percent_significant"]),
                    mean(summary_trait_tbl[which(summary_trait_tbl$category=="hr"),"percent_significant"]),
                    min(summary_trait_tbl[which(summary_trait_tbl$category=="hr"),"percent_significant"]),
                    max(summary_trait_tbl[which(summary_trait_tbl$category=="hr"),"percent_significant"]),
                    "num_significant/768")
summary_table[6,]=c("model_1,both,sdb",
                    median(summary_trait_tbl[which(summary_trait_tbl$category=="sdb"),"num_significant"]),
                    mean(summary_trait_tbl[which(summary_trait_tbl$category=="sdb"),"num_significant"]),
                    min(summary_trait_tbl[which(summary_trait_tbl$category=="sdb"),"num_significant"]),
                    max(summary_trait_tbl[which(summary_trait_tbl$category=="sdb"),"num_significant"]),
                    median(summary_trait_tbl[which(summary_trait_tbl$category=="sdb"),"percent_significant"]),
                    mean(summary_trait_tbl[which(summary_trait_tbl$category=="sdb"),"percent_significant"]),
                    min(summary_trait_tbl[which(summary_trait_tbl$category=="sdb"),"percent_significant"]),
                    max(summary_trait_tbl[which(summary_trait_tbl$category=="sdb"),"percent_significant"]),
                    "num_significant/768")
write.csv(summary_table,file="output/supplementary_tables3.csv",row.names=F)

# Supplementary Table S4
# Compare model 1 and 2 results among only non-timing sleep traits
# check if the model 2 effect estiamtes is within the 95% CI of model 1 effect estimates
atlas_association_both_model1<-comb_df_annot[comb_df_annot$strata=="both"&comb_df_annot$model=="model_1",]
atlas_association_both_model2<-comb_df_annot[comb_df_annot$strata=="both"&comb_df_annot$model=="model_2",]
atlas_association_both_model1$CI_lower <- atlas_association_both_model1$beta - 1.96 * atlas_association_both_model1$se
atlas_association_both_model1$CI_upper <- atlas_association_both_model1$beta + 1.96 * atlas_association_both_model1$se
merged_atlas_association_both <- merge(atlas_association_both_model1, atlas_association_both_model2, by = c("trait", "metabolite"), suffixes = c("_m1", "_m2"))
merged_atlas_association_both$within_CI <- with(merged_atlas_association_both, beta_m2 >= CI_lower & beta_m2 <= CI_upper)
write.csv(merged_atlas_association_both[which(merged_atlas_association_both$within_CI==FALSE),],file="output/supplementary_tables4.csv",row.names=F)


# Supplementary Table S6
# create a summary table based on super pathway
summary_metabolite_tbl<-comb_df_annot[which(comb_df_annot$strata=="both"&comb_df_annot$model=="model_1"),] %>%
  dplyr::group_by(model, strata, super_pathway,metabolite) %>%
  dplyr::summarise(
    num_significant = sum(p_val_fdr < 0.05, na.rm = TRUE),
    percent_significant = mean(p_val_fdr < 0.05, na.rm = TRUE) * 100/40
  )
  summary_superpathway_tbl <- summary_metabolite_tbl %>%
    dplyr::group_by(super_pathway) %>%
    dplyr::summarise(
      sum_num_significant = sum(num_significant),
      num_unique_metabolites = n_distinct(metabolite),
      min_num_significant = min(num_significant),
      max_num_significant = max(num_significant),
      avg_num_significant=sum_num_significant/num_unique_metabolites
    )
write.csv(summary_superpathway_tbl,file="output/supplementary_tables6.csv",row.names=F)

 # Supplementary Table S7 
 # make a summary table for sub_pathway and sleep trait category
 summary_subpathway_category_tbl <- comb_df_annot[which(comb_df_annot$model=="model_1"&comb_df_annot$strata=="both"),] %>%
   merge(.,ref_table,by="trait",all.x = T)%>%
   dplyr::group_by(super_pathway,sub_pathway, category) %>%
   dplyr::summarise(
     num_total = n(),
     num_significant = sum(p_val_fdr < 0.05),
     percent_significant = (num_significant / num_total) * 100
   )
 write.csv(summary_subpathway_category_tbl,file="output/supplementary_tables7.csv",row.names=F)
 
# Supplementary Table S8
 # make a summary table for super_pathway and sleep trait category
  summary_superpathway_category_tbl <- comb_df_annot[which(comb_df_annot$model=="model_1"&comb_df_annot$strata=="both"),] %>%
    merge(.,ref_table,by="trait",all.x = T)%>%
    dplyr::group_by(super_pathway, category) %>%
    dplyr::summarise(
      num_total = n(),
      num_significant = sum(p_val_fdr < 0.05),
      percent_significant = (num_significant / num_total) * 100
    )
  write.csv(summary_superpathway_category_tbl,file="output/supplementary_tables8.csv",row.names=F) 
 

# Supplementary Table s9
# Bipartite network property metrics

both_md1_sig<-merge(both_md1_sig,hchs_annot[,c("chemical_name","metabolite")],by="metabolite",all.x=T)
both_md1_web<-frame2webs(both_md1_sig,varnames = c("chemical_name","label","super_pathway"))
both_md1_sig$freq_char="one"
both_md1_uni<-frame2webs(both_md1_sig,varnames = c("chemical_name","label","freq_char"))

both_md1_array<-frame2webs(both_md1_sig,varnames = c("chemical_name","label","super_pathway"),type.out="array")
both_md1_web.full<-frame2webs(both_md1_sig,varnames = c("chemical_name","label","super_pathway"),emptylist = F)
both_md1_array.full<-frame2webs(both_md1_sig,varnames = c("chemical_name","label","super_pathway"),type.out="array",emptylist = F)
# calculate indices describing nework topography
network_ttl<-networklevel(both_md1_uni$one)
network_lipid<-networklevel(both_md1_web$`Lipid`)
network_aa<-networklevel(both_md1_web$`Amino Acid`)
network_carbohydrate<-networklevel(both_md1_web$`Carbohydrate`)
network_vitamin<-networklevel(both_md1_web$`Cofactors and Vitamins`)
network_energy<-networklevel(both_md1_web$`Energy`)
network_nucleotide<-networklevel(both_md1_web$`Nucleotide`)
network_partial<-networklevel(both_md1_web$`Partially Characterized Molecules`)
network_peptide<-networklevel(both_md1_web$`Peptide`)
network_xenobiotics<-networklevel(both_md1_web$`Xenobiotics`)
V<-networklevel(both_md1_web$`Unnamed`)
network_output<-cbind(as.data.frame(network_lipid),as.data.frame(network_aa),as.data.frame(network_carbohydrate),
                      as.data.frame(network_vitamin),as.data.frame(network_energy),as.data.frame(network_nucleotide),
                      as.data.frame(network_partial),as.data.frame(network_peptide),as.data.frame(network_xenobiotics),
                      as.data.frame(network_unnamed),as.data.frame(network_ttl))
write.csv(network_output,file="output/supplementary_tables9.csv")

############################
# visualize the bipartite network measures
# drop the partially characterized metabolites
network_output<-network_output[,!colnames(network_output)%in%c("network_partial")]
network_t<-as.data.frame(t(network_output))
network_t$super_pathway<-c("lipid","amino_acid","carbohydrate","cofactors_vitamins","energy","nucleotide","peptide","xenobiotics","unnamed","total")
# convert wide table to long table for selected cols
network_long1<-tidyr::pivot_longer(network_t[,c("super_pathway","connectance","cluster coefficient","modularity Q","weighted nestedness","number.of.species.LL")],
                                     # cols=c("connectance","cluster coefficient","modularity Q","weighted nestedness"),
                                     cols=c("connectance","cluster coefficient","weighted nestedness"),
                                     names_to="measure",
                                     values_to="value")%>%
  arrange(number.of.species.LL) # sort the table based on the number of metabolites
# create a custom sorting order for the super_pathway 
pathway_order=unique(network_long1$super_pathway)

network_long2<-tidyr::pivot_longer(network_t[,c("super_pathway","NODF","Shannon diversity","linkage density","links per species","number.of.species.LL")],
                                   cols=c("NODF","Shannon diversity","linkage density","links per species"),
                                   names_to="measure",
                                   values_to="value")%>%
  arrange(number.of.species.LL) # sort the table based on the number of metabolites

# Supplementary Figure S9
network_properties_plot1<-ggplot(network_long1,aes(x=factor(super_pathway,levels=pathway_order),y=value,color=measure, group = measure))+
  theme_bw()+geom_line(lwd=1)+xlab("Super Pathway")+
  # scale_y_continuous(name = "range_1", limits = c(0, 1), sec.axis = sec_axis(~ . * 100, name = "range_2")) +
  theme(plot.title = element_text(hjust = 1),axis.text.x = element_text(size = 10))+
  scale_color_brewer(palette="Set1")
network_properties_plot2<-ggplot(network_long2,aes(x=factor(super_pathway,levels=pathway_order),y=value,color=measure, group = measure))+
  theme_bw()+geom_line(lwd=1)+xlab("Super Pathway")+
  # scale_y_continuous(name = "range_1", limits = c(0, 1), sec.axis = sec_axis(~ . * 100, name = "range_2")) +
  theme(plot.title = element_text(hjust = 1),axis.text.x = element_text(size = 10))+
  scale_color_brewer(palette="Set2")

network_properties_plot<-ggpubr::ggarrange(network_properties_plot1, network_properties_plot2, ncol = 1)
# ggsave(network_properties_plot, file = "output/atlas/imputed/combine_batch/unknown/network_properties_plot.png", width = 12, height = 6, units = "in") # Adjust dimensions as needed
ggsave(network_properties_plot, file = "output/supplementary_figures9.png", width = 12, height = 6, units = "in") 