###########################################################################################
#########                                                                         #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             13/09/2018                                  #########
#########                              EMBL-EBI                                   #########
#########                           Thornton group                                #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####              Biological insights into the epigenetic ageing clock           #########
###########################################################################################
##### Perform the screening for the developmental disorders (i.e. check whether they  #####
##### have an accelerated epigenetic age).                                            #####
##### NOTE: 'delta' is another name for 'Epigenetic Age Acceleration' (EAA).          #####
#####       'ccc' stands for cell composition correction.                             #####
#####       'Extrinsic + intrinsic delta' is equivalent to 'EAA without ccc'.         #####
#####       'Intrinsic delta' is equivalent to 'EAA with ccc'.                        #####
###########################################################################################
##### USAGE: manual                                                                   #####
###########################################################################################

###########################################################
##################### Dependencies ########################
###########################################################

library(data.table);
library(ggplot2);
library(ggpubr);
library(ggthemes);
library(gridExtra);

setwd('~/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/syndromes_screen/');


################################################################
################## Run the pipeline ############################
################################################################

##### 1. Read the curated data for cases and controls. ##### 

raw_controls <- as.data.frame(fread('final_control_data.tsv'));
raw_cases <- as.data.frame(fread('final_cases_data.tsv'));


##### 2. Process the data. ##### 

## Merge relevant cases.

raw_cases$Disease_status[raw_cases$Disease_status=='LEOPARD_PTPN11' | raw_cases$Disease_status=='Noonan_PTPN11'] <- 'Noonan_PTPN11';

## Filter the data for the cases. 

n_thr1 <- 5; # Minimum number of samples for a developmental disorder
age_thr2 <- 20; 
n_age_thr2 <- 2; # Minimum number of samples for a developmental disorder with age >= age_thr2

diseases_thr1 <- names(table(raw_cases$Disease_status))[table(raw_cases$Disease_status) >= n_thr1];
print(paste0('Discarding ', names(table(raw_cases$Disease_status))[table(raw_cases$Disease_status) < n_thr1], '...'));
raw_cases <- raw_cases[which(raw_cases$Disease_status %in% diseases_thr1),];

diseases_thr2 <- c();

for(d in unique(raw_cases$Disease_status)){
  
  temp_data <- raw_cases[which(raw_cases$Disease_status==d),];
  temp_n_age <- sum(temp_data$Age_years >= age_thr2);
  if(temp_n_age >= n_age_thr2){
    diseases_thr2 <- c(diseases_thr2,d);
  }else{
    print(paste0('Discarding ', d, '...'));
  }
}

raw_cases <- raw_cases[which(raw_cases$Disease_status %in% diseases_thr2),];

#additional <- which(raw_cases$Disease_status=='ASD' | raw_cases$Disease_status=='FFS' | raw_cases$Disease_status=='Undiagnosed_DD/ID'); # Include in the analysis these disease groups even though there genes are not clearly defined (and have NA for Pathogenic column)
additional <- which(raw_cases$Disease_status=='ASD');
cases_data_final_all <- raw_cases[c(additional, which(raw_cases$Pathogenic=='YES' | raw_cases$Pathogenic=='YES_predicted')), ];
#cases_data_final_all <- raw_cases[c(additional, which(raw_cases$Pathogenic=='YES')), ];
cases_data_final_all$delta_ext_int <- NA;
cases_data_final_all$delta_int <- NA;

## Select age distribution for controls. 

controls_data_final <- raw_controls;
control_max_age <- max(cases_data_final_all$Age_years); # Maximum age in cases: 55 years
controls_data_final <- controls_data_final[which(controls_data_final$Age_years<=control_max_age),]; # Take only controls in 0-55 years range, which is the age range for cases. 

## Fit linear models to controls.

lm_formula_ext_int <- paste0('DNAmAge_noob~Age_years+Sex+', paste(colnames(controls_data_final)[grep('PC',colnames(controls_data_final))], collapse='+'));
lm_formula_int <- paste0('DNAmAge_noob~Age_years+Sex+Gran+CD4T+CD8T+B+Mono+NK+', paste(colnames(controls_data_final)[grep('PC',colnames(controls_data_final))], collapse='+'));

lm_control_ext_int <- lm(lm_formula_ext_int, data=controls_data_final);
lm_control_int <- lm(lm_formula_int, data=controls_data_final);
controls_data_final$delta_ext_int <- lm_control_ext_int$residuals;
controls_data_final$delta_int <- lm_control_int$residuals;


##### 3. Perform the main screening. ##### 

all_diseases <- unique(cases_data_final_all$Disease_status)[order(unique(cases_data_final_all$Disease_status))];
my_palette <- c('orange', 'grey');
plots_delta_ext_int <- list();
plots_delta_int <- list();
plots_scatterplot_DNAmAge <- list();
plots_scatterplot_delta_ext_int <- list();
plots_scatterplot_delta_int <- list();
p_values_ext_int <- c();
p_values_int <- c();
deltas_ext_int_cases <- c();
deltas_int_cases <- c();
disease_vector_deltas <- c();

for(d in all_diseases){
  
  print(d);
  
  # Calculate the deltas for the disease cases and store them.
  
  disease_vector <- c(controls_data_final$Disease_status, cases_data_final_all$Disease_status[which(cases_data_final_all$Disease_status==d)]);
  temp_data <- data.frame(Disease_status=disease_vector,
                          Age_years=c(controls_data_final$Age_years, cases_data_final_all$Age_years[which(cases_data_final_all$Disease_status==d)]),
                          DNAmAge_noob = c(controls_data_final$DNAmAge_noob, cases_data_final_all$DNAmAge_noob[which(cases_data_final_all$Disease_status==d)]),
                          delta_ext_int=c(controls_data_final$delta_ext_int, rep(NA, length(cases_data_final_all$Disease_status[which(cases_data_final_all$Disease_status==d)]))),
                          delta_int=c(controls_data_final$delta_int, rep(NA, length(cases_data_final_all$Disease_status[which(cases_data_final_all$Disease_status==d)]))));
  temp_data$Disease_status <- factor(temp_data$Disease_status, levels=c(d, "Control")); # Right order for plotting
  temp_data$delta_ext_int[temp_data$Disease_status==d] <- cases_data_final_all$delta_ext_int[which(cases_data_final_all$Disease_status==d)] <-
    cases_data_final_all$DNAmAge_noob[cases_data_final_all$Disease_status==d] - as.numeric(predict(lm_control_ext_int, newdata=cases_data_final_all[cases_data_final_all$Disease_status==d,]));
  deltas_ext_int_cases <- c(deltas_ext_int_cases, temp_data$delta_ext_int[temp_data$Disease_status==d]);
  temp_data$delta_int[temp_data$Disease_status==d] <- cases_data_final_all$delta_int[which(cases_data_final_all$Disease_status==d)] <-
    cases_data_final_all$DNAmAge_noob[cases_data_final_all$Disease_status==d] - as.numeric(predict(lm_control_int, newdata=cases_data_final_all[cases_data_final_all$Disease_status==d,]));
  deltas_int_cases <- c(deltas_int_cases, temp_data$delta_int[temp_data$Disease_status==d]);
  disease_vector_deltas <- c(disease_vector_deltas,cases_data_final_all$Disease_status[which(cases_data_final_all$Disease_status==d)]);
  
 # Compare the deltas for cases and controls and store the p-values. 
  
  mc_1 <- list(c(d, "Control"));

  plots_delta_ext_int[[d]] <- ggboxplot(temp_data, x = "Disease_status", y = "delta_ext_int", color = "Disease_status", palette = my_palette, fill="Disease_status", alpha=0.5) +
    theme_classic() +
    theme(axis.text=element_text(size=12, angle=90),
          axis.title=element_text(size=14,face="bold")) +
    xlab("") + ylab("EAA without CCC (years)") +
    ylim(c(-50,50)) +
    geom_abline(slope=0, intercept=0, linetype=2) +
    stat_compare_means(comparisons = mc_1) +
    guides(colour=FALSE, fill=FALSE);
  
  p_values_ext_int <- c(p_values_ext_int, 
                        wilcox.test(x=temp_data$delta_ext_int[temp_data$Disease_status=='Control'], y=temp_data$delta_ext_int[temp_data$Disease_status==d])$p.value);
  
  plots_delta_int[[d]] <- ggboxplot(temp_data, x = "Disease_status", y = "delta_int", color = "Disease_status", palette = my_palette, fill="Disease_status", alpha=0.5) +
    theme_classic() +
    theme(axis.text=element_text(size=12, angle=90),
          axis.title=element_text(size=14,face="bold")) +
    xlab("") + ylab("EAA with CCC (years)") +
    ylim(c(-50,50)) +
    geom_abline(slope=0, intercept=0, linetype=2) +
    stat_compare_means(comparisons = mc_1) +
    guides(colour=FALSE, fill=FALSE);
  
  p_values_int <- c(p_values_int,
                    wilcox.test(x=temp_data$delta_int[temp_data$Disease_status=='Control'], y=temp_data$delta_int[temp_data$Disease_status==d])$p.value);
  
  # Create scatterplots for DNAmAge and EAAs vs chronological age. 
  
  plots_scatterplot_DNAmAge[[d]] <- ggplot(data=temp_data, aes(x=Age_years, y=DNAmAge_noob, col=Disease_status)) + 
    geom_point() + scale_color_manual(values=my_palette) +
    theme_classic() +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) +
    xlab("Chronological age (years)") + ylab("DNAmAge (years)") +
    labs(title = paste0("Control: N=", sum(temp_data$Disease_status=='Control'), "\n",
                        d,": N=", sum(temp_data$Disease_status==d)), 
         subtitle=paste0()) + 
    guides(col=guide_legend(title="Disease status")) + 
    xlim(c(-1,55)) + ylim(c(-1,80)) + geom_abline(slope=1, intercept=0, linetype=2);

  plots_scatterplot_delta_ext_int[[d]] <- ggplot(data=temp_data[temp_data$Disease_status==d,], aes(x=Age_years, y=delta_ext_int)) + 
    geom_point(data=temp_data[temp_data$Disease_status=='Control',], aes(x=Age_years, y=delta_ext_int), col='grey') + 
    geom_abline(slope=0, intercept=0, linetype=2) +
    geom_smooth(method="lm", formula=y~x, show.legend=F, fill='khaki', alpha = 0.3, col='gold') + geom_point(col='orange') + 
    theme_classic() +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) +
    xlab("Chronological age (years)") + ylab("EAA without CCC (years)") +
    labs(title = paste0("Control: N=", sum(temp_data$Disease_status=='Control'), "\n",
                        d,": N=", sum(temp_data$Disease_status==d)), 
         subtitle=paste0()) +  
    xlim(c(-1,55)) + ylim(c(-30,40));
  
  plots_scatterplot_delta_int[[d]] <- ggplot(data=temp_data[temp_data$Disease_status==d,], aes(x=Age_years, y=delta_int)) + 
    geom_point(data=temp_data[temp_data$Disease_status=='Control',], aes(x=Age_years, y=delta_int), col='grey') + 
    geom_abline(slope=0, intercept=0, linetype=2) +
    geom_smooth(method="lm", formula=y~x, show.legend=F, fill='khaki', alpha = 0.3, col='gold') + geom_point(col='orange') + 
    theme_classic() +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) +
    xlab("Chronological age (years)") + ylab("EAA with CCC (years)") +
    labs(title = paste0("Control: N=", sum(temp_data$Disease_status=='Control'), "\n",
                        d,": N=", sum(temp_data$Disease_status==d)), 
         subtitle=paste0()) +  
    xlim(c(-1,55)) + ylim(c(-30,40));
  
  if(d=='Sotos'){Sotos_data <- temp_data}; # Store data for Sotos
  
}


## Export the results of the screening for downstream analysis. ##

write.table(x=controls_data_final, file='controls_data_downstream.tsv', quote=F, sep='\t', row.names=F);
write.table(x=cases_data_final_all, file='cases_data_downstream.tsv', quote=F, sep='\t', row.names=F);

## Export the results for the supplementary information in the paper. ##

controls_supp <- controls_data_final;
controls_supp <- controls_supp[,c(1:9, 17, 36, 37, 18, 10, 11:16, 19:35)];
colnames(controls_supp)[c(1,10:13)] <- c('Sample_ID','DNAmAge_years', 'EAA_without_CCC', 'EAA_with_CCC',
                                         'meanMethBySample');
controls_supp <- controls_supp[order(controls_supp$Disease_status, controls_supp$Batch, controls_supp$Sample_ID),];
write.table(x=controls_supp, file='Additional_file_3.tsv', 
            quote=F, sep='\t', row.names=F);

cases_supp <- cases_data_final_all;
cases_supp <- cases_supp[,c(1:9, 24, 43, 44, 25, 10:23, 26:42)];
colnames(cases_supp)[c(1,10:13)] <- c('Sample_ID','DNAmAge_years', 'EAA_without_CCC', 'EAA_with_CCC',
                                      'meanMethBySample');
cases_supp <- cases_supp[order(cases_supp$Disease_status, cases_supp$Batch, cases_supp$Sample_ID),];
write.table(x=cases_supp, file='Additional_file_2.tsv', 
            quote=F, sep='\t', row.names=F);


##### 4. Create plots to summarise the main screening. ##### 

## Plot the p-values for all the diseases.

overall_p_values_df <- data.frame(Disease_status=rep(all_diseases, 2), 
                                  type=rep(c('Without CCC', 'With CCC'), each=length(all_diseases)),
                                  log10=c(-log10(p_values_ext_int), -log10(p_values_int)));
alpha <- -log10(0.01/(length(all_diseases)));

plot_overall_p_values <- ggplot(data=overall_p_values_df, aes(x=Disease_status, y=log10, fill=type)) + 
  geom_bar(position="dodge", stat="identity") + scale_fill_manual(values=c('red', 'blue')) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.text.x=element_text(angle = -90, hjust = 0, size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5, face="bold"),
        plot.subtitle = element_text(hjust = 0.5),
        plot.margin=unit(c(1,1,-0.5,1), "cm")) +
  xlab("") + ylab(expression(bold(paste(-log[10](P-value))))) +
  labs(title = "", 
       subtitle=paste0("Age range in control: 0-", control_max_age, 
                       " years\nMedian age in control: ", median(controls_data_final$Age_years), 
                       " years \nNumber of samples in control: ", nrow(controls_data_final))) + 
  guides(fill=guide_legend(title="EAA model")) + 
  geom_abline(slope=0, intercept=alpha, linetype=2, col='forestgreen', size=0.8);
ggsave("plots/plot_overall_p_values.pdf", height=6, width=5);


## Plot the deltas for all the diseases.

overall_deltas_df <- data.frame(Disease_status=c(disease_vector_deltas, disease_vector_deltas),
                                type=rep(c('Without CCC', 'With CCC'), each=length(disease_vector_deltas)),
                                deltas=c(deltas_ext_int_cases, deltas_int_cases));

plot_overall_deltas <- ggplot(data=overall_deltas_df, aes(x=Disease_status, y=deltas, fill=type)) + 
  geom_boxplot() + scale_fill_manual(values=c('red', 'blue')) +
  theme_classic() +
  scale_y_continuous(breaks=seq(-40,40,10)) +
  scale_x_discrete(position = "top") + 
  theme(axis.text=element_text(size=12),
        axis.text.x=element_blank(), 
        #axis.text.x=element_text(angle = 90, hjust = 0),
        axis.title=element_text(size=14,face="bold"),
        plot.margin=unit(c(-0.5,1,1,1), "cm")) + 
  xlab("") + ylab("Epigenetic age acceleration (years)") +
  guides(fill=guide_legend(title="EAA model")) + 
  geom_abline(slope=0, intercept=0, linetype=2, size=0.25)+
  #geom_abline(slope=0, intercept=10, linetype=2, size=0.25)+
  #geom_abline(slope=0, intercept=-10, linetype=2, size=0.25);
ggsave("plots/plot_overall_EEA.pdf", height=5, width=5);

pdf("plots/plot_overall_p_values_and_EEA.pdf", height=8, width=6);
grid.arrange(plot_overall_p_values, plot_overall_deltas, ncol=1);
dev.off();

## Scatterplots for Sotos. 

pdf("plots/plot_Sotos_DNAmAge.pdf", height=4, width=4.5);
plots_scatterplot_DNAmAge[['Sotos']];
dev.off();
pdf("plots/plot_Sotos_DNAmAge_no_legend.pdf", height=4, width=4);
plots_scatterplot_DNAmAge[['Sotos']] + theme(legend.position="none");
dev.off();
pdf("plots/plot_Sotos_EEA_without_CCC.pdf", height=4, width=4);
plots_scatterplot_delta_ext_int[['Sotos']];
dev.off();
pdf("plots/plot_Sotos_EEA_with_CCC.pdf", height=4, width=4);
plots_scatterplot_delta_int[['Sotos']];
dev.off();

## Check how EAA behaves for Sotos patients. 

only_Sotos_data <- Sotos_data[Sotos_data$Disease_status=='Sotos',];
summary(lm(delta_ext_int~Age_years, data=only_Sotos_data)); # EAA without CCC, all data, p-value=0.00514
summary(lm(delta_ext_int~Age_years, data=only_Sotos_data[only_Sotos_data$Age_years!=41,])); # EAA without CCC, remove outlier, p-value=0.1087
summary(lm(delta_int~Age_years, data=only_Sotos_data)); # EAA with CCC, all data, p-value=0.00569
summary(lm(delta_int~Age_years, data=only_Sotos_data[only_Sotos_data$Age_years!=41,])); # EAA with CCC, remove outlier, p-value=0.1785

## Scatterplots for all the disorders in screening (supplementary).

supp_scatterplots_list <- list();

z <- 1;

for(i in 1:length(all_diseases)){
  supp_scatterplots_list[[z]] <- plots_scatterplot_DNAmAge[[i]];
  supp_scatterplots_list[[z+1]] <- plots_scatterplot_delta_ext_int[[i]];
  supp_scatterplots_list[[z+2]] <- plots_scatterplot_delta_int[[i]];
  z <- z+3;
}

glist <- lapply(supp_scatterplots_list, ggplotGrob)
ggsave("plots/all_scatterplots_disorders.pdf", marrangeGrob(glist, nrow=5, ncol=3, 
                                layout_matrix=matrix(1:15, nrow=5, ncol=3, byrow=T), top=NULL),
       height=16, width=12);


#### 5. Check the influence of the age range in the controls on the result of the screening. 

max_ages <- seq(5,100, by=1); # This can be changed

diseases_vector <- c();
median_d_age_vector <- c();
n_d_vector <- c();
max_control_ages_vector <- c();
median_control_age_vector <- c();
n_control_vector <- c();
p_values_ext_int_vector <- c(); 
p_values_int_vector <- c();

for (d in all_diseases){
  
  print(paste0('Checking disease ', d, '...'));
  
  # Get the info for the disease.
  
  temp_d <- cases_data_final_all[which(cases_data_final_all$Disease_status==d),];
  diseases_vector <- c(diseases_vector, rep(d, length(max_ages)));
  median_d_age_vector <- c(median_d_age_vector, rep(median(temp_d$Age_years), length(max_ages)));
  n_d_vector <- c(n_d_vector, rep(nrow(temp_d), length(max_ages)));
  
  # Loop over the different age ranges for control. 
  
  for(max_age in max_ages){
    
    #print(paste0('Max. age in control: ', max_age));
    max_control_ages_vector <- c(max_control_ages_vector, max_age);
    
    temp_control <- raw_controls[raw_controls$Age_years <= max_age,];
    median_control_age_vector <- c(median_control_age_vector, median(temp_control$Age_years));
    n_control_vector <- c(n_control_vector, nrow(temp_control));
    
    # Fit linear model to control subset. 
    
    lm_temp_ext_int <- lm(lm_formula_ext_int, data=temp_control);
    lm_temp_int <- lm(lm_formula_int, data=temp_control);
    
    # Obtain deltas for the cases (EAAs).
    
    temp_d_ext_int <- temp_d$DNAmAge_noob - as.numeric(predict(lm_temp_ext_int, newdata=temp_d));
    temp_d_int <- temp_d$DNAmAge_noob - as.numeric(predict(lm_temp_int, newdata=temp_d));
    
    # Compare against control and store the p-values.
    
    p_values_ext_int_vector <- c(p_values_ext_int_vector, 
                                 wilcox.test(x=lm_temp_ext_int$residuals, y=temp_d_ext_int)$p.value);
    p_values_int_vector <- c(p_values_int_vector, 
                             wilcox.test(x=lm_temp_int$residuals, y=temp_d_int)$p.value);
  }
}

## Make the plots. 

ar_df <- data.frame(cbind(diseases_vector, median_d_age_vector, n_d_vector, max_control_ages_vector,
                          median_control_age_vector, n_control_vector, p_values_ext_int_vector, 
                          p_values_int_vector), stringsAsFactors=F);
ar_plots <- list();

for(d in all_diseases){
  
  temp_ar <- data.frame(median_age_controls=rep(as.numeric(ar_df$median_control_age_vector[ar_df$diseases_vector==d]),2),
                        p_values=c(as.numeric(ar_df$p_values_ext_int_vector[ar_df$diseases_vector==d]),
                                   as.numeric(ar_df$p_values_int_vector[ar_df$diseases_vector==d])),
                        model_type=c(rep(c('Without CCC', 'With CCC'), each=length(max_ages)))); 
  md <- as.numeric(unique(ar_df$median_d_age_vector[ar_df$diseases_vector==d]));
  n <- as.numeric(unique(ar_df$n_d_vector[ar_df$diseases_vector==d]));
  
  ar_plots[[d]] <- ggplot(data=temp_ar, aes(x=median_age_controls, y=-log10(p_values), col=model_type)) + 
    geom_point() + geom_line() + scale_color_manual(values=c('red', 'blue')) +
    theme_classic() +
    theme(axis.text=element_text(size=12),
          axis.text.x=element_text(hjust = 0),
          axis.title=element_text(size=14,face="bold"),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) +
    xlab("Median age in control (years)") + ylab(expression(bold(paste(-log[10](P-value))))) +
    labs(title = paste0(d, ' (N=', n, ')')) + 
    guides(col=guide_legend(title="EAA model")) +
    ylim(c(0,13)) +
    geom_vline(xintercept=md, linetype=2, col='orange', size=0.8) +
    geom_abline(slope=0, intercept=-log10(0.01/(length(all_diseases))), linetype=2, col='forestgreen', size=0.8);
}

glist <- lapply(ar_plots, ggplotGrob)
ggsave("plots/effect_median_age_control.pdf", marrangeGrob(glist, nrow=5, ncol=3, 
                                                            layout_matrix=matrix(1:15, nrow=5, ncol=3, byrow=T), top=NULL),
       height=18, width=13.5);


################################################################
################## End of the script ###########################
################################################################
