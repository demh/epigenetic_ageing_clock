###########################################################################################
#########                                                                         #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             22/05/2018                                  #########
#########                              EMBL-EBI                                   #########
#########                           Thornton group                                #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####              Biological insights into the epigenetic ageing clock           #########
###########################################################################################
##### Study whether global hypomethylation in Sotos patients affects epigenetic age   #####
##### acceleration according to Horvath's clock.                                      #####
###########################################################################################
##### USAGE: manual                                                                   #####
###########################################################################################

###########################################################
##################### Dependencies ########################
###########################################################

library(data.table);
library(ggplot2);
library(ggpubr);

setwd('~/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_ageing_clock/paper/additional_analyses/hypomethylation/');


###########################################################
##################### Functions ###########################
###########################################################

## Function: inverse of the transforming function.
#
# linear_part: linear part of the regression model (i.e.
#             weighted average methylation based on sample data)
# a: adult age (for humans, 20 years)

F_inverse_transf <- function(linear_part, a=20){
  
  if(linear_part <= 0){
    
    DNAm_age <- exp(linear_part) * (a+1) - 1;
    
  }
  
  if(linear_part > 0) {
    
    DNAm_age <- linear_part * (a+1) + a;
    
  }
  
  return(DNAm_age);
  
}


############################################################
################## Running the pipeline ####################
############################################################

##### 1. Read all the data. #####

## Read the control data.

control_metadata <- as.data.frame(fread('../../syndromes_screen/controls_data_downstream.tsv')); # Small control
control_files <- list.files(path='../../syndromes_screen/horvath_clock_sites/blood_control/', full.names = TRUE);

for(i in 1:length(control_files)){
  
  temp_betas <- as.data.frame(fread(control_files[i]));
  rownames(temp_betas) <- temp_betas[,1];
  temp_betas <- temp_betas[,-1];
  
  if(i==1){
    control_betas <- temp_betas;
  }else{
    control_betas <- cbind(control_betas, temp_betas);
  }
}

control_betas <- as.data.frame(t(control_betas[,which(colnames(control_betas) %in% control_metadata$GEO_sample)])); # Keep only the betas for those samples that were selected after QC.
control_betas$GEO_sample <- rownames(control_betas); rownames(control_betas) <- NULL;
control_final <- merge(control_metadata, control_betas, by='GEO_sample');
rm(control_metadata, temp_betas, control_files);

## Read the information regarding Sotos DMPs. 

sotos_dmps <- as.data.frame(fread('~/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/differential_methylation/Sotos_DMPs_final.csv')); # 428266 probes
sotos_dmps_sig_hypo <- sotos_dmps[sotos_dmps$pval < (0.01/nrow(sotos_dmps)) & sotos_dmps$beta<0,]; # 15062 significant 'Hypo' Sotos DMPs
prop_hypo <- nrow(sotos_dmps_sig_hypo)/nrow(sotos_dmps); # 3.52% of all probes
n_site <- round(prop_hypo*353); # We will change 12 sites of the Horvath clock to make the 'Sotos-like' hypomethylation signature
median_hypo <- median(sotos_dmps_sig_hypo$beta);


##### 2. Calculate DNA methylation age (Horvath's clock) for the controls (raw). #####

horvath_coeffs <- as.data.frame(fread('../../../utils/AdditionalFile3.csv'))[,1:2];
control_betas_final <- control_betas[,-ncol(control_betas)];
rownames(control_betas_final) <- control_betas$GEO_sample;
control_betas_final <- control_betas_final[,match(horvath_coeffs$CpGmarker[-1], colnames(control_betas_final))];
all(colnames(control_betas_final) == horvath_coeffs$CpGmarker[-1]);
horvath_coeffs_m <- matrix(horvath_coeffs$Coefficient[-1], nrow=353, ncol=1);
control_predictions <- as.data.frame(as.matrix(control_betas_final) %*% horvath_coeffs_m + horvath_coeffs[1,2]);
control_predictions <- apply(control_predictions, 1, F_inverse_transf);
control_predictions <- data.frame(GEO_sample=names(control_predictions), DNAmAge_new=as.numeric(control_predictions));
control_final <- merge(control_final, control_predictions, by='GEO_sample');
plot(control_final$DNAmAge_noob, control_final$DNAmAge_new); # Check that the DNAmAge calculated without internal normalisation is similar --> OK


##### 3. Calculate the effects of 'Sotos-like' hypomethylation on the DNAmAge of the controls. #####

n_iterations <- 1000;
medians_EAA_with_CCC_controls_Sotos_like <- rep(NA, n_iterations); # Effect size: -0.18; 3.52% sites (12)  

for(i in 1:n_iterations){
  
  print(i);
  
  # Hypomethylate a subset of clock sites.
  
  control_betas_Sotos_like <- control_betas_final;
  set.seed(1+i); # Allow for randomness to happen between iterations
  control_betas_Sotos_like_after <- t(apply(control_betas_Sotos_like, 1, function(x){
    to_change <- sample(1:353,n_site);
    new_x <- x; 
    new_x[to_change] <- new_x[to_change] + median_hypo;
    return(new_x)}));
  control_betas_Sotos_like_after[control_betas_Sotos_like_after<0] <- 0; # Do not allow negative beta-values
  #print(head(control_betas_Sotos_like_after)[1:40]);
  
  # Calculate DNAmAge for the changed controls. 
  
  after_predictions <- as.data.frame(control_betas_Sotos_like_after %*% horvath_coeffs_m + horvath_coeffs[1,2]);
  after_predictions <- apply(after_predictions, 1, F_inverse_transf);
  after_predictions <- data.frame(GEO_sample=paste0(names(after_predictions), '_2'), DNAmAge_new=as.numeric(after_predictions));
  #print(head(after_predictions$DNAmAge_new));
  control_final_i <- control_final[,-ncol(control_final)]; control_final_i$Disease_status <- 'Sotos-like control';
  control_final_i$GEO_sample <- paste0(control_final_i$GEO_sample, '_2');
  control_final_i$Sample_name <- paste0(control_final_i$Sample_name, '_2');
  control_final_i <- merge(control_final_i, after_predictions, by='GEO_sample');
  control_final_merged_i <- rbind(control_final, control_final_i);
  
  # Calculate EAA when compared with unchanged control.
  
  lm_formula_int_new <- paste0('DNAmAge_new~Age_years+Sex+Gran+CD4T+CD8T+B+Mono+NK+',
                              paste(colnames(control_final_merged_i)[grep('PC',colnames(control_final_merged_i))], collapse='+'));
  lm_control_int_new <- lm(lm_formula_int_new, data=control_final_merged_i[which(control_final_merged_i$Disease_status=='Control'),]);
  control_final_merged_i$new_delta_int[control_final_merged_i$Disease_status=='Control'] <- lm_control_int_new$residuals;
  control_final_merged_i$new_delta_int[control_final_merged_i$Disease_status=='Sotos-like control'] <-
    control_final_merged_i$DNAmAge_new[control_final_merged_i$Disease_status=='Sotos-like control'] - as.numeric(predict(lm_control_int_new, 
                                                                                                            newdata=control_final_merged_i[control_final_merged_i$Disease_status=='Sotos-like control',]));
  #print(head(control_final_merged_i$new_delta_int[control_final_merged_i$Disease_status=='Sotos-like control']));
  medians_EAA_with_CCC_controls_Sotos_like[i] <- median(control_final_merged_i$new_delta_int[control_final_merged_i$Disease_status=='Sotos-like control']);
  
  # # Plot the results. 
  # 
  # my_palette <- c('grey', 'orange');
  # plots_scatterplot_DNAmAge_new <- ggplot(data=control_final_merged_i, aes(x=Age_years, y=DNAmAge_new, col=Disease_status)) +
  #   geom_point() + scale_color_manual(values=my_palette) +
  #   theme_classic() +
  #   theme(axis.text=element_text(size=12),
  #         axis.title=element_text(size=14,face="bold"),
  #         plot.title = element_text(hjust = 0.5),
  #         plot.subtitle = element_text(hjust = 0.5)) +
  #   xlab("Chronological age (years)") + ylab("DNAmAge (years)") +
  #   labs(title = paste0("Control: N=", sum(control_final_merged_i$Disease_status=='Control'), "\n",
  #                       "'Sotos-like' control: N=", sum(control_final_merged_i$Disease_status=='Sotos-like control')),
  #        subtitle=paste0()) +
  #   guides(col=guide_legend(title="Disease status")) +
  #   xlim(c(-1,55)) + ylim(c(-1,85)) + geom_abline(slope=1, intercept=0, linetype=2);
  # plots_scatterplot_DNAmAge_new;
  # 
  # mc_1 <- list(c("Sotos-like control", "Control"));
  # plots_delta_int_new <- ggboxplot(control_final_merged_i, x = "Disease_status", y = "new_delta_int", color = "Disease_status", palette = my_palette, fill="Disease_status", alpha=0.5) +
  #   theme_classic() +
  #   theme(axis.text=element_text(size=12, angle=90),
  #         axis.title=element_text(size=14,face="bold")) +
  #   xlab("") + ylab("EAA with CCC (years)") +
  #   ylim(c(-50,50)) +
  #   geom_abline(slope=0, intercept=0, linetype=2) +
  #   stat_compare_means(comparisons = mc_1) +
  #   guides(colour=FALSE, fill=FALSE);
  # plots_delta_int_new;

}


##### 4. Calculate the effects of a more extreme case of hypomethylation on the DNAmAge of the controls. #####

n_iterations <- 1000;
medians_EAA_with_CCC_controls_extreme <- rep(NA, n_iterations); # Effect size: -0.5; 50% sites (176)

for(i in 1:n_iterations){
  
  print(i);
  
  # Hypomethylate a subset of clock sites.
  
  control_betas_Sotos_like <- control_betas_final;
  set.seed(1+i); # Allow for randomness to happen between iterations
  control_betas_Sotos_like_after <- t(apply(control_betas_Sotos_like, 1, function(x){
    to_change <- sample(1:353,176);
    new_x <- x; 
    new_x[to_change] <- new_x[to_change] - 0.5;
    return(new_x)}));
  control_betas_Sotos_like_after[control_betas_Sotos_like_after<0] <- 0; # Do not allow negative beta-values
  #print(head(control_betas_Sotos_like_after)[1:40]);
  
  # Calculate DNAmAge for the changed controls. 
  
  after_predictions <- as.data.frame(control_betas_Sotos_like_after %*% horvath_coeffs_m + horvath_coeffs[1,2]);
  after_predictions <- apply(after_predictions, 1, F_inverse_transf);
  after_predictions <- data.frame(GEO_sample=paste0(names(after_predictions), '_2'), DNAmAge_new=as.numeric(after_predictions));
  #print(head(after_predictions$DNAmAge_new));
  control_final_i <- control_final[,-ncol(control_final)]; control_final_i$Disease_status <- 'Sotos-like control';
  control_final_i$GEO_sample <- paste0(control_final_i$GEO_sample, '_2');
  control_final_i$Sample_name <- paste0(control_final_i$Sample_name, '_2');
  control_final_i <- merge(control_final_i, after_predictions, by='GEO_sample');
  control_final_merged_i <- rbind(control_final, control_final_i);
  
  # Calculate EAA when compared with unchanged control.
  
  lm_formula_int_new <- paste0('DNAmAge_new~Age_years+Sex+Gran+CD4T+CD8T+B+Mono+NK+',
                               paste(colnames(control_final_merged_i)[grep('PC',colnames(control_final_merged_i))], collapse='+'));
  lm_control_int_new <- lm(lm_formula_int_new, data=control_final_merged_i[which(control_final_merged_i$Disease_status=='Control'),]);
  control_final_merged_i$new_delta_int[control_final_merged_i$Disease_status=='Control'] <- lm_control_int_new$residuals;
  control_final_merged_i$new_delta_int[control_final_merged_i$Disease_status=='Sotos-like control'] <-
    control_final_merged_i$DNAmAge_new[control_final_merged_i$Disease_status=='Sotos-like control'] - as.numeric(predict(lm_control_int_new, 
                                                                                                                         newdata=control_final_merged_i[control_final_merged_i$Disease_status=='Sotos-like control',]));
  #print(head(control_final_merged_i$new_delta_int[control_final_merged_i$Disease_status=='Sotos-like control']));
  medians_EAA_with_CCC_controls_extreme[i] <- median(control_final_merged_i$new_delta_int[control_final_merged_i$Disease_status=='Sotos-like control']);
  
}


##### 5. Compare the results against the EAA observed in Sotos syndrome patients. #####

cases_metadata <- as.data.frame(fread('../../syndromes_screen/cases_data_downstream.tsv'));
cases_Sotos <- cases_metadata[cases_metadata$Disease_status=='Sotos',];
print(median(cases_Sotos$delta_int));

all_median_EAA <- data.frame(Median=c(medians_EAA_with_CCC_controls_Sotos_like, medians_EAA_with_CCC_controls_extreme),
                             Group=rep(c("'Sotos-like' hypomethylation \nEffect size: -0.18; 3.52% CpG sites",
                                         'Extreme hypomethylation \nEffect size: -0.50; 50.0% CpG sites'), 
                                       each=n_iterations),
                             Median_median=rep(c(median(medians_EAA_with_CCC_controls_Sotos_like), 
                                               median(medians_EAA_with_CCC_controls_extreme)),
                                               each=n_iterations));

effects_hypo <- ggplot(data=all_median_EAA, aes(x=Median, fill=Group, col=Group))+
  geom_density(alpha=0.7) +   
  geom_vline(data=all_median_EAA, aes(xintercept=Median_median, color=Group), linetype="dashed") +
  geom_vline(xintercept=median(cases_Sotos$delta_int), linetype="dashed", col='orange') +
  scale_color_manual(values=c("red", "blue", "orange")) + 
  scale_fill_manual(values=c("red", "blue")) +
  ylab('Density') + xlab('Median EAA with CCC (years)') +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5));
ggsave('plots/effects_hypomethylation_on_EAA.pdf', height=6, width=9);

################################################################
################## End of the script ###########################
################################################################
