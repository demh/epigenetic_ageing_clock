###########################################################################################
#########                                                                         #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             06/05/2019                                  #########
#########                              EMBL-EBI                                   #########
#########                           Thornton group                                #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####              Biological insights into the epigenetic ageing clock           #########
###########################################################################################
##### Check whether Sotos patients present epigenetic age acceleration according to   #####
##### other epigenetic clocks (Hannum, SkinBlood, Lin).                               #####
###########################################################################################
##### USAGE: manual                                                                   #####
###########################################################################################

###########################################################
##################### Dependencies ########################
###########################################################

library(data.table);
library(ggplot2);
library(ggpubr);

setwd('~/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_ageing_clock/paper/additional_analyses/other_clocks');


###########################################################
##################### Functions ###########################
###########################################################

## Function: inverse of the transforming function
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


################################################################
################## Run the pipeline ############################
################################################################

##### 1. Hannum's clock. #####

## Read metadata for Sotos cases and controls. ##

metadata_cases <- as.data.frame(fread('cases_data_downstream_GSE74432.tsv'));
metadata_controls <- as.data.frame(fread('controls_data_downstream.tsv'));
all_metadata <- rbind(metadata_cases[,which(colnames(metadata_cases) %in% colnames(metadata_controls))], metadata_controls);
all_metadata <- all_metadata[which(all_metadata$Disease_status != 'Weaver'),]; # Remove Weaver samples

## Read beta-values for Sotos cases and controls. ##

hannum_betas_cases <- as.data.frame(fread('betas_for_probes_hannum_GSE74432.csv'));
hannum_betas_controls <- as.data.frame(fread('betas_for_probes_hannum_controls.csv'));
all(rownames(hannum_betas_cases) == rownames(hannum_betas_controls));
hannum_betas_all <- cbind(hannum_betas_cases, hannum_betas_controls[,-1]);
na_hannum <- which(is.na(hannum_betas_all$ProbeID) | hannum_betas_all$ProbeID=='NA.1' | hannum_betas_all$ProbeID=='NA.2');
hannum_betas_all <- hannum_betas_all[-na_hannum,];

## Calculate epigenetic ages according to Hannum's clock. ##

hannum_coeffs <- as.data.frame(fread('hannum_clock_coefs.tsv'));
hannum_coeffs <- hannum_coeffs[-na_hannum,];
if(!all(hannum_betas_all$ProbeID==hannum_coeffs$ProbeID)){stop('The probeIDs do not match')};

hannum_betas_f <- t(hannum_betas_all[,-1]);
colnames(hannum_betas_f) <- hannum_betas_all$ProbeID;
hannum_coeffs_f <- matrix(hannum_coeffs$Coefficient, nrow=nrow(hannum_coeffs), ncol=1);
rownames(hannum_coeffs_f) <- hannum_coeffs$ProbeID;
hannum_predictions <- as.data.frame(hannum_betas_f %*% hannum_coeffs_f);
hannum_predictions$GEO_sample <- rownames(hannum_predictions);
colnames(hannum_predictions)[1] <- 'hannum_age';
all_metadata <- merge(all_metadata, hannum_predictions, by='GEO_sample');
all_metadata$Disease_status <- as.factor(all_metadata$Disease_status);

## Scatterplot. ##

my_palette <- c('grey', 'orange');
scatterplot_hannum <- ggplot(data=all_metadata[all_metadata$Disease_status=='Control',], aes(x=Age_years, y=hannum_age)) +
  geom_point(aes(col='Control')) +
  geom_point(data=all_metadata[all_metadata$Disease_status=='Sotos',], aes(x=Age_years, y=hannum_age, col='Sotos')) +
  scale_colour_manual(name="Disease status",values=c("Control"='grey', "Sotos"='orange')) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  xlab("Chronological age (years)") + ylab("HannumAge (years)") +
  labs(title = paste0("Control: N=", sum(all_metadata$Disease_status=='Control'), "\n",
                      "Sotos: N=", sum(all_metadata$Disease_status=='Sotos')),
       subtitle=paste0()) +
  xlim(c(-1,55)) + ylim(c(-20,85)) + geom_abline(slope=1, intercept=0, linetype=2);
ggsave("plots/plot_Sotos_HannumAge.pdf", height=5, width=5.5);

## Calculate EAA according to Hannum's clock. ##

lm_formula_ext_int_hannum <- paste0('hannum_age~Age_years+Sex+', 
                                    paste(colnames(all_metadata)[grep('PC',colnames(all_metadata))], collapse='+'));
lm_formula_int_hannum <- paste0('hannum_age~Age_years+Sex+Gran+CD4T+CD8T+B+Mono+NK+', 
                                paste(colnames(all_metadata)[grep('PC',colnames(all_metadata))], collapse='+'));

lm_control_ext_int_hannum <- lm(lm_formula_ext_int_hannum, data=all_metadata[which(all_metadata$Disease_status=='Control'),]);
lm_control_int_hannum <- lm(lm_formula_int_hannum, data=all_metadata[which(all_metadata$Disease_status=='Control'),]);
all_metadata$hannum_delta_ext_int[all_metadata$Disease_status=='Control'] <- lm_control_ext_int_hannum$residuals;
all_metadata$hannum_delta_int[all_metadata$Disease_status=='Control'] <- lm_control_int_hannum$residuals;

all_metadata$hannum_delta_ext_int[all_metadata$Disease_status=='Sotos'] <-
  all_metadata$hannum_age[all_metadata$Disease_status=='Sotos'] - as.numeric(predict(lm_control_ext_int_hannum, newdata=all_metadata[all_metadata$Disease_status=='Sotos',]));
all_metadata$hannum_delta_int[all_metadata$Disease_status=='Sotos'] <-
  all_metadata$hannum_age[all_metadata$Disease_status=='Sotos'] - as.numeric(predict(lm_control_int_hannum, newdata=all_metadata[all_metadata$Disease_status=='Sotos',]));

## EAA boxplots. ##

mc_1 <- list(c("Sotos", "Control"));

plot_delta_ext_int_hannum <- ggboxplot(all_metadata, x = "Disease_status", y = "hannum_delta_ext_int", color = "Disease_status", palette = my_palette, 
                                       fill="Disease_status", alpha=0.5) +
  theme_classic() +
  theme(axis.text=element_text(size=12, angle=90),
        axis.title=element_text(size=14,face="bold")) +
  xlab("") + ylab("Hannum EAA without CCC (years)") +
  ylim(c(-50,50)) +
  geom_abline(slope=0, intercept=0, linetype=2) +
  stat_compare_means(comparisons = mc_1) +
  guides(colour=FALSE, fill=FALSE);
ggsave("plots/plot_Sotos_HannumAge_EAA_without_CCC.pdf", height=5, width=4.5);

wilcox.test(x=all_metadata$hannum_delta_ext_int[all_metadata$Disease_status=='Control'], 
            y=all_metadata$hannum_delta_ext_int[all_metadata$Disease_status=='Sotos'])$p.value;

plot_delta_int_hannum <- ggboxplot(all_metadata, x = "Disease_status", y = "hannum_delta_int", color = "Disease_status", palette = my_palette, 
                                       fill="Disease_status", alpha=0.5) +
  theme_classic() +
  theme(axis.text=element_text(size=12, angle=90),
        axis.title=element_text(size=14,face="bold")) +
  xlab("") + ylab("Hannum EAA with CCC (years)") +
  ylim(c(-50,50)) +
  geom_abline(slope=0, intercept=0, linetype=2) +
  stat_compare_means(comparisons = mc_1) +
  guides(colour=FALSE, fill=FALSE);
ggsave("plots/plot_Sotos_HannumAge_EAA_with_CCC.pdf", height=5, width=4.5);


wilcox.test(x=all_metadata$hannum_delta_int[all_metadata$Disease_status=='Control'], 
            y=all_metadata$hannum_delta_int[all_metadata$Disease_status=='Sotos'])$p.value;


##### 2. Skin-blood clock. #####

## Read beta-values for Sotos cases and controls. ##

sb_betas_cases <- as.data.frame(fread('betas_for_probes_skin_blood_GSE74432.csv'));
sb_betas_controls <- as.data.frame(fread('betas_for_probes_skin_blood_controls.csv'));
all(rownames(sb_betas_cases) == rownames(sb_betas_controls));
sb_betas_all <- cbind(sb_betas_cases, sb_betas_controls[,-1]);
na_sb <- which(is.na(sb_betas_all$ProbeID) | sb_betas_all$ProbeID=='NA.1' | sb_betas_all$ProbeID=='NA.2' |
                 sb_betas_all$ProbeID=='NA.3' | sb_betas_all$ProbeID=='NA.4' | sb_betas_all$ProbeID=='NA.5');
sb_betas_all <- sb_betas_all[-na_sb,];

## Calculate epigenetic ages according to skin-blood clock. ##

sb_coeffs <- as.data.frame(fread('skin_blood_clock_coeffs.csv'));
sb_intercept <- sb_coeffs[1,2]; sb_coeffs <- sb_coeffs[-1,];
sb_coeffs <- sb_coeffs[-na_sb,];
if(!all(sb_betas_all$ProbeID==sb_coeffs$ProbeID)){stop('The probeIDs do not match')};

sb_betas_f <- t(sb_betas_all[,-1]);
colnames(sb_betas_f) <- sb_betas_all$ProbeID;
sb_coeffs_f <- matrix(sb_coeffs$Coef, nrow=nrow(sb_coeffs), ncol=1);
rownames(sb_coeffs_f) <- sb_coeffs$ProbeID;
sb_predictions <- as.data.frame(sb_betas_f %*% sb_coeffs_f + sb_intercept);
sb_predictions <- apply(sb_predictions, 1, F_inverse_transf);
sb_predictions <- data.frame(GEO_sample=names(sb_predictions), SkinBloodAge=as.numeric(sb_predictions));
all_metadata <- merge(all_metadata, sb_predictions, by='GEO_sample');

## Scatterplot. ##

scatterplot_sb <- ggplot(data=all_metadata[all_metadata$Disease_status=='Control',], aes(x=Age_years, y=SkinBloodAge)) +
  geom_point(aes(col='Control')) +
  geom_point(data=all_metadata[all_metadata$Disease_status=='Sotos',], aes(x=Age_years, y=SkinBloodAge, col='Sotos')) +
  scale_colour_manual(name="Disease status",values=c("Control"='grey', "Sotos"='orange')) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  xlab("Chronological age (years)") + ylab("SkinBloodAge (years)") +
  labs(title = paste0("Control: N=", sum(all_metadata$Disease_status=='Control'), "\n",
                      "Sotos: N=", sum(all_metadata$Disease_status=='Sotos')),
       subtitle=paste0()) +
  xlim(c(-1,55)) + ylim(c(-5,85)) + geom_abline(slope=1, intercept=0, linetype=2);
ggsave("plots/plot_Sotos_SkinBloodAge.pdf", height=5, width=5.5);

## Calculate EAA according to skin-blood clock. ##

lm_formula_ext_int_sb <- paste0('SkinBloodAge~Age_years+Sex+', 
                                    paste(colnames(all_metadata)[grep('PC',colnames(all_metadata))], collapse='+'));
lm_formula_int_sb <- paste0('SkinBloodAge~Age_years+Sex+Gran+CD4T+CD8T+B+Mono+NK+', 
                                paste(colnames(all_metadata)[grep('PC',colnames(all_metadata))], collapse='+'));

lm_control_ext_int_sb <- lm(lm_formula_ext_int_sb, data=all_metadata[which(all_metadata$Disease_status=='Control'),]);
lm_control_int_sb <- lm(lm_formula_int_sb, data=all_metadata[which(all_metadata$Disease_status=='Control'),]);
all_metadata$sb_delta_ext_int[all_metadata$Disease_status=='Control'] <- lm_control_ext_int_sb$residuals;
all_metadata$sb_delta_int[all_metadata$Disease_status=='Control'] <- lm_control_int_sb$residuals;

all_metadata$sb_delta_ext_int[all_metadata$Disease_status=='Sotos'] <-
  all_metadata$SkinBloodAge[all_metadata$Disease_status=='Sotos'] - as.numeric(predict(lm_control_ext_int_sb, newdata=all_metadata[all_metadata$Disease_status=='Sotos',]));
all_metadata$sb_delta_int[all_metadata$Disease_status=='Sotos'] <-
  all_metadata$SkinBloodAge[all_metadata$Disease_status=='Sotos'] - as.numeric(predict(lm_control_int_sb, newdata=all_metadata[all_metadata$Disease_status=='Sotos',]));

## EAA boxplots. ##

plot_delta_ext_int_sb <- ggboxplot(all_metadata, x = "Disease_status", y = "sb_delta_ext_int", color = "Disease_status", palette = my_palette, 
                                       fill="Disease_status", alpha=0.5) +
  theme_classic() +
  theme(axis.text=element_text(size=12, angle=90),
        axis.title=element_text(size=14,face="bold")) +
  xlab("") + ylab("SkinBlood EAA without CCC (years)") +
  ylim(c(-50,50)) +
  geom_abline(slope=0, intercept=0, linetype=2) +
  stat_compare_means(comparisons = mc_1) +
  guides(colour=FALSE, fill=FALSE);
ggsave("plots/plot_Sotos_SkinBloodAge_EAA_without_CCC.pdf", height=5, width=4.5);

wilcox.test(x=all_metadata$sb_delta_ext_int[all_metadata$Disease_status=='Control'], 
            y=all_metadata$sb_delta_ext_int[all_metadata$Disease_status=='Sotos'])$p.value;

plot_delta_int_sb <- ggboxplot(all_metadata, x = "Disease_status", y = "sb_delta_int", color = "Disease_status", palette = my_palette, 
                                   fill="Disease_status", alpha=0.5) +
  theme_classic() +
  theme(axis.text=element_text(size=12, angle=90),
        axis.title=element_text(size=14,face="bold")) +
  xlab("") + ylab("SkinBlood EAA with CCC (years)") +
  ylim(c(-50,50)) +
  geom_abline(slope=0, intercept=0, linetype=2) +
  stat_compare_means(comparisons = mc_1) +
  guides(colour=FALSE, fill=FALSE);
ggsave("plots/plot_Sotos_SkinBloodAge_EAA_with_CCC.pdf", height=5, width=4.5);


wilcox.test(x=all_metadata$sb_delta_int[all_metadata$Disease_status=='Control'], 
            y=all_metadata$sb_delta_int[all_metadata$Disease_status=='Sotos'])$p.value;


##### 3. Lin clock. #####

## Read beta-values for Sotos cases and controls. ##

lin_betas_cases <- as.data.frame(fread('betas_for_probes_lin_GSE74432.csv'));
lin_betas_controls <- as.data.frame(fread('betas_for_probes_lin_controls.csv'));
all(rownames(lin_betas_cases) == rownames(lin_betas_controls));
lin_betas_all <- cbind(lin_betas_cases, lin_betas_controls[,-1]);
na_lin <- which(is.na(lin_betas_all$ProbeID) | lin_betas_all$ProbeID=='NA.1');
lin_betas_all <- lin_betas_all[-na_lin,];

## Calculate epigenetic ages according to Lin clock. ##

lin_coeffs <- as.data.frame(fread('lin_clock_coefs.tsv'));
lin_intercept <- lin_coeffs[1,2]; lin_coeffs <- lin_coeffs[-1,];
lin_coeffs <- lin_coeffs[-na_lin,];
if(!all(lin_betas_all$ProbeID==lin_coeffs$ProbeID)){stop('The probeIDs do not match')};

lin_betas_f <- t(lin_betas_all[,-1]);
colnames(lin_betas_f) <- lin_betas_all$ProbeID;
lin_coeffs_f <- matrix(lin_coeffs$Coef, nrow=nrow(lin_coeffs), ncol=1);
rownames(lin_coeffs_f) <- lin_coeffs$ProbeID;
lin_predictions <- as.data.frame(lin_betas_f %*% lin_coeffs_f + lin_intercept);
lin_predictions <- data.frame(GEO_sample=rownames(lin_predictions), LinAge=as.numeric(lin_predictions$V1));
all_metadata <- merge(all_metadata, lin_predictions, by='GEO_sample');

## Scatterplot. ##

scatterplot_lin <- ggplot(data=all_metadata[all_metadata$Disease_status=='Control',], aes(x=Age_years, y=LinAge)) +
  geom_point(aes(col='Control')) +
  geom_point(data=all_metadata[all_metadata$Disease_status=='Sotos',], aes(x=Age_years, y=LinAge, col='Sotos')) +
  scale_colour_manual(name="Disease status",values=c("Control"='grey', "Sotos"='orange')) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  xlab("Chronological age (years)") + ylab("LinAge (years)") +
  labs(title = paste0("Control: N=", sum(all_metadata$Disease_status=='Control'), "\n",
                      "Sotos: N=", sum(all_metadata$Disease_status=='Sotos')),
       subtitle=paste0()) +
  xlim(c(-1,55)) + ylim(c(-60,71)) + geom_abline(slope=1, intercept=0, linetype=2);
ggsave("plots/plot_Sotos_LinAge.pdf", height=5, width=5.5);

## Calculate EAA according to Lin clock. ##

lm_formula_ext_int_lin <- paste0('LinAge~Age_years+Sex+', 
                                paste(colnames(all_metadata)[grep('PC',colnames(all_metadata))], collapse='+'));
lm_formula_int_lin <- paste0('LinAge~Age_years+Sex+Gran+CD4T+CD8T+B+Mono+NK+', 
                            paste(colnames(all_metadata)[grep('PC',colnames(all_metadata))], collapse='+'));

lm_control_ext_int_lin <- lm(lm_formula_ext_int_lin, data=all_metadata[which(all_metadata$Disease_status=='Control'),]);
lm_control_int_lin <- lm(lm_formula_int_lin, data=all_metadata[which(all_metadata$Disease_status=='Control'),]);
all_metadata$lin_delta_ext_int[all_metadata$Disease_status=='Control'] <- lm_control_ext_int_lin$residuals;
all_metadata$lin_delta_int[all_metadata$Disease_status=='Control'] <- lm_control_int_lin$residuals;

all_metadata$lin_delta_ext_int[all_metadata$Disease_status=='Sotos'] <-
  all_metadata$LinAge[all_metadata$Disease_status=='Sotos'] - as.numeric(predict(lm_control_ext_int_lin, newdata=all_metadata[all_metadata$Disease_status=='Sotos',]));
all_metadata$lin_delta_int[all_metadata$Disease_status=='Sotos'] <-
  all_metadata$LinAge[all_metadata$Disease_status=='Sotos'] - as.numeric(predict(lm_control_int_lin, newdata=all_metadata[all_metadata$Disease_status=='Sotos',]));

## EAA boxplots. ##

plot_delta_ext_int_lin <- ggboxplot(all_metadata, x = "Disease_status", y = "lin_delta_ext_int", color = "Disease_status", palette = my_palette, 
                                   fill="Disease_status", alpha=0.5) +
  theme_classic() +
  theme(axis.text=element_text(size=12, angle=90),
        axis.title=element_text(size=14,face="bold")) +
  xlab("") + ylab("Lin EAA without CCC (years)") +
  ylim(c(-70,70)) +
  geom_abline(slope=0, intercept=0, linetype=2) +
  stat_compare_means(comparisons = mc_1) +
  guides(colour=FALSE, fill=FALSE);
ggsave("plots/plot_Sotos_LinAge_EAA_without_CCC.pdf", height=5, width=4.5);

wilcox.test(x=all_metadata$lin_delta_ext_int[all_metadata$Disease_status=='Control'], 
            y=all_metadata$lin_delta_ext_int[all_metadata$Disease_status=='Sotos'])$p.value;

plot_delta_int_lin <- ggboxplot(all_metadata, x = "Disease_status", y = "lin_delta_int", color = "Disease_status", palette = my_palette, 
                               fill="Disease_status", alpha=0.5) +
  theme_classic() +
  theme(axis.text=element_text(size=12, angle=90),
        axis.title=element_text(size=14,face="bold")) +
  xlab("") + ylab("Lin EAA with CCC (years)") +
  ylim(c(-70,70)) +
  geom_abline(slope=0, intercept=0, linetype=2) +
  stat_compare_means(comparisons = mc_1) +
  guides(colour=FALSE, fill=FALSE);
ggsave("plots/plot_Sotos_LinAge_EAA_with_CCC.pdf", height=5, width=4.5);


wilcox.test(x=all_metadata$lin_delta_int[all_metadata$Disease_status=='Control'], 
            y=all_metadata$lin_delta_int[all_metadata$Disease_status=='Sotos'])$p.value;

################################################################
################## End of the script ###########################
################################################################
