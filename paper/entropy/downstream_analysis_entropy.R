###########################################################################################
#########                                                                         #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             10/12/2018                                  #########
#########                              EMBL-EBI                                   #########
#########                           Thornton group                                #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####              Biological insights into the epigenetic ageing clock           #########
###########################################################################################
##### Analyse Shannon entropy results for the controls and Sotos syndrome.            #####
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

setwd('/Users/dem44/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/entropy');


################################################################
################## Functions ###################################
################################################################

#### Function: calculate Shannon entropy given the beta-values for a sample.
# The formula used is the same as the one reported by Hannum et al. 2013
#
# betas: numeric character containing the beta-values for a sample. 

calculate_entropy <- function(betas){
  betas[betas==0] <- 0.000000000000001;
  betas[betas==1] <- 0.999999999999999;
  return(sum(betas*log2(betas) + (1-betas)*log2(1-betas))/(length(betas)*log2(1/2)));
}


############################################################
################## Running the pipeline ####################
############################################################

#### 1. Read the raw data. ####

entropy_raw <- as.data.frame(fread('entropy_results.csv'));
controls_metadata <- as.data.frame(fread('../syndromes_screen/controls_data_downstream.tsv'));
cases_metadata <- as.data.frame(fread('../syndromes_screen/cases_data_downstream.tsv'));
horvath_sites_control_path <- '/Users/dem44/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/syndromes_screen/horvath_clock_sites/blood_control/';
horvath_sites_cases_path <- '/Users/dem44/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/syndromes_screen/horvath_clock_sites/cases/';


#### 2. Results for genome-wide entropy. ####

## Merge the entropy results with the metadata.

controls_merge <- merge(controls_metadata[,colnames(controls_metadata) %in% colnames(cases_metadata)], entropy_raw, by='GEO_sample');
sotos_merge <- merge(cases_metadata[cases_metadata$Disease_status=='Sotos',], entropy_raw, by='GEO_sample');
all_merge <- rbind(controls_merge, sotos_merge[,colnames(sotos_merge) %in% colnames(controls_merge)]);
cor(controls_merge$Age_years, controls_merge$entropy, method="spearman");
cor.test(controls_merge$Age_years, controls_merge$entropy, method="spearman")$p.value;
mean(controls_merge$entropy);

# Correlation result when removing 'GSE41273' batch --> same conclusion.

cor(controls_merge$Age_years[controls_merge$Batch!='GSE41273'], controls_merge$entropy[controls_merge$Batch!='GSE41273'], method="spearman");
cor.test(controls_merge$Age_years[controls_merge$Batch!='GSE41273'], controls_merge$entropy[controls_merge$Batch!='GSE41273'], method="spearman")$p.value;

# Correlation result when removing 'GSE59065' batch --> same conclusion.

cor(controls_merge$Age_years[controls_merge$Batch!='GSE59065'], controls_merge$entropy[controls_merge$Batch!='GSE59065'], method="spearman");
cor.test(controls_merge$Age_years[controls_merge$Batch!='GSE59065'], controls_merge$entropy[controls_merge$Batch!='GSE59065'], method="spearman")$p.value;


## Does genome-wide entropy from Sotos patients differ from the healthy controls? --> No

lm_formula_int <- paste0('entropy~Age_years+Sex+Gran+CD4T+CD8T+B+Mono+NK+', paste(colnames(all_merge)[grep('PC',colnames(all_merge))], collapse='+'));
lm_entropy_control <- lm(lm_formula_int, data=all_merge[all_merge$Disease_status=='Control',]);
all_merge$entropy_diff <- NA;
all_merge$entropy_diff[all_merge$Disease_status=='Control'] <- lm_entropy_control$residuals;
all_merge$entropy_diff[all_merge$Disease_status=='Sotos'] <- all_merge$entropy[all_merge$Disease_status=='Sotos'] - 
  as.numeric(predict(lm_entropy_control, newdata=all_merge[all_merge$Disease_status=='Sotos',]));
wilcox.test(all_merge$entropy_diff[all_merge$Disease_status=='Control'], all_merge$entropy_diff[all_merge$Disease_status=='Sotos'])$p.value;

## Plots. 

my_palette <- c('grey', 'orange');

entropy_Sotos_scatterplot <- ggplot(data=all_merge, aes(x=Age_years, y=entropy, col=Disease_status)) + 
  geom_point() + scale_color_manual(values=my_palette) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  xlab("Chronological age (years)") + ylab("Genome-wide Shannon entropy") +
  labs(title = paste0("Control: N=", sum(all_merge$Disease_status=='Control'), "\n",
                      "Sotos: N=", sum(all_merge$Disease_status=='Sotos')), 
       subtitle=paste0()) + 
  guides(col=guide_legend(title="Disease status")) + 
  xlim(c(-1,55)) + ylim(c(0.30, 0.60));
ggsave("plots/entropy_Sotos_scatterplot.pdf", height=5, width=5);

entropy_batch_scatterplot <- ggplot(data=all_merge, aes(x=Age_years, y=entropy, col=Batch)) + 
  geom_point() + scale_colour_stata() +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  xlab("Chronological age (years)") + ylab("Genome-wide Shannon entropy") +
  labs(title = paste0("Control: N=", sum(all_merge$Disease_status=='Control'), "\n",
                      "Sotos: N=", sum(all_merge$Disease_status=='Sotos')), 
       subtitle=paste0()) + 
  guides(col=guide_legend(title="Batch")) + 
  xlim(c(-1,55)) + ylim(c(0.30, 0.60));
ggsave("plots/entropy_batch_scatterplot.pdf", height=5, width=5);

mc_1 <- list(c("Control", "Sotos"));

plots_entropy_diff <- ggboxplot(all_merge, x = "Disease_status", y = "entropy_diff", color = "Disease_status", palette = my_palette, fill="Disease_status", alpha=0.5) +
  theme_classic() +
  theme(axis.text=element_text(size=12, angle=90),
        axis.title=element_text(size=13,face="bold")) +
  xlab("") + ylab("Genome-wide Shannon entropy acceleration") +
  geom_abline(slope=0, intercept=0, linetype=2) +
  stat_compare_means(comparisons = mc_1) +
  guides(colour=FALSE, fill=FALSE);
ggsave("plots/entropy_diff_comparison.pdf", height=5, width=5);

plots_scatterplot_entropy_diff_Sotos <- ggplot(data=all_merge[all_merge$Disease_status=='Sotos',], aes(x=Age_years, y=entropy_diff)) + 
  geom_point(data=all_merge[all_merge$Disease_status=='Control',], aes(x=Age_years, y=entropy_diff), col='grey') + 
  geom_abline(slope=0, intercept=0, linetype=2) +
  geom_smooth(method="lm", formula=y~x, show.legend=F, fill='khaki', alpha = 0.3, col='gold') + geom_point(col='orange') + 
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  xlab("Chronological age (years)") + ylab("Genome-wide Shannon entropy acceleration") +
  labs(title = paste0("Control: N=", sum(all_merge$Disease_status=='Control'), "\n",
                      "Sotos: N=", sum(all_merge$Disease_status=='Sotos')), 
       subtitle=paste0()) +  
  xlim(c(-1,55)) + ylim(c(-0.06, 0.06));
ggsave("plots/entropy_diff_Sotos_scatterplot.pdf", height=5, width=5);

# plots_scatterplot_entropy_diff_batch <- ggplot(data=all_merge, aes(x=Age_years, y=entropy_diff, col=Batch)) +
#   geom_point() + scale_colour_stata() +
#   geom_abline(slope=0, intercept=0, linetype=2) +
#   theme_classic() +
#   theme(axis.text=element_text(size=12),
#         axis.title=element_text(size=14,face="bold"),
#         plot.title = element_text(hjust = 0.5),
#         plot.subtitle = element_text(hjust = 0.5)) +
#   xlab("Chronological age (years)") + ylab("Genome-wide Shannon entropy acceleration") +
#   labs(title = paste0("Control: N=", sum(all_merge$Disease_status=='Control'), "\n",
#                       "Sotos: N=", sum(all_merge$Disease_status=='Sotos')),
#        subtitle=paste0()) +
#   guides(col=guide_legend(title="Batch")) +
#   xlim(c(-1,55)) + ylim(c(-0.06, 0.06));
# ggsave("plots/entropy_diff_batch_scatterplot.pdf", height=5, width=5);


#### 3. Results for entropy in the 353 clock CpG sites. ####

## Read the control data. 

control_files <- list.files(path=horvath_sites_control_path, full.names = TRUE);

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

control_betas <- as.data.frame(t(control_betas[,which(colnames(control_betas) %in% controls_metadata$GEO_sample)])); # Keep only the betas for those samples that were selected after QC.
control_betas$GEO_sample <- rownames(control_betas); rownames(control_betas) <- NULL;
controls_merge <- merge(controls_merge, control_betas, by='GEO_sample');
rm(control_betas, temp_betas, control_files);

## Read the cases data. 

cases_files <- list.files(path=horvath_sites_cases_path, full.names = TRUE);

for(i in 1:length(cases_files)){
  
  temp_betas <- as.data.frame(fread(cases_files[i]));
  rownames(temp_betas) <- temp_betas[,1];
  temp_betas <- temp_betas[,-1];
  
  if(i==1){
    cases_betas <- temp_betas;
  }else{
    cases_betas <- cbind(cases_betas, temp_betas);
  }
}

sotos_betas <- as.data.frame(t(cases_betas[,which(colnames(cases_betas) %in% sotos_merge$GEO_sample)])); # Keep only the betas for those samples that were selected after QC.
sotos_betas$GEO_sample <- rownames(sotos_betas); rownames(sotos_betas) <- NULL;
sotos_merge <- merge(sotos_merge, sotos_betas, by='GEO_sample');
rm(cases_betas, temp_betas, cases_files);

## Merge all the data.

all_merge <- rbind(controls_merge, sotos_merge[,colnames(sotos_merge) %in% colnames(controls_merge)]);

## Calculate clock sites entropy. 

all_merge$clock_entropy <- as.numeric(apply(all_merge[,grep('cg', colnames(all_merge))], 1, calculate_entropy));
cor(all_merge$Age_years[all_merge$Disease_status=='Control'], all_merge$clock_entropy[all_merge$Disease_status=='Control'], method="spearman");
cor.test(all_merge$Age_years[all_merge$Disease_status=='Control'], all_merge$clock_entropy[all_merge$Disease_status=='Control'], method="spearman")$p.value;
mean(all_merge$clock_entropy[all_merge$Disease_status=='Control']);

# Correlation result when removing 'Europe' batch --> entropy increases with age this time (different conclusion).

cor(all_merge$Age_years[all_merge$Disease_status=='Control' & all_merge$Batch != 'Europe'], 
    all_merge$clock_entropy[all_merge$Disease_status=='Control' & all_merge$Batch != 'Europe'], method="spearman");
cor.test(all_merge$Age_years[all_merge$Disease_status=='Control' & all_merge$Batch != 'Europe'], 
    all_merge$clock_entropy[all_merge$Disease_status=='Control' & all_merge$Batch != 'Europe'], method="spearman")$p.value;


## Does the clock sites entropy from Sotos patients differ from the healthy controls? --> Yes

lm_formula_int <- paste0('clock_entropy~Age_years+Sex+Gran+CD4T+CD8T+B+Mono+NK+', paste(colnames(all_merge)[grep('PC',colnames(all_merge))], collapse='+'));
lm_entropy_control <- lm(lm_formula_int, data=all_merge[all_merge$Disease_status=='Control',]);
all_merge$clock_entropy_diff <- NA;
all_merge$clock_entropy_diff[all_merge$Disease_status=='Control'] <- lm_entropy_control$residuals;
all_merge$clock_entropy_diff[all_merge$Disease_status=='Sotos'] <- all_merge$clock_entropy[all_merge$Disease_status=='Sotos'] - 
  as.numeric(predict(lm_entropy_control, newdata=all_merge[all_merge$Disease_status=='Sotos',]));
wilcox.test(all_merge$clock_entropy_diff[all_merge$Disease_status=='Control'], all_merge$clock_entropy_diff[all_merge$Disease_status=='Sotos'])$p.value;


## Plots. 

clock_entropy_Sotos_scatterplot <- ggplot(data=all_merge, aes(x=Age_years, y=clock_entropy, col=Disease_status)) + 
  geom_point() + scale_color_manual(values=my_palette) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  xlab("Chronological age (years)") + ylab("Shannon entropy for the clock sites") +
  labs(title = paste0("Control: N=", sum(all_merge$Disease_status=='Control'), "\n",
                      "Sotos: N=", sum(all_merge$Disease_status=='Sotos')), 
       subtitle=paste0()) + 
  guides(col=guide_legend(title="Disease status")) + 
  xlim(c(-1,55)) + ylim(c(0.30, 0.60));
ggsave("plots/clock_entropy_Sotos_scatterplot.pdf", height=5, width=5);

clock_entropy_batch_scatterplot <- ggplot(data=all_merge, aes(x=Age_years, y=clock_entropy, col=Batch)) + 
  geom_point() + scale_colour_stata() +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  xlab("Chronological age (years)") + ylab("Shannon entropy for the clock sites") +
  labs(title = paste0("Control: N=", sum(all_merge$Disease_status=='Control'), "\n",
                      "Sotos: N=", sum(all_merge$Disease_status=='Sotos')), 
       subtitle=paste0()) + 
  guides(col=guide_legend(title="Batch")) + 
  xlim(c(-1,55)) + ylim(c(0.30, 0.60));
ggsave("plots/clock_entropy_batch_scatterplot.pdf", height=5, width=5);

plots_clock_entropy_diff <- ggboxplot(all_merge, x = "Disease_status", y = "clock_entropy_diff", color = "Disease_status", palette = my_palette, fill="Disease_status", alpha=0.5) +
  theme_classic() +
  theme(axis.text=element_text(size=12, angle=90),
        axis.title=element_text(size=13,face="bold")) +
  xlab("") + ylab("Shannon entropy acceleration for the clock sites") +
  geom_abline(slope=0, intercept=0, linetype=2) +
  stat_compare_means(comparisons = mc_1) +
  guides(colour=FALSE, fill=FALSE);
ggsave("plots/clock_entropy_diff_comparison.pdf", height=5, width=5);

plots_scatterplot_clock_entropy_diff_Sotos <- ggplot(data=all_merge[all_merge$Disease_status=='Sotos',], aes(x=Age_years, y=clock_entropy_diff)) + 
  geom_point(data=all_merge[all_merge$Disease_status=='Control',], aes(x=Age_years, y=clock_entropy_diff), col='grey') + 
  geom_abline(slope=0, intercept=0, linetype=2) +
  geom_smooth(method="lm", formula=y~x, show.legend=F, fill='khaki', alpha = 0.3, col='gold') + geom_point(col='orange') + 
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  xlab("Chronological age (years)") + ylab("Shannon entropy acceleration for the clock sites") +
  labs(title = paste0("Control: N=", sum(all_merge$Disease_status=='Control'), "\n",
                      "Sotos: N=", sum(all_merge$Disease_status=='Sotos')), 
       subtitle=paste0()) +  
  xlim(c(-1,55)) + ylim(c(-0.06, 0.06));
ggsave("plots/clock_entropy_diff_Sotos_scatterplot.pdf", height=5, width=5);

# plots_scatterplot_clock_entropy_diff_batch <- ggplot(data=all_merge, aes(x=Age_years, y=clock_entropy_diff, col=Batch)) +
#   geom_point() + scale_colour_stata() +
#   geom_abline(slope=0, intercept=0, linetype=2) +
#   theme_classic() +
#   theme(axis.text=element_text(size=12),
#         axis.title=element_text(size=14,face="bold"),
#         plot.title = element_text(hjust = 0.5),
#         plot.subtitle = element_text(hjust = 0.5)) +
#   xlab("Chronological age (years)") + ylab("Shannon entropy acceleration for the clock sites") +
#   labs(title = paste0("Control: N=", sum(all_merge$Disease_status=='Control'), "\n",
#                       "Sotos: N=", sum(all_merge$Disease_status=='Sotos')),
#        subtitle=paste0()) +
#   guides(col=guide_legend(title="Batch")) +
#   xlim(c(-1,55)) + ylim(c(-0.06, 0.06));
# ggsave("plots/clock_entropy_diff_batch_scatterplot.pdf", height=5, width=5);

################################################################
################## End of the script ###########################
################################################################
