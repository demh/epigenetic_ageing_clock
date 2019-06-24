###########################################################################################
#########                                                                         #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             09/12/2018                                  #########
#########                              EMBL-EBI                                   #########
#########                           Thornton group                                #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####              Biological insights into the epigenetic ageing clock           #########
###########################################################################################
##### Downstream analysis of the results for the epigenetic mitotic clock (pcgtAge)   #####
##### for Sotos compared with healthy samples.                                        #####
###########################################################################################
##### USAGE: manual                                                                   #####
###########################################################################################

###########################################################
##################### Dependencies ########################
###########################################################

library(data.table);
library(ggplot2);
library(ggpubr);

setwd('/Users/dem44/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/pcgtAge');


################################################################
################## Run the pipeline ############################
################################################################

#### 1. Read the raw data. ####

pcgtAge_raw <- as.data.frame(fread('pcgtAge_results.csv'));
controls_metadata <- as.data.frame(fread('../syndromes_screen/controls_data_downstream.tsv'));
cases_metadata <- as.data.frame(fread('../syndromes_screen/final_cases_data.tsv'));


#### 2. Merge the pcgtAge results with the metadata. ####

controls_merge <- merge(controls_metadata[,colnames(controls_metadata) %in% colnames(cases_metadata)], pcgtAge_raw, by='GEO_sample');
sotos_merge <- merge(cases_metadata[cases_metadata$Disease_status=='Sotos',], pcgtAge_raw, by='GEO_sample');
all_merge <- rbind(controls_merge, sotos_merge[,colnames(sotos_merge) %in% colnames(controls_merge)]);


#### 3. Obtain the results. ####

## Does pcgtAge from Sotos patients differ from the healthy controls? --> Yes?

lm_formula_int <- paste0('pcgtAge~Age_years+Sex+Gran+CD4T+CD8T+B+Mono+NK+', paste(colnames(all_merge)[grep('PC',colnames(all_merge))], collapse='+'));
lm_pcgtAge_control <- lm(lm_formula_int, data=all_merge[all_merge$Disease_status=='Control',]);
all_merge$pcgtAge_diff <- NA;
all_merge$pcgtAge_diff[all_merge$Disease_status=='Control'] <- lm_pcgtAge_control$residuals;
all_merge$pcgtAge_diff[all_merge$Disease_status=='Sotos'] <- all_merge$pcgtAge[all_merge$Disease_status=='Sotos'] - 
  as.numeric(predict(lm_pcgtAge_control, newdata=all_merge[all_merge$Disease_status=='Sotos',]));
wilcox.test(all_merge$pcgtAge_diff[all_merge$Disease_status=='Control'], all_merge$pcgtAge_diff[all_merge$Disease_status=='Sotos'])$p.value; # All data, p-value=0.01122896
wilcox.test(all_merge$pcgtAge_diff[all_merge$Disease_status=='Control'], all_merge$pcgtAge_diff[all_merge$Disease_status=='Sotos'&all_merge$Age_years!=41])$p.value; # Without Sotos outlier, p-value=0.02454965

## Plots. 

my_palette <- c('grey', 'orange');

pcgtAge_scatterplot <- ggplot(data=all_merge, aes(x=Age_years, y=pcgtAge, col=Disease_status)) + 
  geom_point() + scale_color_manual(values=my_palette) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  xlab("Chronological age (years)") + ylab("pcgtAge") +
  labs(title = paste0("Control: N=", sum(all_merge$Disease_status=='Control'), "\n",
                      "Sotos: N=", sum(all_merge$Disease_status=='Sotos')), 
       subtitle=paste0()) + 
  guides(col=guide_legend(title="Disease status")) + 
  xlim(c(-1,55));
ggsave("plots/pcgtAge_Sotos_scatterplot.pdf", height=4, width=4.5);

pcgtAge_scatterplot_no_legend <- ggplot(data=all_merge, aes(x=Age_years, y=pcgtAge, col=Disease_status)) + 
  geom_point() + scale_color_manual(values=my_palette) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  xlab("Chronological age (years)") + ylab("pcgtAge") +
  labs(title = paste0("Control: N=", sum(all_merge$Disease_status=='Control'), "\n",
                      "Sotos: N=", sum(all_merge$Disease_status=='Sotos')), 
       subtitle=paste0()) + 
  guides(col=FALSE) + 
  xlim(c(-1,55));
ggsave("plots/pcgtAge_Sotos_scatterplot_no_legend.pdf", height=4, width=4);

mc_1 <- list(c("Control", "Sotos"));

plots_pcgtAge_diff <- ggboxplot(all_merge, x = "Disease_status", y = "pcgtAge_diff", color = "Disease_status", palette = my_palette, fill="Disease_status", alpha=0.5) +
  theme_classic() +
  theme(axis.text=element_text(size=12, angle=90),
        axis.title=element_text(size=14,face="bold")) +
  xlab("") + ylab("pcgtAge acceleration") +
  geom_abline(slope=0, intercept=0, linetype=2) +
  stat_compare_means(comparisons = mc_1) +
  guides(colour=FALSE, fill=FALSE);
ggsave("plots/pcgtAge_diff_Sotos_comparison.pdf", height=4, width=4);

plots_scatterplot_pcgtAge_diff <- ggplot(data=all_merge[all_merge$Disease_status=='Sotos',], aes(x=Age_years, y=pcgtAge_diff)) + 
  geom_point(data=all_merge[all_merge$Disease_status=='Control',], aes(x=Age_years, y=pcgtAge_diff), col='grey') + 
  geom_abline(slope=0, intercept=0, linetype=2) +
  geom_smooth(method="lm", formula=y~x, show.legend=F, fill='khaki', alpha = 0.3, col='gold') + geom_point(col='orange') + 
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  xlab("Chronological age (years)") + ylab("pcgtAge acceleration") +
  labs(title = paste0("Control: N=", sum(all_merge$Disease_status=='Control'), "\n",
                      "Sotos: N=", sum(all_merge$Disease_status=='Sotos')), 
       subtitle=paste0()) +  
  xlim(c(-1,55));
ggsave("plots/pcgtAge_diff_Sotos_scatterplot.pdf", height=4, width=4);

################################################################
################## End of the script ###########################
################################################################