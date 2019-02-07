###########################################################################################
#########                                                                         #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             12/12/2018                                  #########
#########                              EMBL-EBI                                   #########
#########                           Thornton group                                #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####              Biological insights into the epigenetic ageing clock           #########
###########################################################################################
##### Look at the beta-values of the 353 clock CpG sites for Sotos and healthy samples. ###
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
library(gridExtra);
library(circlize);
library(viridis);
library(colorRamps);
library(ComplexHeatmap);

setwd('/Users/dem44/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/differential_methylation/');


################################################################
################## Run the pipeline ############################
################################################################

##### 1. Read all the data. ####

sotos_dmps_path <- '/Users/dem44/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/differential_methylation/Sotos_DMPs_final.csv';
admps_path <- '/Users/dem44/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/differential_methylation/aDMPs_final.csv';

## Read the control data. 

control_metadata <- as.data.frame(fread('../syndromes_screen/controls_data_downstream.tsv')); # Small control

control_files <- list.files(path='../syndromes_screen/horvath_clock_sites/blood_control/', full.names = TRUE);

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
rm(control_betas, control_metadata, temp_betas, control_files);

## Read the cases data. 

cases_metadata <- as.data.frame(fread('../syndromes_screen/cases_data_downstream.tsv'));

# Read the beta-values. 

cases_files <- list.files(path='../syndromes_screen/horvath_clock_sites/cases/', full.names = TRUE);

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

cases_betas <- as.data.frame(t(cases_betas[,which(colnames(cases_betas) %in% cases_metadata$GEO_sample)])); # Keep only the betas for those samples that were selected after QC.
cases_betas$GEO_sample <- rownames(cases_betas); rownames(cases_betas) <- NULL;
cases_final <- merge(cases_metadata, cases_betas, by='GEO_sample');
rm(cases_betas, cases_metadata, temp_betas, cases_files);

## Merge all the data.

cases_final_edited <- cases_final[,match(colnames(control_final), colnames(cases_final))];
all_final <- rbind(control_final, cases_final_edited);


##### 2. Create matrices that contain deviations from the expected beta-values in the clock CpGs for the cases of interest.  
## Two matrices are created: one without cell composition correction and one with cell composition correction. 

select_dis <- c('Sotos'); 
selected_cases <- cases_final[cases_final$Disease_status %in% select_dis,];
cg_ids <- colnames(selected_cases)[grep('cg', colnames(selected_cases))];
cg_df <- data.frame(CpGmarker=cg_ids, sp_ours=rep(NA, length(cg_ids)));
selected_cases$hm_names <- paste0(selected_cases$GEO_sample, '_', selected_cases$Disease_status);
diff_betas_ei_df <- matrix(NA, ncol=nrow(selected_cases), nrow=length(cg_ids)); # Without cell composition correction
diff_betas_i_df <- matrix(NA, ncol=nrow(selected_cases), nrow=length(cg_ids)); # With cell composition correction
colnames(diff_betas_ei_df) <- colnames(diff_betas_i_df) <- selected_cases$hm_names;
rownames(diff_betas_ei_df) <- rownames(diff_betas_i_df) <- cg_ids;

for(i in 1:length(cg_ids)){
  
  cg <- cg_ids[i];
  
  # Fit linear models to the controls
  
  formula_ei <- paste0(cg, '~Age_years+I(Age_years^2)+Sex+', paste(colnames(control_final)[grep('PC',colnames(control_final))], collapse='+'));
  formula_i <- paste0(cg, '~Age_years+I(Age_years^2)+Sex+Gran+CD4T+CD8T+B+Mono+NK+', paste(colnames(control_final)[grep('PC',colnames(control_final))], collapse='+'));
  model_ei <- lm(formula_ei, data=control_final);
  model_i <- lm(formula_i, data=control_final);
  
  # Store residuals for the cases.
  
  diff_betas_ei_df[i,] <- selected_cases[,cg] - predict(model_ei, newdata=selected_cases);
  diff_betas_i_df[i,] <- selected_cases[,cg] - predict(model_i, newdata=selected_cases);
  
  # Calculate additional information for CpGs.
  
  cg_df$sp_ours[i] <- cor(control_final$Age_years, control_final[,cg], method='spearman');
  cg_df$Sotos_ei_median[i] <- median(selected_cases[,cg] - predict(model_ei, newdata=selected_cases)); 
  cg_df$Sotos_i_median[i] <- median(selected_cases[,cg] - predict(model_i, newdata=selected_cases));
  
}


#### 3. Plot the distribution of hypermethylated and hypomethylated clock sites in Sotos. ####

## Obtain final annotation for the 353 CpG probes. 

clock_sites_ann <- as.data.frame(fread('clock_CpGs_extended_annotation.csv'));
clock_sites_ann$`Marginal Age Relationship` <- ifelse(clock_sites_ann$`Marginal Age Relationship`=='positive', 'Hyper', 'Hypo');
colnames(clock_sites_ann)[colnames(clock_sites_ann)=='Marginal Age Relationship'] <- 'hyperhypo_horvath';

final_cg_ann <- merge(clock_sites_ann, cg_df, by='CpGmarker');
sotos_dmps <- as.data.frame(fread(sotos_dmps_path));
sotos_dmps <- sotos_dmps[sotos_dmps$pval<(0.01/nrow(sotos_dmps)),]; # Filter significant DMPs
colnames(sotos_dmps)[1] <- 'CpGmarker'; 
final_cg_ann <- merge(final_cg_ann, sotos_dmps[,c(1,3,5)], all.x=TRUE, by='CpGmarker');
#final_cg_ann$Sotos_diff_i_discrete <- ifelse(final_cg_ann$Sotos_i_median < -0.05, 'Hypomethylated', 
#                                             ifelse(final_cg_ann$Sotos_i_median > 0.05, 'Hypermethylated', NA));
final_cg_ann$Sotos_diff_i_discrete <- ifelse(final_cg_ann$beta < 0, 'Hypomethylated',
                                             ifelse(final_cg_ann$beta > 0, 'Hypermethylated', NA));

admps <- as.data.frame(fread(admps_path));
admps <- admps[admps$pval<(0.01/nrow(admps)),]; # Filter significant DMPs
colnames(admps)[1] <- 'CpGmarker'; 
final_cg_ann <- merge(final_cg_ann, admps[,c(1,3,5)], all.x=TRUE, by='CpGmarker');
final_cg_ann$aDMP <- ifelse(final_cg_ann$beta.y < 0, 'Hypomethylated',
                                             ifelse(final_cg_ann$beta.y > 0, 'Hypermethylated', NA));


final_cg_ann <- final_cg_ann[match(rownames(diff_betas_ei_df), final_cg_ann$CpGmarker),]; # Order
write.table(x=final_cg_ann[,-c(46,47,49,50)], file='/Users/dem44/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/differential_methylation/clock_CpGs_extended_annotation_final.csv', 
            quote=F, sep=',', row.names=F);

## Plots. 

# plot_density_hyperhypo_Sotos <- ggplot(data=final_cg_ann, aes(x=Sotos_i_median)) + 
#   geom_density(fill='orange', col='orange') +
#   theme_classic() +
#   theme(axis.text=element_text(size=12),
#         axis.title=element_text(size=14,face="bold"),
#         plot.title = element_text(hjust = 0.5),
#         plot.subtitle = element_text(hjust = 0.5)) +
#   xlab("Median beta-value difference (Sotos-control)") + ylab("Density") +
#   geom_vline(xintercept=0, linetype=2);
# ggsave('plots/distribution_diff_betas_Sotos.pdf', height=6, width=6);


counts_all <- c(sum(final_cg_ann$Sotos_i_median < 0), # Hypomethylated
                sum(final_cg_ann$Sotos_i_median >0)); # Hypermethylated
#counts_filtered <- c(sum(final_cg_ann$Sotos_i_median[abs(final_cg_ann$Sotos_i_median) > 0.05] < 0), # Hypomethylated
#                     sum(final_cg_ann$Sotos_i_median[abs(final_cg_ann$Sotos_i_median) > 0.05] > 0)); # Hypermethylated
counts_filtered <- c(sum(final_cg_ann$Sotos_diff_i_discrete=='Hypomethylated', na.rm=TRUE),
                     sum(final_cg_ann$Sotos_diff_i_discrete=='Hypermethylated', na.rm=TRUE));
#hyperhypo_Sotos_df <- data.frame(counts=c(counts_all, counts_filtered),
#                                 filtering=c('All sites', 'All sites', 'Sites with \neffect size > 5%', 'Sites with \neffect size > 5%'),
#                                 sites=c('Hypomethylated', 'Hypermethylated', 'Hypomethylated', 'Hypermethylated')); 
hyperhypo_Sotos_df <- data.frame(counts=c(counts_all, counts_filtered),
                                 filtering=c('All sites', 'All sites', 'Sotos DMPs', 'Sotos DMPs'),
                                 sites=c('Hypomethylated', 'Hypermethylated', 'Hypomethylated', 'Hypermethylated')); 
# plot_hist_hyperhypo_Sotos <- ggplot(data=hyperhypo_Sotos_df, aes(x=filtering, y=counts, fill=sites)) +
#   geom_bar(stat="identity") + scale_fill_manual(values=c('Red', 'Blue')) +
#   theme_classic() +
#   theme(axis.text=element_text(size=12),
#         axis.title=element_text(size=14,face="bold"),
#         plot.title = element_text(hjust = 0.5),
#         plot.subtitle = element_text(hjust = 0.5)) +
#   xlab("") + ylab("Number of clock CpG sites") + 
#   guides(fill=guide_legend(title="Methylation change in Sotos"));
# ggsave('plots/hist_diff_betas_Sotos.pdf', height=6, width=6);


##### 4. Plot the beta-values for two Sotos DMPs.

#cgs <- c(final_cg_ann$CpGmarker[final_cg_ann$Sotos_i_median < -0.3 & final_cg_ann$sp_ours < -0.4],
#         final_cg_ann$CpGmarker[final_cg_ann$Sotos_i_median == max(final_cg_ann$Sotos_i_median)]);
cgs <- final_cg_ann$CpGmarker[order(final_cg_ann$pval.x)[1:2]];

cg_plots <- list();

for(cg in cgs){
  temp_df <- all_final[,which(colnames(all_final) %in% c('Age_years', cg, 'Sex', 'Disease_status'))];
  temp_df <- temp_df[temp_df$Disease_status %in% c('Control', 'Sotos'),];
  colnames(temp_df)[which(colnames(temp_df) == cg)] <- 'Beta';
  cg_plots[[cg]] <- ggplot(data=temp_df[temp_df$Disease_status=='Control',], aes(x=Age_years, y=Beta, col=Disease_status)) + geom_point() + 
    theme_classic() + scale_colour_manual(values=c('grey', 'orange')) +
    theme(axis.text=element_text(size=12, angle=90),
          axis.title=element_text(size=14,face="bold"),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) +
    labs(title=cg, subtitle=paste0(""),
         colour="Disease status") + 
    xlab('Chronological age (years)') + ylab('Beta-value') +
    geom_smooth(method='lm',formula=y~x+I(x^2), show.legend = F, col='gray40') + 
    geom_point(data=temp_df[temp_df$Disease_status=='Sotos',], aes(x=Age_years, y=Beta, col=Disease_status)) +
    ylim(c(0, 1));
}

glist <- lapply(cg_plots, ggplotGrob)
ggsave("plots/cg_beta_plots.pdf", marrangeGrob(glist, nrow=1, ncol=2, 
                                               layout_matrix=matrix(1:2, nrow=1, ncol=2, byrow=T), top=NULL),
       height=5, width=10);


#### 5. Create heatmaps to represent the matrices of the beta-value differences for the clock CpGs. ####

## Heatmap annotation. 

col_sex <- c('Male'='Pink', 'Female'='Deepskyblue2');
col_delta_ext_int <- colorRamp2(c(min(selected_cases$delta_ext_int), max(selected_cases$delta_ext_int)), inferno(2)[c(2,1)]);
col_delta_int <- colorRamp2(c(min(selected_cases$delta_int), max(selected_cases$delta_int)), inferno(2)[c(2,1)]);
col_age <- colorRamp2(c(0,50), c('White', 'Black'));

col_weights <- colorRamp2(c(-1, 0, 1), magenta2green(3)[c(3,2,1)]); #c("#d8bfd8", "White", "#ffe200"));
#col_horvath <- c("Hyper"="#cd9d0b", "Hypo"="#614051");
#col_sp <- colorRamp2(c(-1, 0, 1), c('#614051', 'White', '#cd9d0b'));
col_aDMP <- c('Hypermethylated'='#cd9d0b', "Hypomethylated"='#614051');
col_diff <- c('Hypermethylated'='Red', "Hypomethylated"='Blue');
col_chrhmm <- c('Active TSS'=rgb(255,0,0, maxColorValue = 255), 'Promoter'=rgb(255,69,0, maxColorValue = 255), 'Transcribed'=rgb(0,128,0, maxColorValue = 255), 'Weakly transcribed'=rgb(0,150,0, maxColorValue = 255),
                'Transcribed/regulatory'=rgb(194,225,5, maxColorValue = 255), 'Active enhancer'=rgb(255,195,77, maxColorValue = 255), 'Weak enhancer'=rgb(255,255,0, maxColorValue = 255),
                'DNase'=rgb(255,255,102, maxColorValue = 255), 'ZNF/repeats'=rgb(102,205,170, maxColorValue = 255), 'Heterochromatin'=rgb(138,145,208, maxColorValue = 255), 'Poised promoter'=rgb(230,184,183, maxColorValue = 255),
                'Bivalent promoter'=rgb(112,48,160, maxColorValue = 255), 'Repressed polycomb'=rgb(128,128,128, maxColorValue = 255), 'Quiescent/low'="lightgrey");
col_rna <- colorRamp2(c(-2, 0, 2), c('#8A91D0', 'White', '#008000'));
col_H3K36me3 <- colorRamp2(c(-2,0,2), c('deeppink3', 'White', 'Blue4'))
col_gb <- c('Yes'='Blue4', 'No'='gray');

main_col_scale <- colorRamp2(c(-0.2,-0.1, 0, 0.1, 0.2),viridis(5));

# ha_samples_ei <- rowAnnotation(df=data.frame(Sex=selected_cases$Sex, Delta_ext_int=selected_cases$delta_ext_int, Age_years=selected_cases$Age_years),
#                                col=list(Sex=col_sex, Delta_ext_int=col_delta_ext_int, Age_years=col_age),
#                                annotation_legend_param=list(Delta_ext_int=list(title='EAA without ccc (years)'),
#                                                             Age_years=list(title='Age (years)')));    # Using info from selected_cases

ha_samples_i <- rowAnnotation(df=data.frame(Sex=selected_cases$Sex, Delta_int=selected_cases$delta_int, Age_years=selected_cases$Age_years),
                              col=list(Sex=col_sex, Delta_int=col_delta_int, Age_years=col_age),
                              annotation_legend_param=list(Delta_int=list(title='EAA with CCC (years)'),
                                                           Age_years=list(title='Chronological age (years)')));

ha_cgs <- HeatmapAnnotation(df=data.frame(Sotos_diff=final_cg_ann$Sotos_diff_i_discrete,
                                          #Spearman=final_cg_ann$sp_ours,
                                          aDMP=final_cg_ann$aDMP,
                                          Weight=final_cg_ann$CoefficientTraining,
                                          ChrHMM=final_cg_ann$ChrHMM_state_reann,
                                          RNA=final_cg_ann$RNA,
                                          H3K36me3=final_cg_ann$H3K36me3,
                                          Gene_body=final_cg_ann$Gene_body),
                            col=list(Sotos_diff=col_diff,
                                     #Spearman=col_sp,
                                     aDMP=col_aDMP,
                                     Weight=col_weights,
                                     ChrHMM=col_chrhmm,
                                     RNA=col_rna,
                                     H3K36me3=col_H3K36me3,
                                     Gene_body=col_gb),
                            annotation_legend_param=list(Sotos_diff=list(title="Sotos DMPs"),
                                                         #Spearman=list(title='Age correlation \n(in controls)',
                                                         #             labels=c('-1 (hypomethylation)', '-0.5', '0', '0.5', '1 (hypermethylation)')),
                                                         aDMP=list(title="aDMPs"),
                                                         Weight=list(title="Weight \nin model"),
                                                         ChrHMM=list(title="ChrHMM state \n(in K562)"),
                                                         RNA=list(title="RNA \n(in PBMC)"),
                                                         H3K36me3=list(title="H3K36me3 \n(in PBMC)"),
                                                         Gene_body=list(title="In gene body"))); # Using info from final_cg_ann

## Plot without cell composition correction.

# final_hm_ei <- Heatmap(t(diff_betas_ei_df), col=main_col_scale,
#                        show_row_names = FALSE, column_names_gp = gpar(fontsize = 2), 
#                        bottom_annotation=ha_cgs, heatmap_legend_param = list(title = "Beta-value difference"),
#                        column_title="Horvath's clock CpGs", column_title_gp = gpar(fontsize = 15, fontface = "bold"), column_title_side = "top",
#                        row_title="Sotos samples", row_title_gp = gpar(fontsize = 15, fontface = "bold")) + ha_samples_ei;
# pdf('plots/heatmap_betadiffs_ei.pdf', width=10, height=8);
# draw(final_hm_ei, annotation_legend_side = "bottom");
# dev.off();

## Plot with cell composition correction.

final_hm_i <- Heatmap(t(diff_betas_i_df), col=main_col_scale,
                      show_row_names = FALSE, column_names_gp = gpar(fontsize = 2), 
                      bottom_annotation=ha_cgs, heatmap_legend_param = list(title = "Beta-value difference"),
                      column_title="Horvath's clock CpGs", column_title_gp = gpar(fontsize = 15, fontface = "bold"), column_title_side = "top",
                      row_title="Sotos samples", row_title_gp = gpar(fontsize = 15, fontface = "bold")) + ha_samples_i;
pdf('plots/heatmap_betadiffs_i.pdf', width=10, height=8);
draw(final_hm_i, annotation_legend_side = "bottom");
dev.off();


################################################################
################## End of the script ###########################
################################################################