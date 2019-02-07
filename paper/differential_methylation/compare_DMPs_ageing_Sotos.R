###########################################################################################
#########                                                                         #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             05/12/2018                                  #########
#########                              EMBL-EBI                                   #########
#########                           Thornton group                                #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####              Biological insights into the epigenetic ageing clock           #########
###########################################################################################
##### Compare the differentially methylated positions (DMPs) found during ageing and  #####
##### in Sotos syndrome. This includes finding genomic features that are characteristic ###
##### of the different subsets.                                                       #####
###########################################################################################
##### USAGE: manual                                                                   #####
###########################################################################################

###########################################################
##################### Dependencies ########################
###########################################################

library(data.table);
library(ggplot2);
library(ggthemes);
library(limma);
library(VennDiagram);
library(gridExtra);
library(ggpubr);
library(viridis);

set.seed(1);
setwd('/Users/dem44/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/differential_methylation');


################################################################
################## Run the pipeline ############################
################################################################

#### 1. Read the DMPs. ####

age_dmps <- as.data.frame(fread('aDMPs_final.csv'));
thr <- 0.01 / nrow(age_dmps);
age_dmps_filtered <- age_dmps[age_dmps$pval < thr,];
sotos_dmps <- as.data.frame(fread('Sotos_DMPs_final.csv'));
sotos_dmps_filtered <- sotos_dmps[sotos_dmps$pval < thr,];


#### 2. Overview of the overlap between aDMPs and Sotos DMPs. ####

## Obtain some useful variables.

all_cpgs <- age_dmps$ProbeID;
ageing_hypo <- nrow(age_dmps_filtered[age_dmps_filtered$beta<0,]);
ageing_hyper <- nrow(age_dmps_filtered[age_dmps_filtered$beta>0,]);
sotos_hypo <- nrow(sotos_dmps_filtered[sotos_dmps_filtered$beta<0,]);
sotos_hyper <- nrow(sotos_dmps_filtered[sotos_dmps_filtered$beta>0,]);
summary_DMPs <- data.frame(Condition=rep(c('Ageing', 'Sotos'), each=2), 
                           Direction=rep(c('Hypomethylated', 'Hypermethylated'), 2),
                           N=c(ageing_hypo, ageing_hyper, sotos_hypo, sotos_hyper));

## Barplots to show the numbers of DMPs.

barplot_DMPs <- ggplot(data=summary_DMPs, aes(x=Condition, y=N, fill=Direction)) +
  geom_bar(stat="identity") + scale_fill_manual(values=c('Red', 'Blue')) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  xlab("") + ylab("Number of DMPs") + 
  guides(fill=guide_legend(title="Methylation change"));
ggsave('plots/barplots_DMPs_ageing_Sotos.pdf', height=4, width=4);


## Histograms to show the direction of the DMPs.

effects <- data.frame(Condition=c(rep('Ageing', nrow(age_dmps_filtered)), rep('Sotos', nrow(sotos_dmps_filtered))),
                      Effect=c(age_dmps_filtered$beta, sotos_dmps_filtered$beta));
hist_ageing_sotos_effects <- ggplot(data=effects, aes(x=Effect)) +
  geom_histogram(color="black", fill="white") + facet_grid(. ~ Condition, scales = "free") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  xlab("Effect size") + ylab("Frequency");
ggsave('plots/hist_effect_ageing_Sotos.pdf', height=4, width=8);
  

## Venn diagram.

venn_df <- data.frame('Hyper aDMPs'=(age_dmps$ProbeID %in% age_dmps_filtered[age_dmps_filtered$beta>0,1]),
                      'Hypo aDMPs'=(age_dmps$ProbeID %in% age_dmps_filtered[age_dmps_filtered$beta<0,1]),
                      'Hyper Sotos DMPs'=(age_dmps$ProbeID %in% sotos_dmps_filtered[sotos_dmps_filtered$beta>0,1]),
                      'Hypo Sotos DMPs'=(age_dmps$ProbeID %in% sotos_dmps_filtered[sotos_dmps_filtered$beta<0,1]));
venn_counts <- vennCounts(venn_df);
pdf('plots/venn_ageing_sotos.pdf', height=6, width=7);
vennDiagram(venn_counts, cex=c(1.2,1,0.7),
            names=c("Hyper aDMPs", "Hypo aDMPs", "Hyper Sotos DMPs", "Hypo Sotos DMPs"));
dev.off();


# #### 3. Process the genomic features for the 450K probes. ####
# 
# ## Read all the features information available for the CpG probes.
# 
# epi_PBMC_info_450K <- as.data.frame(fread('../epigenetic_annotation/ENCODE_data/output_all_450k_CpGs_ENCODE_PBMC_FC_QCd_200bp.csv')); # Histone marks in PBMC
# EZH2_Bcell_info_450K <- as.data.frame(fread('../epigenetic_annotation/ENCODE_data/output_all_450k_CpGs_ENCODE_Bcell_FC_200bp.csv')); # EZH2 in B cells
# RNF2_K562_info_450K <- as.data.frame(fread('../epigenetic_annotation/ENCODE_data/output_all_450k_CpGs_ENCODE_K562_polycomb_200bp.csv')); # RNF2 in K562 cells
# RNA_PBMC_info_450K <- as.data.frame(fread('../epigenetic_annotation/ENCODE_data/output_all_450k_CpGs_ENCODE_PBMC_RNA_200bp.csv')); # Total RNA (rRNA depleted) in PBMC
# rep_lamin_IMR90_info_450K <- as.data.frame(fread('../epigenetic_annotation/replication_lamin/output_all_450k_CpGs_replication_lamin_200bp.csv')); # Replication timing and lamin B1 in IMR90 cells
# genomic_info_450K <- as.data.frame(fread('../epigenetic_annotation/genomic_annotation/genomic_annotation_450k.csv')); # Gene body, CGI, shore, shelf, ChrHMM state
# 
# ## Process the ENCODE Chip-Seq datasets.
# 
# # Make NAs equivalent to 0. 
# 
# epi_PBMC_info_450K[is.na(epi_PBMC_info_450K)] <- 0;
# EZH2_Bcell_info_450K[is.na(EZH2_Bcell_info_450K)] <- 0;
# RNF2_K562_info_450K[is.na(RNF2_K562_info_450K)] <- 0;
# RNA_PBMC_info_450K[is.na(RNA_PBMC_info_450K)] <- 0;
# rep_lamin_IMR90_info_450K[is.na(rep_lamin_IMR90_info_450K)] <- 0;
# genomic_info_450K[is.na(genomic_info_450K)] <- 0;
# 
# # Substitute fold change values of 0 for the minimum value.
# 
# subs_zeros <- function(df){
#   df_new <- df;
#   if(ncol(df)>3){
#     df_new[,3:ncol(df)] <- apply(df[,3:ncol(df)], 2, function(x){y <- x; y[y==0] <- min(y[y!=0]); return(y)});
#   }else{
#     df_new[df_new[,3]==0,3] <- min(df_new[df_new[,3]!=0,3]);
#   }
#   return(df_new);
# }
# 
# epi_PBMC_info_450K <- subs_zeros(epi_PBMC_info_450K);
# EZH2_Bcell_info_450K <- subs_zeros(EZH2_Bcell_info_450K);
# RNF2_K562_info_450K <- subs_zeros(RNF2_K562_info_450K);
# 
# # Scale them.
# 
# epi_PBMC_info_450K[,3:ncol(epi_PBMC_info_450K)] <- scale(epi_PBMC_info_450K[,3:ncol(epi_PBMC_info_450K)]);
# EZH2_Bcell_info_450K[,3:ncol(EZH2_Bcell_info_450K)] <- scale(EZH2_Bcell_info_450K[,3:ncol(EZH2_Bcell_info_450K)]);
# RNF2_K562_info_450K[,3:ncol(RNF2_K562_info_450K)] <- scale(RNF2_K562_info_450K[,3:ncol(RNF2_K562_info_450K)]);
# 
# ## Process the ENCODE RNA tracks.
# 
# RNA_PBMC_info_450K$RNA_seq_agg <- log2(1+RNA_PBMC_info_450K$RNA_seq_ENCFF754LBN+RNA_PBMC_info_450K$RNA_seq_ENCFF398HDS); # Aggregate both strands
# RNA_PBMC_info_450K$RNA_seq_agg <- as.numeric(scale(RNA_PBMC_info_450K$RNA_seq_agg)); # Calculate Z-scores
# 
# ## Reannotate ChromHMM states as in https://egg2.wustl.edu/roadmap/web_portal/imputed.html#chr_imp
# 
# genomic_info_450K$ChrHMM_state_reann <- genomic_info_450K$ChrHMM_state;
# genomic_info_450K$ChrHMM_state_reann <- ifelse(genomic_info_450K$ChrHMM_state_reann == 1, 'Active TSS', ifelse(
#   genomic_info_450K$ChrHMM_state_reann %in% c(2,3,4), 'Promoter', ifelse(
#     genomic_info_450K$ChrHMM_state_reann %in% c(5,6,7), 'Transcribed', ifelse(
#       genomic_info_450K$ChrHMM_state_reann == 8, 'Weakly transcribed', ifelse(
#         genomic_info_450K$ChrHMM_state_reann %in% c(9,10,11,12), 'Transcribed/regulatory', ifelse(
#           genomic_info_450K$ChrHMM_state_reann %in% c(13,14,15), 'Active enhancer', ifelse(
#             genomic_info_450K$ChrHMM_state_reann %in% c(16,17,18), 'Weak enhancer', ifelse(
#               genomic_info_450K$ChrHMM_state_reann == 19, 'DNase', ifelse(
#                 genomic_info_450K$ChrHMM_state_reann == 20, 'ZNF/repeats', ifelse(
#                   genomic_info_450K$ChrHMM_state_reann == 21, 'Heterochromatin', ifelse(
#                     genomic_info_450K$ChrHMM_state_reann == 22, 'Poised promoter', ifelse(
#                       genomic_info_450K$ChrHMM_state_reann == 23, 'Bivalent promoter', ifelse(
#                         genomic_info_450K$ChrHMM_state_reann == 24, 'Repressed polycomb', 'Quiescent/low')))))))))))));
# genomic_info_450K$ChrHMM_state_reann <- factor(genomic_info_450K$ChrHMM_state_reann);
# 
# ## Reannotate other features.
# 
# genomic_info_450K$Gene_body <- ifelse(genomic_info_450K$Gene_body, 'Yes', 'No');
# genomic_info_450K$CGI <- ifelse(genomic_info_450K$CGI, 'Yes', 'No');
# genomic_info_450K$Shore <- ifelse(genomic_info_450K$Shore, 'Yes', 'No');
# genomic_info_450K$Shelf <- ifelse(genomic_info_450K$Shelf, 'Yes', 'No');
# 
# ## Create a dataframe with all the individual features.
# 
# all_features_indiv <- merge(epi_PBMC_info_450K, EZH2_Bcell_info_450K, by=c("CpGmarker","chr:coord_hg19"));
# all_features_indiv <- merge(all_features_indiv, RNF2_K562_info_450K, by=c("CpGmarker","chr:coord_hg19"));
# all_features_indiv <- merge(all_features_indiv, RNA_PBMC_info_450K, by=c("CpGmarker","chr:coord_hg19"));
# all_features_indiv <- merge(all_features_indiv, rep_lamin_IMR90_info_450K, by=c("CpGmarker","chr:coord_hg19"));
# all_features_indiv <- merge(all_features_indiv, genomic_info_450K, by=c("CpGmarker","chr:coord_hg19"));
# 
# ## Create a dataframe with the aggregated features (mean between biological replicates).
# 
# multi_feats <- c("H3K27ac", "H3K4me3", "H3K36me3", "H3K27me3", "H3K9ac", "H3K4me1", "H3K9me3", "RNF2");
# all_features_agg <- data.frame(CpGmarker=all_features_indiv$CpGmarker,
#                                'chr:coord_hg19'=all_features_indiv$`chr:coord_hg19`,
#                                H3K27ac=apply(all_features_indiv[,grep('H3K27ac', colnames(all_features_indiv))], 1, mean),
#                                H3K4me3=apply(all_features_indiv[,grep('H3K4me3', colnames(all_features_indiv))], 1, mean),
#                                H3K36me3=apply(all_features_indiv[,grep('H3K36me3', colnames(all_features_indiv))], 1, mean),
#                                H3K27me3=apply(all_features_indiv[,grep('H3K27me3', colnames(all_features_indiv))], 1, mean),
#                                H3K9ac=apply(all_features_indiv[,grep('H3K9ac', colnames(all_features_indiv))], 1, mean),
#                                H3K4me1=apply(all_features_indiv[,grep('H3K4me1', colnames(all_features_indiv))], 1, mean),
#                                H3K9me3=apply(all_features_indiv[,grep('H3K9me3', colnames(all_features_indiv))], 1, mean),
#                                RNF2=apply(all_features_indiv[,grep('RNF2', colnames(all_features_indiv))], 1, mean),
#                                EZH2=all_features_indiv$EZH2_ENCFF516PTT,
#                                RNA=all_features_indiv$RNA_seq_agg,
#                                Replication_timing=all_features_indiv$Replication_timing_GSM923447,
#                                LaminB1=all_features_indiv$LaminB1_GSM1289416,
#                                Gene_body=all_features_indiv$Gene_body,
#                                CGI=all_features_indiv$CGI,
#                                Shore=all_features_indiv$Shore,
#                                Shelf=all_features_indiv$Shelf,
#                                ChrHMM_state=all_features_indiv$ChrHMM_state,
#                                ChrHMM_state_reann=all_features_indiv$ChrHMM_state_reann);
# write.table(all_features_agg, file='../epigenetic_annotation/all_features_agg_450k.csv', quote=F,
#             sep=',', row.names=F);

all_features_agg <- as.data.frame(fread('../epigenetic_annotation/all_features_agg_450k.csv'));


#### 4. Compare the continuous genomic features accross different subsets of CpGs. ####

## Define the subsets.

subsets <- c("Hyper aDMPs", "Hypo aDMPs", "Hyper Sotos DMPs", "Hypo Sotos DMPs", "Hypo-Hypo DMPs", "Hyper-Hypo DMPs");
all_cpgs <- all_features_agg$CpGmarker;
hyper_admps <- all_cpgs[all_cpgs %in% age_dmps_filtered$ProbeID[which(age_dmps_filtered$beta>0)]];
hypo_admps <- all_cpgs[all_cpgs %in% age_dmps_filtered$ProbeID[which(age_dmps_filtered$beta<0)]];
hyper_sotos_dmps <- all_cpgs[all_cpgs %in% sotos_dmps_filtered$ProbeID[which(sotos_dmps_filtered$beta>0)]];
hypo_sotos_dmps <- all_cpgs[all_cpgs %in% sotos_dmps_filtered$ProbeID[which(sotos_dmps_filtered$beta<0)]];
hypo_hypo_dmps <- intersect(hypo_admps, hypo_sotos_dmps);
hyper_hypo_dmps <- intersect(hyper_admps, hypo_sotos_dmps);
subsets_list <- list(hyper_admps, hypo_admps, hyper_sotos_dmps, hypo_sotos_dmps, hypo_hypo_dmps, hyper_hypo_dmps);
names(subsets_list) <- subsets;
  

## Create dataframes for the continuous features. 

comp_cont_features <- all_features_agg[,3:14];
rownames(comp_cont_features) <- all_features_agg$CpGmarker;
final_cont_df <- data.frame(matrix(NA, ncol=4, nrow=ncol(comp_cont_features)*length(subsets)*nrow(comp_cont_features)));
colnames(final_cont_df) <- c('Feature', 'Subset', 'Control_or_subset', 'Score');
final_cont_pv <- data.frame(matrix(NA, ncol=3, nrow=ncol(comp_cont_features)*length(subsets)));
colnames(final_cont_pv) <- c('Feature', 'Subset', 'log10');
i <- 1;
j <- 1;

for(f in colnames(comp_cont_features)){
  print(f);
  for(s in subsets){
    print(s);
    
    final_cont_df[i:(i+nrow(comp_cont_features)-1),1] <- rep(f, nrow(comp_cont_features));
    final_cont_df[i:(i+nrow(comp_cont_features)-1),2] <- rep(s, nrow(comp_cont_features));
    final_cont_df[i:(i+nrow(comp_cont_features)-1),3] <- ifelse(rownames(comp_cont_features) %in% subsets_list[[s]], 'In subset', 'Control');
    final_cont_df[i:(i+nrow(comp_cont_features)-1),4] <- comp_cont_features[, which(colnames(comp_cont_features)==f)];
    i <- i + nrow(comp_cont_features);
    
    in_subset <- (rownames(comp_cont_features) %in% subsets_list[[s]]);
    temp_log <- -log10(wilcox.test(comp_cont_features[!in_subset, which(colnames(comp_cont_features)==f)], 
                                   comp_cont_features[in_subset, which(colnames(comp_cont_features)==f)])$p.value);
    final_cont_pv[j,1] <- f; final_cont_pv[j,2] <-s; final_cont_pv[j,3] <- temp_log;
    j <- j + 1;
  }
}


## Create plots for comparisons. 

my_palette <- c('grey', 'yellow');
mc_1 <- list(c("In subset", "Control"));

give.n <- function(x){  # Adapted from https://stackoverflow.com/questions/15660829/how-to-add-a-number-of-observations-per-group-and-use-group-mean-in-ggplot2-boxp
  ref <- ifelse(median(x)>50, 0, -3);
  step <- ifelse(median(x)>50, 16, 1.5);
  return(c(y = ref+step, label=length(x)))
}
give.median <- function(x){
  ref <- ifelse(median(x)>50, 0, -3);
  step <- ifelse(median(x)>50, 8, 0.5);
  return(c(y = ref+step, label=round(median(x), digits=3)))
}

boxplot_subset_comparison_list <- list();

for(f in colnames(comp_cont_features)){
  
  ylim_min <- ifelse(f=='Replication_timing', 0, -3);
  ylim_max <- ifelse(f=='Replication_timing', 100, 4);
  y_label_edit <- ifelse(f=="RNA", "NRE", ifelse(f=="Replication_timing", "WTS", ifelse(f=="LaminB1", "NRC", "NFC")));
  boxplot_subset_comparison_list[[f]] <- ggboxplot(final_cont_df[final_cont_df$Feature==f,], x = "Control_or_subset", y = "Score", 
                                                   fill = "Control_or_subset", palette = my_palette, alpha=0.5, 
                                                   title=paste0('Feature: ', f), outlier.shape = NA) +
    facet_wrap(~Subset, nrow=1)+
    theme_classic() +
    theme(axis.text=element_text(size=12, angle=90),
          axis.title=element_text(size=14,face="bold"),
          plot.title = element_text(hjust = 0.5, face="bold"),
          strip.text = element_text(size=6)) +
    xlab("") + ylab(y_label_edit) + 
    stat_compare_means(comparisons = mc_1, label.y=ifelse(f=='Replication_timing', 90,3), tip.length = 0) +
    guides(colour=FALSE, fill=FALSE) + coord_cartesian(ylim=c(ylim_min, ylim_max)) +
    stat_summary(fun.data = give.n, geom = "text", fun.y = median, col='darkgreen', size=3) +
    stat_summary(fun.data = give.median, geom = "text", fun.y = mean, col='darkred', size=3);
}

# Plot all the features.

glist <- lapply(boxplot_subset_comparison_list, ggplotGrob);
ggsave("plots/ageing_vs_sotos_subset_boxplot_comparisons_cont.pdf", marrangeGrob(glist, nrow=4, ncol=3, 
                                                                       layout_matrix=matrix(1:12, nrow=4, ncol=3, byrow=T), top=NULL),
       height=20, width=19);

# barplot_subset_comparison_p_values <- ggplot(data=final_cont_pv, aes(x=Feature, y=log10, fill=Subset)) + 
#   geom_bar(position="dodge", stat="identity") + scale_fill_stata() +
#   theme_classic() +
#   theme(axis.text=element_text(size=12),
#         axis.text.x=element_text(angle = -90, hjust = 0, size=10),
#         axis.title=element_text(size=14,face="bold"),
#         plot.title = element_text(hjust = 0.5, face="bold"),
#         plot.subtitle = element_text(hjust = 0.5),
#         plot.margin=unit(c(1,1,-0.5,1), "cm")) +
#   xlab("Feature") + ylab(expression(bold(paste(-log[10](P-value))))) +
#   guides(fill=guide_legend(title="Subset")) + 
#   geom_abline(slope=0, intercept=alpha, linetype=2, col='black', size=0.8);
# ggsave("plots/ageing_vs_sotos_subset_barplots_p_values_cont.pdf", height=6, width=10);

# Plot only RNA.

give.n.RNA <- function(x){  # Adapted from https://stackoverflow.com/questions/15660829/how-to-add-a-number-of-observations-per-group-and-use-group-mean-in-ggplot2-boxp
  return(c(y = -0.45, label=length(x)))
}
give.median.RNA <- function(x){
  return(c(y = -0.5, label=round(median(x), digits=3)))
}

RNA_plot <- ggboxplot(final_cont_df[final_cont_df$Feature=='RNA' & final_cont_df$Subset %in% c('Hypo aDMPs', 'Hypo Sotos DMPs', 'Hypo-Hypo DMPs'),], 
                      x = "Control_or_subset", y = "Score", 
                      fill = "Control_or_subset", palette = my_palette, alpha=0.5, 
                      title='Feature: RNA', outlier.shape = NA) +
  facet_wrap(~Subset, nrow=1)+
  theme_classic() +
  theme(axis.text=element_text(size=12, angle=90),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5, face="bold")) +
  xlab("") + ylab("NRE") + 
  stat_compare_means(comparisons = mc_1, label.y=0.5, tip.length = 0) +
  guides(colour=FALSE, fill=FALSE) + coord_cartesian(ylim=c(-0.6, 0.6)) +
  stat_summary(fun.data = give.n.RNA, geom = "text", fun.y = median, col='darkgreen') +
  stat_summary(fun.data = give.median.RNA, geom = "text", fun.y = mean, col='darkred');
ggsave('plots/RNA_comparison_ageing_Sotos.pdf', height=6, width=6);

# Plot only H3K36me3.

give.n.H3K36me3 <- function(x){  # Adapted from https://stackoverflow.com/questions/15660829/how-to-add-a-number-of-observations-per-group-and-use-group-mean-in-ggplot2-boxp
  return(c(y = -0.50, label=length(x)))
}
give.median.H3K36me3 <- function(x){
  return(c(y = -0.6, label=round(median(x), digits=3)))
}

H3K36me3_plot <- ggboxplot(final_cont_df[final_cont_df$Feature=='H3K36me3' & final_cont_df$Subset %in% c('Hypo aDMPs', 'Hypo Sotos DMPs', 'Hypo-Hypo DMPs'),], 
                      x = "Control_or_subset", y = "Score", 
                      fill = "Control_or_subset", palette = my_palette, alpha=0.5, 
                      title='Feature: H3K36me3', outlier.shape = NA) +
  facet_wrap(~Subset, nrow=1)+
  theme_classic() +
  theme(axis.text=element_text(size=12, angle=90),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5, face="bold")) +
  xlab("") + ylab("NFC") + 
  stat_compare_means(comparisons = mc_1, label.y=0.5, tip.length = 0) +
  guides(colour=FALSE, fill=FALSE) + coord_cartesian(ylim=c(-0.6, 0.6)) +
  stat_summary(fun.data = give.n.H3K36me3, geom = "text", fun.y = median, col='darkgreen') +
  stat_summary(fun.data = give.median.H3K36me3, geom = "text", fun.y = mean, col='darkred');
ggsave('plots/H3K36me3_comparison_ageing_Sotos.pdf', height=6, width=6);


#### 5. Compare the categorical genomic features accross different subsets of CpGs. ####

## Create dataframe for the categorical features. 

comp_cat_features <- all_features_agg[,15:19];
rownames(comp_cat_features) <- all_features_agg$CpGmarker;

# Expand the ChromHMM states.

comp_cat_features <- cbind(comp_cat_features, as.data.frame(matrix(NA, nrow=nrow(comp_cat_features), ncol=25)));
for(state in 1:25){
  comp_cat_features[,state+5] <- ifelse(comp_cat_features$ChrHMM_state==state, 'Yes', 'No');
}
colnames(comp_cat_features)[6:ncol(comp_cat_features)] <- c('Active TSS', 'Promoter Upstream TSS', 'Promoter Downstream TSS 1', 'Promoter Downstream TSS 2', 
                                                            "Transcribed - 5' preferential",'Strong transcription', "Transcribed - 3' preferential", 'Weak transcription', 
                                                            'Transcribed & regulatory (Prom/Enh)', "Transcribed 5' preferential and Enh", "Transcribed 3' preferential and Enh",
                                                            "Transcribed and Weak Enhancer", 'Active Enhancer 1', 'Active Enhancer 2', 'Active Enhancer Flank', 'Weak Enhancer 1',
                                                            'Weak Enhancer 2', 'Primary H3K27ac possible Enhancer', 'Primary DNase', 'ZNF genes & repeats', 'Heterochromatin', 
                                                            'Poised promoter', 'Bivalent promoter', 'Repressed polycomb', 'Quiescent/low');
comp_cat_features <- comp_cat_features[,-5];

final_cat_df <- data.frame(matrix(NA, ncol=6, 
                                  nrow=(ncol(comp_cat_features)*length(subsets))));
colnames(final_cat_df) <- c('Feature', 'Subset', 'OR', "log10Pvalue", "LowerOR", "UpperOR");
i <- 1;

for(f in colnames(comp_cat_features)){
  print(f);
  for(s in subsets){
    print(s);
    final_cat_df[i,1:2] <- c(f,s);
    insubset <- comp_cat_features[rownames(comp_cat_features) %in% subsets_list[[s]],colnames(comp_cat_features)==f];
    outsubset <- comp_cat_features[!(rownames(comp_cat_features) %in% subsets_list[[s]]),colnames(comp_cat_features)==f];
    #random_s <- rownames(comp_cat_features)[sample(1:nrow(comp_cat_features), 353, replace=F)];
    #insubset <- comp_cat_features[rownames(comp_cat_features) %in% random_s,colnames(comp_cat_features)==f];
    #outsubset <- comp_cat_features[!(rownames(comp_cat_features) %in% random_s),colnames(comp_cat_features)==f];
    cont_table <- matrix(c(sum(insubset=="Yes"), sum(insubset=="No"), sum(outsubset=="Yes"), sum(outsubset=="No")),
                         byrow = TRUE, ncol=2);
    fisher_test <- fisher.test(cont_table);
    final_cat_df[i,3:4] <- c(as.numeric(fisher_test$estimate), -log10(as.numeric(fisher_test$p.value)));
    final_cat_df[i,5:6] <- as.numeric(fisher_test$conf.int);
    i <- i + 1;
  }
}

## Create plots for comparisons. 

alpha <- -log10(0.01/(nrow(final_cat_df)));


# Put cutoffs in the data displayed (in this case: 0.005 <= OR <= 10). 
final_cat_df_mod <- final_cat_df;
max_lim_log10 <- 100;
final_cat_df_mod$log10Pvalue <- ifelse(is.infinite(final_cat_df_mod$log10Pvalue) | final_cat_df_mod$log10Pvalue>max_lim_log10,
                                       max_lim_log10, final_cat_df_mod$log10Pvalue);
max_lim_OR <- 10;
final_cat_df_mod$OR <- ifelse(final_cat_df_mod$OR<0.005, 0.005, final_cat_df_mod$OR);
final_cat_df_mod$LowerOR <- ifelse(final_cat_df_mod$LowerOR<0.005, 0.005, final_cat_df_mod$LowerOR);
final_cat_df_mod$UpperOR <- ifelse(final_cat_df_mod$UpperOR>max_lim_OR, max_lim_OR, final_cat_df_mod$UpperOR);

final_cat_df_mod$log10Pvaluecol <- final_cat_df_mod$log10Pvalue;
final_cat_df_mod$log10Pvaluecol <- ifelse(final_cat_df_mod$log10Pvaluecol<alpha, NA, final_cat_df_mod$log10Pvaluecol);

# barplot_subset_comparison_p_values_cat <- ggplot(data=final_cat_df_mod, aes(x=Feature, y=log10Pvalue, fill=Feature)) + 
#   geom_bar(position="dodge", stat="identity") + scale_fill_viridis(discrete=T) +
#   facet_wrap(~Subset, ncol=2, nrow=3) +
#   theme_classic() +
#   theme(axis.text=element_text(size=12),
#         axis.text.x=element_text(angle = -90, hjust = 0, size=10),
#         axis.title=element_text(size=14,face="bold"),
#         plot.title = element_text(hjust = 0.5, face="bold"),
#         plot.subtitle = element_text(hjust = 0.5),
#         plot.margin=unit(c(1,1,-0.5,1), "cm")) +
#   xlab("Feature") + ylab(expression(bold(paste(-log[10](P-value))))) +
#   guides(fill=FALSE) + 
#   geom_abline(slope=0, intercept=alpha, linetype=2, col='black', size=0.6);
# ggsave("plots/ageing_vs_sotos_subset_barplots_p_values_cat.pdf", height=10, width=10);

OR_plot_subset_comparison_cat <- ggplot(final_cat_df_mod, aes(x = Feature, y = OR)) + 
  facet_wrap(~Subset, ncol=2, nrow=3) +
  geom_point(aes(col=log10Pvaluecol), size = 2) + scale_color_viridis(discrete=F, na.value='grey') +
  geom_vline(aes(xintercept = 0), size = .25, linetype = "dashed") + 
  geom_errorbar(aes(ymax = UpperOR, ymin = LowerOR), size = .35, color = "grey66") +
  theme_classic()+
  theme(axis.text=element_text(size=12),
        axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.5, size=10),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5, face="bold"),
        plot.subtitle = element_text(hjust = 0.5),
        plot.margin=unit(c(1,1,-0.5,1), "cm")) +
  scale_y_log10(limits=c(0.005,max_lim_OR)) +
  geom_abline(slope=0, intercept=0, linetype=2, col='black', size=0.6) +
  ylab("Odds ratio") + xlab("Feature") +labs(colour=expression(bold(paste(-log[10](P-value)))));
ggsave("plots/ageing_vs_sotos_subset_OR_plot_cat.pdf", height=10, width=10);

# Plot only some of the subsets.

OR_plot_selected <- ggplot(final_cat_df_mod[final_cat_df_mod$Subset %in% c('Hypo aDMPs', 'Hypo Sotos DMPs', 'Hypo-Hypo DMPs'),],
                                        aes(x = Feature, y = OR)) + 
  facet_wrap(~Subset, ncol=1, nrow=3) +
  geom_point(aes(col=log10Pvaluecol), size = 2) + scale_color_viridis(discrete=F, na.value='grey') +
  geom_vline(aes(xintercept = 0), size = .25, linetype = "dashed") + 
  geom_errorbar(aes(ymax = UpperOR, ymin = LowerOR), size = .35, color = "grey66") +
  theme_classic()+
  theme(axis.text=element_text(size=12),
        axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.5, size=10),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5, face="bold"),
        plot.subtitle = element_text(hjust = 0.5),
        plot.margin=unit(c(1,1,-0.5,1), "cm"),
        legend.text=element_text(size=8),
        legend.title=element_text(size=10)) +
  scale_y_log10(limits=c(0.005,max_lim_OR)) +
  geom_abline(slope=0, intercept=0, linetype=2, col='black', size=0.6) +
  ylab("Odds ratio") + xlab("Feature") +labs(colour=expression(bold(paste(-log[10](P-value)))));
ggsave("plots/ageing_vs_sotos_OR_selected.pdf", height=10, width=6.5);


################################################################
################## End of the script ###########################
################################################################