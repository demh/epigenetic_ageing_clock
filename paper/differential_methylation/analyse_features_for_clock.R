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
##### Compare the Horvath clock CpG sites with the differentially methylated positions #### 
##### (DMPs) found during ageing and in Sotos syndrome. This includes finding genomic #####
##### features that are characteristic of the different subsets.                      #####
###########################################################################################
##### USAGE: manual                                                                   #####
###########################################################################################

###########################################################
##################### Dependencies ########################
###########################################################

library(data.table);
library(ggthemes);
library(ComplexHeatmap);
library(circlize);
library(viridis);
library(colorRamps);
library(ggpubr);
library(gridExtra);

setwd('/Users/dem44/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/differential_methylation/');


################################################################
################## Run the pipeline ############################
################################################################

#### 1. Read input data. #### 

## Export the probes (21368) that were used to build the Horvath model (background).

# all_21k_probes <- as.data.frame(fread('../epigenetic_annotation/probeAnnotation21kdatMethUsed.csv'));
# 
# data("IlluminaHumanMethylation450kanno.ilmn12.hg19");
# data("Locations");
# hg19_coords <- Locations;
# all_21k_hg19_coords <- as.data.frame(hg19_coords[which(rownames(hg19_coords) %in% all_21k_probes$Name),]);
# all_21k_hg19_coords$CpGmarker <- rownames(all_21k_hg19_coords);
# rownames(all_21k_hg19_coords) <- NULL;
# write.table(x=all_21k_hg19_coords[,c(4,1,2)], file='all_21k_CpGs_coords_hg19.csv', quote=F,
#             sep=',', row.names=F);

## Horvath's clock sites annotation.

clock_sites_ann <- as.data.frame(fread('../epigenetic_annotation/AdditionalFile3.csv'))[-1,];

## Read all the features information available for the CpG probes.

epi_PBMC_info_21K <- as.data.frame(fread('../epigenetic_annotation/ENCODE_data/output_all_21k_CpGs_ENCODE_PBMC_FC_QCd_200bp.csv')); # Histone marks in PBMC
EZH2_Bcell_info_21K <- as.data.frame(fread('../epigenetic_annotation/ENCODE_data/output_all_21k_CpGs_ENCODE_Bcell_FC_200bp.csv')); # EZH2 in B cells
RNF2_K562_info_21K <- as.data.frame(fread('../epigenetic_annotation/ENCODE_data/output_all_21k_CpGs_ENCODE_K562_polycomb_200bp.csv')); # RNF2 in K562 cells
RNA_PBMC_info_21K <- as.data.frame(fread('../epigenetic_annotation/ENCODE_data/output_all_21k_CpGs_ENCODE_PBMC_RNA_200bp.csv')); # Total RNA (rRNA depleted) in PBMC
rep_lamin_IMR90_info_21K <- as.data.frame(fread('../epigenetic_annotation/replication_lamin/output_all_21k_CpGs_replication_lamin_200bp.csv')); # Replication timing and lamin B1 in IMR90 cells
genomic_info_21K <- as.data.frame(fread('../epigenetic_annotation/genomic_annotation/genomic_annotation_21k.csv')); # Gene body, CGI, shore, shelf, ChrHMM state


#### 2. Process the genomic features. ####

## Process the ENCODE Chip-Seq datasets.

# Make NAs equivalent to 0.

epi_PBMC_info_21K[is.na(epi_PBMC_info_21K)] <- 0;
EZH2_Bcell_info_21K[is.na(EZH2_Bcell_info_21K)] <- 0;
RNF2_K562_info_21K[is.na(RNF2_K562_info_21K)] <- 0;
RNA_PBMC_info_21K[is.na(RNA_PBMC_info_21K)] <- 0;
rep_lamin_IMR90_info_21K[is.na(rep_lamin_IMR90_info_21K)] <- 0;
genomic_info_21K[is.na(genomic_info_21K)] <- 0;

# Substitute fold change values of 0 for the minimum value.

subs_zeros <- function(df){
  df_new <- df;
  if(ncol(df)>3){
    df_new[,3:ncol(df)] <- apply(df[,3:ncol(df)], 2, function(x){y <- x; y[y==0] <- min(y[y!=0]); return(y)});
  }else{
    df_new[df_new[,3]==0,3] <- min(df_new[df_new[,3]!=0,3]);
  }
  return(df_new);
}

epi_PBMC_info_21K <- subs_zeros(epi_PBMC_info_21K);
EZH2_Bcell_info_21K <- subs_zeros(EZH2_Bcell_info_21K);
RNF2_K562_info_21K <- subs_zeros(RNF2_K562_info_21K);

# Scale them.

epi_PBMC_info_21K[,3:ncol(epi_PBMC_info_21K)] <- scale(epi_PBMC_info_21K[,3:ncol(epi_PBMC_info_21K)]);
EZH2_Bcell_info_21K[,3:ncol(EZH2_Bcell_info_21K)] <- scale(EZH2_Bcell_info_21K[,3:ncol(EZH2_Bcell_info_21K)]);
RNF2_K562_info_21K[,3:ncol(RNF2_K562_info_21K)] <- scale(RNF2_K562_info_21K[,3:ncol(RNF2_K562_info_21K)]);

## Process the ENCODE RNA tracks.

RNA_PBMC_info_21K$RNA_seq_agg <- log2(1+RNA_PBMC_info_21K$RNA_seq_ENCFF754LBN+RNA_PBMC_info_21K$RNA_seq_ENCFF398HDS); # Aggregate both strands
RNA_PBMC_info_21K$RNA_seq_agg <- as.numeric(scale(RNA_PBMC_info_21K$RNA_seq_agg)); # Calculate Z-scores

## Reannotate ChromHMM states as in https://egg2.wustl.edu/roadmap/web_portal/imputed.html#chr_imp

genomic_info_21K$ChrHMM_state_reann <- genomic_info_21K$ChrHMM_state;
genomic_info_21K$ChrHMM_state_reann <- ifelse(genomic_info_21K$ChrHMM_state_reann == 1, 'Active TSS', ifelse(
  genomic_info_21K$ChrHMM_state_reann %in% c(2,3,4), 'Promoter', ifelse(
    genomic_info_21K$ChrHMM_state_reann %in% c(5,6,7), 'Transcribed', ifelse(
      genomic_info_21K$ChrHMM_state_reann == 8, 'Weakly transcribed', ifelse(
        genomic_info_21K$ChrHMM_state_reann %in% c(9,10,11,12), 'Transcribed/regulatory', ifelse(
          genomic_info_21K$ChrHMM_state_reann %in% c(13,14,15), 'Active enhancer', ifelse(
            genomic_info_21K$ChrHMM_state_reann %in% c(16,17,18), 'Weak enhancer', ifelse(
              genomic_info_21K$ChrHMM_state_reann == 19, 'DNase', ifelse(
                genomic_info_21K$ChrHMM_state_reann == 20, 'ZNF/repeats', ifelse(
                  genomic_info_21K$ChrHMM_state_reann == 21, 'Heterochromatin', ifelse(
                    genomic_info_21K$ChrHMM_state_reann == 22, 'Poised promoter', ifelse(
                      genomic_info_21K$ChrHMM_state_reann == 23, 'Bivalent promoter', ifelse(
                        genomic_info_21K$ChrHMM_state_reann == 24, 'Repressed polycomb', 'Quiescent/low')))))))))))));
genomic_info_21K$ChrHMM_state_reann <- factor(genomic_info_21K$ChrHMM_state_reann);

## Reannotate other features. 

genomic_info_21K$Gene_body <- ifelse(genomic_info_21K$Gene_body, 'Yes', 'No');
genomic_info_21K$CGI <- ifelse(genomic_info_21K$CGI, 'Yes', 'No');
genomic_info_21K$Shore <- ifelse(genomic_info_21K$Shore, 'Yes', 'No');
genomic_info_21K$Shelf <- ifelse(genomic_info_21K$Shelf, 'Yes', 'No');

## Create a dataframe with all the individual features.

all_features_indiv <- merge(epi_PBMC_info_21K, EZH2_Bcell_info_21K, by=c("CpGmarker","chr:coord_hg19"));
all_features_indiv <- merge(all_features_indiv, RNF2_K562_info_21K, by=c("CpGmarker","chr:coord_hg19"));
all_features_indiv <- merge(all_features_indiv, RNA_PBMC_info_21K, by=c("CpGmarker","chr:coord_hg19"));
all_features_indiv <- merge(all_features_indiv, rep_lamin_IMR90_info_21K, by=c("CpGmarker","chr:coord_hg19"));
all_features_indiv <- merge(all_features_indiv, genomic_info_21K, by=c("CpGmarker","chr:coord_hg19"));

## Create a dataframe with the aggregated features (mean between biological replicates).

multi_feats <- c("H3K27ac", "H3K4me3", "H3K36me3", "H3K27me3", "H3K9ac", "H3K4me1", "H3K9me3", "RNF2");
all_features_agg <- data.frame(CpGmarker=all_features_indiv$CpGmarker, 
                               'chr:coord_hg19'=all_features_indiv$`chr:coord_hg19`,
                               H3K27ac=apply(all_features_indiv[,grep('H3K27ac', colnames(all_features_indiv))], 1, mean),
                               H3K4me3=apply(all_features_indiv[,grep('H3K4me3', colnames(all_features_indiv))], 1, mean),
                               H3K36me3=apply(all_features_indiv[,grep('H3K36me3', colnames(all_features_indiv))], 1, mean),
                               H3K27me3=apply(all_features_indiv[,grep('H3K27me3', colnames(all_features_indiv))], 1, mean),
                               H3K9ac=apply(all_features_indiv[,grep('H3K9ac', colnames(all_features_indiv))], 1, mean),
                               H3K4me1=apply(all_features_indiv[,grep('H3K4me1', colnames(all_features_indiv))], 1, mean),
                               H3K9me3=apply(all_features_indiv[,grep('H3K9me3', colnames(all_features_indiv))], 1, mean),
                               RNF2=apply(all_features_indiv[,grep('RNF2', colnames(all_features_indiv))], 1, mean),
                               EZH2=all_features_indiv$EZH2_ENCFF516PTT,
                               RNA=all_features_indiv$RNA_seq_agg,
                               Replication_timing=all_features_indiv$Replication_timing_GSM923447,
                               LaminB1=all_features_indiv$LaminB1_GSM1289416,
                               Gene_body=all_features_indiv$Gene_body,
                               CGI=all_features_indiv$CGI,
                               Shore=all_features_indiv$Shore,
                               Shelf=all_features_indiv$Shelf,
                               ChrHMM_state=all_features_indiv$ChrHMM_state,
                               ChrHMM_state_reann=all_features_indiv$ChrHMM_state_reann);
clock_sites_ann_with_agg_features <- merge(clock_sites_ann, all_features_agg, by='CpGmarker');
write.table(clock_sites_ann_with_agg_features, file='clock_CpGs_extended_annotation.csv', quote=F,
            sep=',', row.names=F); # Create the extended annotation to be used by analyse_horvath_clock_sites_only_Sotos.R


#### 3. Create heatmap / hierarchical clustering with the ChIP-seq features for the clock CpG sites. ####

## Read clock sites extended annotation after adding info with analyse_horvath_clock_sites_only_Sotos.R

clock_sites_ann_with_agg_features <- as.data.frame(fread('clock_CpGs_extended_annotation_final.csv'));

## Format data. 

continuous_data_all <- all_features_indiv[all_features_indiv$CpGmarker %in% clock_sites_ann_with_agg_features$CpGmarker,3:25];
rownames(continuous_data_all) <- all_features_indiv$CpGmarker[all_features_indiv$CpGmarker %in% clock_sites_ann_with_agg_features$CpGmarker];
clock_sites_ann_with_agg_features <- clock_sites_ann_with_agg_features[match(rownames(continuous_data_all), clock_sites_ann_with_agg_features$CpGmarker),]; # Order

## Annotation. 

main_col_scale <- colorRamp2(c(-3,-1.5, 0, 1.5, 3),viridis(5));
col_Sotos <- c('Hypermethylated'='Red', "Hypomethylated"='Blue');
#col_sp <- colorRamp2(c(-1, 0, 1), c('#614051', 'White', '#cd9d0b'));
col_aDMP <- c('Hypermethylated'='#cd9d0b', "Hypomethylated"='#614051');
#col_horvath <- c("Hyper"="#cd9d0b", "Hypo"="#614051");
col_weights <- colorRamp2(c(-1, 0, 1), magenta2green(3)[c(3,2,1)]); #c("#d8bfd8", "White", "#ffe200"));
col_chrhmm <- c('Active TSS'=rgb(255,0,0, maxColorValue = 255), 'Promoter'=rgb(255,69,0, maxColorValue = 255), 'Transcribed'=rgb(0,128,0, maxColorValue = 255), 'Weakly transcribed'=rgb(0,150,0, maxColorValue = 255),
                'Transcribed/regulatory'=rgb(194,225,5, maxColorValue = 255), 'Active enhancer'=rgb(255,195,77, maxColorValue = 255), 'Weak enhancer'=rgb(255,255,0, maxColorValue = 255),
                'DNase'=rgb(255,255,102, maxColorValue = 255), 'ZNF/repeats'=rgb(102,205,170, maxColorValue = 255), 'Heterochromatin'=rgb(138,145,208, maxColorValue = 255), 'Poised promoter'=rgb(230,184,183, maxColorValue = 255),
                'Bivalent promoter'=rgb(112,48,160, maxColorValue = 255), 'Repressed polycomb'=rgb(128,128,128, maxColorValue = 255), 'Quiescent/low'="lightgrey");
col_RNA <- colorRamp2(c(-2, 0, 2), c('#8A91D0', 'White', '#008000'));
col_gb <- c('Yes'='Blue4', 'No'='gray');

ha_cgs <- HeatmapAnnotation(df=data.frame(Sotos_diff=clock_sites_ann_with_agg_features$Sotos_diff_i_discrete,
                                          #Spearman=clock_sites_ann_with_agg_features$sp_ours,
                                          aDMP=clock_sites_ann_with_agg_features$aDMP,
                                          #Horvath_hyperhypo=ifelse(clock_sites_ann_with_agg_features$hyperhypo_horvath=='Hyper', "Hypermethylated", "Hypomethylated"),
                                          Weight=clock_sites_ann_with_agg_features$CoefficientTraining,
                                          #CGI=clock_sites_ann_with_agg_features$CGI,
                                          ChrHMM=clock_sites_ann_with_agg_features$ChrHMM_state_reann,
                                          #rep=clock_sites_ann_with_agg_features$Replication_timing_GSM923447,
                                          #lamin=clock_sites_ann_with_agg_features$LaminB1_GSM1289416,
                                          RNA=clock_sites_ann_with_agg_features$RNA,
                                          Gene_bodies=clock_sites_ann_with_agg_features$Gene_body),
                            col=list(Sotos_diff=col_Sotos,
                                     #Spearman=col_sp,
                                     aDMP=col_aDMP,
                                     #Horvath_hyperhypo=col_horvath,
                                     Weight=col_weights,
                                     ChrHMM=col_chrhmm,
                                     RNA=col_RNA,
                                     Gene_bodies=col_gb),
                            annotation_legend_param=list(Sotos_diff=list(title="Sotos DMPs"),
                                                         aDMP=list(title="aDMPs"),
                                                         #Spearman=list(title='Age correlation \n(in controls)',
                                                         #             labels=c('-1 (hypomethylation)', '-0.5', '0', '0.5', '1 (hypermethylation)')),
                                                         Weight=list(title="Weight \nin model"),
                                                         ChrHMM=list(title="ChrHMM state \n(in K562)"),
                                                         RNA=list(title="RNA \n(in PBMC)"),
                                                         Gene_bodies=list(title="In gene body"))); 


feature_df <- data.frame(Feature=colnames(continuous_data_all), Cell_type=ifelse(colnames(continuous_data_all)=='EZH2_ENCFF516PTT', 'B cell', ifelse(
  grepl('RNF2', colnames(continuous_data_all)), 'K562', 'PBMC')));

ha_features <- rowAnnotation(df=data.frame(Cell_type=feature_df$Cell_type),
                             col=list(Cell_type=c('B cell'= stata_pal()(3)[1], 'K562'=stata_pal()(3)[2], 'PBMC'=stata_pal()(3)[3])),
                             annotation_legend_param=list(Cell_type=list(title='Cell type')));


# Plot. 

heatmap_clock_sites <- Heatmap(t(continuous_data_all), col=main_col_scale,
                               row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 2), row_names_side="left", 
                               bottom_annotation=ha_cgs, heatmap_legend_param = list(title = "Z-score \n(in PBMC)"),
                               column_title="Horvath's clock CpGs", column_title_gp = gpar(fontsize = 15, fontface = "bold"), column_title_side = "top",
                               row_title="Features", row_title_gp = gpar(fontsize = 15, fontface = "bold")) + ha_features;

pdf('plots/heatmap_clock_sites_annotated.pdf', width=10, height=8);
draw(heatmap_clock_sites, annotation_legend_side = "bottom");
dev.off();


#### 4. Compare the continuous genomic features accross different subsets of CpGs (21K, 353 clock sites, Hypo Sotos, ...). ####

## Define the subsets (change for 21K / 450K). 

subsets <- c("All Horvath", "Hyper aDMPs", "Hypo aDMPs", "Hypo Sotos DMPs");
cgs_allclock <- clock_sites_ann_with_agg_features$CpGmarker;
cgs_hyperclock <- clock_sites_ann_with_agg_features$CpGmarker[which(clock_sites_ann_with_agg_features$aDMP=="Hypermethylated")];
cgs_hypoclock <- clock_sites_ann_with_agg_features$CpGmarker[which(clock_sites_ann_with_agg_features$aDMP=="Hypomethylated")];
cgs_hyposotos <- clock_sites_ann_with_agg_features$CpGmarker[which(clock_sites_ann_with_agg_features$Sotos_diff_i_discrete=="Hypomethylated")];
subsets_list <- list(cgs_allclock, cgs_hyperclock, cgs_hypoclock, cgs_hyposotos);
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
          plot.title = element_text(hjust = 0.5, face="bold")) +
    xlab("") + ylab(y_label_edit) + 
    stat_compare_means(comparisons = mc_1, label.y=ifelse(f=='Replication_timing', 90,3), tip.length = 0) +
    guides(colour=FALSE, fill=FALSE) + coord_cartesian(ylim=c(ylim_min, ylim_max)) +
    stat_summary(fun.data = give.n, geom = "text", fun.y = median, col='darkgreen') +
    stat_summary(fun.data = give.median, geom = "text", fun.y = mean, col='darkred');
}

glist <- lapply(boxplot_subset_comparison_list, ggplotGrob)
ggsave("plots/clock_subset_boxplot_comparisons_cont.pdf", marrangeGrob(glist, nrow=4, ncol=3, 
                                                            layout_matrix=matrix(1:12, nrow=4, ncol=3, byrow=T), top=NULL),
       height=20, width=18);

alpha <- -log10(0.01/(nrow(final_cont_pv)));

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
# ggsave("plots/clock_subset_barplots_p_values_cont.pdf", height=6, width=10);
 

#### 5. Compare the categorical genomic features accross different subsets of CpGs (21K, 353 clock sites, Hypo Sotos, ...). ####

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

# barplot_subset_comparison_p_values_cat <- ggplot(data=final_cat_df, aes(x=Feature, y=log10Pvalue, fill=Feature)) + 
#   geom_bar(position="dodge", stat="identity") + scale_fill_viridis(discrete=T) +
#   facet_wrap(~Subset, ncol=2, nrow=2) +
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
# ggsave("plots/clock_subset_barplots_p_values_cat.pdf", height=6, width=10);


# Put cutoffs in the data displayed (in this case: 0.005 <= OR <= 10). 
final_cat_df_mod <- final_cat_df;
max_lim_log10 <- 100;
final_cat_df_mod$log10Pvalue <- ifelse(is.infinite(final_cat_df_mod$log10Pvalue) | final_cat_df_mod$log10Pvalue>max_lim_log10,
                                       max_lim_log10, final_cat_df_mod$log10Pvalue);
max_lim_OR <- 15;
final_cat_df_mod$OR <- ifelse(final_cat_df_mod$OR<0.005, 0.005, final_cat_df_mod$OR);
final_cat_df_mod$LowerOR <- ifelse(final_cat_df_mod$LowerOR<0.005, 0.005, final_cat_df_mod$LowerOR);
final_cat_df_mod$UpperOR <- ifelse(final_cat_df_mod$UpperOR>max_lim_OR, max_lim_OR, final_cat_df_mod$UpperOR);

final_cat_df_mod$log10Pvaluecol <- final_cat_df_mod$log10Pvalue;
final_cat_df_mod$log10Pvaluecol <- ifelse(final_cat_df_mod$log10Pvaluecol<alpha, NA, final_cat_df_mod$log10Pvaluecol);

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
ggsave("plots/clock_subset_OR_plot_cat.pdf", height=8, width=10);

################################################################
################## End of the script ###########################
################################################################
