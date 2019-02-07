###########################################################################################
#########                                                                         #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             28/08/2018                                  #########
#########                              EMBL-EBI                                   #########
#########                           Thornton group                                #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####              Biological insights into the epigenetic ageing clock           #########
###########################################################################################
##### Benchmark the background correction in the context of the Horvath clock.        #####
##### i.e. should we apply NOOB before calculating DNAmAge or not?.                   #####
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
library(ggthemes);
library(tidyr);
library(ComplexHeatmap);
library(circlize);

setwd('~/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/blood_control/benchmark_background_correction/');


################################################################
################## Functions ###################################
################################################################

#### Function: find the elbow in a curve f(x). 
#              Taken from https://github.com/BIMSBbioinfo/AmpliconBiSeq/blob/master/R/AmpliconViews_function.R
#
# vars: vector with the f(x) values.

findElbow<-function(vars){
  
  # if there are too few points return 1
  if( length(vars) <= 2){
    return(1)
  }
  
  nPoints = length(vars)
  allCoord <- cbind(1:nPoints,vars)              
  
  # pull out first point
  firstPoint = allCoord[1,];
  
  # get vector between first and last point - this is the line
  lineVec = allCoord[nrow(allCoord),] - firstPoint;
  
  # normalize the line vector
  lineVecN = lineVec / sqrt(sum(lineVec**2));
  
  # find the distance from each point to the line:
  # vector between all points and first point
  vecFromFirst =sweep(allCoord,2,firstPoint,FUN= "-")
  scalarProduct = vecFromFirst %*% lineVecN
  vecFromFirstParallel = do.call("rbind",lapply(scalarProduct,function(x) x*lineVecN))
  vecToLine = vecFromFirst - vecFromFirstParallel
  distToLine = sqrt(rowSums(vecToLine**2))
  distToLine = (distToLine-min(distToLine))/(max(distToLine)-min(distToLine))
  dists=abs(distToLine-max(distToLine))
  
  # this is to get the last component if there is no
  # elbow, where everything is on a straight line
  if( all(is.na(dists)) ){
    return(length(dists))
  }
  
  # this bit is to remove componets that are close
  # to elbow. If their distance to the line
  # is reasonably close to max dist and if they explain more
  cutoff=0.05
  if( any(dists<cutoff & dists != 0) ){
    
    my.order=order(dists)
    return(my.order[my.order %in%  which(dists<cutoff & dists != 0)][1] )
  }
  
  return(which.max(distToLine))
}


#### Function: create a plot (for a given background correction technique) which finds the optimal number of PCs to 
#              minimise the MAE in the controls. 
#
# merged_data: input dataframe with the columns from blood_control_merged_data_benchmark_bc.csv and the calculated PCs. 
# bc_type: background correction ('raw', 'noob').

plot_optimal_strategy <- function(merged_data, bc_type){
  
  ## Initialise variables.
  
  r_F_F <- c();
  r_F_T <- c();
  r_T_F <- c();
  r_T_T <- c();
  n_PCs <- length(grep('PC', colnames(merged_data)));
  
  ## Model fitting with different number of PCs.
  
  for(cc_correction in c(FALSE,TRUE)){
    
    for (batch_correction in c(FALSE,TRUE)){
      
      print(paste(cc_correction, batch_correction));
      temp_mae <- c();
      
      for(i in 0:n_PCs){ # Iterate over number of PCs
        
        print(i);
        temp_data <- merged_data[,1:(i+32)];
        if(cc_correction==FALSE & i == 0){lm_formula <- paste0('DNAmAge_', bc_type, '~Age_years+Sex')};
        if(cc_correction==TRUE & i == 0){lm_formula <- paste0('DNAmAge_', bc_type,'~Age_years+Sex+Gran+CD4T+CD8T+B+Mono+NK')};
        if(cc_correction==FALSE & batch_correction==FALSE & i > 0){lm_formula <- paste0('DNAmAge_', bc_type, '~Age_years+Sex')};
        if(cc_correction==FALSE & batch_correction==TRUE & i > 0){lm_formula <- paste0('DNAmAge_', bc_type, '~Age_years+Sex+', paste(colnames(temp_data)[grep('PC',colnames(temp_data))], collapse='+'))};
        if(cc_correction==TRUE & batch_correction==FALSE & i > 0){lm_formula <- paste0('DNAmAge_', bc_type,'~Age_years+Sex+Gran+CD4T+CD8T+B+Mono+NK')};
        if(cc_correction==TRUE & batch_correction==TRUE & i > 0){lm_formula <- paste0('DNAmAge_', bc_type, '~Age_years+Sex+Gran+CD4T+CD8T+B+Mono+NK+', paste(colnames(temp_data)[grep('PC',colnames(temp_data))], collapse='+'))};
        lm_DNAmAge <- lm(lm_formula, data=temp_data);
        temp_data$delta <- lm_DNAmAge$residuals;
        mae <- median(abs(temp_data$delta)); # Median absolute error
        temp_mae <- c(temp_mae, mae);
  
      }
      
      if(cc_correction==FALSE & batch_correction==FALSE){r_F_F <- temp_mae};
      if(cc_correction==FALSE & batch_correction==TRUE){r_F_T <- temp_mae};
      if(cc_correction==TRUE & batch_correction==FALSE){r_T_F <- temp_mae};
      if(cc_correction==TRUE & batch_correction==TRUE){r_T_T <- temp_mae};
    }
  }
  
  
  ## Calculate optimal number of PCs. 
  
  PC_vector <- 0:n_PCs;
  average_vector <- (r_F_T + r_T_T)/2;
  opt_index <- findElbow(average_vector);
  opt_mae <- average_vector[opt_index];
  opt_PC <- PC_vector[opt_index];

  ## Make the plot.
  
  PC_plot_df <- as.data.frame(cbind(r_F_F, r_F_T, r_T_F, r_T_T));
  PC_plot_df$PCs <- 0:n_PCs;
  colnames(PC_plot_df)[1:4] <- c('CCC: No | Batch: No', 'CCC: No | Batch: Yes',
                                 'CCC: Yes | Batch: No', 'CCC: Yes | Batch: Yes');
  PC_plot_df_gather <- gather(PC_plot_df, key=Corrections, value=MAE, -PCs);
  PC_benchmarking <- ggplot(data=PC_plot_df_gather, aes(x=PCs, y=MAE, col=Corrections)) + 
    geom_line() + theme_classic() +
    theme(axis.text=element_text(size=12, angle=90),
          axis.title=element_text(size=14,face="bold"),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) +
    xlab("Number of PCs") + ylab('MAE in control') +
    labs(title = paste0("Background correction: ", ifelse(bc_type=='raw', 'None', 'noob')),
         subtitle = paste0("Optimal number of PCs: ", opt_PC, '\nOptimal mean MAE: ', round(opt_mae, digits=4))) +
    ylim(c(2.5, 4)) + geom_vline(xintercept=opt_PC, linetype=2);
  
}


#### Function: create a plot which contains the boxplots of the deltas split by batch. In this version of the function,
#              different DNAmAge coming from different background correction techniques can be assessed. 
#
# merged_data: dataframe with the columns from blood_control_merged_data_benchmark_bc.csv and the calculated PCs. 
# bc_type: background correction ('raw', 'noob').
# cc_correction: apply cell composition correction? (TRUE, FALSE).
# batch_correction: apply batch effect correction using the PCs from the control probes? (TRUE, FALSE).

plot_deltas_by_batch_with_bc <- function(merged_data, bc_type, batch_correction, cc_correction){
  
  ## Fit linear model and calculate deltas and MAE. 
  
  if(cc_correction==FALSE & batch_correction==FALSE){lm_formula <- paste0('DNAmAge_', bc_type, '~Age_years+Sex')};
  if(cc_correction==FALSE & batch_correction==TRUE){lm_formula <- paste0('DNAmAge_', bc_type, '~Age_years+Sex+', paste(colnames(merged_data)[grep('PC',colnames(merged_data))], collapse='+'))};
  if(cc_correction==TRUE & batch_correction==FALSE){lm_formula <- paste0('DNAmAge_', bc_type,'~Age_years+Sex+Gran+CD4T+CD8T+B+Mono+NK')};
  if(cc_correction==TRUE & batch_correction==TRUE){lm_formula <- paste0('DNAmAge_', bc_type, '~Age_years+Sex+Gran+CD4T+CD8T+B+Mono+NK+', paste(colnames(merged_data)[grep('PC',colnames(merged_data))], collapse='+'))};
  
  lm_DNAmAge <- lm(lm_formula, data=merged_data);
  merged_data$delta <- lm_DNAmAge$residuals;
  mae <- median(abs(merged_data$delta)); # Median absolute error
  
  ## Make the plot.
  
  delta_by_batch <- ggplot(data=merged_data, aes(x=Batch, y=delta, col=Batch)) +
    geom_boxplot() + scale_colour_stata() +
    theme_classic() +
    theme(axis.text=element_text(size=12, angle=90),
          axis.title=element_text(size=14,face="bold"),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) +
    xlab("Batch") + ylab(ifelse(cc_correction, "EAA with CCC (years)", "EAA without CCC (years)")) +
    labs(subtitle=paste0(#"Background correction: ", ifelse(bc_type=='raw', 'None', 'noob'), 
                         "\nBatch effect correction: ", as.character(batch_correction), 
                         #"\nCell composition correction (CCC): ", as.character(cc_correction),
                         "\nMAE: ", round(mae, digits=4))) +
    ylim(c(-40,40)) + geom_abline(slope=0, intercept=0, linetype=2) + 
    guides(colour=FALSE);
  
}


################################################################
################## Run the pipeline ############################
################################################################

##### 1. Read the raw data. #####.

## Metadata + DNAmAge calculations for blood control. 

merged_data_controls <- as.data.frame(fread('blood_control_merged_data_benchmark_bc.tsv'));
merged_data_controls$Slide_ID <- as.character(merged_data_controls$Slide_ID);

## Control probe intensities for blood control. 

int_files <- list.files(path='../control_intensities_control/', full.names = TRUE);

for(i in 1:length(int_files)){
  if(i==1){
    ints <- fread(int_files[i]); 
  }else{
    ints <- merge(ints, fread(int_files[i]), by=c('ControlAddress', 'ControlChannel', 'ControlType'));
  }
}

ints <- as.data.frame(ints);
diff_sample_names <- which(sapply(strsplit(colnames(ints), '_'), function(x){length(x)<3}));
colnames(ints)[-diff_sample_names] <- sapply(strsplit(colnames(ints)[-diff_sample_names], '_'), function(x){x[1]}); # Edit all names except for Aref and GSE40279 datasets
rownames(ints) <- paste(ints[,1], ints[,2], ints[,3], sep='-');
ints <- ints[,-c(1:3)];
ints <- t(ints);
ints_control <- ints[rownames(ints) %in% merged_data_controls$GEO_sample,]; # Keep only those samples that were selected in control metadata.

## Metadata + DNAmAge calculations for cases.

merged_data_cases <- as.data.frame(fread('cases_merged_data_benchmark_bc.tsv'));
merged_data_cases$Slide_ID <- as.character(merged_data_cases$Slide_ID);

## Control probe intensities for cases. 

int_files <- list.files(path='../../syndromes_screen/control_intensities_cases/', full.names = TRUE);

for(i in 1:length(int_files)){
  if(i==1){
    ints <- fread(int_files[i]); 
  }else{
    ints <- merge(ints, fread(int_files[i]), by=c('ControlAddress', 'ControlChannel', 'ControlType'));
  }
}

ints <- as.data.frame(ints);
diff_sample_names <- which(sapply(strsplit(colnames(ints), '_'), function(x){length(x)<3}));
colnames(ints)[-diff_sample_names] <- sapply(strsplit(colnames(ints)[-diff_sample_names], '_'), function(x){x[1]}); # Edit all names except for Aref, GSE41273 and GSE62298 datasets.
rownames(ints) <- paste(ints[,1], ints[,2], ints[,3], sep='-');
ints <- ints[,-c(1:3)];
ints <- t(ints);
ints_cases <- ints[rownames(ints) %in% merged_data_cases$GEO_sample,]; # Keep only those samples that were selected in cases metadata.

## Merge all the control probe data. 

all_ints <- rbind(ints_control, ints_cases);


##### 2. PCA on control probes intensities. #####

## Perform PCA.

ints_PCA <- prcomp(all_ints);

## Merge with the data for controls + cases. 

all_metadata <- rbind(merged_data_controls, merged_data_cases[,colnames(merged_data_cases) %in% colnames(merged_data_controls)]);
PCs_df <- as.data.frame(cbind(rownames(ints_PCA$x), ints_PCA$x[,1:100]), stringsAsFactors=FALSE); # Restrict analysis to first 100 PCs
colnames(PCs_df)[1] <- 'GEO_sample';
PCs_df[,2:ncol(PCs_df)] <- apply(PCs_df[,2:ncol(PCs_df)], 2, function(x) {as.numeric(x)});
all_metadata_with_PCs <- merge(all_metadata, PCs_df, by='GEO_sample');
control_metadata_with_PCs <- all_metadata_with_PCs[all_metadata_with_PCs$Disease_status=='Control',];

## Find the optimal strategy and number of PCs to minimise the MAE. 

opt_plot_raw <- plot_optimal_strategy(merged_data=control_metadata_with_PCs, bc_type='raw');
ggsave("plots/plot_benchmark_PCs_raw.pdf", height=5, width=6);
opt_plot_noob <- plot_optimal_strategy(merged_data=control_metadata_with_PCs, bc_type='noob');
ggsave("plots/plot_benchmark_PCs_noob.pdf", height=5, width=6);

## We select 'noob' with 17 PCs as the best strategy and we export this for the final screening.

final_control_data <- control_metadata_with_PCs[,c(1:16, 25, 28, 33:49)];
final_control_data <- final_control_data[order(final_control_data$Batch, final_control_data$GEO_sample),]; # Order by batch and GEO sample name
rownames(final_control_data) <- NULL;
write.table(x=final_control_data, file='../../syndromes_screen/final_control_data.tsv', quote=F, 
            sep='\t', row.names=F);

case_metadata_with_PCs <- merge(merged_data_cases, PCs_df, by='GEO_sample');
final_cases_data <- case_metadata_with_PCs[,c(1:23, 32, 35, 40:56)];  
final_cases_data <- final_cases_data[order(final_cases_data$Disease_status, final_cases_data$Batch, 
                                           final_cases_data$GEO_sample),]; # Order by disease status, batch and GEO sample name
rownames(final_cases_data) <- NULL;
write.table(x=final_cases_data, file='../../syndromes_screen/final_cases_data.tsv', quote=F, 
            sep='\t', row.names=F);


##### 3. Some additional plots. #####

## PC1 vs PC2 plots.

var_summary <- summary(ints_PCA)$importance;
PC1_v <- var_summary[2,1]*100;
PC2_v <- var_summary[2,2]*100;

PC1_vs_PC2_controls <- ggplot(data=all_metadata_with_PCs[all_metadata_with_PCs$Disease_status=='Control',], aes(x=PC1, y=PC2, col=Batch)) + 
  geom_point(size=1) + 
  theme_classic() + scale_colour_stata() + 
  theme(axis.text=element_text(size=12, angle=90),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  labs(title = "Controls") +
  #labs(title = paste0("Controls (N=", nrow(all_metadata_with_PCs[all_metadata_with_PCs$Disease_status=='Control',]), ")")) +
  xlab(paste0('PC1 (', round(PC1_v, digits=2), '%)')) + ylab(paste0('PC2 (', round(PC2_v, digits=2), '%)')) +
  xlim(c(-91000, 120000)) + ylim(c(-35000, 40000));
ggsave("plots/plot_PC1_vs_PC2_controls.pdf", height=5, width=6);

PC1_vs_PC2_cases <- ggplot(data=all_metadata_with_PCs[all_metadata_with_PCs$Disease_status!='Control',], aes(x=PC1, y=PC2, col=Batch)) + 
  geom_point(size=1) + 
  theme_classic() + scale_colour_stata() + 
  theme(axis.text=element_text(size=12, angle=90),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  labs(title='Cases') +
  #labs(title = paste0("Cases (N=", nrow(all_metadata_with_PCs[all_metadata_with_PCs$Disease_status!='Control',]), ")")) +
  xlab(paste0('PC1 (', round(PC1_v, digits=2), '%)')) + ylab(paste0('PC2 (', round(PC2_v, digits=2), '%)')) +
  xlim(c(-91000, 60000)) + ylim(c(-55000, 40000));
ggsave("plots/plot_PC1_vs_PC2_cases.pdf", height=5, width=6);

## Plot showing % variance for each PC in PCA. 

pervar_df <- data.frame(PC_int=1:25, pervar=var_summary[2,1:25]*100, 
                        cumper=var_summary[3,1:25]*100);
pervar_df <- gather(pervar_df, key=type, value=perc, -PC_int);

pervar_plot <- ggplot(data=pervar_df, aes(x=PC_int, y=perc, col=type)) + geom_point() + geom_line() +  
  theme_classic() + scale_colour_manual(values=c('blue', 'red'), 
                                        labels=c("Cumulative % variance","% variance for PC")) + 
  theme(axis.text=element_text(size=12, angle=90),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title=element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.position = c(0.81, 0.5)) +
  labs(title = "") +
  xlab('Principal components (PCs)') + ylab('% variance') +
  ylim(c(0, 100)) + geom_vline(xintercept=17, linetype=2);
ggsave("plots/plot_PCA_variance.pdf", height=5, width=5);


# ## Contribution of control probes to first PCs (1-17). 
# 
# loadings_sq <- ints_PCA$rotation[,1:17]^2;
# contributions <- loadings_sq*100; # Equivalent to var$contrib in http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/
# 
# # Filter out those control probes that almost don't contribute to the first 17 PCs. Cutoff: 1% contribution in at least one of the 17 PCs.
# 
# var_rm <- c();
# 
# for(i in 1:nrow(contributions)){
#   if(all(contributions[i,] < 5)){var_rm <- c(var_rm,i)}; 
# }
# 
# contributions_filtered <- contributions[-var_rm,];
# 
# # Plot.
# 
# col_scale <- colorRamp2(c(0, max(contributions_filtered)),c("darkblue","yellow2"));
# pdf('plots/variable_contribution_PCA_heatmap.pdf', width=10, height=8);
# Heatmap(contributions_filtered, col=col_scale, 
#                      row_names_side = "left", row_dend_side = "right",
#                      heatmap_legend_param = list(title = "Variable\ncontribution (%)"),
#                      row_names_gp = gpar(fontsize = 8),
#                      row_title="Control probe variables", row_title_gp = gpar(fontsize = 15, fontface = "bold"),
#                      column_title="Principal Components", column_title_gp = gpar(fontsize = 15, fontface = "bold"), column_title_side='bottom');
# dev.off();


## By batch plots. 

control_metadata_with_17PCs <- control_metadata_with_PCs[,1:49];
control_metadata_with_17PCs <- control_metadata_with_17PCs[order(control_metadata_with_17PCs$Batch, control_metadata_with_17PCs$GEO_sample),]; # Order by batch and GEO sample name

plot_noob_no_no <- plot_deltas_by_batch_with_bc(merged_data=control_metadata_with_17PCs, bc_type='noob', batch_correction=FALSE, cc_correction=FALSE);
ggsave("plots/plot_by_batch_noob_no_no.pdf", height=5, width=8);
plot_noob_no_yes <- plot_deltas_by_batch_with_bc(merged_data=control_metadata_with_17PCs, bc_type='noob', batch_correction=FALSE, cc_correction=TRUE);
ggsave("plots/plot_by_batch_noob_no_yes.pdf", height=5, width=8);
plot_noob_yes_no <- plot_deltas_by_batch_with_bc(merged_data=control_metadata_with_17PCs, bc_type='noob', batch_correction=TRUE, cc_correction=FALSE);
ggsave("plots/plot_by_batch_noob_yes_no_17PCs.pdf", height=5, width=8);
plot_noob_yes_yes <- plot_deltas_by_batch_with_bc(merged_data=control_metadata_with_17PCs, bc_type='noob', batch_correction=TRUE, cc_correction=TRUE);
ggsave("plots/plot_by_batch_noob_yes_yes_17PCs.pdf", height=5, width=8);


## Plot to identify reason for deviations from 0 in median delta per control batch.

lm_formula <- paste0('DNAmAge_noob~Age_years+Sex+Gran+CD4T+CD8T+B+Mono+NK+', paste(colnames(control_metadata_with_17PCs)[grep('PC',colnames(control_metadata_with_17PCs))], collapse='+'));
lm_DNAmAge <- lm(lm_formula, data=control_metadata_with_17PCs);
control_metadata_with_17PCs$delta <- lm_DNAmAge$residuals;

reasons_delta_df <- aggregate(control_metadata_with_17PCs[,c(9,50)], list(control_metadata_with_17PCs$Batch), median);
colnames(reasons_delta_df)[1] <- 'Batch';
reasons_delta_df <- merge(reasons_delta_df,
                          data.frame(Batch=names(table(control_metadata_with_17PCs$Batch)), N=as.numeric(table(control_metadata_with_17PCs$Batch))), by='Batch');

reasons_plot <- ggplot(data=reasons_delta_df, aes(x=Age_years, y=delta, size=N, col=Batch)) + geom_point() +
  scale_size_area(max_size = 10) + scale_colour_stata() +
  theme_classic() +
  theme(axis.text=element_text(size=12, angle=90),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  xlab("Median age (years)") + ylab("Median EAA with CCC (years)") +
  ylim(c(-4.5,4)) + geom_abline(slope=0, intercept=0, linetype=2) + geom_abline(slope=-1/15, intercept=-1, col='grey') +
  guides(size=guide_legend(title="Number of samples in batch"));
ggsave("plots/reasons_deviations_deltas_control_batch_plot.pdf", height=5, width=6);


################################################################
################## End of the script ###########################
################################################################
