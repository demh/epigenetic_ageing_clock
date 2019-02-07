###########################################################################################
#########                                                                         #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             14/12/2018                                  #########
#########                              EMBL-EBI                                   #########
#########                           Thornton group                                #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####              Biological insights into the epigenetic ageing clock           #########
###########################################################################################
##### Calculate Sotos DMPs (differentially methylated positions in Sotos samples),    #####
##### while accounting for covariates.                                                #####
###########################################################################################
##### USAGE: manual                                                                   #####
###########################################################################################

###########################################################
##################### Dependencies ########################
###########################################################

library(data.table);
library(minfi);
library(limma);
library(siggenes);


################################################################
################## Arguments ###################################
################################################################

input_path <- "/Users/dem44/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/differential_methylation/raw_data_control/GSE74432_betas_matrix.csv"; # Path where all the beta-values matrices are located
output_path <- "/Users/dem44/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/differential_methylation/";
metadata_control_path <- "/Users/dem44/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/syndromes_screen/controls_data_downstream.tsv";
metadata_sotos_path <- "/Users/dem44/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/syndromes_screen/cases_data_downstream.tsv";


################################################################
################## Functions ###################################
################################################################

#### Function: calculate Sotos DMPs accounting for the required covariates. ####
# Adapted from minfi::dmpFinder
# In our case with fit the following linear model: Beta~Disease_status+Age_years+Sex+Gran+CD4T+CD8T+B+Mono+NK+PC1+...+PC17

# betas: dataframe with the beta values (rows: probes, columns: samples in the same order as metadata).
# metadata: dataframe with all the metadata (covariates) information (rows:samples, columns:covariates).

calculate_Sotos_DMPs <- function(betas, metadata){
  
  pheno <- metadata;
  design <- model.matrix(~Disease_status+Age_years+Sex+Gran+CD4T+CD8T+B+Mono+NK+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17, 
                         data=pheno);
  print('Fitting linear models ...');
  fit <- lmFit(betas, design);
  sigma <- fit$sigma;
  df <- fit$df.residual;
  coef <- fit$coefficients;
  stdev <- fit$stdev.unscaled;
  t <- coef[, 2]/(stdev[, 2] * sigma);
  pval <- 2 * pt(abs(t), df = df, lower.tail = FALSE);
  tab <- data.frame(intercept = coef[, 1], beta = coef[, 2], t = t, pval = pval);
  p0 <- pi0.est(tab$pval[!is.na(tab$pval)])$p0;
  tab$qval <- qvalue.cal(tab$pval, p0);
  o <- order(tab$pval);
  tab <- tab[o, ];
  return(tab);
  
}


################################################################
################## Run the pipeline ############################
################################################################

#### 1. Read all the input data. ####

final_metadata_control <- as.data.frame(fread(metadata_control_path));
final_metadata_control <- final_metadata_control[final_metadata_control$Batch=='GSE74432',];
final_metadata_sotos <- as.data.frame(fread(metadata_sotos_path));
final_metadata_sotos <- final_metadata_sotos[final_metadata_sotos$Disease_status=='Sotos',which(colnames(final_metadata_sotos) %in% colnames(final_metadata_control))];
final_metadata <- rbind(final_metadata_control, final_metadata_sotos);
print(paste0('There are ', nrow(final_metadata), ' samples in the original metadata.'));

temp_betas <- as.data.frame(fread(input_path));
rownames(temp_betas) <- temp_betas$ProbeID;
temp_betas <- temp_betas[,-1];
colnames(temp_betas) <- sapply(strsplit(colnames(temp_betas), '_'), function(x){x[1]});
final_betas <- temp_betas[,which(colnames(temp_betas) %in% final_metadata$GEO_sample)]; # Select only those samples that are in metadata (and therefore passed QC)
print(paste0('There are ', ncol(final_betas), ' samples in the final beta-values matrix.'));

final_metadata <- final_metadata[match(colnames(final_betas), final_metadata$GEO_sample),];


#### 2. Calculate the Sotos DMPs. ####

print('Calculating Sotos DMPs ...');
dmps <- calculate_Sotos_DMPs(betas=final_betas, metadata=final_metadata);
print('Sotos DMP calculation finished.');
dmps <- cbind(rownames(dmps), dmps);
colnames(dmps)[1] <- 'ProbeID'; rownames(dmps) <- NULL;
print('Exporting Sotos DMPs ...');
write.table(x=dmps, file=paste0(output_path, '/Sotos_DMPs_final.csv'), quote=F, sep=',', row.names=F);

print('The script finished correctly.');

################################################################
################## End of the script ###########################
################################################################
