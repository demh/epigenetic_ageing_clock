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
##### Calculate aDMPs (differentially methylated positions with age) in our control   #####
##### blood samples, while accounting for covariates.                                 #####
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

input_path <- "/nfs/nobackup/thornton/dem44/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/blood_control/beta_values/"; # Path where all the beta-values matrices are located
output_path <- "/nfs/nobackup/thornton/dem44/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/blood_control/aDMPs/";
metadata_path <- "/nfs/nobackup/thornton/dem44/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/blood_control/aDMPs/controls_data_downstream.tsv";


################################################################
################## Functions ###################################
################################################################

#### Function: calculate aDMPs accounting for the required covariates. ####
# Adapted from minfi::dmpFinder
# In our case with fit the following linear model: Beta~Age_years+Sex+Gran+CD4T+CD8T+B+Mono+NK+PC1+...+PC17

# betas: dataframe with the beta values (rows: probes, columns: samples in the same order as metadata).
# metadata: dataframe with all the metadata (covariates) information (rows:samples, columns:covariates).

calculate_aDMPs <- function(betas, metadata){
  
  pheno <- metadata;
  design <- model.matrix(~Age_years+Sex+Gran+CD4T+CD8T+B+Mono+NK+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17, 
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

final_metadata <- as.data.frame(fread(metadata_path));
print(paste0('There are ', nrow(final_metadata), ' control samples in the original metadata.'));

beta_files <- list.files(input_path, full.names = T);

for(i in 1:length(beta_files)){
  
  batch <- strsplit(strsplit(beta_files[i], '/')[[1]][length(strsplit(beta_files[i], '/')[[1]])], '_')[[1]][1];
  print(paste0('Processing batch ', batch, ' ...'));
  temp_betas <- as.data.frame(fread(beta_files[i]));
  rownames(temp_betas) <- temp_betas$ProbeID;
  cgs <- rownames(temp_betas);
  temp_betas <- temp_betas[,-1];
  if(!(batch %in% c('GSE40279', 'GSE41273', 'Europe'))){colnames(temp_betas) <- sapply(strsplit(colnames(temp_betas), '_'), function(x){x[1]})};
  temp_betas <- temp_betas[,which(colnames(temp_betas) %in% final_metadata$GEO_sample)]; # Select only those samples that are in metadata (and therefore passed QC)
  
  if(i==1){
    
    final_betas <- temp_betas;
    
  }else{
    
    if(!all(cgs==rownames(final_betas))){stop('The probe IDs in this matrix is not in the same order as the previous ones.')};
    final_betas <- cbind(final_betas, temp_betas);
  }
}

print(paste0('There are ', ncol(final_betas), ' control samples in the final beta-values matrix.'));

# Edit metadata so it matches the samples we have. 

final_metadata <- final_metadata[match(colnames(final_betas), final_metadata$GEO_sample),];


#### 2. Calculate the aDMPs. ####

print('Calculating aDMPs ...');
dmps <- calculate_aDMPs(betas=final_betas, metadata=final_metadata);
print('aDMP calculation finished.');
dmps <- cbind(rownames(dmps), dmps);
colnames(dmps)[1] <- 'ProbeID'; rownames(dmps) <- NULL;
print('Exporting aDMPs ...');
write.table(x=dmps, file=paste0(output_path, '/aDMPs_final.csv'), quote=F, sep=',', row.names=F);

print('The script finished correctly.');

#### Extra. Create file with the coordinates for the 450K probes that survived preprocessing (428266). ####

# all_450k_probes <- readLines('probes_450K_after_preprocessing.txt');
# 
# data("IlluminaHumanMethylation450kanno.ilmn12.hg19");
# data("Locations");
# hg19_coords <- Locations;
# all_450k_hg19_coords <- as.data.frame(hg19_coords[which(rownames(hg19_coords) %in% all_450k_probes),]);
# all_450k_hg19_coords$CpGmarker <- rownames(all_450k_hg19_coords);
# rownames(all_450k_hg19_coords) <- NULL;
# write.table(x=all_450k_hg19_coords[,c(4,1,2)], file='all_450k_CpGs_coords_hg19.csv', quote=F,
#             sep=',', row.names=F);

################################################################
################## End of the script ###########################
################################################################
