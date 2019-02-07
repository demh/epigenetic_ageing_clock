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
##### Calculate Shannon entropy from the beta-values for each of our samples.          ####
###########################################################################################
##### USAGE: manual                                                                    ####
###########################################################################################

###########################################################
##################### Dependencies ########################
###########################################################

library(data.table);


###########################################################
#####################  Arguments ##########################
###########################################################

input_path <- "/nfs/nobackup/thornton/dem44/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/blood_control/beta_values/"; # Path where all the beta-values matrices are located
output_path <- "/nfs/nobackup/thornton/dem44/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/blood_control/entropy/";


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

#### 1. Read all the input data. ####

beta_files <- list.files(input_path, full.names = T);

for(i in 1:length(beta_files)){
  
  batch <- strsplit(strsplit(beta_files[i], '/')[[1]][length(strsplit(beta_files[i], '/')[[1]])], '_')[[1]][1];
  print(paste0('Processing batch ', batch, ' ...'));
  temp_betas <- as.data.frame(fread(beta_files[i]));
  temp_betas <- temp_betas[,-1];
  if(!(batch %in% c('GSE40279', 'GSE41273', 'Europe'))){colnames(temp_betas) <- sapply(strsplit(colnames(temp_betas), '_'), function(x){x[1]})};
  
  if(i==1){
    
    final_results <- data.frame(GEO_sample=colnames(temp_betas),
                                entropy=as.numeric(apply(temp_betas,2,calculate_entropy)));
    
  }else{
    
    temp_results <- data.frame(GEO_sample=colnames(temp_betas),
                               entropy=as.numeric(apply(temp_betas,2,calculate_entropy)));
    final_results <- rbind(final_results, temp_results);
  }
}


#### 2. Export the results. ####

print('Exporting the results ...');
write.table(x=final_results, file=paste0(output_path, '/entropy_results.csv'), quote=F, sep=',', row.names=F);

print('The script finished correctly.');

#########################################################
########## End of the script ############################
#########################################################
