###########################################################################################
#########                                                                         #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             07/12/2018                                  #########
#########                              EMBL-EBI                                   #########
#########                           Thornton group                                #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####              Biological insights into the epigenetic ageing clock           #########
###########################################################################################
##### Calculate pcgtAge (mitotic clock age) from the beta-values of our samples.       ####
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
output_path <- "/nfs/nobackup/thornton/dem44/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/blood_control/pcgtAge/";
probes_path <- "/nfs/research1/thornton/dem44/methylation_clock/raw_data_human/epigenetic_syndromes/utils/pcgt_Age_probes.txt";


############################################################
################## Running the pipeline ####################
############################################################

#### 1. Read all the input data. ####

pcgtAge_probes <- readLines(probes_path);
beta_files <- list.files(input_path, full.names = T);

for(i in 1:length(beta_files)){
  
  batch <- strsplit(strsplit(beta_files[i], '/')[[1]][length(strsplit(beta_files[i], '/')[[1]])], '_')[[1]][1];
  print(paste0('Processing batch ', batch, ' ...'));
  temp_betas <- as.data.frame(fread(beta_files[i]));
  rownames(temp_betas) <- temp_betas$ProbeID;
  cgs <- rownames(temp_betas);
  temp_betas <- temp_betas[,-1];
  if(!(batch %in% c('GSE40279', 'GSE41273', 'Europe'))){colnames(temp_betas) <- sapply(strsplit(colnames(temp_betas), '_'), function(x){x[1]})};
  temp_betas <- temp_betas[rownames(temp_betas) %in% pcgtAge_probes,];
  if(nrow(temp_betas) == 378){
    print(paste0('Number of pcgtAge probes retrieved: ', nrow(temp_betas)));
  }else{
    stop('The number of pcgtAge is not 378.');
  }

  if(i==1){
    
    final_results <- data.frame(GEO_sample=colnames(temp_betas),
                                pcgtAge=as.numeric(apply(temp_betas,2,mean)));
    
  }else{
    
    temp_results <- data.frame(GEO_sample=colnames(temp_betas),
                               pcgtAge=as.numeric(apply(temp_betas,2,mean)));
    final_results <- rbind(final_results, temp_results);
  }
}


#### 2. Export the results. ####

print('Exporting the results ...');
write.table(x=final_results, file=paste0(output_path, '/pcgtAge_results.csv'), quote=F, sep=',', row.names=F);

print('The script finished correctly.');

#########################################################
########## End of the script ############################
#########################################################
