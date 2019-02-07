###########################################################################################
#########                                                                         #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             01/05/2018                                  #########
#########                              EMBL-EBI                                   #########
#########                           Thornton group                                #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####              Biological insights into the epigenetic ageing clock           #########
###########################################################################################
##### Edit the names of the samples in GSE41273 to allow merging.                     #####
###########################################################################################
##### USAGE: manual                                                                   #####
###########################################################################################

library(data.table);
setwd('/Users/dem44/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/blood_control/');

#### 1. Obtain metadata with all the different names. ####

GSE41273_raw <- as.data.frame(fread('metadata_files/metadata_GSE41273.tsv', header=T, sep='\t'));
GSE41273_ages <- as.data.frame(fread('metadata_files/FXS_samples.csv', header=T, sep=','));
GSE41273_files <- readLines('metadata_files/GSE41273_file_names.txt');
GSE41273_bridge <- data.frame(sample_number=sapply(strsplit(sapply(strsplit(GSE41273_files, '/'), function(x){x[3]}), '_'), function(x){x[1]}),
                              ID450k=sapply(strsplit(sapply(strsplit(GSE41273_files, '/'), function(x){x[3]}), '_'), function(x){paste0(x[2],'_',x[3])}));
GSE41273_bridge <- GSE41273_bridge[!duplicated(GSE41273_bridge),];
GSE41273_m1 <- merge(GSE41273_ages, GSE41273_bridge, by='ID450k');
colnames(GSE41273_m1)[2] <- 'title';
GSE41273_m2 <- merge(GSE41273_raw, GSE41273_m1, by='title');
GSE41273_equiv <- data.frame(SampleID=GSE41273_m2$sample_number, GEO=GSE41273_m2$geo_accession);

#### 2. Edit names in DNAmAge file (raw). #### 

GSE41273_DNAmAge <- fread('~/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/old_stuff/GSE41273/DNAm_age_of_DNAmAge_methylation_matrix_from_idat_GSE41273.csv');
GSE41273_DNAmAge_final <- merge(GSE41273_DNAmAge, GSE41273_equiv, by='SampleID');
GSE41273_DNAmAge_final$SampleID <- GSE41273_DNAmAge_final$GEO;
GSE41273_DNAmAge_final <- GSE41273_DNAmAge_final[,-10];
write.table(GSE41273_DNAmAge_final, 'DNAmAge/DNAm_age_of_DNAmAge_methylation_matrix_from_idat_GSE41273_edited.csv', sep=',', row.names=F, quote=F);

#### 3. Edit names in DNAmAge file (with NOOB). #### 

GSE41273_DNAmAge <- fread('~/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/old_stuff/GSE41273/DNAm_age_of_DNAmAge_methylation_matrix_from_idat_GSE41273_noob.csv');
GSE41273_DNAmAge_final <- merge(GSE41273_DNAmAge, GSE41273_equiv, by='SampleID');
GSE41273_DNAmAge_final$SampleID <- GSE41273_DNAmAge_final$GEO;
GSE41273_DNAmAge_final <- GSE41273_DNAmAge_final[,-10];
write.table(GSE41273_DNAmAge_final, 'DNAmAge_noob/DNAm_age_of_DNAmAge_methylation_matrix_from_idat_GSE41273_edited.csv', sep=',', row.names=F, quote=F);

#### 4. Edit names in cc file. ####

GSE41273_cc <- fread('~/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/old_stuff/GSE41273/GSE41273_cc_predictions.csv');
GSE41273_cc$SampleID <- sapply(strsplit(GSE41273_cc$SampleID, '_'), function(x){x[1]});
GSE41273_cc_final <- merge(GSE41273_cc, GSE41273_equiv, by='SampleID');
GSE41273_cc_final$SampleID <- GSE41273_cc_final$GEO;
GSE41273_cc_final <- GSE41273_cc_final[,-8];
write.table(GSE41273_cc_final, 'cc_files/GSE41273_cc_predictions_edited.csv', sep=',', row.names=F, quote=F);

#### 5. Edit names in QC matrix. ####

GSE41273_QC <- fread('~/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/old_stuff/GSE41273/GSE41273_QC_matrix.csv');
GSE41273_QC_final <- merge(GSE41273_QC, GSE41273_equiv, by='SampleID');
GSE41273_QC_final$SampleID <- GSE41273_QC_final$GEO;
GSE41273_QC_final <- GSE41273_QC_final[,-8];
write.table(GSE41273_QC_final, 'QC_files/GSE41273_QC_matrix_edited.csv', sep=',', row.names=F, quote=F);

#### 6. Edit names in control probe intensities file. ####

GSE41273_control <- as.data.frame(fread('~/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/old_stuff/GSE41273/GSE41273_control_intensities.csv'));
GSE41273_control_old_sample_names <- sapply(strsplit(colnames(GSE41273_control)[-c(1:3)], '_'), function(x){x[1]});
rm_samples <- c(which(!GSE41273_control_old_sample_names %in% GSE41273_equiv$SampleID)+3);
GSE41273_control_with_age <- GSE41273_control[,-rm_samples]; # Remove those samples (columns) for which we have no ages
GSE41273_control_old_sample_names <- sapply(strsplit(colnames(GSE41273_control_with_age)[-c(1:3)], '_'), function(x){x[1]});
GSE41273_control_new_sample_names <- GSE41273_equiv$GEO[match(GSE41273_control_old_sample_names, as.character(GSE41273_equiv$SampleID))];
colnames(GSE41273_control_with_age)[4:length(colnames(GSE41273_control_with_age))] <- as.character(GSE41273_control_new_sample_names);
write.table(GSE41273_control_with_age, 'control_intensities_control/GSE41273_control_intensities_edited.csv', sep=',', row.names=F, quote=F);

#### 7. Edit names in horvath clock sites file. ####

GSE41273_sites <- as.data.frame(fread('~/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/old_stuff/GSE41273/GSE41273_horvath_clock_beta_values.csv'));
GSE41273_sites_old_sample_names <- sapply(strsplit(colnames(GSE41273_sites)[-1], '_'), function(x){x[1]});
rm_samples <- c(which(!GSE41273_sites_old_sample_names %in% GSE41273_equiv$SampleID)+1);
GSE41273_sites_with_age <- GSE41273_sites[,-rm_samples]; # Remove those samples (columns) for which we have no ages
GSE41273_sites_old_sample_names <- sapply(strsplit(colnames(GSE41273_sites_with_age)[-1], '_'), function(x){x[1]});
GSE41273_sites_new_sample_names <- GSE41273_equiv$GEO[match(GSE41273_sites_old_sample_names, as.character(GSE41273_equiv$SampleID))];
colnames(GSE41273_sites_with_age)[2:length(colnames(GSE41273_sites_with_age))] <- as.character(GSE41273_sites_new_sample_names);
write.table(GSE41273_sites_with_age, 'horvath_clock_sites/data/GSE41273_horvath_clock_beta_values_edited.csv', sep=',', row.names=F, quote=F);

#### 8. Edit the names in the beta values file. ####

GSE41273_sites <- as.data.frame(fread('/Users/dem44/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/differential_methylation/raw_data_control/GSE41273_betas_matrix.csv'));
GSE41273_sites_old_sample_names <- sapply(strsplit(colnames(GSE41273_sites)[-1], '_'), function(x){x[1]});
rm_samples <- c(which(!GSE41273_sites_old_sample_names %in% GSE41273_equiv$SampleID)+1);
GSE41273_sites_with_age <- GSE41273_sites[,-rm_samples]; # Remove those samples (columns) for which we have no ages
GSE41273_sites_old_sample_names <- sapply(strsplit(colnames(GSE41273_sites_with_age)[-1], '_'), function(x){x[1]});
GSE41273_sites_new_sample_names <- GSE41273_equiv$GEO[match(GSE41273_sites_old_sample_names, as.character(GSE41273_equiv$SampleID))];
colnames(GSE41273_sites_with_age)[2:length(colnames(GSE41273_sites_with_age))] <- as.character(GSE41273_sites_new_sample_names);
write.table(GSE41273_sites_with_age, '/Users/dem44/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/differential_methylation/raw_data_control/GSE41273_betas_matrix_edited.csv', sep=',', row.names=F, quote=F);


#### End of the script. ####
