###########################################################################################
#########                                                                         #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             17/08/2018                                  #########
#########                              EMBL-EBI                                   #########
#########                           Thornton group                                #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####              Biological insights into the epigenetic ageing clock           #########
###########################################################################################
##### Build blood control using DNA methylation data from a healthy population of     #####
##### individuals from different studies.                                             #####
##### Do it with DNAmAge calculations that come from data without ('raw') and with    #####
##### ('noob') background correction.                                                 #####
###########################################################################################
##### USAGE: manual                                                                   #####
###########################################################################################

###########################################################
##################### Dependencies ########################
###########################################################

library(data.table);
library(ggplot2);

setwd('/Users/dem44/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/blood_control/');

################################################################
################## Run the pipeline ############################
################################################################

##### 1. Read and curate the metadata information. #####

col_names <- c('Sample_name','GEO_sample','Batch','Platform','Tissue','Slide_ID','Array_ID','Gender','Age_years'); # Gender has now been changed for (biological) sex.

## GSE59065.

GSE59065_raw <- as.data.frame(fread('metadata_files/metadata_GSE59065.tsv', header=T, sep='\t'));
GSE59065_metadata <- GSE59065_raw[GSE59065_raw$`cell/tissue type:ch1` == "Peripheral blood", # Select only peripheral blood samples
                                  c('title', 'geo_accession', 'platform_id', 'cell/tissue type:ch1', 'slide_id:ch1', 'array_id:ch1', 'gender:ch1', 'age (yr):ch1')];
GSE59065_metadata <- cbind(GSE59065_metadata[,c(1,2)], 'GSE59065', GSE59065_metadata[,c(3,4,5,6,7)],GSE59065_metadata[,8]);
colnames(GSE59065_metadata) <- col_names;
GSE59065_metadata$Gender <- ifelse(GSE59065_metadata$Gender=='male', 'Male', 'Female');
GSE59065_metadata$Batch <- as.character(GSE59065_metadata$Batch);
GSE59065_metadata$Slide_ID <- as.character(GSE59065_metadata$Slide_ID);
GSE59065_metadata$Age_years <- as.numeric(GSE59065_metadata$Age_years);
rm(GSE59065_raw);

## GSE111629.

GSE111629_raw <- as.data.frame(fread('metadata_files/metadata_GSE111629.tsv', header=T, sep='\t'));
GSE111629_metadata <- GSE111629_raw[GSE111629_raw$`disease state:ch1`=='PD-free control', # Select only control samples
                                    c('title', 'geo_accession', 'platform_id', 'tissue:ch1', 'gender:ch1', 'age:ch1')];
GSE111629_metadata <- cbind(GSE111629_metadata[,c(1,2)], 'GSE111629', GSE111629_metadata[,c(3,4)],
                            sapply(sapply(strsplit(as.character(GSE111629_metadata[,1]), ' '), function(x){strsplit(x[4], '_')}), function(x){x[1]}),
                            sapply(sapply(strsplit(as.character(GSE111629_metadata[,1]), ' '), function(x){strsplit(x[4], '_')}), function(x){x[2]}),
                            GSE111629_metadata[,c(5,6)]);
colnames(GSE111629_metadata) <- col_names;
GSE111629_metadata$Tissue <- "Whole blood";
GSE111629_metadata$Batch <- as.character(GSE111629_metadata$Batch);
GSE111629_metadata$Slide_ID <- as.character(GSE111629_metadata$Slide_ID);
GSE111629_metadata$Array_ID <- as.character(GSE111629_metadata$Array_ID);
GSE111629_metadata$Age_years <- as.numeric(GSE111629_metadata$Age_years);
rm(GSE111629_raw);

## GSE104812.

GSE104812_raw <- as.data.frame(fread('metadata_files/metadata_GSE104812.tsv', header=T, sep='\t'));
GSE104812_metadata <- GSE104812_raw[, c('title', 'geo_accession', 'platform_id', 'gender:ch1', 'age (y):ch1')];
GSE104812_metadata <- cbind(GSE104812_metadata[,c(1,2)], 'GSE104812', GSE104812_metadata[,3], 'Whole blood', NA, NA, GSE104812_metadata[,4], GSE104812_metadata[,5]);
colnames(GSE104812_metadata) <- col_names;
GSE104812_metadata$Batch <- as.character(GSE104812_metadata$Batch);
GSE104812_metadata$Platform <- as.character(GSE104812_metadata$Platform);
GSE104812_metadata$Tissue <- as.character(GSE104812_metadata$Tissue);
GSE104812_metadata$Slide_ID <- as.character(GSE104812_metadata$Slide_ID);
GSE104812_metadata$Array_ID <- as.character(GSE104812_metadata$Array_ID);
GSE104812_metadata$Gender <- as.character(GSE104812_metadata$Gender);
rm(GSE104812_raw);

## GSE81961.

GSE81961_raw <- as.data.frame(fread('metadata_files/metadata_GSE81961.tsv', header=T, sep='\t'));
GSE81961_metadata <- GSE81961_raw[GSE81961_raw$`characteristics_ch1.1`=='disease state: Control', # Select only control samples
                                  c('title', 'geo_accession', 'platform_id', 'tissue:ch1', 'gender:ch1', 'age (yr):ch1' ,'supplementary_file')];
GSE81961_metadata <- cbind(GSE81961_metadata[,c(1,2)], 'GSE81961', GSE81961_metadata[,c(3,4)],
                           sapply(strsplit(GSE81961_metadata$supplementary_file, '_'), function(x){x[2]}),
                           sapply(strsplit(GSE81961_metadata$supplementary_file, '_'), function(x){x[3]}),
                           GSE81961_metadata[,5], GSE81961_metadata[,6]);
colnames(GSE81961_metadata) <- col_names;
GSE81961_metadata$Tissue <- "Peripheral blood";
GSE81961_metadata$Batch <- as.character(GSE81961_metadata$Batch);
GSE81961_metadata$Slide_ID <- as.character(GSE81961_metadata$Slide_ID);
GSE81961_metadata$Array_ID <- as.character(GSE81961_metadata$Array_ID);
GSE81961_metadata$Gender <- as.character(GSE81961_metadata$Gender);
rm(GSE81961_raw);

## GSE61496.

GSE61496_raw <- as.data.frame(fread('metadata_files/metadata_GSE61496.tsv', header=T, sep='\t'));
GSE61496_metadata <- GSE61496_raw[!duplicated(GSE61496_raw$`pair id:ch1`), # Select only one member of each twins pair
                                  c('title', 'geo_accession', 'platform_id', 'tissue:ch1', 'sex, 1=m, 2=f:ch1', 'age:ch1', 'supplementary_file')];
GSE61496_metadata <- GSE61496_metadata[!is.na(GSE61496_metadata$`age:ch1`),]; # Remove samples with NA in age and gender
GSE61496_metadata <- cbind(GSE61496_metadata[,c(1,2)], 'GSE61496', GSE61496_metadata[,c(3,4)], 
                           sapply(strsplit(GSE61496_metadata$supplementary_file, '_'), function(x){x[2]}),
                           sapply(strsplit(GSE61496_metadata$supplementary_file, '_'), function(x){x[3]}),
                           ifelse(GSE61496_metadata[,5]==1, 'Male', 'Female'), GSE61496_metadata[,6]);
colnames(GSE61496_metadata) <- col_names;
GSE61496_metadata$Tissue <- 'Whole blood';
GSE61496_metadata$Batch <- as.character(GSE61496_metadata$Batch);
GSE61496_metadata$Slide_ID <- as.character(GSE61496_metadata$Slide_ID);
GSE61496_metadata$Array_ID <- as.character(GSE61496_metadata$Array_ID);
GSE61496_metadata$Gender <- as.character(GSE61496_metadata$Gender);
GSE61496_metadata$Age_years <- as.numeric(GSE61496_metadata$Age_years);
rm(GSE61496_raw);

## Aref-Eshghi dataset.

AE_raw <- as.data.frame(fread('metadata_files/results_Erfan.csv', header=T, sep=','));
AE_disease <- as.data.frame(fread('metadata_files/Aref_Eshghi_disease_metadata_curated.csv'));
AE_all <- merge(AE_raw, AE_disease, by='SampleID');
AE_controls <- AE_all[AE_all$Disease_status=='Control',]; # Select only controls
AE_controls <- AE_controls[!is.na(AE_controls$age),]; # Remove samples without ages
AE_metadata <- as.data.frame(cbind(AE_controls$id,
                     AE_controls$SampleID,
                     AE_controls$batch,
                     'GPL13534',
                     'Peripheral blood',
                     sapply(strsplit(AE_controls$idat, '_'), function(x){x[1]}),
                     sapply(strsplit(AE_controls$idat, '_'), function(x){x[2]}),
                     ifelse(AE_controls$sex=='f', 'Female', 'Male'),
                     AE_controls$age));
colnames(AE_metadata) <- col_names;
rm(AE_raw, AE_disease, AE_all, AE_controls);
AE_metadata$Sample_name <- as.character(AE_metadata$Sample_name);
AE_metadata$GEO_sample <- as.character(AE_metadata$GEO_sample);
AE_metadata$Batch <- as.character(AE_metadata$Batch);
AE_metadata$Platform <- as.character(AE_metadata$Platform);
AE_metadata$Tissue <- as.character(AE_metadata$Tissue);
AE_metadata$Slide_ID <- as.character(AE_metadata$Slide_ID);
AE_metadata$Array_ID <- as.character(AE_metadata$Array_ID);
AE_metadata$Gender <- as.character(AE_metadata$Gender);
AE_metadata$Age_years <- as.numeric(as.character(AE_metadata$Age_years));

## GSE55491.

GSE55491_raw <- as.data.frame(fread('metadata_files/metadata_GSE55491.tsv', header=T, sep='\t'));
GSE55491_add <- as.data.frame(fread('metadata_files/samples4Daniel.csv', header=T, sep=','));
colnames(GSE55491_add)[1] <- 'source_name_ch1';
GSE55491_raw <- merge(GSE55491_raw, GSE55491_add, by='source_name_ch1'); # Add ages to metadata
GSE55491_raw <- GSE55491_raw[GSE55491_raw$SampleGroup=='Ctrl', ]; # Select only controls
GSE55491_metadata <- as.data.frame(cbind(GSE55491_raw$title, GSE55491_raw$geo_accession, 'GSE55491', GSE55491_raw$platform_id,
                                         'Peripheral blood', as.character(GSE55491_raw$Sentrix_Barcode), GSE55491_raw$Sentrix_Position, 
                                         ifelse(GSE55491_raw$Sex=='M', 'Male', 'Female'), as.character(GSE55491_raw$Age)));
rm(GSE55491_raw, GSE55491_add);
colnames(GSE55491_metadata) <- col_names;
GSE55491_metadata$Sample_name <- as.character(GSE55491_metadata$Sample_name);
GSE55491_metadata$GEO_sample <- as.character(GSE55491_metadata$GEO_sample);
GSE55491_metadata$Batch <- as.character(GSE55491_metadata$Batch);
GSE55491_metadata$Platform <- as.character(GSE55491_metadata$Platform);
GSE55491_metadata$Tissue <- as.character(GSE55491_metadata$Tissue);
GSE55491_metadata$Slide_ID <- as.character(GSE55491_metadata$Slide_ID);
GSE55491_metadata$Array_ID <- as.character(GSE55491_metadata$Array_ID);
GSE55491_metadata$Gender <- as.character(GSE55491_metadata$Gender);
GSE55491_metadata$Age_years <- as.numeric(as.character(GSE55491_metadata$Age_years));

## GSE97362.

GSE97362_raw <- as.data.frame(fread('metadata_files/metadata_GSE97362.tsv', header=T, sep='\t'));
GSE97362_raw <- GSE97362_raw[GSE97362_raw$characteristics_ch1.3=="disease state: Control",]; # Select only controls
GSE97362_metadata <- as.data.frame(cbind(GSE97362_raw$title, GSE97362_raw$geo_accession, 'GSE97362', GSE97362_raw$platform_id,'Whole blood',
                                         sapply(strsplit(sapply(strsplit(GSE97362_raw$supplementary_file, '/'), function(x){x[9]}), '_'), function(x){x[2]}),
                                         sapply(strsplit(sapply(strsplit(GSE97362_raw$supplementary_file, '/'), function(x){x[9]}), '_'), function(x){x[3]}),
                                         ifelse(GSE97362_raw$characteristics_ch1=='gender: female', 'Female', 'Male'),
                                         sapply(strsplit(GSE97362_raw$characteristics_ch1.1, ': '), function(x){x[2]})));
colnames(GSE97362_metadata) <- col_names;
GSE97362_metadata <- GSE97362_metadata[!(GSE97362_metadata$Age_years=='-'),]; # Remove samples without ages
rm(GSE97362_raw);
GSE97362_metadata$Sample_name <- as.character(GSE97362_metadata$Sample_name);
GSE97362_metadata$GEO_sample <- as.character(GSE97362_metadata$GEO_sample);
GSE97362_metadata$Batch <- as.character(GSE97362_metadata$Batch);
GSE97362_metadata$Platform <- as.character(GSE97362_metadata$Platform);
GSE97362_metadata$Tissue <- as.character(GSE97362_metadata$Tissue);
GSE97362_metadata$Slide_ID <- as.character(GSE97362_metadata$Slide_ID);
GSE97362_metadata$Array_ID <- as.character(GSE97362_metadata$Array_ID);
GSE97362_metadata$Gender <- as.character(GSE97362_metadata$Gender);
GSE97362_metadata$Age_years <- as.numeric(as.character(GSE97362_metadata$Age_years));


## GSE41273.

GSE41273_raw <- as.data.frame(fread('metadata_files/metadata_GSE41273.tsv', header=T, sep='\t'));
GSE41273_ages <- as.data.frame(fread('metadata_files/FXS_samples.csv', header=T, sep=','));
GSE41273_files <- readLines('metadata_files/GSE41273_file_names.txt');
GSE41273_bridge <- data.frame(sample_number=sapply(strsplit(sapply(strsplit(GSE41273_files, '/'), function(x){x[3]}), '_'), function(x){x[1]}),
                              ID450k=sapply(strsplit(sapply(strsplit(GSE41273_files, '/'), function(x){x[3]}), '_'), function(x){paste0(x[2],'_',x[3])}));
GSE41273_bridge <- GSE41273_bridge[!duplicated(GSE41273_bridge),];
GSE41273_m1 <- merge(GSE41273_ages, GSE41273_bridge, by='ID450k');
colnames(GSE41273_m1)[2] <- 'title';
GSE41273_m2 <- merge(GSE41273_raw, GSE41273_m1, by='title');
GSE41273_m2 <- GSE41273_m2[GSE41273_m2$STATUS=='Control',]; # Select only controls
GSE41273_metadata <- as.data.frame(cbind(as.character(GSE41273_m2$title), GSE41273_m2$geo_accession, 'GSE41273', GSE41273_m2$platform_id, 'Peripheral blood',
                                         sapply(strsplit(GSE41273_m2$ID450k, '_'), function(x){x[1]}), sapply(strsplit(GSE41273_m2$ID450k, '_'), function(x){x[2]}),
                                         ifelse(GSE41273_m2$SEX=='male', 'Male', 'Female'), GSE41273_m2$Ados_Age_months/12));
colnames(GSE41273_metadata) <- col_names;
rm(GSE41273_raw, GSE41273_ages, GSE41273_files, GSE41273_bridge, GSE41273_m1, GSE41273_m2);
GSE41273_metadata$Sample_name <- as.character(GSE41273_metadata$Sample_name);
GSE41273_metadata$GEO_sample <- as.character(GSE41273_metadata$GEO_sample);
GSE41273_metadata$Batch <- as.character(GSE41273_metadata$Batch);
GSE41273_metadata$Platform <- as.character(GSE41273_metadata$Platform);
GSE41273_metadata$Tissue <- as.character(GSE41273_metadata$Tissue);
GSE41273_metadata$Slide_ID <- as.character(GSE41273_metadata$Slide_ID);
GSE41273_metadata$Array_ID <- as.character(GSE41273_metadata$Array_ID);
GSE41273_metadata$Gender <- as.character(GSE41273_metadata$Gender);
GSE41273_metadata$Age_years <- as.numeric(as.character(GSE41273_metadata$Age_years));

## GSE42861.

GSE42861_raw <- as.data.frame(fread('metadata_files/metadata_GSE42861.tsv', header=T, sep='\t'));
GSE42861_raw <- GSE42861_raw[GSE42861_raw$`disease state:ch1`=="Normal",]; # Select only controls
GSE42861_metadata <- as.data.frame(cbind(GSE42861_raw$title, GSE42861_raw$geo_accession, 'GSE42861', GSE42861_raw$platform_id,'Peripheral blood',
                                         sapply(strsplit(sapply(strsplit(GSE42861_raw$supplementary_file, '/'), function(x){x[9]}), '_'), function(x){x[2]}),
                                         sapply(strsplit(sapply(strsplit(GSE42861_raw$supplementary_file, '/'), function(x){x[9]}), '_'), function(x){x[3]}),
                                         ifelse(GSE42861_raw$`gender:ch1`=='m', 'Male', 'Female'),
                                         GSE42861_raw$`age:ch1`));
colnames(GSE42861_metadata) <- col_names;
rm(GSE42861_raw);
GSE42861_metadata$Sample_name <- as.character(GSE42861_metadata$Sample_name);
GSE42861_metadata$GEO_sample <- as.character(GSE42861_metadata$GEO_sample);
GSE42861_metadata$Batch <- as.character(GSE42861_metadata$Batch);
GSE42861_metadata$Platform <- as.character(GSE42861_metadata$Platform);
GSE42861_metadata$Tissue <- as.character(GSE42861_metadata$Tissue);
GSE42861_metadata$Slide_ID <- as.character(GSE42861_metadata$Slide_ID);
GSE42861_metadata$Array_ID <- as.character(GSE42861_metadata$Array_ID);
GSE42861_metadata$Gender <- as.character(GSE42861_metadata$Gender);
GSE42861_metadata$Age_years <- as.numeric(as.character(GSE42861_metadata$Age_years));

## GSE40279.

GSE40279_raw <- as.data.frame(fread('metadata_files/metadata_GSE40279.tsv', header=T, sep='\t'));
GSE40279_raw$title <- sapply(strsplit(GSE40279_raw$title, ' '), function(x){x[3]});
GSE40279_keys <- as.data.frame(fread('metadata_files/GSE40279_sample_key.txt', header=F, sep='\t'));
colnames(GSE40279_keys) <- c('col1', 'title', 'GEO_sample');
GSE40279_merged <- merge(GSE40279_raw, GSE40279_keys, by='title');
GSE40279_metadata <- as.data.frame(cbind(GSE40279_merged$title, GSE40279_merged$GEO_sample, 'GSE40279', GSE40279_merged$platform_id, 'Whole blood',
                                         sapply(strsplit(GSE40279_merged$GEO_sample, '_'), function(x){x[1]}),
                                         sapply(strsplit(GSE40279_merged$GEO_sample, '_'), function(x){x[2]}),
                                         ifelse(GSE40279_merged$`gender:ch1`=='F', 'Female', 'Male'),
                                         GSE40279_merged$`age (y):ch1`));
colnames(GSE40279_metadata) <- col_names;
GSE40279_na_samples <- c('5901393011_R05C01', '5901393011_R01C02', '5901393027_R01C01'); # Remove the 3 samples that were not present in raw data
GSE40279_metadata <- GSE40279_metadata[!(GSE40279_metadata$GEO_sample %in% GSE40279_na_samples),];
rm(GSE40279_raw, GSE40279_keys, GSE40279_merged);
GSE40279_metadata$Sample_name <- as.character(GSE40279_metadata$Sample_name);
GSE40279_metadata$GEO_sample <- as.character(GSE40279_metadata$GEO_sample);
GSE40279_metadata$Batch <- as.character(GSE40279_metadata$Batch);
GSE40279_metadata$Platform <- as.character(GSE40279_metadata$Platform);
GSE40279_metadata$Tissue <- as.character(GSE40279_metadata$Tissue);
GSE40279_metadata$Slide_ID <- as.character(GSE40279_metadata$Slide_ID);
GSE40279_metadata$Array_ID <- as.character(GSE40279_metadata$Array_ID);
GSE40279_metadata$Gender <- as.character(GSE40279_metadata$Gender);
GSE40279_metadata$Age_years <- as.numeric(as.character(GSE40279_metadata$Age_years));

## GSE51032.

GSE51032_raw <- as.data.frame(fread('metadata_files/metadata_GSE51032.tsv', header=T, sep='\t'));
GSE51032_controls <- GSE51032_raw[is.na(GSE51032_raw$`time to diagnosis:ch1`),]; # Remove people diagnosed with cancer
GSE51032_metadata <- as.data.frame(cbind(GSE51032_controls$title, GSE51032_controls$geo_accession, 'GSE51032', GSE51032_controls$platform_id, 'Peripheral blood', 
                                         sapply(strsplit(GSE51032_controls$title, '_'), function(x){x[1]}),
                                         sapply(strsplit(GSE51032_controls$title, '_'), function(x){x[2]}),
                                         ifelse(GSE51032_controls$`gender:ch1`=="F", "Female", "Male"), GSE51032_controls$`age:ch1`));
colnames(GSE51032_metadata) <- col_names;
rm(GSE51032_raw, GSE51032_controls);
GSE51032_metadata$Sample_name <- as.character(GSE51032_metadata$Sample_name);
GSE51032_metadata$GEO_sample <- as.character(GSE51032_metadata$GEO_sample);
GSE51032_metadata$Batch <- as.character(GSE51032_metadata$Batch);
GSE51032_metadata$Platform <- as.character(GSE51032_metadata$Platform);
GSE51032_metadata$Tissue <- as.character(GSE51032_metadata$Tissue);
GSE51032_metadata$Slide_ID <- as.character(GSE51032_metadata$Slide_ID);
GSE51032_metadata$Array_ID <- as.character(GSE51032_metadata$Array_ID);
GSE51032_metadata$Gender <- as.character(GSE51032_metadata$Gender);
GSE51032_metadata$Age_years <- as.numeric(as.character(GSE51032_metadata$Age_years));

## GSE74432. 

GSE74432_metadata <- as.data.frame(fread('metadata_files/metadata_GSE74432_controls.tsv'));
GSE74432_metadata$Slide_ID <- as.character(GSE74432_metadata$Slide_ID);
GSE74432_metadata$Age_years <- as.numeric(GSE74432_metadata$Age_years);

## Merge all the metadata.

metadata_all <- rbind(GSE59065_metadata, GSE111629_metadata, GSE104812_metadata, GSE81961_metadata, GSE61496_metadata, AE_metadata,
                      GSE55491_metadata, GSE97362_metadata, GSE41273_metadata, GSE42861_metadata, GSE40279_metadata, GSE51032_metadata, GSE74432_metadata);
metadata_all$Disease_status <- 'Control';
colnames(metadata_all)[8] <- 'Sex';


##### 2. Perform QC. #####

## Read the QC information.

qc_files <- list.files(path='QC_files/', full.names = TRUE);

for(i in 1:length(qc_files)){
  if(i==1){
    QC_matrix <- fread(qc_files[i]); 
  }else{
    QC_matrix <- rbind(QC_matrix, fread(qc_files[i]));
  }
}

QC_matrix <- as.data.frame(QC_matrix);
colnames(QC_matrix)[1] <- 'GEO_sample'; 
QC_merge <- merge(metadata_all, QC_matrix, by='GEO_sample'); # 1 sample less than metadata_all (GSM1235539 from GSE51032). It was removed from IDAT files

## Apply QC criteria, keep those samples that satisfy: (median(log2(M)) + median(log2(U)))/2 >= 10.5; predicted sex with minfi == reported sex in metadata; sample finished BMIQ normalisation.

QC_remove1 <- as.character(QC_merge[((QC_merge$mMed + QC_merge$uMed)/2)<10.5, 1]);
QC_remove2 <- as.character(QC_merge[ifelse(QC_merge$predictedSexQC=='F', 'Female', 'Male')!=QC_merge$Sex,1]);
QC_remove3 <- as.character(QC_merge[QC_merge$BMIQ=='FAIL',1]);
QC_discard <- QC_merge[QC_merge$GEO_sample %in% unique(c(QC_remove1,QC_remove2,QC_remove3)),];
metadata_all <- metadata_all[!metadata_all$GEO_sample%in%QC_discard$GEO_sample,];


##### 3. Cell composition. #####

## Read cell composition information and merge. 

cc_files <- list.files(path='cc_files/', full.names = TRUE);

for(i in 1:length(cc_files)){
  if(i==1){
    ccs <- fread(cc_files[i]); 
  }else{
    ccs <- rbind(ccs, fread(cc_files[i]));
  }
}

colnames(ccs)[1] <- 'GEO_sample';
diff_sample_names <- which(sapply(strsplit(ccs$GEO_sample, '_'), function(x){length(x)<3}));
ccs$GEO_sample[-diff_sample_names] <- sapply(strsplit(ccs$GEO_sample[-diff_sample_names], '_'), function(x){x[1]}); # Edit all names except for Aref and GSE40279 datasets
metadata_with_cc <- merge(metadata_all, ccs, by='GEO_sample');


##### 4. DNAmAge raw information (i.e. no background correction).

DNAm_raw_files <- list.files(path='DNAmAge_raw/', full.names = TRUE);

for(i in 1:length(DNAm_raw_files)){
  if(i==1){
    DNAm_raw_ages <- fread(DNAm_raw_files[i]); 
  }else{
    DNAm_raw_ages <- rbind(DNAm_raw_ages, fread(DNAm_raw_files[i]));
  }
}

colnames(DNAm_raw_ages)[1] <- 'GEO_sample';
merged_data_1 <- as.data.frame(merge(metadata_with_cc, DNAm_raw_ages, by='GEO_sample'));
colnames(merged_data_1)[17:24] <- paste0(colnames(merged_data_1)[17:24], '_raw');


##### 5. DNAmAge noob information (i.e. noob background correction).

DNAm_noob_files <- list.files(path='DNAmAge_noob/', full.names = TRUE);

for(i in 1:length(DNAm_noob_files)){
  if(i==1){
    DNAm_noob_ages <- fread(DNAm_noob_files[i]); 
  }else{
    DNAm_noob_ages <- rbind(DNAm_noob_ages, fread(DNAm_noob_files[i]));
  }
}

colnames(DNAm_noob_ages)[1] <- 'GEO_sample';

merged_data_2 <- as.data.frame(merge(merged_data_1, DNAm_noob_ages, by='GEO_sample'));
colnames(merged_data_2)[25:32] <- paste0(colnames(merged_data_2)[25:32], '_noob');
merged_data_2$Comment_raw[merged_data_2$Comment_raw==""] <- NA;
merged_data_2$Comment_noob[merged_data_2$Comment_noob==""] <- NA;

write.table(x=merged_data_2, file='benchmark_background_correction/blood_control_merged_data_benchmark_bc.tsv', quote=F, 
            sep='\t', row.names=F); # Export data so far to benchmark background correction method

################################################################
################## End of the script ###########################
################################################################
