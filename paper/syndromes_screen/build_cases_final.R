###########################################################################################
#########                                                                         #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             03/07/2018                                  #########
#########                              EMBL-EBI                                   #########
#########                           Thornton group                                #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####              Biological insights into the epigenetic ageing clock           #########
###########################################################################################
##### Build a dataframe that contains all the information for the cases (samples with #####
##### developmental disorders).                                                       #####
##### Do it with DNAmAge calculations that come from data without ('raw') and with    #####
##### ('noob') background correction.                                                 #####
###########################################################################################
##### USAGE: manual                                                                   #####
###########################################################################################

###########################################################
##################### Dependencies ########################
###########################################################

library(data.table);

setwd('/Users/dem44/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/syndromes_screen/');


################################################################
################## Run the pipeline ############################
################################################################

##### 1. Read and curate the metadata information. #####

col_names <- c('Sample_name','GEO_sample','Batch','Platform','Tissue','Slide_ID','Array_ID','Gender','Age_years', 
               'Disease_status', 'Gene', 'Genome_assembly', 'Mutation', 'Protein_change', 'Mutation_effect', 'Pathogenic', 'Other_comments');

## GSE74432.

GSE74432_metadata <- as.data.frame(fread('metadata_files/metadata_GSE74432_cases.tsv'));
GSE74432_metadata$Slide_ID <- as.character(GSE74432_metadata$Slide_ID);

## GSE55491.

GSE55491_raw <- as.data.frame(fread('metadata_files/GSE55491/metadata_GSE55491.tsv', header=T, sep='\t'));
GSE55491_ages <- as.data.frame(fread('metadata_files/GSE55491/samples4Daniel_edited.csv', header=T, sep=','));
colnames(GSE55491_ages)[1] <- 'source_name_ch1';
GSE55491_merge <- merge(GSE55491_raw, GSE55491_ages, by='source_name_ch1'); 
GSE55491_merge <- GSE55491_merge[GSE55491_merge$SampleGroup=='SRS',]; # Remove 6 controls
GSE55491_metadata <- data.frame(cbind(GSE55491_merge$title, GSE55491_merge$geo_accession, 'GSE55491',
                                         GSE55491_merge$platform_id, 'Whole blood', 
                                         as.character(GSE55491_merge$Sentrix_Barcode), as.character(GSE55491_merge$Sentrix_Position),
                                         ifelse(GSE55491_merge$Sex=='F', 'Female', 'Male'),
                                         GSE55491_merge$Age_years, 'Silver_Russell', NA, NA,
                                         NA,NA,NA,'YES', NA), stringsAsFactors=FALSE);
colnames(GSE55491_metadata) <- col_names;
GSE55491_metadata$Age_years <- as.numeric(GSE55491_metadata$Age_years);
rm(GSE55491_raw,GSE55491_ages,GSE55491_merge);

## GSE97362.

GSE97362_raw <- as.data.frame(fread('metadata_files/GSE97362/metadata_GSE97362.tsv', header=T, sep='\t'));
GSE97362_raw$SampleID <- sapply(strsplit(GSE97362_raw$title, ' '), function(x){x[1]});
GSE97362_raw <- GSE97362_raw[-grep('control', GSE97362_raw$SampleID),]; # Remove controls
GSE97362_raw <- GSE97362_raw[GSE97362_raw$characteristics_ch1.2 != 'sample type: Validation cohort',]; # Remove 'Validation cohort' batch
GSE97362_charge_muts <- as.data.frame(fread('metadata_files/GSE97362/CHD7_info_final.tsv', header=T, sep='\t'));
GSE97362_kabuki_muts <- as.data.frame(fread('metadata_files/GSE97362/KMT2D_info_final.tsv', header=T, sep='\t'));
GSE97362_muts_all <- rbind(GSE97362_charge_muts,GSE97362_kabuki_muts);
GSE97362_merge <- merge(GSE97362_raw, GSE97362_muts_all, by='SampleID');
GSE97362_metadata <- data.frame(cbind(GSE97362_merge$title, GSE97362_merge$geo_accession, 'GSE97362', GSE97362_merge$platform_id, 'Whole blood',
                                      sapply(strsplit(GSE97362_merge$supplementary_file, '_'), function(x){x[2]}),
                                      sapply(strsplit(GSE97362_merge$supplementary_file, '_'), function(x){x[3]}),
                                      ifelse(GSE97362_merge$characteristics_ch1 == 'gender: female', 'Female', 'Male'),
                                      as.numeric(sapply(strsplit(GSE97362_merge$characteristics_ch1.1, ' '), function(x){x[3]})),
                                      GSE97362_merge$Disease_status, GSE97362_merge$Gene, GSE97362_merge$Genome_assembly, GSE97362_merge$Mutation, GSE97362_merge$Protein_change,
                                      GSE97362_merge$Mutation_effect, GSE97362_merge$Pathogenic, GSE97362_merge$Other_comments), stringsAsFactors=FALSE);
colnames(GSE97362_metadata) <- col_names;
GSE97362_metadata$Age_years <- as.numeric(GSE97362_metadata$Age_years);
rm(GSE97362_raw,GSE97362_charge_muts,GSE97362_kabuki_muts,GSE97362_muts_all,GSE97362_merge);
  
## GSE41273. 

GSE41273_raw <- as.data.frame(fread('metadata_files/GSE41273/metadata_GSE41273.tsv', header=T, sep='\t'));
GSE41273_ages <- as.data.frame(fread('metadata_files/GSE41273/FXS_samples.csv', header=T, sep=','));
colnames(GSE41273_ages)[1] <- 'title';
GSE41273_merge <- merge(GSE41273_raw,GSE41273_ages,by='title');
GSE41273_merge <- GSE41273_merge[GSE41273_merge$STATUS=='Case',]; # Select only cases
GSE41273_metadata <- data.frame(cbind(GSE41273_merge$title, GSE41273_merge$geo_accession, 'GSE41273', GSE41273_merge$platform_id, 'Peripheral blood',
                                      sapply(strsplit(GSE41273_merge$ID450k, '_'), function(x){x[1]}),
                                      sapply(strsplit(GSE41273_merge$ID450k, '_'), function(x){x[2]}),
                                      'Male', as.numeric(GSE41273_merge$Ados_Age_months)/12, 'FXS', 'FMR1', NA, NA, NA, 'CGG repeat expansion', 
                                      'YES', ''), stringsAsFactors=FALSE);
colnames(GSE41273_metadata) <- col_names;
GSE41273_metadata$Age_years <- as.numeric(GSE41273_metadata$Age_years);
rm(GSE41273_raw, GSE41273_ages, GSE41273_merge);

## Aref-Eshghi dataset.

AE_raw <- as.data.frame(fread('metadata_files/Aref_Eshghi_dataset/results_Erfan.csv', header=T, sep=','));
AE_disease <- as.data.frame(fread('metadata_files/Aref_Eshghi_dataset/Aref_Eshghi_disease_metadata_curated.csv'));
AE_all <- merge(AE_raw, AE_disease, by='SampleID');
AE_remove <- c('ASD_replicate', 'ATR-X_bad_annotation', 'Atypical_Rett', 'Beckwith_Wiedemann_bad', 'Coffin_Lowry_bad_annotation', 
               'Control', 'Control_replicate', 'FXS_47XXY', 'FXS_bad_sample', 'FXS_grey', 'FXS_normal', 'FXS_replicate', 
               'MECP2_non_Rett', 'Noonan_MAP2K1_bad_annotation', 'Rett_Angelman??', 'Rett_bad_annotation', 'Saethre_Chotzen_bad_annotation',
               'TWIST_sex_mismatch');
AE_cases <- AE_all[!AE_all$Disease_status%in%AE_remove,]; # Select only cases with correct annotation
AE_cases <- AE_cases[(!(is.na(AE_cases$age)) | AE_cases$Gene=='DNMT1'),]; # Remove samples without ages, but keep DNMT1 mutants
AE_metadata <- as.data.frame(cbind(AE_cases$id, AE_cases$SampleID, AE_cases$batch, 'GPL13534', 'Peripheral blood',
                                   sapply(strsplit(AE_cases$idat, '_'), function(x){x[1]}), sapply(strsplit(AE_cases$idat, '_'), function(x){x[2]}),
                                   ifelse(AE_cases$sex=='f', 'Female', 'Male'), AE_cases$age, AE_cases$Disease_status, AE_cases$Gene, AE_cases$Genome_assembly,
                                   AE_cases$Mutation, AE_cases$Protein_change, AE_cases$Mutation_effect, AE_cases$Pathogenic, AE_cases$Other_comments), stringsAsFactors=FALSE);
colnames(AE_metadata) <- col_names;
AE_metadata$Age_years <- as.numeric(AE_metadata$Age_years);
rm(AE_raw, AE_disease, AE_all, AE_cases);

## GSE116300.

GSE116300_metadata <- as.data.frame(fread('metadata_files/GSE116300/metadata_GSE116300_curated.csv'));
GSE116300_metadata$Slide_ID <- as.character(GSE116300_metadata$Slide_ID);
GSE116300_metadata <- GSE116300_metadata[!is.na(GSE116300_metadata$Age_years),]; 


## Merge all the metadata. 

metadata_all <- rbind(GSE74432_metadata, GSE55491_metadata, GSE97362_metadata, 
                      GSE41273_metadata, AE_metadata, GSE116300_metadata);
metadata_all$Gene[metadata_all$Gene==""] <- NA;
metadata_all$Mutation[metadata_all$Mutation==""] <- NA;
metadata_all$Protein_change[metadata_all$Protein_change==""] <- NA;
metadata_all$Mutation_effect[metadata_all$Mutation_effect==""] <- NA;
metadata_all$Pathogenic[metadata_all$Pathogenic==""] <- NA;
metadata_all$Other_comments[metadata_all$Other_comments==""] <- NA;
metadata_all$Genome_assembly[is.na(metadata_all$Gene)] <- NA;
metadata_all$Genome_assembly[!is.na(metadata_all$Gene)] <- 'hg19';
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
QC_merge <- merge(metadata_all, QC_matrix, by='GEO_sample');

## Apply QC criteria, keep those samples that satisfy: (median(log2(M)) + median(log2(U)))/2 >= 10.5; predicted sex with minfi == reported sex in metadata; sample finished BMIQ normalisation.

QC_remove1 <- as.character(QC_merge[((QC_merge$mMed + QC_merge$uMed)/2)<10.5, 1]);
QC_remove2 <- as.character(QC_merge[ifelse(QC_merge$predictedSexQC=='F', 'Female', 'Male')!=QC_merge$Sex,1]);
QC_remove3 <- as.character(QC_merge[QC_merge$BMIQ=='FAIL',1]);
QC_discard <- QC_merge[QC_merge$GEO_sample %in% unique(c(QC_remove1,QC_remove2,QC_remove3)),];
metadata_all <- metadata_all[!metadata_all$GEO_sample%in%QC_discard$GEO_sample,]; # We remove all the DNMT1 mutants


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
ccs$GEO_sample[-diff_sample_names] <- sapply(strsplit(ccs$GEO_sample[-diff_sample_names], '_'), function(x){x[1]}); # Edit all names needed
metadata_with_cc <- merge(metadata_all, ccs, by='GEO_sample');


##### 4. DNAmAge raw information (i.e. no background correction). #####

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
colnames(merged_data_1)[24:31] <- paste0(colnames(merged_data_1)[24:31], '_raw');


##### 5. DNAmAge noob information (i.e. noob background correction). #####

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
colnames(merged_data_2)[32:39] <- paste0(colnames(merged_data_2)[32:39], '_noob');
merged_data_2$Comment_raw[merged_data_2$Comment_raw==""] <- NA;
merged_data_2$Comment_noob[merged_data_2$Comment_noob==""] <- NA;

write.table(x=merged_data_2, file='../blood_control/benchmark_background_correction/cases_merged_data_benchmark_bc.tsv', quote=F, 
            sep='\t', row.names=F); # Export data so far to benchmark background correction method


################################################################
################## End of the script ###########################
################################################################