###########################################################################################
#########                                                                         #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             15/08/2018                                  #########
#########                              EMBL-EBI                                   #########
#########                           Thornton group                                #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####              Biological insights into the epigenetic ageing clock           #########
###########################################################################################
##### Create the metadata for GSE74432, including Sotos, Weaver and control samples.   ####
###########################################################################################
##### USAGE: manual                                                                    ####
###########################################################################################

###########################################################
##################### Dependencies ########################
###########################################################

library(data.table);
setwd('~/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/blood_control/metadata_files/metadata_GSE74432/');


################################################################
################## Run the pipeline ############################
################################################################

##### 1. Read all the sources of information and edit as necessary. #####

# Sotos discovery: 19   -- Supp. Data 1 (we remove 1 without age: 18 left)
# Controls: 57 -- Supp. Data 2  (we remove 4 fibroblasts and 2 without age: 51 left)
# Sotos validation (Hong Kong cohort): 19   -- Supp. Data 5   MISSING AGES
# Sotos fibroblast: 3 -- Supp. Data 1 (not included in this analysis)
# NSD1 VUS: 16 (9/16 classified as Sotos according to classifier) MISSING AGES. 3 of them in 'Copy of Supplementary Table 6_EZH2 vairants.xlsx' (i.e. with ages). 
#            2 of them classified as Sotos and with the age available. 
# Weaver: 8  -- Supp. Data 6 + Ages = Copy of Supplementary Table 6_EZH2 vairants.xlsx 

GEO_metadata <- fread('metadata_GSE74432.tsv');
GEO_metadata <- GEO_metadata[GEO_metadata$`tissue:ch1`=='whole blood',]; # Remove fibroblasts
GEO_metadata$title <- sapply(strsplit(GEO_metadata$title, ' '), function(x){x[1]}); # Edit sample names

supp_data_1 <- fread('ncomms10207-s2_edited.tsv'); # Sotos discovery patients, we do not include fibroblast samples
supp_data_1 <- cbind(supp_data_1[,1:10], 'Sotos'); # Select useful columns
colnames(supp_data_1)[1] <- 'title';
colnames(supp_data_1)[11] <- 'Disease_status';
supp_data_1 <- supp_data_1[!is.na(supp_data_1$Age_years),]; # Remove samples with NA in age
supp_data_1[1,1] <- '11D/0326'; # Edit sample name to make it match with GEO name
supp_data_1[2,1] <- '11D/0328'; # Edit sample name to make it match with GEO name

supp_data_7 <- fread('ncomms10207-s8_edited.tsv'); # Sotos VUS patients, with positive signature and ages. 
supp_data_7 <- cbind(supp_data_7[,1:10], 'Sotos'); # Select useful columns
colnames(supp_data_7)[1] <- 'title';
colnames(supp_data_7)[11] <- 'Disease_status';

weaver_ages <- fread('Copy_of_Supplementary_Table_6_EZH2_variants_edited.csv'); # Weaver patients
weaver_ages <- cbind(weaver_ages, 'Weaver');
colnames(weaver_ages)[1] <- 'title';
colnames(weaver_ages)[11] <- 'Disease_status';
weaver_ages[1,1] <- 'A123W'; # Edit sample name to make it match with GEO name
weaver_ages[2,1] <- 'A134W'; # Edit sample name to make it match with GEO name
weaver_ages[4,1] <- 'A120W'; # Edit sample name to make it match with GEO name

all_cases <- rbind(supp_data_1, supp_data_7, weaver_ages);

supp_data_2 <- fread('ncomms10207-s3_edited.csv'); # Controls
supp_data_2 <- supp_data_2[supp_data_2$Sample_Group=='Control-blood',]; # Remove fibroblasts
supp_data_2 <- supp_data_2[,c(1,3,4)]; # Select useful columns
colnames(supp_data_2) <- c('title', 'Gender', 'Age_years');
supp_data_2$Age_years <- as.numeric(supp_data_2$Age_years);
supp_data_2 <- supp_data_2[!is.na(supp_data_2$Age_years),]; # Remove samples with NA in age


##### 2. Merge all the data, select useful columns and export. #####

## Cases. This will be used in the screening of the syndromes. 

col_names_cases <- c('Sample_name','GEO_sample','Batch','Platform','Tissue','Slide_ID','Array_ID','Gender','Age_years', 
               'Disease_status', 'Gene', 'Genome_assembly', 'Mutation', 'Protein_change', 'Mutation_effect', 'Pathogenic', 'Other_comments');
merged_cases <- merge(GEO_metadata, all_cases, by=c('title'));
final_cases <- as.data.frame(cbind(merged_cases$title, merged_cases$geo_accession, 'GSE74432', merged_cases$platform_id, 'Whole blood',
                                   sapply(strsplit(merged_cases$supplementary_file, '_'), function(x){x[2]}),
                                   sapply(strsplit(merged_cases$supplementary_file, '_'), function(x){x[3]}),
                                   ifelse(merged_cases$Gender=='F', 'Female', 'Male'), merged_cases$Age_years,
                                   merged_cases$Disease_status, merged_cases$Gene, merged_cases$Genome_assembly, merged_cases$Mutation,
                                   merged_cases$Protein_change, merged_cases$Mutation_effect, merged_cases$Pathogenic, merged_cases$Other_comments));
colnames(final_cases) <- col_names_cases;
write.table(x=final_cases, file='~/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/syndromes_screen/metadata_files/metadata_GSE74432_cases.tsv', quote=F, sep='\t', row.names=F, col.names=T);

## Controls. This will be used in to build the blood control.

col_names_controls <- c('Sample_name','GEO_sample','Batch','Platform','Tissue','Slide_ID','Array_ID','Gender','Age_years');
merged_controls <- merge(GEO_metadata, supp_data_2, by='title');
final_controls <- as.data.frame(cbind(merged_controls$title, merged_controls$geo_accession, 'GSE74432', merged_controls$platform_id, 'Whole blood',
                                      sapply(strsplit(merged_controls$supplementary_file, '_'), function(x){x[2]}),
                                      sapply(strsplit(merged_controls$supplementary_file, '_'), function(x){x[3]}),
                                      ifelse(merged_controls$Gender=='F', 'Female', 'Male'), merged_controls$Age_years));
colnames(final_controls) <- col_names_controls;
write.table(x=final_controls, file='~/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/blood_control/metadata_files/metadata_GSE74432/metadata_GSE74432_controls.tsv', quote=F, sep='\t', row.names=F, col.names=T);

#########################################################
########## End of the script ############################
#########################################################
