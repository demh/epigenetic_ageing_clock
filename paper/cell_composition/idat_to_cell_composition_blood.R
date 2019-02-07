###########################################################################################
#########                                                                         #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             04/06/2018                                  #########
#########                              EMBL-EBI                                   #########
#########                           Thornton group                                #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####              Biological insights into the epigenetic ageing clock           #########
###########################################################################################
##### Estimate the cell-type composition of a blood sample (whole blood, PBMC) given   ####
##### the IDAT files in a folder.                                                      ####
##### The idol_NFB_houseman strategy will be used, which is only applicable to 450K data. #
###########################################################################################
##### USAGE: Rscript path/to/idats /path/to/output /path/to/annotation                 ####
###########################################################################################
##### NOTE: The annotation folder must contain:                                        ####
#####      - Cross-reactive probes from Chen et al.: cross_reactive_probes_Chen_2013.txt  #
#####      - Beta values for the processed IDOL reference: idol_NFB_reference.csv      ####
###########################################################################################

###########################################################
##################### Dependencies ########################
###########################################################

library(minfi);
library(wateRmelon);
library(IlluminaHumanMethylation27kmanifest);
library(EpiDISH);
library(FlowSorted.Blood.450k);


###########################################################
#####################  Arguments ##########################
###########################################################

print('Getting the input arguments ...');

#args <- commandArgs(trailingOnly=TRUE);

## Input, output and annotation paths

path_to_raw_idat <- "~/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/cell_composition/gold_standard/GSE77797/raw_idat/";
#path_to_raw_idat <- args[1];
output_path <- "~/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/cell_composition/gold_standard/";
#output_path <- args[2];
ann_path <- '~/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/utils/'
#ann_path <- args[3];

folder_name <- strsplit(path_to_raw_idat, '/')[[1]][length(strsplit(path_to_raw_idat, '/')[[1]])-1];


################################################################
################## Functions ###################################
################################################################

#### Function: reorganise paths of IDAT files so they can be used by the minfi package.

# path_to_raw: path to the input folder.

move_idat_files <- function(path_to_raw){
  
  all_idat_files_paths <- list.files(path_to_raw, recursive = TRUE, full.names = TRUE);
  
  for(f in all_idat_files_paths){
    
    file_name <- strsplit(f, '/')[[1]][length(strsplit(f, '/')[[1]])];
    slide_f <- strsplit(file_name, '_')[[1]][2];
    array_f <- strsplit(file_name, '_')[[1]][3];
    dir.create(file.path(path_to_raw, slide_f), showWarnings = FALSE);
    
    new_path <- file.path(path_to_raw,slide_f,file_name);
    file.rename(from=f, to=new_path);
    
  }
}


##### Function: pipeline to preprocess the raw data (RG object) to beta-values in a customised way.

## input_rg: RG input object
## bgcor: perform background correction and dye-bias equalization with Noob? (TRUE/FALSE)
## filter: filter out array probes with SNPs (interrogation site, single nucleotide extension), 
#          cross-reactive (Chen et al. 2013) or belonging to sex chromosome? (TRUE/FALSE)
## BMIQ: perform BMIQ normalisation? (TRUE/FALSE).
## path_cr: path to the file that contains the probe IDs of the cross-reactive probes from Chen et al. 2013

preprocess_RG <- function(input_rg, bgcor, filter, BMIQ, path_cr=NA){
  
  ## Background correction and dye-bias normalisation.
  
  if(bgcor){
    
    print('Performing Noob background correction ...');
    input_ms <- preprocessNoob(input_rg); 
    
  }else{
    
    print('No background correction will be performed ...');
    input_ms <- preprocessRaw(input_rg);
    
  }
  
  print(paste0('At the start, there are ', length(input_ms@NAMES), ' probes.'));
  
  ## Filtering steps. 
  
  if(filter){
    
    # Filter out the probes associated with SNPs.
    # See https://www.bioconductor.org/help/course-materials/2014/BioC2014/minfi_BioC2014.pdf
    
    input_ms <- dropLociWithSnps(mapToGenome(input_ms));
    
    # Obtain beta-values.
    
    input_betas <- as.data.frame(getBeta(input_ms, type="Illumina"));
    print(paste0('After removing probes associated with SNPs, there are ', nrow(input_betas), ' probes left.'));
    
    # Filter out cross-reactive probes.
    
    cr_probes <- readLines(path_cr);
    input_betas <- input_betas[!(rownames(input_betas) %in% cr_probes),];
    print(paste0('After removing cross-reactive probes, there are ', nrow(input_betas),
                 ' probes left.'));
    
    # Filter out probes from sex chromosomes.
    
    ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19);
    sex_probes <- ann450k$Name[ann450k$chr %in% c("chrX","chrY")];
    input_betas <- input_betas[!(rownames(input_betas) %in% sex_probes),];
    print(paste0('After removing probes in sex chromosomes, there are ', nrow(input_betas), ' probes left.'));
    
    
  }else{
    
    # Obtain beta-values.
    
    input_betas <- as.data.frame(getBeta(input_ms, type="Illumina"));
    print(paste0('After converting to beta-values, there are ', nrow(input_betas), ' probes left.'));
    
  }
  
  
  ## BMIQ normalisation.
  
  if(BMIQ){
    
    print('Performing BMIQ normalisation ...');
    info_450K_II <- getProbeInfo(IlluminaHumanMethylation450kmanifest, type = c("II"));
    type_II <- which(rownames(input_betas) %in% info_450K_II$Name);
    design_vector <- rep(1,nrow(input_betas));
    design_vector[type_II] <- 2;
    input_norm <- apply(input_betas, 2, BMIQ, design.v=design_vector, plots=F, pri=T);
    input_final_df <- as.data.frame(sapply(input_norm, function(x){
      return(x[[1]]);
    }));
    return(input_final_df);
    
  }else{
    
    return(input_betas);
    
  }
}


############################################################
################## Running the pipeline ####################
############################################################

#### 1. Process the input dataset.

## Check if the folder has any IDAT files.

if(!file.exists(path_to_raw_idat)){stop(paste0('The project ', folder_name, ' has no IDAT files.'))};

## Move the IDAT files to a correct directory tree. ##

print('Rearranging IDAT files in a new directory tree ...');
move_idat_files(path_to_raw_idat);

## Classify the files into 27K / 450K. ##

print('Classifying the array platforms ...');
new_files_paths <- list.files(path_to_raw_idat, full.names = TRUE, recursive=TRUE);
array_class <- sapply(new_files_paths, function(f){
  ifelse(file.size(f) < 800000, "27K", ifelse(file.size(f) < 9000000, "450K", NA));
});

if(sum(is.na(array_class)) > 0){
  stop("The array platform could not be identified in some samples !!");
}

if(sum(array_class=="27K")){
  stop('There are 27K samples among the input IDAT files. This version of the cell-type composition estimation can only handle 450K data.');
}

## Create sample annotation files. ##

print('Creating MINFI annotation for samples ...');
sample_ann <- unique(data.frame(
  Array=as.character(sapply(names(which(array_class == "450K")), function(x){
    a <- strsplit(strsplit(x, '/')[[1]][length(strsplit(names(which(array_class == "450K")), '/')[[1]])], '_')[[1]][3];
    return(a);
  })),
  Slide=as.character(sapply(names(which(array_class == "450K")), function(x){
    s <- strsplit(x, '/')[[1]][length(strsplit(names(which(array_class == "450K")), '/')[[1]])-1];
    return(s);
  })),
  Basename=gsub('_Red.idat', '', gsub('_Grn.idat', '', names(which(array_class == "450K"))))));

print(paste0('The folder ', folder_name, ' contains ', length(list.files(path_to_raw_idat, recursive = TRUE, full.names = TRUE)),
             ' IDAT files.'));

## Process the input dataset. 

# Read the IDAT files.

print('Reading input IDATs ...');
input_rg <-read.metharray.exp(targets = sample_ann);
names(colData(input_rg))[4] <- 'Sample_Name'; # Fix bug in minfi::estimateCellCounts

# Process.

if(!file.exists(paste0(ann_path, '/cross_reactive_probes_Chen_2013.txt'))){
  stop('The file specifying the cross-reactive probes (cross_reactive_probes_Chen_2013.txt) could not be found in the annotation path.');
}

print('Preprocessing RG data ...');
input_NFB_df <- preprocess_RG(input_rg=input_rg, bgcor=T, filter=T, BMIQ=T, path_cr=paste0(ann_path, '/cross_reactive_probes_Chen_2013.txt')); # With filtering


#### 3. Obtain the cell-type composition estimations. 

## Predictions for IDOL with Houseman algorithm. 

# Load the reference.  

idol_our_NFB <- read.csv(paste0(ann_path, '/idol_NFB_reference.csv'));
rownames(idol_our_NFB) <- idol_our_NFB$ProbeID;
idol_our_NFB <- idol_our_NFB[,-1];
idol_our_NFB <- as.matrix(idol_our_NFB);

# Select our input data.

idol_input_NFB <- input_NFB_df[rownames(input_NFB_df) %in% rownames(idol_our_NFB),];
idol_input_NFB <- as.matrix(idol_input_NFB[order(rownames(idol_input_NFB)),]);

# Predictions.

idol_pred <- epidish(avdata.m=idol_input_NFB, ref.m=idol_our_NFB, method="CP");
idol_pred_df <- as.data.frame(idol_pred$estF);
idol_pred_df <- cbind(rownames(idol_pred_df), idol_pred_df);
colnames(idol_pred_df)[1] <- 'SampleID';
rownames(idol_pred_df) <- NULL;

## Export the predictions.

write.table(idol_pred_df, file=paste0(output_path, '/', folder_name, '_cc_predictions.csv'),
            sep=',', row.names=F, quote=F);


################################################################
######################## Extra code ############################
################################################################
########## Create the IDOL NFB reference from scratch. #########
################################################################

# #### 1. Process the Reinius et al. 2012 dataset.
# 
# ## Load as a RG object.
# 
# reinius_rg <- get('FlowSorted.Blood.450k');
# 
# ## Preprocess the data.
# 
# reinius_final_df_NFB <- preprocess_RG(input_rg=reinius_rg, bgcor=T, filter=T, BMIQ=T, path_cr=path_crossr); # With filtering
# 
# ## Create our own references (all the probes).
# 
# # Remove the samples with whole or peripheral blood. 
# 
# reinius_cells_NFB <- as.matrix(reinius_final_df_NFB[,-(c(grep('WB', colnames(reinius_final_df_NFB)),grep('PBMC', colnames(reinius_final_df_NFB))))]);
# 
# # Average the beta-values per cell-type.
# 
# cell_types <- unique(sapply(strsplit(colnames(reinius_cells_NFB), '_'), function(x){x[1]}));
# reinius_averages_NFB <- matrix(NA, ncol=length(cell_types), nrow=nrow(reinius_cells_NFB));
# colnames(reinius_averages_NFB) <- cell_types;
# rownames(reinius_averages_NFB) <- rownames(reinius_cells_NFB);
# 
# i <- 1;
# 
# for(ct in cell_types){
#   c_name <- paste0('mean_', ct);
#   reinius_averages_NFB[,i] <- rowMeans2(reinius_cells_NFB, cols=grep(ct,colnames(reinius_cells_NFB)));
#   i <- i+1;
# }
# 
# # Select only appropiate columns. 
# 
# reinius_averages_NFB <- as.data.frame(reinius_averages_NFB[,1:6]);
# colnames(reinius_averages_NFB) <- c('Gran', 'CD4T', 'CD8T', 'B', 'Mono', 'NK');
# 
# 
# #### 2. Select the IDOL probes from the Reinius data and export the reference. 
# 
# idol_probes <- as.character(read.csv('~/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/cell_composition/IDOL_Table_S3.csv')[,1]);
# idol_our_NFB <- reinius_averages_NFB[rownames(reinius_averages_NFB)%in%idol_probes,];
# idol_our_NFB <- as.matrix(idol_our_NFB[order(rownames(idol_our_NFB)),]); # Order rows
# 
# idol_our_NFB_df <- data.frame(ProbeID=rownames(idol_our_NFB), Gran=idol_our_NFB[,1], CD4T=idol_our_NFB[,2], 
#                               CD8T=idol_our_NFB[,3], B=idol_our_NFB[,4], Mono=idol_our_NFB[,5], NK=idol_our_NFB[,6]);
# rownames(idol_our_NFB_df) <- NULL;
# write.table(idol_our_NFB_df, file='~/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/cell_composition/idol_NFB_reference.csv',
#             quote=F, sep=',', row.names=F);

## The reference can oly be used with 450K and posterior versions of the array.

#########################################################
########## End of the script ############################
#########################################################