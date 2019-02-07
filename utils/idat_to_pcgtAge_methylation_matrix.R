###########################################################################################
#########                                                                         #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             07/11/2017                                  #########
#########                              EMBL-EBI                                   #########
#########                           Thornton group                                #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####              Biological insights into the epigenetic ageing clock           #########
###########################################################################################
##### Create the DNA methylation matrices (needed to calculate pcgtAge) given the      ####
##### IDAT files in a folder (e.g. GSE). Since the probes needed are only found in the ####
##### 450K data, we will only include these files. Furthermore, BMIQ normalisation will be
##### applied, as described in the pcgtAge paper.                                      ####
###########################################################################################
##### USAGE: Rscript idat_to_pcgtAge_methylation_matrix_with_BMIQ.R  path/to/idats path/to/output path/to/annotation
##### ARGS: path/to/idats: absolute path to the folder containing the raw IDAT data.   ####
#####                      The name of the batch should be in the folder previous to the final directory.
#####                      e.g. /data/GSE12345/raw_idat/ for the GSE12345 batch.       ####
#####       path/to/output: absolute path to the folder where the matrix will be outputted.
#####       path/to/annotation: absolute path to the folder which contains the mitotic clock probes file.
###########################################################################################
##### NOTE1: the IDAT files must have a name with the following format: samplename_slide_array_channel.idat
#####        e.g. GSM3101870_9985131140_R01C01_Grn.idat                                ####
##### NOTE2: The annotation folder must contain:                                       ####
#####      - Epigenetic mitotic clock probes: pcgt_Age_probes.txt                      #### 
##### WARNING: the IDAT files will e organised inside the /data/batch/raw_idat/ folder ####
#####          into subfolders (see move_idat_files function).                         ####
###########################################################################################
###########################################################################################

###########################################################
##################### Dependencies ########################
###########################################################

library(minfi);
library(wateRmelon);
library(IlluminaHumanMethylation450kmanifest);


###########################################################
#####################  Arguments ##########################
###########################################################

print('Getting the input arguments ...');

## Fix the paths.

args <- commandArgs(trailingOnly=TRUE);

path_to_raw_idat <- args[1];
#path_to_raw_idat <- "~/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/GSE97362/raw_idat/";
output_path <- args[2];
#output_path <- "/Users/dem44/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/pcgtAge";
pcgtAge_probes_path <- paste0(as.character(args[3]), '/pcgt_Age_probes.txt');
#pcgtAge_probes_path <- "/nfs/research2/thornton/dem44/methylation_clock/raw_data_human/epigenetic_syndromes/utils/pcgt_Age_probes.txt";

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


################################################################
################## Running the pipeline ########################
################################################################

### 1. Preliminary steps. ###

## Check if the project has any IDAT files. 

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

if(all(array_class == "27K")){
  stop("All the IDAT files are from the 27K array platform and no methylation matrix will be generated.");
}

## Create sample annotation files. ##

print('Creating MINFI annotation for samples ...');

sample_ann_450K <- unique(data.frame(
  Array=as.character(sapply(names(which(array_class == "450K")), function(x){
    a <- strsplit(strsplit(x, '/')[[1]][length(strsplit(names(which(array_class == "450K")), '/')[[1]])], '_')[[1]][3];
    return(a);
  })),
  Slide=as.character(sapply(names(which(array_class == "450K")), function(x){
    s <- strsplit(x, '/')[[1]][length(strsplit(names(which(array_class == "450K")), '/')[[1]])-1];
    return(s);
  })),
  Basename=gsub('_Red.idat', '', gsub('_Grn.idat', '', names(which(array_class == "450K"))))));


### 2. Obtain some useful metadata. ###

pcgtAge_probes <- readLines(pcgtAge_probes_path);

print(paste0('The folder ', folder_name, ' contains ', length(list.files(path_to_raw_idat, recursive = TRUE, full.names = TRUE)),
             ' IDAT files.'));
print(paste0('Of these, ', sum(array_class=="450K"), " correspond to 450K data."));


### 3. Process 450K data. ###

print('Reading 450K data ...');
RGset_450K <-read.metharray.exp(targets = sample_ann_450K);

print('Converting 450K data to MethylSet ...');
MsetEx_450K <- preprocessRaw(RGset_450K);

print('Obtaining Beta-values for 450K data ...');
betas_450K <- as.data.frame(getBeta(MsetEx_450K, type="Illumina"));

print('Normalising Beta-values ...');
info_450K_II <- getProbeInfo(IlluminaHumanMethylation450kmanifest, type = c("II"));
type_II <- which(rownames(betas_450K) %in% info_450K_II$Name);
design_vector <- rep(1,nrow(betas_450K));
design_vector[type_II] <- 2;
betas_450K_BMIQ_norm <- apply(betas_450K, 2, BMIQ, design.v=design_vector, plots=FALSE, pri=FALSE);
betas_final_df <- as.data.frame(sapply(betas_450K_BMIQ_norm, function(x){
  return(x[[1]]);
}));

print('Editing sample names for 450K data ...');
sample_new_names_450K <- sapply(colnames(betas_final_df), function(x){
  s <- strsplit(x, '_')[[1]][1];
  return(s);
});
colnames(betas_final_df) <- sample_new_names_450K;

print('Selecting the pcgtAge probes for the 450K data ...');
pcgtAge_450K <- betas_final_df[match(pcgtAge_probes, rownames(betas_final_df)),];
if(nrow(pcgtAge_450K) != length(pcgtAge_probes)){stop('We could not find all the pcgtAge probes in the 450K data.')};
rownames(pcgtAge_450K) <- NULL;
  

### 4. Create the final methylation matrix with the following columns: ProbeID,Sample1,Sample2,...,SampleN ###

meth_matrix <- cbind(data.frame(ProbeID=pcgtAge_probes), pcgtAge_450K);

# Export the data.

print('Creating DNA methylation matrix ...');

write.csv(x=meth_matrix, file=paste0(output_path, '/pcgtAge_methylation_matrix_from_idat_', gsub('-', '_', folder_name), '.csv'),
          quote=FALSE, row.names=FALSE);

print('The script finished correctly.');

#########################################################
########## End of the script ############################
#########################################################
