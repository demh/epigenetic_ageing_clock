###########################################################################################
#########                                                                         #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             21/06/2018                                  #########
#########                              EMBL-EBI                                   #########
#########                           Thornton group                                #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####              Biological insights into the epigenetic ageing clock           #########
###########################################################################################
##### Create the DNA methylation matrices (needed to calculate DNAmAge) given the      ####
##### IDAT files in a folder.                                                          ####
###########################################################################################
##### USAGE: Rscript idat_to_DNAmAge_methylation_matrix.R path/to/idats path/to/output path/to/annotation  preprocessing_option
##### ARGS: path/to/idats: absolute path to the folder containing the raw IDAT data.   ####
#####                      The name of the batch should be in the folder previous to the final directory.
#####                      e.g. /data/GSE12345/raw_idat/ for the GSE12345 batch.       ####
#####       path/to/output: absolute path to the folder where the matrix will be outputted.
#####       path/to/annotation: absolute path to the folder which contains the common_probes_27K_450K.txt file.
#####       proprocessing_option: one of the following strings: 'raw' (no pre-processing), 'noob' (apply background correction, recommended)
###########################################################################################
##### NOTE1: the IDAT files must have a name with the following format: samplename_slide_array_channel.idat
#####        e.g. GSM3101870_9985131140_R01C01_Grn.idat                                ####
##### NOTE2: The annotation folder must contain:                                       ####
#####      - Common probes between 27K and 450K arrays: common_probes_27K_450K.txt     ####
##### WARNING: the IDAT files will e organised inside the /data/batch/raw_idat/ folder ####
#####          into subfolders (see move_idat_files function).                         ####
###########################################################################################

###########################################################
##################### Dependencies ########################
###########################################################

library(minfi);
library(IlluminaHumanMethylation27kmanifest);


###########################################################
#####################  Arguments ##########################
###########################################################

print('Getting the input arguments ...');

## Fix the paths.

args <- commandArgs(trailingOnly=TRUE);

#path_to_raw_idat <- "~/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/GSE97362/raw_idat/";
path_to_raw_idat <- args[1];
#output_path <- "/Users/dem44/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/DNAmAge/prepare_matrices/";
output_path <- args[2];
#common_probes_path <- "/nfs/research1/thornton/dem44/methylation_clock/raw_data_human/epigenetic_syndromes/utils/common_probes_27K_450K.txt";
common_probes_path <- paste0(as.character(args[3], '/common_probes_27K_450K.txt'));
p_option <- as.character(args[4]);

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

## Create sample annotation files (one for each type of platform). ##

print('Creating MINFI annotation for samples ...');

sample_ann_27K <- unique(data.frame(
  Array=as.character(sapply(names(which(array_class == "27K")), function(x){
    a <- strsplit(strsplit(x, '/')[[1]][length(strsplit(names(which(array_class == "27K")), '/')[[1]])], '_')[[1]][3];
    return(a);
  })),
  Slide=as.character(sapply(names(which(array_class == "27K")), function(x){
    s <- strsplit(x, '/')[[1]][length(strsplit(names(which(array_class == "27K")), '/')[[1]])-1];
    return(s);
  })),
  Basename=gsub('_Red.idat', '', gsub('_Grn.idat', '', names(which(array_class == "27K"))))));

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

common_probes <- readLines(common_probes_path);

print(paste0('The folder ', folder_name, ' contains ', length(list.files(path_to_raw_idat, recursive = TRUE, full.names = TRUE)),
             ' IDAT files.'));


### 3. Process 27K data. ###

if(nrow(sample_ann_27K) > 0){

  print('Reading 27K data ...');
  RGset_27K <-read.metharray.exp(targets = sample_ann_27K);

  print('Converting 27K data to MethylSet ...');
  if(p_option=='raw'){
	MsetEx_27K <- preprocessRaw(RGset_27K);
  }
  if(p_option=='noob'){
  	print('Performing NOOB background correction ...');
 	MsetEx_27K <- preprocessNoob(RGset_27K);
  }

  print('Obtaining Beta-values for 27K data ...');
  betas_27K <- as.data.frame(getBeta(MsetEx_27K, type="Illumina"));

  print('Editing sample names for 27K data ...');
  sample_new_names_27K <- sapply(colnames(betas_27K), function(x){
    s <- strsplit(x, '_')[[1]][1];
    return(s);
  });
  colnames(betas_27K) <- sample_new_names_27K;

  print('Selecting the common probes for the 27K data ...');
  common_27K <- betas_27K[match(common_probes, rownames(betas_27K)),];
  if(nrow(common_27K) != length(common_probes)){stop('We could not find all the common probes in the 27K data.')};
  rownames(common_27K) <- NULL;


}else{

  print('No data was found for the 27K platform.')
  common_27K <- NULL;

}


### 4. Process 450K data. ###

if(nrow(sample_ann_450K) > 0){

  print('Reading 450K data ...');
  RGset_450K <-read.metharray.exp(targets = sample_ann_450K);

  print('Converting 450K data to MethylSet ...');
  if(p_option=='raw'){
        MsetEx_450K <- preprocessRaw(RGset_450K);
  }
  if(p_option=='noob'){
    	print('Performing NOOB background correction ...');
        MsetEx_450K <- preprocessNoob(RGset_450K);
  }

  print('Obtaining Beta-values for 450K data ...');
  betas_450K <- as.data.frame(getBeta(MsetEx_450K, type="Illumina"));

  print('Editing sample names for 450K data ...');
  sample_new_names_450K <- sapply(colnames(betas_450K), function(x){
    s <- strsplit(x, '_')[[1]][1];
    return(s);
  });
  colnames(betas_450K) <- sample_new_names_450K;

  print('Selecting the common probes for the 450K data ...');
  common_450K <- betas_450K[match(common_probes, rownames(betas_450K)),];
  if(nrow(common_450K) != length(common_probes)){stop('We could not find all the common probes in the 450K data.')};
  rownames(common_450K) <- NULL;

}else{

  print('No data was found for the 450K platform.')
  common_450K <- NULL;

}


### 5. Create the final matrix. ###

# Merge the data if necessary.

if(length(common_27K)>0 & length(common_450K)>0){

  print('Merging 27K and 450K data ...');
  meth_matrix <- cbind(data.frame(ProbeID=common_probes), common_27K, common_450K);

}

if(length(common_27K)>0 & length(common_450K)==0){

  print('Only considering 27K data ...');
  meth_matrix <- cbind(data.frame(ProbeID=common_probes), common_27K);

}

if(length(common_27K)==0 & length(common_450K)>0){

  print('Only considering 450K data ...');
  meth_matrix <- cbind(data.frame(ProbeID=common_probes), common_450K);

}


# Export the data.

print('Creating DNA methylation matrix ...');

write.csv(x=meth_matrix, file=paste0(output_path, '/DNAmAge_methylation_matrix_from_idat_', gsub('-', '_', folder_name), '.csv'),
          quote=FALSE, row.names=FALSE);

print('The script finished correctly.');

#########################################################
########## End of the script ############################
#########################################################
