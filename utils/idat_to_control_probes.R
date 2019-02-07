###########################################################################################
#########                                                                         #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             11/06/2018                                  #########
#########                              EMBL-EBI                                   #########
#########                           Thornton group                                #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####              Biological insights into the epigenetic ageing clock           #########
###########################################################################################
##### Extract the control probes from IDAT files. The intensities from the control probes #
##### will be used to estimate the technical variation and correct for batch effects   ####
##### later on. At the moment it only works for 450K arrays.                           ####
###########################################################################################
##### USAGE: Rscript idat_to_control_probes.R path/to/idats /path/to/output            ####
##### ARGS: path/to/idats: absolute path to the folder containing the raw IDAT data.   ####
#####                      The name of the batch should be in the folder previous to the final directory.
#####                      e.g. /data/GSE12345/raw_idat/ for the GSE12345 batch.       ####
#####       path/to/output: absolute path to the folder where the matrix will be outputted.
###########################################################################################
##### NOTE: the IDAT files must have a name with the following format: samplename_slide_array_channel.idat
#####        e.g. GSM3101870_9985131140_R01C01_Grn.idat                                ####
##### WARNING: the IDAT files will be organised inside the /data/batch/raw_idat/ folder ###
#####          into subfolders (see move_idat_files function).                         ####
###########################################################################################
###########################################################################################

###########################################################
##################### Dependencies ########################
###########################################################

library(IlluminaHumanMethylation450kmanifest);
library(minfi);


###########################################################
#####################  Arguments ##########################
###########################################################

print('Getting the input arguments ...');

args <- commandArgs(trailingOnly=TRUE);

## Input, output and annotation paths

#path_to_raw_idat <- "~/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/cell_composition/debug/GSE51032/raw_idat_quick//";
path_to_raw_idat <- args[1];
#output_path <- "~/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/cell_composition/debug/GSE51032/";
output_path <- args[2];

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


############################################################
################## Running the pipeline ####################
############################################################

#### 1. Read the input dataset.

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
  stop('There are 27K samples among the input IDAT files. This version of the script can only handle 450K data.');
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

## Read the input dataset. 

print('Reading input IDATs ...');
input_rg <-read.metharray.exp(targets = sample_ann);
names(colData(input_rg))[4] <- 'Sample_Name'; # Fix bug in minfi::estimateCellCounts


#### 2. Extract the control probes intensities and export them. 

## Create a dataframe with all the control probes information.
# It contains the following columns: ControlAddress,ControlChannel,ControlType,Sample1,Sample2,...,SampleN
# This code is inspired by shinySummarize function in shinyMethyl package.

print('Separating Red and Green intensities ...');
r <- minfi::getRed(input_rg);
g <- minfi::getGreen(input_rg);
if(!all(rownames(r) == rownames(g)) |
   !all(colnames(r) == colnames(g))){stop('The Red and Green matrices are not in the same order')};

controlType <- c("RESTORATION", "BISULFITE CONVERSION I", "BISULFITE CONVERSION II", 
                 "EXTENSION", "HYBRIDIZATION", "NEGATIVE", "NON-POLYMORPHIC", 
                 "NORM_A", "NORM_C", "NORM_G", "NORM_T", "SPECIFICITY I", 
                 "SPECIFICITY II", "TARGET REMOVAL", "STAINING");

print('Extracting control probe information ...');
for (i in 1:length(controlType)) { 
  
  ctrlAddress <- minfi::getControlAddress(input_rg, controlType = controlType[i]);
  redControls <- r[ctrlAddress, ];
  greenControls <- g[ctrlAddress, ];

  if(i==1){ # For the restoration probes
   control_df <- data.frame(ControlAddress=rep(ctrlAddress, 2),
                            ControlChannel=c('Grn', 'Red'),
                            ControlType=rep(controlType[i], 2));
   controls1 <- rbind(greenControls, redControls);
   rownames(controls1) <- NULL;
   control_df <- cbind(control_df, controls1);
   
  }else{
    if(!all(rownames(redControls) == rownames(greenControls))){stop('The Red and Green channels do not have the same control probes.')};
    current_df <- data.frame(ControlAddress=rep(rownames(redControls), 2),
                             ControlChannel=rep(c('Grn', 'Red'), each=length(rownames(redControls))),
                             ControlType=rep(controlType[i], nrow(redControls)*2));
    current_df <- cbind(current_df, rbind(greenControls, redControls));
    control_df <- rbind(control_df, current_df);
  }
  
}

## Edit and export the controls dataframe. 

control_df <- control_df[order(control_df[,1], control_df[,2]),];
rownames(control_df) <- NULL;
print(paste0('Number of control intensities found: ', nrow(control_df), ' i.e. ', nrow(control_df)/2, ' control probes.'));
print(paste0('Number of samples: ', ncol(control_df)-3));

write.table(control_df, file=paste0(output_path, '/', folder_name, '_control_intensities.csv'),
            sep=',', row.names=F, quote=F);

print('The script finished correctly.');

#########################################################
########## End of the script ############################
#########################################################
