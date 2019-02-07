###########################################################################################
#########                                                                         #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             01/06/2018                                  #########
#########                              EMBL-EBI                                   #########
#########                           Thornton group                                #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####              Biological insights into the epigenetic ageing clock           #########
###########################################################################################
##### Create the DNA methylation matrices (needed to calculate DNAmAge) given the      ####
##### raw signals from GSE40279 (Hannum clock dataset).                                ####
###########################################################################################
##### USAGE: Rsript signals_GSE40279_to_DNAmAge_methylation_matrix_with_prepoc_option.R path/to/output path/to/annotation  preprocessing_option
##### ARGS: path/to/output: absolute path to the folder where the matrix will be outputted.
#####       path/to/annotation: absolute path to the folder which contains the common_probes_27K_450K.txt file.
#####       proprocessing_option: one of the following strings: 'raw' (no pre-processing), 'noob' (apply background correction, recommended)
###########################################################################################
##### NOTE: The annotation folder must contain:                                       ####
#####      - Common probes between 27K and 450K arrays: common_probes_27K_450K.txt    ####
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

path_to_raw_red <- "/nfs/research1/thornton/dem44/methylation_clock/raw_data_human/epigenetic_syndromes/GSE40279/raw_idat/Hannum_Red_Channel_noRS_paperSamples.csv";
path_to_raw_green <- "/nfs/research1/thornton/dem44/methylation_clock/raw_data_human/epigenetic_syndromes/GSE40279/raw_idat/Hannum_Green_Channel_noRS_paperSamples.csv";
#output_path <- "/nfs/nobackup/thornton/dem44/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/blood_control/DNAmAge/";
output_path <- args[1];
#common_probes_path <- "/nfs/research1/thornton/dem44/methylation_clock/raw_data_human/epigenetic_syndromes/utils/common_probes_27K_450K.txt";
common_probes_path <- paste0(as.character(args[2], '/common_probes_27K_450K.txt'));
p_option <- as.character(args[3]);


################################################################
################## Running the pipeline ########################
################################################################

### 1. Preliminary steps. ###

## Read the raw signals. ##

print('Reading the raw data ...');

raw_red <- read.csv(path_to_raw_red, header=TRUE,row.names=1);
mat_red <- as.matrix(raw_red);
colnames(mat_red) <- gsub('X', '', colnames(raw_red));
rm(raw_red);
raw_green <- read.csv(path_to_raw_green, header=TRUE,row.names=1);
mat_green <- as.matrix(raw_green);
colnames(mat_green) <- gsub('X', '', colnames(raw_green));
rm(raw_green);


### 2. Obtain some useful metadata. ###

common_probes <- readLines(common_probes_path);


### 3. Process 450K data. ###

print('Creating the RG object ...');
RGset_450K <- RGChannelSet(Green= mat_green,Red= mat_red,annotation = c(array="IlluminaHumanMethylation450k", annotation="ilmn12.hg19"));

print('Converting 450K data to MethylSet ...');
if(p_option=='raw'){
        MsetEx_450K <- preprocessRaw(RGset_450K);
}
if(p_option=='noob'){
        print('Performing NOOB backaground correction ...');
        MsetEx_450K <- preprocessNoob(RGset_450K);
}
rm(mat_red, mat_green, RGset_450K);

print('Obtaining Beta-values for 450K data ...');
betas_450K <- as.data.frame(getBeta(MsetEx_450K, type="Illumina"));
rm(MsetEx_450K);

print('Selecting the common probes for the 450K data ...');
common_450K <- betas_450K[match(common_probes, rownames(betas_450K)),];
if(nrow(common_450K) != length(common_probes)){stop('We could not find all the common probes in the 450K data.')};
rownames(common_450K) <- NULL;


### 4. Create the final matrix. ###

print('Creating the final matrix ...');
meth_matrix <- cbind(data.frame(ProbeID=common_probes), common_450K);

# Export the data.

print('Creating DNA methylation matrix ...');
write.csv(x=meth_matrix, file=paste0(output_path, '/DNAmAge_methylation_matrix_from_idat_GSE40279.csv'),
          quote=FALSE, row.names=FALSE);

print('The script finished correctly.');

#########################################################
########## End of the script ############################
#########################################################
