###########################################################################################
#########                                                                         #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             13/04/2018                                  #########
#########                              EMBL-EBI                                   #########
#########                           Thornton group                                #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####              Biological insights into the epigenetic ageing clock           #########
###########################################################################################
##### Download the supplementary data (IDAT files, ...) and the metadata from a GEO    ####
##### series given a GSE identifier.                                                   ####
###########################################################################################
##### USAGE: Rscript get_data_from_GEO.R 'GSE12345' output_path                        ####
###########################################################################################

###########################################################
##################### Dependencies ########################
###########################################################

library(GEOquery);


###########################################################
#####################  Arguments ##########################
###########################################################

print('Getting the input arguments ...');

## Input, output and annotation paths

args <- commandArgs(trailingOnly=TRUE)
GSE_ID <- as.character(args[1]);
output_path <- as.character(args[2]);

############################################################
################## Running the pipeline ####################
############################################################

print('Setting working folder ...');
setwd(output_path);

print('Downloading supplementary data ...');
getGEOSuppFiles(GSE_ID);

print('Obtaining metadata ...');
gse <- getGEO(GSE_ID, getGPL=FALSE);
metadata <- as.data.frame(pData(phenoData(gse[[1]])));
write.table(x=metadata, file=paste0('metadata_', GSE_ID, '.tsv'), quote=FALSE, sep='\t', row.names=FALSE);

## Afterwards if needed:
# tar -xvf GSE43976_RAW.tar
# gzip -d *.gz

print('The script finished correctly.');

#########################################################
########## End of the script ############################
#########################################################
