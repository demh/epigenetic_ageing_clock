###########################################################################################
#########                                                                         #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             13/06/2018                                  #########
#########                              EMBL-EBI                                   #########
#########                           Thornton group                                #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####              Biological insights into the epigenetic ageing clock           #########
###########################################################################################
##### Extract the control probes from the raw signals in GSE40279.                     ####
##### The intensities from the control probes will be used to estimate the technical   ####
##### variation and correct for batch effects later on.                                ####
###########################################################################################
##### USAGE: manual                                                                    ####
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


## Input, output and annotation paths

path_to_raw_red <- "/nfs/research1/thornton/dem44/methylation_clock/raw_data_human/epigenetic_syndromes/GSE40279/raw_idat/Hannum_Red_Channel_noRS_paperSamples.csv";
path_to_raw_green <- "/nfs/research1/thornton/dem44/methylation_clock/raw_data_human/epigenetic_syndromes/GSE40279/raw_idat/Hannum_Green_Channel_noRS_paperSamples.csv";
output_path <- "/nfs/nobackup/thornton/dem44/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/blood_control/control_intensities/";

folder_name <- 'GSE40279';



############################################################
################## Running the pipeline ####################
############################################################

#### 1. Process the input dataset.

# Read the raw signals.

print('Reading the raw data ...');

raw_red <- read.csv(path_to_raw_red, header=TRUE,row.names=1);
mat_red <- as.matrix(raw_red);
colnames(mat_red) <- gsub('X', '', colnames(raw_red));
rm(raw_red);
raw_green <- read.csv(path_to_raw_green, header=TRUE,row.names=1);
mat_green <- as.matrix(raw_green);
colnames(mat_green) <- gsub('X', '', colnames(raw_green));
rm(raw_green);

# Creating the RG object.

print('Creating the RG object ...');
input_rg <- RGChannelSet(Green= mat_green,Red= mat_red,annotation = c(array="IlluminaHumanMethylation450k", annotation="ilmn12.hg19"));


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
