###########################################################################################
#########                                                                         #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             17/05/2018                                  #########
#########                              EMBL-EBI                                   #########
#########                           Thornton group                                #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####              Biological insights into the epigenetic ageing clock           #########
###########################################################################################
##### Comparison for the different reference-based cell-type deconvolution strategies. ####
##### Obtain the predictions for the different strategies.                             ####
###########################################################################################
##### USAGE: manual                                                                    ####
###########################################################################################

###########################################################
##################### Dependencies ########################
###########################################################

library(minfi);
library(FlowSorted.Blood.450k);
library(wateRmelon);
library(EpiDISH);
library(data.table);
set.seed(1);


############################################################
################## Functions ###############################
############################################################

##### Function: pipeline to preprocess the raw data (RG object) to beta-values in a customised way.

## input_rg: RG input object
## bgcor: perform background correction and dye-bias equalization with Noob? (TRUE/FALSE)
## filter: filter out array probes with SNPs (interrogation site, single nucleotide extension), 
#          cross-reactive (Chen et al. 2013) or belonging to sex chromosome? (TRUE/FALSE)
## BMIQ: perform BMIQ normalisation? (TRUE/FALSE).

preprocess_RG <- function(input_rg, bgcor, filter, BMIQ){
  
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
    
    cr_probes <- readLines('~/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/utils/cross_reactive_probes_Chen_2013.txt');
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

#### 1. Process the Reinius et al. 2012 dataset [REFERENCE]  

## Load as a RG object.

reinius_rg <- get('FlowSorted.Blood.450k');

## Preprocess the data.

reinius_final_df_NB <- preprocess_RG(input_rg=reinius_rg, bgcor=T, filter=F, BMIQ=T);  # Without filtering 
reinius_final_df_NFB <- preprocess_RG(input_rg=reinius_rg, bgcor=T, filter=T, BMIQ=T); # With filtering

## Create our own references (all the probes).

# Remove the samples with whole or peripheral blood. 

reinius_cells_NB <- as.matrix(reinius_final_df_NB[,-(c(grep('WB', colnames(reinius_final_df_NB)),grep('PBMC', colnames(reinius_final_df_NB))))]);
reinius_cells_NFB <- as.matrix(reinius_final_df_NFB[,-(c(grep('WB', colnames(reinius_final_df_NFB)),grep('PBMC', colnames(reinius_final_df_NFB))))]);

# Average the beta-values per cell-type.

cell_types <- unique(sapply(strsplit(colnames(reinius_cells_NFB), '_'), function(x){x[1]}));
reinius_averages_NB <- matrix(NA, ncol=length(cell_types), nrow=nrow(reinius_cells_NB));
reinius_averages_NFB <- matrix(NA, ncol=length(cell_types), nrow=nrow(reinius_cells_NFB));
colnames(reinius_averages_NB) <- colnames(reinius_averages_NFB) <- cell_types;

rownames(reinius_averages_NB) <- rownames(reinius_cells_NB);
rownames(reinius_averages_NFB) <- rownames(reinius_cells_NFB);

i <- 1;

for(ct in cell_types){
  c_name <- paste0('mean_', ct);
  reinius_averages_NB[,i] <- rowMeans2(reinius_cells_NB, cols=grep(ct,colnames(reinius_cells_NB)));
  reinius_averages_NFB[,i] <- rowMeans2(reinius_cells_NFB, cols=grep(ct,colnames(reinius_cells_NFB)));
  i <- i+1;
}

# Select only appropiate columns. 

reinius_averages_NB <- as.data.frame(reinius_averages_NB[,1:6]);
reinius_averages_NFB <- as.data.frame(reinius_averages_NFB[,1:6]);

colnames(reinius_averages_NB) <- colnames(reinius_averages_NFB) <- c('Gran', 'CD4T', 'CD8T', 'B', 'Mono', 'NK');


#### 2. Process dataset GSE77797 [GOLD-STANDARD]

## Classify the files into 27K / 450K. ##

path_to_gs <- '~/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/cell_composition/gold_standard/GSE77797/raw_idat/';

new_files_paths <- list.files(path_to_gs, full.names = TRUE, recursive=TRUE);
array_class <- sapply(new_files_paths, function(f){
  ifelse(file.size(f) < 800000, "27K", ifelse(file.size(f) < 9000000, "450K", NA));
});

if(sum(is.na(array_class)) > 0){
  stop("The array platform could not be identified in some samples !!");
}

## Create sample annotation files. ##

sample_ann_gs <- unique(data.frame(
  Array=as.character(sapply(names(which(array_class == "450K")), function(x){
    a <- strsplit(strsplit(x, '/')[[1]][length(strsplit(names(which(array_class == "450K")), '/')[[1]])], '_')[[1]][3];
    return(a);
  })),
  Slide=as.character(sapply(names(which(array_class == "450K")), function(x){
    s <- strsplit(x, '/')[[1]][length(strsplit(names(which(array_class == "450K")), '/')[[1]])-1];
    return(s);
  })),
  Basename=gsub('_Red.idat', '', gsub('_Grn.idat', '', names(which(array_class == "450K"))))));


## Process the goldstandard. 

# Read the IDAT files.

gs_rg <-read.metharray.exp(targets = sample_ann_gs);
names(colData(gs_rg))[4] <- 'Sample_Name'; # Fix bug in minfi::estimateCellCounts

# Process.

gs_NB_df <- preprocess_RG(input_rg=gs_rg, bgcor=T, filter=F, BMIQ=T);  # Without filtering 
gs_NFB_df <- preprocess_RG(input_rg=gs_rg, bgcor=T, filter=T, BMIQ=T); # With filtering


### Load the real cell-type values in the goldstandard. 

metadata_path <- '~/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/cell_composition/gold_standard/GSE77797/metadata/metadata_GSE77797_edited.txt';
raw_metadata <- as.data.frame(fread(file=metadata_path, verbose=F, showProgress=F));

real_cc_gs <- data.frame(Sample=as.character(raw_metadata$geo_accession), Gran=as.numeric(raw_metadata$`granulocyte (%):ch1`),
                         CD4T=as.numeric(raw_metadata$`cd4+ t cell (%):ch1`), CD8T=as.numeric(raw_metadata$`cd8+ t cell (%):ch1`),
                         B=as.numeric(raw_metadata$`b cell (%):ch1`), Mono=as.numeric(raw_metadata$`monocyte (%):ch1`), NK=as.numeric(raw_metadata$`natural killer cell (%):ch1`));


#### 3. Obtain the cell-type composition estimations using the different methodologies. 

## Predictions for minfi::estimateCellCounts function (default parameters).

minfi_pred <- estimateCellCounts(rgSet=gs_rg, referencePlatform = "IlluminaHumanMethylation450k");
minfi_pred <- as.data.frame(minfi_pred);
colnames(minfi_pred)[4] <- 'B';
minfi_pred <- minfi_pred[,colnames(reinius_averages_NB)];

## Create a dataframe to store all the results from the predictions. 

all_predictions <- data.frame(Name=rep('minfi', nrow(minfi_pred)),
                              Sample=rownames(minfi_pred));
all_predictions <- cbind(all_predictions, minfi_pred);
rownames(all_predictions) <- NULL;
all_predictions$Name <- as.character(all_predictions$Name);


## Predictions for centDHSbloodDMC.m (EpiDISH).

# Default reference.

data(centDHSbloodDMC.m);
cent_default <- as.data.frame(centDHSbloodDMC.m[,1:6]);

cent_default_NB <- cent_default[rownames(centDHSbloodDMC.m) %in% rownames(gs_NB_df),]; 
cent_default_NB <- cent_default_NB[,colnames(reinius_averages_NB)]; # Order columns as in reinius_averages 
cent_default_NB <- as.matrix(cent_default_NB[order(rownames(cent_default_NB)),]); # Order rows

cent_default_NFB <- cent_default[rownames(centDHSbloodDMC.m) %in% rownames(gs_NFB_df),]; # Take only the probes that survived our processing 
cent_default_NFB <- cent_default_NFB[,colnames(reinius_averages_NFB)]; # Order columns as in reinius_averages 
cent_default_NFB <- as.matrix(cent_default_NFB[order(rownames(cent_default_NFB)),]); # Order rows

# Our references (Noob+BMIQ and Noob+Filtering+BMIQ). 

cent_our_NB <- reinius_averages_NB[rownames(reinius_averages_NB)%in%rownames(centDHSbloodDMC.m),];
cent_our_NB <- as.matrix(cent_our_NB[order(rownames(cent_our_NB)),]); # Order rows

cent_our_NFB <- reinius_averages_NFB[rownames(reinius_averages_NFB)%in%rownames(centDHSbloodDMC.m),];
cent_our_NFB <- as.matrix(cent_our_NFB[order(rownames(cent_our_NFB)),]); # Order rows

# Goldstandard data (Noob+BMIQ and Noob+Filtering+BMIQ).

cent_gs_NB <- gs_NB_df[rownames(gs_NB_df) %in% rownames(cent_our_NB),];
cent_gs_NB <- as.matrix(cent_gs_NB[order(rownames(cent_gs_NB)),]);

cent_gs_NFB <- gs_NFB_df[rownames(gs_NFB_df) %in% rownames(cent_our_NFB),];
cent_gs_NFB <- as.matrix(cent_gs_NFB[order(rownames(cent_gs_NFB)),]);

# Predictions

epi_1 <- epidish(avdata.m=cent_gs_NB, ref.m=cent_default_NB, method="CP");
epi_2 <- epidish(avdata.m=cent_gs_NB, ref.m=cent_our_NB, method="CP");
epi_3 <- epidish(avdata.m=cent_gs_NFB, ref.m=cent_default_NFB, method="CP");
epi_4 <- epidish(avdata.m=cent_gs_NFB, ref.m=cent_our_NFB, method="CP");

epi_5 <- epidish(avdata.m=cent_gs_NB, ref.m=cent_default_NB, method="CBS");
epi_6 <- epidish(avdata.m=cent_gs_NB, ref.m=cent_our_NB, method="CBS");
epi_7 <- epidish(avdata.m=cent_gs_NFB, ref.m=cent_default_NFB, method="CBS");
epi_8 <- epidish(avdata.m=cent_gs_NFB, ref.m=cent_our_NFB, method="CBS");

epi_9 <- epidish(avdata.m=cent_gs_NB, ref.m=cent_default_NB, method="RPC");
epi_10 <- epidish(avdata.m=cent_gs_NB, ref.m=cent_our_NB, method="RPC");
epi_11 <- epidish(avdata.m=cent_gs_NFB, ref.m=cent_default_NFB, method="RPC");
epi_12 <- epidish(avdata.m=cent_gs_NFB, ref.m=cent_our_NFB, method="RPC");

all_epi_df <- data.frame(Name=rep(c('dhs_dif1_houseman', 'dhs_NB_houseman', 'dhs_dif2_houseman', 'dhs_NFB_houseman',	
                                    'dhs_dif1_cibersort', 'dhs_NB_cibersort', 'dhs_dif2_cibersort' ,'dhs_NFB_cibersort',
                                    'dhs_dif1_rpc', 'dhs_NB_rpc', 'dhs_dif2_rpc', 'dhs_NFB_rpc'), each=ncol(cent_gs_NB)),
                         Sample=c(rownames(epi_1$estF), rownames(epi_2$estF), rownames(epi_3$estF), rownames(epi_4$estF),
                                  rownames(epi_5$estF), rownames(epi_6$estF), rownames(epi_7$estF), rownames(epi_8$estF),
                                  rownames(epi_9$estF), rownames(epi_10$estF), rownames(epi_11$estF), rownames(epi_12$estF)));
all_epi_df <- cbind(all_epi_df, rbind(epi_1$estF, epi_2$estF, epi_3$estF, epi_4$estF,
                                      epi_5$estF, epi_6$estF, epi_7$estF, epi_8$estF,
                                      epi_9$estF, epi_10$estF, epi_11$estF, epi_12$estF));
rownames(all_epi_df) <- NULL;
all_predictions <- rbind(all_predictions, all_epi_df);


## Predictions for IDOL.

# Default reference: not available.

# Our references (Noob+BMIQ and Noob+Filtering+BMIQ). 

idol_probes <- as.character(read.csv('~/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/cell_composition/IDOL_Table_S3.csv')[,1]);

idol_our_NB <- reinius_averages_NB[rownames(reinius_averages_NB)%in%idol_probes,];
idol_our_NB <- as.matrix(idol_our_NB[order(rownames(idol_our_NB)),]); # Order rows

idol_our_NFB <- reinius_averages_NFB[rownames(reinius_averages_NFB)%in%idol_probes,];
idol_our_NFB <- as.matrix(idol_our_NFB[order(rownames(idol_our_NFB)),]); # Order rows

# Goldstandard data (Noob+BMIQ and Noob+Filtering+BMIQ).

idol_gs_NB <- gs_NB_df[rownames(gs_NB_df) %in% rownames(idol_our_NB),];
idol_gs_NB <- as.matrix(idol_gs_NB[order(rownames(idol_gs_NB)),]);

idol_gs_NFB <- gs_NFB_df[rownames(gs_NFB_df) %in% rownames(idol_our_NFB),];
idol_gs_NFB <- as.matrix(idol_gs_NFB[order(rownames(idol_gs_NFB)),]);

# Predictions

idol_1 <- epidish(avdata.m=idol_gs_NB, ref.m=idol_our_NB, method="CP");
idol_2 <- epidish(avdata.m=idol_gs_NFB, ref.m=idol_our_NFB, method="CP");

idol_3 <- epidish(avdata.m=idol_gs_NB, ref.m=idol_our_NB, method="CBS");
idol_4 <- epidish(avdata.m=idol_gs_NFB, ref.m=idol_our_NFB, method="CBS");

idol_5 <- epidish(avdata.m=idol_gs_NB, ref.m=idol_our_NB, method="RPC");
idol_6 <- epidish(avdata.m=idol_gs_NFB, ref.m=idol_our_NFB, method="RPC");

all_idol_df <- data.frame(Name=rep(c('idol_NB_houseman', 'idol_NFB_houseman', 
                                     'idol_NB_cibersort', 'idol_NFB_cibersort', 
                                     'idol_NB_rpc', 'idol_NFB_rpc'), each=ncol(idol_gs_NB)),
                         Sample=c(rownames(idol_1$estF), rownames(idol_2$estF), rownames(idol_3$estF),
                                  rownames(idol_4$estF), rownames(idol_5$estF), rownames(idol_6$estF)));
all_idol_df <- cbind(all_idol_df, rbind(idol_1$estF, idol_2$estF, idol_3$estF, idol_4$estF, idol_5$estF, idol_6$estF));
rownames(all_idol_df) <- NULL;

all_predictions <- rbind(all_predictions, all_idol_df);


## Edit the dataframe with all predictions.

all_predictions$Sample <- sapply(strsplit(as.character(all_predictions$Sample), '_'), function(x){x[1]});
all_predictions[,-c(1,2)] <- all_predictions[,-c(1,2)] * 100; # Display values as percentages


#### 4. Export the predictions and the actual gold-standard values.

write.table(x=all_predictions,
            file='~/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/cell_composition/compare_cc_strategies_predictions.csv',
            quote=F, sep=',', row.names=F, col.names=T);
write.table(x=real_cc_gs,
            file='~/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/cell_composition/compare_cc_strategies_real_values.csv',
            quote=F, sep=',', row.names=F, col.names=T);


############################################################
################## End of the script  ######################
############################################################
