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
##### In this version of the script, if the BMIQ normalisation fails for a sample, this ###
##### sample is removed from the final output but the script does not fail.            ####
###########################################################################################
##### USAGE: Rscript idat_to_cell_composition_blood.R path/to/idats /path/to/output /path/to/annotation
##### ARGS: path/to/idats: absolute path to the folder containing the raw IDAT data.   ####
#####                      The name of the batch should be in the folder previous to the final directory.
#####                      e.g. /data/GSE12345/raw_idat/ for the GSE12345 batch.       ####
#####       path/to/output: absolute path to the folder where the matrix will be outputted.
#####       path/to/annotation: absolute path to the folder which contains the cross-reactive
#####                           probes file and the reference.                         ####
###########################################################################################
##### NOTE1: the IDAT files must have a name with the following format: samplename_slide_array_channel.idat
#####        e.g. GSM3101870_9985131140_R01C01_Grn.idat                                ####
##### NOTE2: The annotation folder must contain:                                       ####
#####      - Cross-reactive probes from Chen et al.: cross_reactive_probes_Chen_2013.txt  #
#####      - Beta values for the processed IDOL reference: idol_NFB_reference.csv      ####
##### WARNING: the IDAT files will be organised inside the /data/batch/raw_idat/ folder ###
#####          into subfolders (see move_idat_files function).                         ####
###########################################################################################

###########################################################
##################### Dependencies ########################
###########################################################

library(minfi);
library(RPMM);
library(IlluminaHumanMethylation27kmanifest);
library(EpiDISH);
library(FlowSorted.Blood.450k);


###########################################################
#####################  Arguments ##########################
###########################################################

print('Getting the input arguments ...');

args <- commandArgs(trailingOnly=TRUE);

## Input, output and annotation paths

#path_to_raw_idat <- "~/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/cell_composition/debug/GSE51032/raw_idat/";
path_to_raw_idat <- args[1];
#output_path <- "~/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/cell_composition/debug/GSE51032/";
output_path <- args[2];
#ann_path <- '~/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/utils/'
ann_path <- args[3];

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


##### Function: BMIQ normalisation function (as implemented in the 'wateRmelon' package).

BMIQ <- function (beta.v, design.v, nL = 3, doH = TRUE, nfit = 50000, 
          th1.v = c(0.2, 0.75), th2.v = NULL, niter = 5, tol = 0.001, 
          plots = TRUE, sampleID = 1, pri = TRUE) 
{
  if (!library(RPMM, logical.return = TRUE, quietly = TRUE)) {
    stop("need RPMM package")
  }
  good <- !is.na(beta.v)
  out <- beta.v
  beta.v <- beta.v[good]
  design.v <- design.v[good]
  print <- function(x) {
    if (pri) 
      base::print(x)
  }
  type1.idx <- which(design.v == 1)
  type2.idx <- which(design.v == 2)
  beta1.v <- beta.v[type1.idx]
  beta2.v <- beta.v[type2.idx]
  if (min(beta1.v) == 0) {
    beta1.v[beta1.v == 0] <- min(setdiff(beta1.v, 0))
  }
  if (min(beta2.v) == 0) {
    beta2.v[beta2.v == 0] <- min(setdiff(beta2.v, 0))
  }
  if (max(beta1.v) == 1) {
    beta1.v[beta1.v == 1] <- max(setdiff(beta1.v, 1))
  }
  if (max(beta2.v) == 1) {
    beta2.v[beta2.v == 1] <- max(setdiff(beta2.v, 1))
  }
  w0.m <- matrix(0, nrow = length(beta1.v), ncol = nL)
  w0.m[which(beta1.v <= th1.v[1]), 1] <- 1
  w0.m[intersect(which(beta1.v > th1.v[1]), which(beta1.v <= 
                                                    th1.v[2])), 2] <- 1
  w0.m[which(beta1.v > th1.v[2]), 3] <- 1
  print("Fitting EM beta mixture to type1 probes")
  rand.idx <- sample(1:length(beta1.v), nfit, replace = FALSE)
  em1.o <- blc(matrix(beta1.v[rand.idx], ncol = 1), w = w0.m[rand.idx, 
                                                             ], maxiter = niter, tol = tol, verbose = pri)
  subsetclass1.v <- apply(em1.o$w, 1, which.max)
  subsetth1.v <- c(mean(max(beta1.v[rand.idx[subsetclass1.v == 
                                               1]]), min(beta1.v[rand.idx[subsetclass1.v == 2]])), mean(max(beta1.v[rand.idx[subsetclass1.v == 
                                                                                                                               2]]), min(beta1.v[rand.idx[subsetclass1.v == 3]])))
  class1.v <- rep(2, length(beta1.v))
  class1.v[which(beta1.v < subsetth1.v[1])] <- 1
  class1.v[which(beta1.v > subsetth1.v[2])] <- 3
  nth1.v <- subsetth1.v
  print("Done")
  if (plots) {
    print("Check")
    tmpL.v <- as.vector(rmultinom(1:nL, length(beta1.v), 
                                  prob = em1.o$eta))
    tmpB.v <- vector()
    for (l in 1:nL) {
      tmpB.v <- c(tmpB.v, rbeta(tmpL.v[l], em1.o$a[l, 1], 
                                em1.o$b[l, 1]))
    }
    pdf(paste("Type1fit-", sampleID, ".pdf", sep = ""), width = 6, 
        height = 4)
    plot(density(beta1.v))
    d.o <- density(tmpB.v)
    points(d.o$x, d.o$y, col = "green", type = "l")
    legend(x = 0.5, y = 3, legend = c("obs", "fit"), fill = c("black", 
                                                              "green"), bty = "n")
    dev.off()
  }
  d1U.o <- density(beta1.v[class1.v == 1])
  d1M.o <- density(beta1.v[class1.v == 3])
  mod1U <- d1U.o$x[which.max(d1U.o$y)]
  mod1M <- d1M.o$x[which.max(d1M.o$y)]
  d2U.o <- density(beta2.v[which(beta2.v < 0.4)])
  d2M.o <- density(beta2.v[which(beta2.v > 0.6)])
  mod2U <- d2U.o$x[which.max(d2U.o$y)]
  mod2M <- d2M.o$x[which.max(d2M.o$y)]
  th2.v <- vector()
  th2.v[1] <- nth1.v[1] + (mod2U - mod1U)
  th2.v[2] <- nth1.v[2] + (mod2M - mod1M)
  w0.m <- matrix(0, nrow = length(beta2.v), ncol = nL)
  w0.m[which(beta2.v <= th2.v[1]), 1] <- 1
  w0.m[intersect(which(beta2.v > th2.v[1]), which(beta2.v <= 
                                                    th2.v[2])), 2] <- 1
  w0.m[which(beta2.v > th2.v[2]), 3] <- 1
  print("Fitting EM beta mixture to type2 probes")
  rand.idx <- sample(1:length(beta1.v), nfit, replace = FALSE)
  em2.o <- blc(matrix(beta2.v[rand.idx], ncol = 1), w = w0.m[rand.idx, 
                                                             ], maxiter = niter, tol = tol, verbose = pri)
  print("Done")
  subsetclass2.v <- apply(em2.o$w, 1, which.max)
  subsetth2.v <- c(mean(max(beta2.v[rand.idx[subsetclass2.v == 
                                               1]]), min(beta2.v[rand.idx[subsetclass2.v == 2]])), mean(max(beta2.v[rand.idx[subsetclass2.v == 
                                                                                                                               2]]), min(beta2.v[rand.idx[subsetclass2.v == 3]])))
  class2.v <- rep(2, length(beta2.v))
  class2.v[which(beta2.v < subsetth2.v[1])] <- 1
  class2.v[which(beta2.v > subsetth2.v[2])] <- 3
  if (plots) {
    tmpL.v <- as.vector(rmultinom(1:nL, length(beta2.v), 
                                  prob = em2.o$eta))
    tmpB.v <- vector()
    for (lt in 1:nL) {
      tmpB.v <- c(tmpB.v, rbeta(tmpL.v[lt], em2.o$a[lt, 
                                                    1], em2.o$b[lt, 1]))
    }
    pdf(paste("Type2fit-", sampleID, ".pdf", sep = ""), width = 6, 
        height = 4)
    plot(density(beta2.v))
    d.o <- density(tmpB.v)
    points(d.o$x, d.o$y, col = "green", type = "l")
    legend(x = 0.5, y = 3, legend = c("obs", "fit"), fill = c("black", 
                                                              "green"), bty = "n")
    dev.off()
  }
  classAV1.v <- vector()
  classAV2.v <- vector()
  for (l in 1:nL) {
    classAV1.v[l] <- em1.o$mu[l, 1]
    classAV2.v[l] <- em2.o$mu[l, 1]
  }
  print("Start normalising type 2 probes")
  nbeta2.v <- beta2.v
  lt <- 1
  selU.idx <- which(class2.v == lt)
  selUR.idx <- selU.idx[which(beta2.v[selU.idx] > classAV2.v[lt])]
  selUL.idx <- selU.idx[which(beta2.v[selU.idx] < classAV2.v[lt])]
  p.v <- pbeta(beta2.v[selUR.idx], em2.o$a[lt, 1], em2.o$b[lt, 
                                                           1], lower.tail = FALSE)
  q.v <- qbeta(p.v, em1.o$a[lt, 1], em1.o$b[lt, 1], lower.tail = FALSE)
  nbeta2.v[selUR.idx] <- q.v
  p.v <- pbeta(beta2.v[selUL.idx], em2.o$a[lt, 1], em2.o$b[lt, 
                                                           1], lower.tail = TRUE)
  q.v <- qbeta(p.v, em1.o$a[lt, 1], em1.o$b[lt, 1], lower.tail = TRUE)
  nbeta2.v[selUL.idx] <- q.v
  lt <- 3
  selM.idx <- which(class2.v == lt)
  selMR.idx <- selM.idx[which(beta2.v[selM.idx] > classAV2.v[lt])]
  selML.idx <- selM.idx[which(beta2.v[selM.idx] < classAV2.v[lt])]
  p.v <- pbeta(beta2.v[selMR.idx], em2.o$a[lt, 1], em2.o$b[lt, 
                                                           1], lower.tail = FALSE)
  q.v <- qbeta(p.v, em1.o$a[lt, 1], em1.o$b[lt, 1], lower.tail = FALSE)
  nbeta2.v[selMR.idx] <- q.v
  if (doH) {
    lt <- 2
    selH.idx <- c(which(class2.v == lt), selML.idx)
    minH <- min(beta2.v[selH.idx])
    maxH <- max(beta2.v[selH.idx])
    deltaH <- maxH - minH
    deltaUH <- -max(beta2.v[selU.idx]) + min(beta2.v[selH.idx])
    deltaHM <- -max(beta2.v[selH.idx]) + min(beta2.v[selMR.idx])
    nmaxH <- min(nbeta2.v[selMR.idx]) - deltaHM
    nminH <- max(nbeta2.v[selU.idx]) + deltaUH
    ndeltaH <- nmaxH - nminH
    hf <- ndeltaH/deltaH
    nbeta2.v[selH.idx] <- nminH + hf * (beta2.v[selH.idx] - 
                                          minH)
  }
  pnbeta.v <- beta.v
  pnbeta.v[type1.idx] <- beta1.v
  pnbeta.v[type2.idx] <- nbeta2.v
  if (plots) {
    print("Generating final plot")
    d1.o <- density(beta1.v)
    d2.o <- density(beta2.v)
    d2n.o <- density(nbeta2.v)
    ymax <- max(d2.o$y, d1.o$y, d2n.o$y)
    pdf(paste("CheckBMIQ-", sampleID, ".pdf", sep = ""), 
        width = 6, height = 4)
    plot(density(beta2.v), type = "l", ylim = c(0, ymax), 
         xlim = c(0, 1))
    points(d1.o$x, d1.o$y, col = "red", type = "l")
    points(d2n.o$x, d2n.o$y, col = "blue", type = "l")
    legend(x = 0.5, y = ymax, legend = c("type1", "type2", 
                                         "type2-BMIQ"), bty = "n", fill = c("red", "black", 
                                                                            "blue"))
    dev.off()
  }
  print(paste("Finished for sample ", sampleID, sep = ""))
  out[good] <- pnbeta.v
  pnbeta.v <- out
  return(list(nbeta = pnbeta.v, class1 = class1.v, class2 = class2.v, 
              av1 = classAV1.v, av2 = classAV2.v, hf = hf, th1 = nth1.v, 
              th2 = th2.v))
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
    failed_BMIQ_samples <- c();
    input_final_df <- as.data.frame(matrix(NA, ncol=ncol(input_betas), nrow=nrow(input_betas)));
    rownames(input_final_df) <- rownames(input_betas);
    colnames(input_final_df) <- colnames(input_betas);
    
    for(i in 1:ncol(input_betas)){
      
      print(paste0('Processing sample ', colnames(input_final_df)[i], ' ...'));
      current_betas <- try(BMIQ(beta.v=input_betas[,i], design.v=design_vector, plots=F, pri=T));
      if(class(current_betas)=="try-error"){
        print(paste0('Sample ', colnames(input_final_df)[i], ' could not be BMIQ normalised.'));
        failed_BMIQ_samples <- c(colnames(input_final_df)[i], failed_BMIQ_samples);
      }else{
        input_final_df[,i] <- current_betas$nbeta;
      }
    }
    
    if(length(failed_BMIQ_samples)>0){ # Remove samples that failed BMIQ
      rm_indexes <- match(failed_BMIQ_samples, colnames(input_final_df))
      input_final_df <- input_final_df[,-rm_indexes];
      print(paste0('Number of samples that did not completed BMIQ: ', length(failed_BMIQ_samples)));
      print(paste0('Sample names: ', paste(failed_BMIQ_samples, collapse=',')));
    }else{
      print('All the samples completed BMIQ.');
    }
    
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

# Filter out those samples with NAs.

nas_in_input <- sum(is.na(idol_input_NFB));
print(paste0('There are ', nas_in_input, ' NAs in the input for the cc calculation.'));
if(nas_in_input > 0){
  cols_with_nas <- c();
  col_names_input <- colnames(idol_input_NFB);
  for(i in 1:ncol(idol_input_NFB)){
    if(sum(is.na(idol_input_NFB[,i]))>0){cols_with_nas <- c(cols_with_nas, i)};
  }
  print(paste0('Number of samples that had NAs for the reference probes: ', length(cols_with_nas)));
  print(paste0('Sample names: ', paste(colnames(idol_input_NFB)[cols_with_nas], collapse=',')));
  idol_input_NFB <- idol_input_NFB[,-cols_with_nas];
  colnames(idol_input_NFB) <- col_names_input[-cols_with_nas];
}
	
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

print('The script finished correctly.');

#########################################################
########## End of the script ############################
#########################################################
