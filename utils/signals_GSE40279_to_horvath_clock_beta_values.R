###########################################################################################
#########                                                                         #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             14/08/2018                                  #########
#########                              EMBL-EBI                                   #########
#########                           Thornton group                                #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####              Biological insights into the epigenetic ageing clock           #########
###########################################################################################
##### Preprocess the raw signals from GSE40279 (NOOB + filtering + BMIQ) to obtain     ####
##### the beta values for the 353 sites of the Horvath's epigenetic clock.             ####
##### In this version of the script, if the BMIQ normalisation fails for a sample, this ###
##### sample is removed from the final output but the script does not fail.            ####
##### A file with QC metrics is also generated and can be used to filter out weird     ####
##### samples in downstream analyses (based on median intensities, predicted sex and   ####
##### failed BMIQ).                                                                    ####
##### This script is only useful for 27K and 450K data.                                ####
###########################################################################################
##### USAGE: manual                                                                    ####
###########################################################################################
##### NOTE: The annotation folder must contain:                                        ####
#####      - Cross-reactive probes from Chen et al.: cross_reactive_probes_Chen_2013.txt  #
#####      - Annotation for the epigenetic clock probes: AdditionalFile3.csv           ####
###########################################################################################


###########################################################
##################### Dependencies ########################
###########################################################

library(minfi);
library(RPMM);
library(IlluminaHumanMethylation27kmanifest);
library(data.table);


###########################################################
#####################  Arguments ##########################
###########################################################

print('Getting the input arguments ...');

args <- commandArgs(trailingOnly=TRUE);

## Input, output and annotation paths

path_to_raw_red <- "/nfs/research1/thornton/dem44/methylation_clock/raw_data_human/epigenetic_syndromes/GSE40279/raw_idat/Hannum_Red_Channel_noRS_paperSamples.csv";
path_to_raw_green <- "/nfs/research1/thornton/dem44/methylation_clock/raw_data_human/epigenetic_syndromes/GSE40279/raw_idat/Hannum_Green_Channel_noRS_paperSamples.csv";
output_path <- "/nfs/nobackup/thornton/dem44/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/blood_control/horvath_clock_beta_values";
ann_path <- '/nfs/research1/thornton/dem44/methylation_clock/raw_data_human/epigenetic_syndromes/utils/'

folder_name <- 'GSE40279';

cross_path <- paste0(ann_path, '/cross_reactive_probes_Chen_2013.txt');
sites_path <- paste0(ann_path, '/AdditionalFile3.csv');

if(!file.exists(cross_path)){
  stop('The file specifying the cross-reactive probes (cross_reactive_probes_Chen_2013.txt) could not be found.');
}

if(!file.exists(sites_path)){
  stop('The file specifying the probes from the epigenetic clock (AdditionalFile3.csv) could not be found.');
}

set.seed(1); # For reproducibility

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
#     It outputs a list where the first element is the matrix of beta-values and the second element the
#     matrix with the QC metrics.

## input_rg: RG input object
## bgcor: perform background correction and dye-bias equalization with Noob? (TRUE/FALSE)
## filter: filter out array probes with SNPs (interrogation site, single nucleotide extension), 
#          cross-reactive (Chen et al. 2013) or belonging to sex chromosome? (TRUE/FALSE)
# a_type: array type ('27K' or '450K').
## BMIQ: perform BMIQ normalisation? (TRUE/FALSE). Only needed for a_type='450K'.
## path_cr: path to the file that contains the probe IDs of the cross-reactive probes from Chen et al. 2013

preprocess_RG <- function(input_rg, bgcor, filter, a_type, BMIQ, path_cr=NA){
  
  ## Background correction and dye-bias normalisation.
  
  if(bgcor){
    
    print('Performing Noob background correction ...');
    input_ms <- preprocessNoob(input_rg); 
    
  }else{
    
    print('No background correction will be performed ...');
    input_ms <- preprocessRaw(input_rg);
    
  }
  
  print(paste0('At the start, there are ', length(input_ms@NAMES), ' probes.'));
  
  ## Perform QC and save QC information to be used in downstream analyses.
  
  qc_output <- as.data.frame(minfiQC(input_ms, fixOutliers=F, verbose=T)$qc);
  qc_output <- cbind(rownames(qc_output), qc_output);
  rownames(qc_output) <- NULL;
  colnames(qc_output)[1] <- 'SampleID';
  colnames(qc_output)[6] <- 'predictedSexQC';
  qc_output$BMIQ <- 'PASS';
  
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
    
    if(a_type=='450K'){
      ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19);
      sex_probes <- ann450k$Name[ann450k$chr %in% c("chrX","chrY")];
    }

    if(a_type=='27K'){
      ann27k <- getAnnotation(IlluminaHumanMethylation27kanno.ilmn12.hg19);
      sex_probes <- ann27k$Name[ann27k$chr %in% c("chrX","chrY")];
    }
    
    input_betas <- input_betas[!(rownames(input_betas) %in% sex_probes),];
    print(paste0('After removing probes in sex chromosomes, there are ', nrow(input_betas), ' probes left.'));
    
    
  }else{
    
    # Obtain beta-values.
    
    input_betas <- as.data.frame(getBeta(input_ms, type="Illumina"));
    print(paste0('After converting to beta-values, there are ', nrow(input_betas), ' probes left.'));
    
  }
  
  
  ## BMIQ normalisation (only 450K data).
  
  if(BMIQ & a_type=='450K'){
    
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

    if(length(failed_BMIQ_samples)>0){ # Remove samples that failed BMIQ and store this information in the QC matrix
      
      rm_indexes <- match(failed_BMIQ_samples, colnames(input_final_df))
      input_final_df <- input_final_df[,-rm_indexes];
      print(paste0('Number of samples that did not completed BMIQ: ', length(failed_BMIQ_samples)));
      print(paste0('Sample names: ', paste(failed_BMIQ_samples, collapse=',')));
      
      qc_output$BMIQ[match(failed_BMIQ_samples, qc_output$SampleID)] <- 'FAIL';
      
    }else{
      print('All the samples completed BMIQ.');
    }
    
    return(list(input_final_df, qc_output));
    
  }else{
    
    return(list(input_final_df, qc_output));
    
  }
}


############################################################
################## Running the pipeline ####################
############################################################

#### 1. Process the input dataset.

a_type<-'450K';

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


# Process.

print('Preprocessing RG data ...');

process_data <- preprocess_RG(input_rg=input_rg, bgcor=T, filter=T, a_type=a_type, BMIQ=T, path_cr=cross_path);

## Export the beta-values for the clock sites and the QC matrix.

clock_annotation <- fread(sites_path);
clock_probes <- clock_annotation$CpGmarker[-1];

NFB_clock_betas <- process_data[[1]][rownames(process_data[[1]]) %in% clock_probes,];
NFB_clock_betas <- cbind(rownames(NFB_clock_betas), NFB_clock_betas);
rownames(NFB_clock_betas) <- NULL;
colnames(NFB_clock_betas)[1] <- 'ProbeID';
QC_matrix <- process_data[[2]];

write.table(x=NFB_clock_betas, file=paste0(output_path, '/', folder_name, '_horvath_clock_beta_values.csv'),
            quote=F, sep=',', row.names=F);
write.table(x=QC_matrix, file=paste0(output_path, '/', folder_name, '_QC_matrix.csv'),
            quote=F, sep=',', row.names=F);

print('The script finished correctly.');

#########################################################
########## End of the script ############################
#########################################################
