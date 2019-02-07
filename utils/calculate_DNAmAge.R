###########################################################################################
#########                                                                         #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             31/07/2017                                  #########
#########                              EMBL-EBI                                   #########
#########                           Thornton group                                #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####              Biological insights into the epigenetic ageing clock           #########
###########################################################################################
##### Calculate the DNA methylation age (DNAm age) using human Illumina methylation array #
##### data. This code is based on the tutorial provided by Steve Horvath in the original ##
##### publication.                                                                     ####
###########################################################################################
##### USAGE: Rscript calculate_DNAm_age.R --help                                       ####
###########################################################################################

###########################################################
##################### Dependencies ########################
###########################################################

suppressWarnings(suppressMessages(library(optparse)));
suppressWarnings(suppressMessages(library(RPMM)));
#suppressWarnings(suppressMessages(library(sqldf)));
suppressWarnings(suppressMessages(library(impute)));
suppressWarnings(suppressMessages(library(data.table)));


###########################################################
############## Command-line arguments #####################
###########################################################

## Create arguments and help in a Pythonic style.

option_list <-  list(
  
  make_option(c('-i', '--input_file'), type='character', default=NULL, 
              help="Absolute path to the input file which contains the DNA methylation information \
                for the different samples. See format in MethylationDataExample55.csv. COMPULSORY",
              metavar='character'),
  
  make_option(c('-o', '--output_path'), type='character', default=NULL,
              help="Absolute path for the output directory.  COMPULSORY",
              metavar='character'),
  
  make_option(c('-a', '--annotation_path'), type='character', default=NULL,
              help="Absolute path to the folder that contains all the necessary annotation files \
                i.e. probeAnnotation21kdatMethUsed.csv, datMiniAnnotation27k.csv and \ 
                AdditionalFile3.csv. COMPULSORY",
              metavar='character')
  
);

## Prepare option list and general description of the script.

opt_parser <- OptionParser(option_list=option_list,
                           description="\nCalculate the DNA methylation age (DNAm age) using human Illumina methylation array \ 
data. This code is based on the tutorial provided by Steve Horvath in the original \
publication.");        

opt <- parse_args(opt_parser);

## Check that all the compulsory command-line arguments are provided.

if(is.null(opt$input_file) | is.null(opt$output_path) | is.null(opt$annotation_path)){
  
  print_help(opt_parser);
  stop("Please, make sure that you have provided all the COMPULSORY arguments.", call.=FALSE);
  
}


## Paths

# Input DNA methylation data. Rows: CpG probes; Cols: ProbeID (same IDs as in annotation for clock sites), Samples 

input_DNA_methylation_path <- as.character(opt$input_file);
#input_DNA_methylation_path <- '~/Desktop/methylation_clock/C_T_hypo/calculate_Horvath_DNAm_age/new_tutorial_files/MethylationDataExample55.csv';
input_name <- strsplit(input_DNA_methylation_path, '/')[[1]][length(strsplit(input_DNA_methylation_path, '/')[[1]])];

# Output path

output_path <- as.character(opt$output_path);
#output_path <- '~/Desktop/methylation_clock/C_T_hypo/calculate_Horvath_DNAm_age/prepare_matrix/';
logfile_name <- paste0(output_path, '/LogFile_', gsub('.csv', '.txt', input_name));

# Annotation information regarding the probes that were used to train the model, 
# including the 'Gold Standard' values.

probeAnnotation21kdatMethUsed_path <- paste0(as.character(opt$annotation_path), "/probeAnnotation21kdatMethUsed.csv");
#probeAnnotation21kdatMethUsed_path <- "~/Desktop/methylation_clock/C_T_hypo/calculate_Horvath_DNAm_age/new_tutorial_files/probeAnnotation21kdatMethUsed.csv";

library(data.table)# Annotation information regarding all the probes in 27K.

probeAnnotation27k_path <- paste0(as.character(opt$annotation_path), "/datMiniAnnotation27k.csv");
#probeAnnotation27k_path <- "~/Desktop/methylation_clock/C_T_hypo/calculate_Horvath_DNAm_age/new_tutorial_files/datMiniAnnotation27k.csv";

# Clock sites annotation data.

datClock_path <- paste0(as.character(opt$annotation_path), "/AdditionalFile3.csv");
#datClock_path <- "~/Desktop/methylation_clock/C_T_hypo/calculate_Horvath_DNAm_age/new_tutorial_files/AdditionalFile3.csv";


## Normalisation

normalizeData <- TRUE; # Should the data be normalized to the goldstandard?

## Imputation
# Slow: using KNN form imputation package
# Fast: using gold standard values for missing CpG probes (WARNING: it comes from blood)

fastImputation <- FALSE; # Should fast imputation be forced independently of the number of NAs? 



################################################################
################## Functions ###################################
################################################################


## Function: correct the beta values from the samples so they are in the interval (0,1).
#
# datM: dataframe/matrix with Illumina beta values to normalise (rows:samples, columns: CpG probes)
# onlyIfOutside: if TRUE, calibration is only performed on those samples that have beta values which are <0 or >1 ?

CalibrateUnitInterval <- function(datM,onlyIfOutside=TRUE){
  
  rangeBySample <- data.frame(lapply(data.frame(t(datM)),range,na.rm=TRUE));
  minBySample <- as.numeric(rangeBySample[1,]);
  maxBySample <- as.numeric(rangeBySample[2,]);
  
  if (onlyIfOutside){ 
    
    indexSamples <- which((minBySample<0 | maxBySample>1) & !is.na(minBySample) & !is.na(maxBySample));
    
  }else{
    
    indexSamples <- 1:length(minBySample);
  }
  
  
  if (length(indexSamples)>=1) {
    
    for (i in indexSamples) {  # Calibrate the required samples using a linear model
      
      y1 <- c(0.001,0.999); 
      x1 <- c(minBySample[i],maxBySample[i]);
      lm1 <- lm( y1 ~ x1 );
      intercept1 <- coef(lm1)[[1]];
      slope1 <- coef(lm1)[[2]];
      datM[i,] <- intercept1+slope1*datM[i,];
    }
  }
  
  return(datM);
}


## Functions: these two functions are modified versions of the blc and betaEst functions from the RPMM package.
#  They were modified by Steve Horvath in order to make the code more robust (e.g. the optimization algorithm is 
#  changed to method="Nelder-Mead"). We have fixed a minor bug in bcl2 to avoid clashing with the WGCNA package.

blc2 <- function (Y, w, maxiter = 25, tol = 1e-06, weights = NULL, verbose = FALSE) 
{
  Ymn = min(Y[Y > 0], na.rm = TRUE)
  Ymx = max(Y[Y < 1], na.rm = TRUE)
  Y = base::pmax(Y, Ymn/2)
  Y = base::pmin(Y, 1 - (1 - Ymx)/2)  # Fix bug here: it was clashing with the pmin function from WGCNA
  Yobs = !is.na(Y)
  J = dim(Y)[2]
  K = dim(w)[2]
  n = dim(w)[1]
  if (n != dim(Y)[1]) 
    stop("Dimensions of w and Y do not agree")
  if (is.null(weights)) 
    weights = rep(1, n)
  mu = a = b = matrix(Inf, K, J)
  crit = Inf
  for (i in 1:maxiter) {
    warn0 = options()$warn
    options(warn = -1)
    eta = apply(weights * w, 2, sum)/sum(weights)
    mu0 = mu
    for (k in 1:K) {
      for (j in 1:J) {
        ab = betaEst2(Y[, j], w[, k], weights)
        a[k, j] = ab[1]
        b[k, j] = ab[2]
        mu[k, j] = ab[1]/sum(ab)
      }
    }
    ww = array(0, dim = c(n, J, K))
    for (k in 1:K) {
      for (j in 1:J) {
        ww[Yobs[, j], j, k] = dbeta(Y[Yobs[, j], j], 
                                    a[k, j], b[k, j], log = TRUE)
      }
    }
    options(warn = warn0)
    w = apply(ww, c(1, 3), sum, na.rm = TRUE)
    wmax = apply(w, 1, max)
    for (k in 1:K) w[, k] = w[, k] - wmax
    w = t(eta * t(exp(w)))
    like = apply(w, 1, sum)
    w = (1/like) * w
    llike = weights * (log(like) + wmax)
    crit = max(abs(mu - mu0))
    if (verbose) 
      print(crit)
    if (crit < tol) 
      break
  }
  return(list(a = a, b = b, eta = eta, mu = mu, w = w, llike = sum(llike)))
}

betaEst2=function (y, w, weights) 
{
  yobs = !is.na(y)
  if (sum(yobs) <= 1) 
    return(c(1, 1))
  y = y[yobs]
  w = w[yobs]
  weights = weights[yobs]
  N = sum(weights * w)
  p = sum(weights * w * y)/N
  v = sum(weights * w * y * y)/N - p * p
  logab = log(c(p, 1 - p)) + log(pmax(1e-06, p * (1 - p)/v - 
                                        1))
  if (sum(yobs) == 2) 
    return(exp(logab))
  opt = try(optim(logab, betaObjf, ydata = y, wdata = w, weights = weights, 
                  method = "Nelder-Mead",control=list(maxit=50) ), silent = TRUE)
  if (inherits(opt, "try-error")) 
    return(c(1, 1))
  exp(opt$par)
}


## Function: normalise the beta values previous to calculate DNAm age. Make distribution of beta values from sample look 
# more similar to distribution of beta values from gold standard. Only consider the 21368 probes used in the modelling.
# In original BMIQ: from typeII distribution to typeI distribution. In BMIQcalibration: from sample distribution to goldstandard distribution.

# datM: dataframe/matrix with Illumina beta values to normalise (rows:samples, columns: CpG probes)
# goldstandard.beta: numeric vector with Illumina beta values that are used as a gold standard for normalising. The length of goldstandard has to equal the number of columns of datM.
# nL: number of components in the beta-mixture model.
# doH: if TRUE, also correct the H (and ML) probes from the samples.
# nfit: maximum number of CpG probes that should be used when fitting the beta-mixture models
# th1.v: numeric vector which contains the two thresholds used to calculate the initial weight matrix (z) for the goldstandard CpG probes (i.e. whether they are U, H or M).
# niter: maximum number of EM iterations used when fitting the beta mixture models.
# tol: covergence tolerance used when fitting the beta mixture models.
# plots: if TRUE, several plots are created that ensure the user that the function is working correctly.
# calibrateUnitInterval: if TRUE, the beta values from the input data are corrected so they are in the interval (0,1).

BMIQcalibration <- function(datM,goldstandard.beta,nL=3,doH=TRUE,nfit=20000,th1.v=c(0.2,0.75),niter=5,tol=0.001,plots=FALSE,calibrateUnitInterval=TRUE){
  
  ### Checks before starting with normalisation.
  
  if (length(goldstandard.beta) != dim(datM)[[2]]){
    
    stop("Error in function BMIQcalibration: length(goldstandard.beta) !=dim(datM)[[2]]. Consider transposing datM.");
    
  }
  
  if (plots){
    
    par(mfrow=c(2,2));
    
  }
  
  ### Initialise some useful variables and potentially calibrate the input data.
  
  beta1.v <- goldstandard.beta;
  
  if (calibrateUnitInterval){
    
    datM <- CalibrateUnitInterval(datM);
    
  }
  
  
  ### A. Fitting the goldstandard probes
  
  ## Estimate initial weight matrix from goldstandard distribution (rows: CpG probes in goldstandard, columns: nL). See z in Ji et al., 2005
  
  w0.m <- matrix(0,nrow=length(beta1.v),ncol=nL);
  w0.m[which(beta1.v <= th1.v[1]),1] <- 1;
  w0.m[intersect(which(beta1.v > th1.v[1]),which(beta1.v <= th1.v[2])),2] <- 1;
  w0.m[which(beta1.v > th1.v[2]),3] <- 1;
  
  
  ## Fit beta mixture to goldstandard probes
  
  print("Fitting EM beta mixture to goldstandard probes ...");
  set.seed(1)
  rand.idx <- sample(1:length(beta1.v),min(c(nfit, length(beta1.v))), replace=FALSE);
  em1.o <- blc(Y=matrix(beta1.v[rand.idx],ncol=1), w=w0.m[rand.idx,], maxiter=niter, tol=tol, verbose=FALSE);
  
  ## Assign classes to goldstandard probes (U, H, M).
  
  subsetclass1.v <- apply(em1.o$w,1,which.max);
  subsetth1.v <- c(mean(max(beta1.v[rand.idx[subsetclass1.v==1]]),
                        min(beta1.v[rand.idx[subsetclass1.v==2]])),
                   mean(max(beta1.v[rand.idx[subsetclass1.v==2]]),
                        min(beta1.v[rand.idx[subsetclass1.v==3]]))); # New thresholds after fitting
  class1.v <- rep(2,length(beta1.v));
  class1.v[which(beta1.v < subsetth1.v[1])] <- 1;
  class1.v[which(beta1.v > subsetth1.v[2])] <- 3;
  nth1.v <- subsetth1.v;
  print("Done.");
  
  ## Generate plot from estimated mixture
  
  if(plots){
    
    print("Check plot !");
    tmpL.v <- as.vector(rmultinom(1,length(beta1.v),prob=em1.o$eta));
    tmpB.v <- vector(); # Sample from fitted model
    for(l in 1:nL){
      tmpB.v <- c(tmpB.v,rbeta(tmpL.v[l],em1.o$a[l,1],em1.o$b[l,1]));
    }
    plot(density(beta1.v),main= "Goldstandard fit"); # Real goldstandard data
    d.o <- density(tmpB.v);
    points(d.o$x,d.o$y,col="green",type="l");
    legend(x=0.5,y=3,legend=c("Observed goldstandard data","Fitted model"),fill=c("black","green"),bty="n");
    
  }
  
  
  ## Estimate modes of U and M goldstandard probes (needed to find th2.v thresholds) 
  
  if (sum(class1.v==1)==1) { mod1U <- beta1.v[class1.v==1]}
  
  if (sum(class1.v==3)==1) { mod1M <- beta1.v[class1.v==3]}
  
  if (sum(class1.v==1)>1){
    
    d1U.o <- density(beta1.v[class1.v==1]);
    mod1U <- d1U.o$x[which.max(d1U.o$y)];
    
  }
  
  if (sum(class1.v==3)>1){
    
    d1M.o <- density(beta1.v[class1.v==3]);
    mod1M <- d1M.o$x[which.max(d1M.o$y)];
    
  }
  
  
  ### B. Fitting the samples probes.
  
  for (ii in 1:dim(datM)[[1]]){
    
    print(paste0("Fitting EM beta mixture to input probes from sample ",ii, " out of ", dim(datM)[[1]], ' ...'));
    sampleID <- ii;
    beta2.v <- as.numeric(datM[ii,]);
    
    ## Estimate modes of U and M probes for the current sample
    
    d2U.o <- density(beta2.v[which(beta2.v<0.4)]);
    d2M.o <-  density(beta2.v[which(beta2.v>0.6)]);
    mod2U <- d2U.o$x[which.max(d2U.o$y)];
    mod2M <- d2M.o$x[which.max(d2M.o$y)];
    
    ## Calculate thresholds for the current sample
    
    th2.v <- vector();
    th2.v[1] <- nth1.v[1] + (mod2U-mod1U);
    th2.v[2] <- nth1.v[2] + (mod2M-mod1M);
    
    ## Estimate initial weight matrix from sample distribution (rows: CpG probes in sample, columns: nL). See z in Ji et al., 2005
    
    w0.m <- matrix(0,nrow=length(beta2.v),ncol=nL);
    w0.m[which(beta2.v <= th2.v[1]),1] <- 1;
    w0.m[intersect(which(beta2.v > th2.v[1]),which(beta2.v <= th2.v[2])),2] <- 1;
    w0.m[which(beta2.v > th2.v[2]),3] <- 1;
    
    ## Fit beta mixture to sample probes
    
    set.seed(1)
    rand.idx <- sample(1:length(beta2.v),min(c(nfit, length(beta2.v)),na.rm=TRUE),replace=FALSE);
    em2.o <- blc2(Y=matrix(beta2.v[rand.idx],ncol=1),w=w0.m[rand.idx,],maxiter=niter,tol=tol);
    print("Done");
    
    
    ## Assign classes to sample probes (U, H, M).
    
    subsetclass2.v <- apply(em2.o$w,1,which.max);
    
    if (sum(subsetclass2.v==2)>0){
      
      subsetth2.v = c(mean(max(beta2.v[rand.idx[subsetclass2.v==1]]),
                           min(beta2.v[rand.idx[subsetclass2.v==2]])),
                      mean(max(beta2.v[rand.idx[subsetclass2.v==2]]),
                           min(beta2.v[rand.idx[subsetclass2.v==3]])));
      
    }else{
      
      subsetth2.v = c(1/2*max(beta2.v[rand.idx[subsetclass2.v==1]])+ 1/2*mean(beta2.v[rand.idx[subsetclass2.v==3]]), 
                      1/3*max(beta2.v[rand.idx[subsetclass2.v==1]])+ 2/3*mean(beta2.v[rand.idx[subsetclass2.v==3]]));
      
    }
    
    class2.v <- rep(2,length(beta2.v));
    class2.v[which(beta2.v <= subsetth2.v[1])] <- 1;
    class2.v[which(beta2.v >= subsetth2.v[2])] <- 3;
    
    
    ## Generate plot from estimated mixture
    
    if(plots){
      
      tmpL.v <- as.vector(rmultinom(1,length(beta2.v),prob=em2.o$eta));
      tmpB.v <- vector(); # Sample from fitted model
      for(lt in 1:nL){
        tmpB.v <- c(tmpB.v,rbeta(tmpL.v[lt],em2.o$a[lt,1],em2.o$b[lt,1]));
      }
      plot(density(beta2.v), main= paste0("Real data fit. Sample: ",sampleID)); # Real data
      d.o <- density(tmpB.v);
      points(d.o$x,d.o$y,col="green",type="l"); 
      legend(x=0.5,y=3,legend=c("Observed real data","Fitted model"),fill=c("black","green"),bty="n");
      
    }
    
    
    ## Store the means (averages) for the 3 fitted beta distributions.
    
    classAV1.v <- vector();
    classAV2.v <- vector();
    
    for(l in 1:nL){
      classAV1.v[l] <-  em1.o$mu[l,1];
      classAV2.v[l] <-  em2.o$mu[l,1];
    }
    
    
    ### C. Normalise the probes from the input samples according to the gold standard
    
    print("Normalising the input probes ...");
    nbeta2.v <- beta2.v;
    
    ## U probes
    
    # Classify in UR and UL based on their position to the U mean
    
    lt <- 1;
    selU.idx <- which(class2.v==lt);
    selUR.idx <- selU.idx[which(beta2.v[selU.idx] > classAV2.v[lt])];
    selUL.idx <- selU.idx[which(beta2.v[selU.idx] < classAV2.v[lt])];
    
    # Transform UR
    
    p.v <- pbeta(beta2.v[selUR.idx],em2.o$a[lt,1],em2.o$b[lt,1],lower.tail=FALSE);
    q.v <- qbeta(p.v,em1.o$a[lt,1],em1.o$b[lt,1],lower.tail=FALSE);
    nbeta2.v[selUR.idx] <- q.v;
    
    # Transform UL
    
    p.v <- pbeta(beta2.v[selUL.idx],em2.o$a[lt,1],em2.o$b[lt,1],lower.tail=TRUE);
    q.v <- qbeta(p.v,em1.o$a[lt,1],em1.o$b[lt,1],lower.tail=TRUE);
    nbeta2.v[selUL.idx] <- q.v;
    
    
    ## M probes
    
    # Classify in MR and ML based on their position to the M mean
    
    lt <- 3;
    selM.idx <- which(class2.v==lt);
    selMR.idx <- selM.idx[which(beta2.v[selM.idx] > classAV2.v[lt])];
    selML.idx <- selM.idx[which(beta2.v[selM.idx] < classAV2.v[lt])];
    
    # Transform MR
    
    p.v <- pbeta(beta2.v[selMR.idx],em2.o$a[lt,1],em2.o$b[lt,1],lower.tail=FALSE);
    q.v <- qbeta(p.v,em1.o$a[lt,1],em1.o$b[lt,1],lower.tail=FALSE);
    nbeta2.v[selMR.idx] <- q.v;
    
    
    ## H probes
    
    if(doH){
      
      # Select H and ML probes (left ML tail is not well described by a beta-distribution).
      
      lt <- 2;
      selH.idx <- c(which(class2.v==lt),selML.idx);
      
      # Calculate necessary variables
      
      minH <- min(beta2.v[selH.idx],na.rm=TRUE);
      maxH <- max(beta2.v[selH.idx],na.rm=TRUE);
      deltaH <- maxH - minH;
      deltaUH <- minH - max(beta2.v[selU.idx],na.rm=TRUE);
      deltaHM <- min(beta2.v[selMR.idx],na.rm=TRUE) - max(beta2.v[selH.idx],na.rm=TRUE); 
      
      nmaxH <- min(nbeta2.v[selMR.idx],na.rm=TRUE) - deltaHM;
      nminH <- max(nbeta2.v[selU.idx],na.rm=TRUE) + deltaUH;
      ndeltaH <- nmaxH - nminH;
      
      # Perform conformal transformation (shift+dilation)
      
      hf <- ndeltaH/deltaH;
      nbeta2.v[selH.idx] <- nminH + hf*(beta2.v[selH.idx]-minH);
      
    }
    
    print('Done.');
    
    ## Generate plot to test normalisation
    
    if(plots){
      
      print("Check plot !");
      d1.o <- density(beta1.v);
      d2.o <- density(beta2.v);
      d2n.o <- density(nbeta2.v);
      ymax <- max(d2.o$y,d1.o$y,d2n.o$y);
      plot(d2.o,type="l",ylim=c(0,ymax),xlim=c(0,1), main=paste0("Check BMIQ. Sample: ",sampleID)); 
      points(d1.o$x,d1.o$y,col="red",type="l");
      points(d2n.o$x,d2n.o$y,col="blue",type="l");
      legend(x=0.5,y=ymax,legend=c("Observed goldstandard","Observed real data","Normalised real data"),bty="n",fill=c("red","black","blue"));
      
    }
    
    ## Update final normalised results
    
    datM[ii,] <- nbeta2.v;
    
  }
  
  return(datM); # Return the normalised data
  
}


## Function: transform chronological age before performing linear regression
#
# c: chronological age of the sample (in years)
# a: adult age (for humans, 20 years)

F_transf <- function(c, a=20){
  
  if(c <= a){
    
    F_c <- log((c+1)/(a+1));
    
  }
  
  if(c > a){
    
    F_c <- (c-a) / (a+1);
    
  }
  
  return(F_c);
  
}


## Function: inverse of the transforming function
#
# linear_part: linear part of the regression model (i.e. 
#             weighted average methylation based on sample data)
# a: adult age (for humans, 20 years)

F_inverse_transf <- function(linear_part, a=20){
  
  if(linear_part <= 0){
    
    DNAm_age <- exp(linear_part) * (a+1) - 1;
    
  }
  
  if(linear_part > 0) {
    
    DNAm_age <- linear_part * (a+1) + a;
    
  }
  
  return(DNAm_age);
  
}


##########################################################################
####### 1. Read necessary information ####################################
##########################################################################


## Read the annotation data.

# Annotation information regarding the probes that were used to train the model, 
# including the 'Gold Standard' values.

probeAnnotation21kdatMethUsed <- read.csv(probeAnnotation21kdatMethUsed_path);

# Annotation information regarding all the probes in 27K.

probeAnnotation27k <- read.csv(probeAnnotation27k_path);

# Clock sites annotation data.

datClock <- read.csv(datClock_path);


## Read the input data (beta values from different samples).
#  Rows: CpG probes
#  Columns: samples

#dat0 <- read.csv.sql(input_DNA_methylation_path);
#dat0[,1] <- gsub(x=dat0[,1], pattern="\"", replacement=""); # Apparently needed when using read.csv.sql function
dat0 <- as.data.frame(fread(input=input_DNA_methylation_path));

# Some useful variables.

nSamples <- dim(dat0)[2] - 1;
nProbes <- dim(dat0)[1];



###############################################################################
###### 2. Create a LOG file to check for errors.  #############################
###############################################################################

### Remove previous LOG file.

if(file.exists(logfile_name)){
  
  file.remove(logfile_name);
  
}


### Create empty LOG file.

file.create(logfile_name);


### Fill in the LOG file.

DoNotProceed=FALSE; # Should we stop the execution of the script?

cat(paste0("The methylation data set contains ", nSamples, " samples (e.g. arrays) and ", 
          nProbes, " probes."), file=logfile_name);

## Several tests.

if(nSamples==0){ 
  
  DoNotProceed=TRUE; 
  cat(paste( "\n ERROR: There must be a data input error since there seem to be no samples.\n",
             "Make sure that you input a comma delimited file (.csv file)\n that can be read",
             "using the R command read.csv.sql . Samples correspond to columns in that file."), 
      file=logfile_name, append=TRUE); 

} 

if(nProbes==0){
  
  DoNotProceed=TRUE; 
  cat(paste( "\n ERROR: There must be a data input error since there seem to be zero probes.\n",
             "Make sure that you input a comma delimited file (.csv file)\n that can be read",
             "using the R command read.csv.sql  CpGs correspond to rows."), 
      file=logfile_name,append=TRUE); 

}

if(nSamples > nProbes){ 
  
  cat(paste( "\n MAJOR WARNING: It worries me a lot that there are more samples than CpG probes.\n",
             "Make sure that probes correspond to rows and samples to columns.\n I wonder whether",
             "you want to first transpose the data and then resubmit them? In any event, I will",
             "proceed with the analysis."),
      file=logfile_name,append=TRUE); 

}

if(is.numeric(dat0[,1])){ 
  
  DoNotProceed=TRUE; 
  cat(paste( "\n ERROR: The first column does not seem to contain probe identifiers (cg numbers from Illumina)",
             "since these entries are numeric values. Make sure that the first column of the file contains probe",
             "identifiers such as cg00000292. Instead it contains ", dat0[1:3,1]),
      file=logfile_name,append=TRUE);  
  
}

if(!is.character(dat0[,1])){  
  
  cat(paste( "\n MAJOR WARNING: The first column does not seem to contain probe identifiers (cg numbers from Illumina)",
             "since these entries are not character values. Make sure that the first column of the file contains CpG probe",
             "identifiers such as cg00000292. Instead it contains ", dat0[1:3,1]  ),
      file=logfile_name,append=TRUE);  
  
}


## Identify those probes, used to train the model, that are missing in our input data. If there are probes
# missing: stop the script.

match1 <- match(probeAnnotation21kdatMethUsed$Name, dat0[,1]);

if(sum(is.na(match1))>0){
  
  missingProbes <- probeAnnotation21kdatMethUsed$Name[!is.element(probeAnnotation21kdatMethUsed$Name, dat0[,1])];    
  DoNotProceed <- TRUE; 
  cat(paste( "\n \n ERROR: You forgot to include the following ", length(missingProbes), 
             " CpG probes (or probe names):\n ", 
             paste( missingProbes, sep="",collapse=", ")),
      file=logfile_name,append=TRUE);
} 


## Identify those samples which contain non-numeric data for the beta-values

nonNumericColumn <- rep(FALSE, nSamples);

for(i in 2:dim(dat0)[2]){ 
  
  nonNumericColumn[i-1] <- !is.numeric(dat0[,i]);
  
}

if(sum(nonNumericColumn) > 0){
  
  cat(paste("\n MAJOR WARNING: Possible input error. The following samples contain",
            "non-numeric beta values: ", colnames(dat0)[-1][nonNumericColumn], 
            "\n Hint: Maybe you use the wrong symbols for missing data. Make sure",
            "to code missing values as NA in the Excel file. To proceed, I will force",
            "the entries into numeric values but make sure this makes sense.\n"),
      file=logfile_name,append=TRUE);
  
}


## Check whether we will be able to predict the gender of the samples.

# Some calculations regarding the probes in X chromosome.

XchromosomalCpGs <- as.character(probeAnnotation27k$Name[probeAnnotation27k$Chr=="X"]);
selectXchromosome <- is.element(dat0[,1], XchromosomalCpGs);
selectXchromosome[is.na(selectXchromosome)] <- FALSE;
meanXchromosome <- rep(NA, nSamples);

if(sum(selectXchromosome) >= 500){
  
  # Obtain the mean beta-value accross X-chromosome probes for each sample.
  
  meanXchromosome <- as.numeric(apply(as.matrix(dat0[selectXchromosome,-1]),2,mean,na.rm=TRUE)); 
  
}

if(sum(is.na(meanXchromosome)) > 0){
  
  cat(paste( "\n \n COMMENT: There are lots of missing values for X chromosomal probes for some",
             "of the samples. This is not a problem when it comes to estimating age but I cannot",
             "predict the gender of these samples.\n " ),
      file=logfile_name,append=TRUE);  
  
}



############################################################################
###### 3. If no errors were found, proceed with the rest of the script #####
############################################################################

if(!DoNotProceed){ # If there were no errors so far
  
  print('Proceed with the analysis ...');
  
  ### 3.1. Extract the 21K probes used in the model from the input data and make them numeric. 
  
  dat1 <- dat0[match1,];
  asnumeric1 <- function(x){as.numeric(as.character(x))};
  dat1[,-1] <- apply(as.matrix(dat1[,-1]),2,asnumeric1);
  

  ### 3.2. Define quality metrics.
  
  # Quality metrics.
  
  meanMethBySample <- as.numeric(apply(as.matrix(dat1[,-1]),2,mean,na.rm=TRUE));
  minMethBySample <- as.numeric(apply(as.matrix(dat1[,-1]),2,min,na.rm=TRUE));
  maxMethBySample  <- as.numeric(apply(as.matrix(dat1[,-1]),2,max,na.rm=TRUE));
  
  # Transpose the array data.
  
  datMethUsed <- t(dat1[,-1]); # Rows: samples; Columns: CpG probes
  colnames(datMethUsed) <- as.character(dat1[,1]);
  
  # Number of missing values (NA) in the array data.
  
  noMissingPerSample <- apply(as.matrix(is.na(datMethUsed)),1,sum);
  table(noMissingPerSample);
  
  
  ### 3.3. Imputation.
  
  set.seed(1); # Make analysis reproducible
  
  ## Slow imputation (using imputation package) - when there are < 3000 missing values
  
  if (!fastImputation & nSamples>1 & max(noMissingPerSample,na.rm=TRUE)<3000){
    
    if (max(noMissingPerSample,na.rm=TRUE)>0){ # If there is at least one value missing
      
      print('Running slow imputation ...');
      
      dimnames1 <- dimnames(datMethUsed);
      datMethUsed <- data.frame(t(impute.knn(t(datMethUsed))$data));
      dimnames(datMethUsed) <- dimnames1;
      
    }else{
      
      print('There are no missing values and no imputation is needed.');
      
    }
  } 
  
  
  ## Fast imputation (taking gold standard value for the CpG probe)
  
  if (max(noMissingPerSample,na.rm=TRUE)>=3000) { # If there are many missing values
    
    fastImputation <- TRUE;
    normalizeData <- FALSE;
    
  }
  
  if(fastImputation | nSamples==1){
    
    if (max(noMissingPerSample,na.rm=TRUE)>0){ # If there is at least one missing
      
      print('Running fast imputation ...')
      
      for (i in which(noMissingPerSample>0) ){
        
        selectMissing1 <- is.na(datMethUsed[i,]);
        datMethUsed[i,selectMissing1] <- as.numeric(probeAnnotation21kdatMethUsed$goldstandard2[selectMissing1]);
        
      }
    }
  } 
  
  
  
  ### 3.4. Data normalization. 

  if(normalizeData){
    
    datMethUsedNormalized <- BMIQcalibration(datM=datMethUsed,
                                             goldstandard.beta= probeAnnotation21kdatMethUsed$goldstandard2,
                                             plots=FALSE);
  }else{ 
    
    datMethUsedNormalized <- datMethUsed;
  
  }
  
  # Remove unused objects
  
  rm(datMethUsed); rm(dat0);

  
  
  ### 3.5. Predict DNAm age 
  
  print('Calculating DNAm age ...');
  
  selectCpGsClock <- is.element(dimnames(datMethUsedNormalized)[[2]], as.character(datClock$CpGmarker[-1]));
  
  if (sum(selectCpGsClock) < dim(datClock)[[1]]-1) {
    
    stop("The CpGs listed in column 1 of the input data did not contain the CpGs needed for calculating DNAm age. Make sure to input cg numbers such as cg00075967.");
  
  }
  
  if (sum(selectCpGsClock) > dim(datClock)[[1]]-1 ){
    
    stop("ERROR: The CpGs listed in column 1 of the input data contain duplicate CpGs. Each row should report only one unique CpG marker (cg number).");
  
  }
  
  
  if(nSamples > 1){
    
    datMethClock0 <- data.frame(datMethUsedNormalized[,selectCpGsClock]);
    datMethClock <- data.frame(datMethClock0[as.character(datClock$CpGmarker[-1])]); # Same order as in datClock

  }
  
  if(nSamples == 1){
    
    datMethClock0 <- t(data.frame(datMethUsedNormalized[,selectCpGsClock]));
    datMethClock <- t(data.frame(datMethClock0[,as.character(datClock$CpGmarker[-1])])); # Same order as in datClock
    
  }
 
  linear_part_results <- as.numeric(datClock$CoefficientTraining[1] + as.matrix(datMethClock) %*% as.numeric(datClock$CoefficientTraining[-1]));
  predictedAge <- as.numeric(sapply(linear_part_results, F_inverse_transf));
  
  
  ### 3.6. Create the output
  
  print('Creating final output ...');
  
  # Add some comments to the samples if needed
  
  Comment <- ifelse(predictedAge < 0, "Negative DNAm age", ifelse (predictedAge > 100, "Old DNAm age", rep("",length(predictedAge))));
  Comment[is.na(predictedAge)] <- "Age prediction was not possible";

  restSamples <- minMethBySample < -0.05 | maxMethBySample > 1.05;
  restSamples[is.na(restSamples)] <- FALSE;
  lab1 <- "MAJOR WARNING: Probably you did not input beta values since either minMethBySample < -0.05 or maxMethBySample > 1.05.";
  Comment[restSamples] <- paste(Comment[restSamples],lab1);
  
  restSamples <- noMissingPerSample > 0 & noMissingPerSample <= 100;
  lab1 <- "WARNING: Some beta values were missing, see noMissingPerSample."; 
  Comment[restSamples] <- paste(Comment[restSamples],lab1);
  restSamples <- noMissingPerSample > 100 & noMissingPerSample <= 3000;
  lab1 <- "MAJOR WARNING: noMissingPerSample>100";
  Comment[restSamples] <- paste(Comment[restSamples],lab1);
  restSamples <- noMissingPerSample > 3000;
  lab1 <- "MAJOR WARNING: More than 3k missing values!!"; 
  Comment[restSamples] <- paste(Comment[restSamples],lab1);
  
  restSamples <- meanMethBySample > 0.35;
  restSamples[is.na(restSamples)] <- FALSE;
  lab1 <- "WARNING: meanMethBySample is > 0.35";
  Comment[restSamples] <- paste(Comment[restSamples],lab1);
  restSamples <- meanMethBySample < 0.25;
  restSamples[is.na(restSamples)] <- FALSE; 
  lab1 <- "WARNING: meanMethBySample is < 0.25";
  Comment[restSamples] <- paste(Comment[restSamples],lab1);
  
  # Store the output so far
  
  datout <- data.frame(SampleID=colnames(dat1)[-1], DNAmAge=predictedAge, Comment, noMissingPerSample, meanMethBySample, minMethBySample, maxMethBySample);
  
  
  # Include information regarding gender prediction based on X chromosome probes
  
  if (!all(is.na(meanXchromosome))){
    
      predictedGender <- ifelse(meanXchromosome > 0.4,"Female",
                             ifelse(meanXchromosome < 0.38, "Male","Unsure"));
      datout <- data.frame(datout,predictedGender=predictedGender,meanXchromosome=meanXchromosome);
  }
  
  
  # Create the output
  
  if(sum(datout$Comment != "")==0) {
    
    cat(paste( "\nThe individual samples appear to be fine. "),
        file=logfile_name,append=TRUE);
  }
  
  if(sum(datout$Comment != "")>0){ 
    
    cat(paste( "\nWarnings were generated for the following samples: \n", 
               datout[,1][datout$Comment != ""], 
               "\nCheck the output file for more details."),
        file=logfile_name,append=TRUE);
  
  }
  
  output_name <- paste0(output_path, "/DNAm_age_of_", input_name);
  write.table(datout, output_name, row.names=F, sep=",", quote = FALSE);
  print('The script finished correctly.');
  
  
}else{
  
  print('Do not proceed with the analysis. Check the LOG file for more information.');
  print('No output file will be created.');
  print('The script did not finish correctly.');
  
}

#### End of the script ####
