# -*- coding: utf-8 -*-
###########################################################################################
#########                                                                         #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             01/11/2018                                  #########
#########                              EMBL-EBI                                   #########
#########                           Thornton group                                #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####              Biological insights into the epigenetic ageing clock           #########
###########################################################################################
##### Extract the scores from different bigWig files for a given set of genomic       #####
##### coordinates (e.g. the clock CpG sites). The different bigWig files can contain  #####
##### information for different epigenetic features (e.g. fold change over control from ###
##### a ChIP-Seq experiment for a given histone mark). A window (+- x bp) can be      #####
##### specified around the genomic coordinates and the averaged scored is reported.   #####
###########################################################################################
##### USAGE: python  extract_scores_from_bigwig.py metadata output                    ##### 
###########################################################################################

## Dependencies

import sys
import numpy
import pyBigWig
import pandas as pd


## Input arguments

# CpG sites to extract the information from. 
# CpGmarker,chr,pos
input_target = "/nfs/research1/thornton/dem44/methylation_clock/raw_data_human/epigenetic_syndromes/utils/all_450k_CpGs_coords_hg19.csv" # 1-based coordinates

# Metadata file with all the information for the epigenetic features. 
# File_ID,Feature_type,Data_type,Genome_assembly,Tissue,Age_years,Sex,File_path 
#input_metadata = "/nfs/research1/thornton/dem44/methylation_clock/raw_data_human/epigenetic_syndromes/epigenetic_annotation/ENCODE_PBMC_data/metadata_ENCODE_PBMC_FC_QCd.csv"
input_metadata = str(sys.argv[1])

# Output file
#output_file = '/nfs/research1/thornton/dem44/methylation_clock/raw_data_human/epigenetic_syndromes/epigenetic_annotation/ENCODE_PBMC_data/output_all_450k_CpGs_ENCODE_PBMC_FC_QCd_200bp.csv'
output_file = str(sys.argv[2])

# Window to consider when calculating the score for the epigenetic feature (+- bp around the CpG coordinate, score is averaged in this region)
window = 200


## Loop over the different epigenetic features to consider.

i = 1
for l in open(input_metadata):
    if(i > 1):        
        l_split = l.strip().split(',')
        input_bigwig = l_split[7]
        temp_feature = l_split[1]+'_'+l_split[0] 
        bw = pyBigWig.open(input_bigwig)
        print('Extracting scores for ' + temp_feature + ' ...')

        # Extract the scores for the current epigenetic feature for all the target sites.

        j = 1
        for line in open(input_target):    
            if(j > 1):    
                cols = line.strip().split(',')
                score = numpy.nanmean(bw.values(cols[1], int(cols[2]) - (1+window), int(cols[2]) + window))
                if(j == 2):
                    temp_df = pd.DataFrame(columns = ['CpGmarker', 'chr:coord_hg19', temp_feature])
                temp_df = temp_df.append({'CpGmarker': cols[0], 'chr:coord_hg19': cols[1]+':'+cols[2], temp_feature:score}, ignore_index=True)
            j+=1
        
        bw.close()

        # Merge with the dataframe containing all the epigenetic features.
            
        print('Merging with final dataframe ...')
        if (i == 2):
            final_df = temp_df
            del temp_df
        else:
            final_df = pd.merge(final_df, temp_df,  how='left', left_on=['CpGmarker', 'chr:coord_hg19'], right_on = ['CpGmarker', 'chr:coord_hg19'])  
    i += 1


## Export the final dataframe.

print('Creating final output file ...')
final_df.to_csv(output_file, index=False, index_label=False)


##### End of the script #####
