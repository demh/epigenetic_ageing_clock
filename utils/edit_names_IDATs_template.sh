#!/bin/bash
###########################################################################################
#########                                                                         #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             29/06/2018                                  #########
#########                              EMBL-EBI                                   #########
#########                           Thornton group                                #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####              Biological insights into the epigenetic ageing clock           #########
###########################################################################################
##### Edit the names from a GSE batch to add an arbitrary sample name, necessary to    ####
##### make the rest of the scripts work.                                               ####
#####  e.g. we go from 200109350018_R01C01_Grn.idat to GSE91375-sample1_200109350018_R01C01_Grn.idat.
###########################################################################################
##### USAGE: manual                                                                    ####
###########################################################################################

i=1

cd /nfs/research1/thornton/dem44/methylation_clock/raw_data_human/epigenetic_syndromes/GSE91375/raw_idat/
root_names=($(ls | awk -F"_Red.idat" '{$0=$1}1' | awk -F"_Grn.idat" '{$0=$1}1' | sort | uniq))

for r in "${root_names[@]}"; do
	mv ${r}_Grn.idat GSE91375-sample${i}_${r}_Grn.idat
	mv ${r}_Red.idat GSE91375-sample${i}_${r}_Red.idat
	i=$((i + 1))
done

