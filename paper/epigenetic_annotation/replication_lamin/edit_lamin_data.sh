#!/bin/bash
###########################################################################################
#########                                                                         #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             05/10/2018                                  #########
#########                              EMBL-EBI                                   #########
#########                           Thornton group                                #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####              Biological insights into the epigenetic ageing clock           #########
###########################################################################################
##### Liftover the lamin B1 data from hg18 to hg19.                                   #####
###########################################################################################
##### USAGE: manual.                                                                  #####
###########################################################################################

bigWigToBedGraph GSM1289416_10032013_LaminB1-Input.OIS.bw GSM1289416_10032013_LaminB1-Input.OIS.bedGraph

liftOver GSM1289416_10032013_LaminB1-Input.OIS.bedGraph hg18ToHg19.over.chain.gz GSM1289416_10032013_LaminB1-Input.OIS_hg19.bedGraph GSM1289416_10032013_LaminB1-Input.OIS_unmapped.bedGraph

bedSort GSM1289416_10032013_LaminB1-Input.OIS_hg19.bedGraph GSM1289416_10032013_LaminB1-Input.OIS_hg19_sorted.bedGraph

bedtools merge -c 4 -o mean -i GSM1289416_10032013_LaminB1-Input.OIS_hg19_sorted.bedGraph > GSM1289416_10032013_LaminB1-Input.OIS_hg19_sorted_merged.bedGraph

bedGraphToBigWig GSM1289416_10032013_LaminB1-Input.OIS_hg19_sorted_merged.bedGraph hg19.chrom.sizes GSM1289416_10032013_LaminB1-Input.OIS_hg19_sorted_merged.bigWig

##### End of the script. #####
