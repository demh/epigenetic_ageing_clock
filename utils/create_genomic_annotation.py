# -*- coding: utf-8 -*-
###########################################################################################
#########                                                                         #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             02/11/2018                                  #########
#########                              EMBL-EBI                                   #########
#########                           Thornton group                                #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####              Biological insights into the epigenetic ageing clock           #########
###########################################################################################
##### Extract categorical genome annotation for a given set of genomic coordinates    ##### 
##### (e.g. the clock CpG sites).                                                     #####
###########################################################################################
##### USAGE: manual                                                                   #####
###########################################################################################

## Dependencies

import pybedtools
from pybedtools import BedTool
import pandas as pd

## Input arguments

path_to_chr_lengths = "/Users/dem44/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/epigenetic_annotation/genomic_annotation/hg19.chrom.sizes"

path_to_gencode_ann = "/Users/dem44/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/epigenetic_annotation/genomic_annotation/gencode.v29lift37.basic.annotation.gtf"
path_to_target = "/Users/dem44/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/epigenetic_annotation/all_21k_CpGs_coords_hg19.csv" # 21K or 450K
path_to_CGI = "/Users/dem44/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/epigenetic_annotation/genomic_annotation/CGI_hg19.gtf"
path_to_ChrHMM = "/Users/dem44/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/epigenetic_annotation/genomic_annotation/hg19_chromHMM_imputed25"

output_file = "/Users/dem44/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/epigenetic_annotation/genomic_annotation/genomic_annotation_21k.csv" # 21K or 450K 


## Create BedTool object for the target CpG sites.

print('Creating BedTool object ...')

i = 1
gff_string = ''
for l in open(path_to_target):
    if(i > 1):
        l_split = l.strip().split(',')
        temp_string = l_split[1] + '\t' + 'na' + '\t' + 'na' + '\t' + l_split[2] + '\t' + l_split[2] + '\t' + '.' + '\t' + '.' + '\t' + '.' + '\t' + l_split[0] + '\n'
        gff_string = gff_string + temp_string
    i += 1
target_gff = BedTool(gff_string, from_string=True)
    


## Annotation files.

print('Calculating annotations ...')

gencode_ann = BedTool(path_to_gencode_ann).sort()
protein_coding_genes_ann = gencode_ann.filter(lambda x: x[2] == 'gene').filter(lambda x: 'gene_type "protein_coding"' in x[8]).sort()

CGI_ann = BedTool(path_to_CGI).sort()
shore_ann = CGI_ann.flank(g=path_to_chr_lengths, b=2000).sort()
shelf_ann = CGI_ann.flank(g=path_to_chr_lengths, b=4000).subtract(shore_ann).sort()

ChrHMM_ann = BedTool(path_to_ChrHMM).sort()


## Intersections

print('Performing gene bodies / CGI intersections ...')

in_gene_bodies_cgs = list(set(list(target_gff.intersect(protein_coding_genes_ann).sort().to_dataframe()['attributes']))) # 15319 / 21368 CpGs are in gene bodies

in_CGI_cgs = list(set(list(target_gff.intersect(CGI_ann).sort().to_dataframe()['attributes']))) # 9319 / 21368 CpGs are in CGIs
in_shore_cgs = list(set(list(target_gff.intersect(shore_ann).sort().to_dataframe()['attributes']))) # 7920 / 21368 CpGs are in shores
in_shelf_cgs = list(set(list(target_gff.intersect(shelf_ann).sort().to_dataframe()['attributes']))) # 1138 / 21368 CpGs are in CGIs
# 2991 / 21368 CpGs are in open sea

# ChrHMM data requires further processing. We will select annotation for the K562 cell line, which is in position 121 (index:120)

print('Performing ChrHMM intersection and building df ...')

ChrHMM_cgs = target_gff.intersect(ChrHMM_ann, wo=True)
ChrHMM_df = pd.DataFrame(columns = ['CpGmarker', 'ChrHMM_state'])

i = 1
for b in ChrHMM_cgs:
    #print(i)
    ChrHMM_df = ChrHMM_df.append({'CpGmarker':str(b[8]),
                                  'ChrHMM_state':str(b[12]).split('[')[1].split(']')[0].split(',')[120]}, ignore_index=True)
    i += 1

del ChrHMM_ann


## Create final dataframe with binary genomic annotation.

print('Building semi-final df ...')

final_df = pd.DataFrame(columns = ['CpGmarker', 'chr:coord_hg19', 'Gene_body', 'CGI', 'Shore', 'Shelf'])

i = 1
for l in target_gff:
    #print(i)
    final_df = final_df.append({'CpGmarker': str(l[8]), 'chr:coord_hg19': str(l[0])+':'+str(l[3]), 
                                'Gene_body':list(pd.Series(str(l[8])).isin(in_gene_bodies_cgs))[0],
                                'CGI':list(pd.Series(str(l[8])).isin(in_CGI_cgs))[0],
                                'Shore':list(pd.Series(str(l[8])).isin(in_shore_cgs))[0],
                                'Shelf':list(pd.Series(str(l[8])).isin(in_shelf_cgs))[0]}, ignore_index=True)
    i += 1

print('Creating final output file ...')

final_final_df = final_df.merge(ChrHMM_df, on='CpGmarker')
final_final_df.to_csv(output_file, index=False, index_label=False)

##### End of the script #####
