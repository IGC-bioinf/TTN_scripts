#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 15:06:23 2022

@author: roma
"""

import pandas as pd
import os
import sys


def readVCFData(file_path, chr): 
    base_name = os.path.basename(file_path)
    print('Reading vcf-file: {}'.format(base_name))
        
    vcf_data = list()
    with open(file_path, 'r') as file:
        for read_line in file:          
            clean_line = read_line.replace('#', '')
            split_line = clean_line.rsplit()
                           
            if '#' not in read_line: vcf_data.append(split_line)            
            if '#CHROM' in read_line: vcf_names = split_line
            
    vcf_frame = pd.DataFrame(vcf_data, columns = vcf_names)
    vcf_frame = vcf_frame.loc[vcf_frame['CHROM'] == str(chr)]
    return vcf_frame

def addPatientsToFrame(df2,vcf):
    df2['HET_PAC']=''
    df2['HOM_PAC']=''
    df2['HET_CON']=''
    df2['HOM_CON']=''
    for i in set(list(df2['SNP'])):
        wf = vcf.loc[vcf['ID'] == i]
        for k  in wf.columns:
            if sys.argv[4] in k:
                if ' 0/1' in str(wf[k]) or ' 0|1' in str(wf[k]) or ' 0/2' in str(wf[k]) :
                  df2["HET_PAC"][int(df2.loc[df2['SNP'] == str(list(wf['ID'][:1])).split("'")[1]].index.values)] = df2["HET_PAC"][int(df2.loc[df2['SNP'] == str(list(wf['ID'][:1])).split("'")[1]].index.values)] + ' ' + k
                elif ' 1/1' in str(wf[k]) or ' 1|1' in str(wf[k]):
                  df2["HOM_PAC"][int(df2.loc[df2['SNP'] == str(list(wf['ID'][:1])).split("'")[1]].index.values)] = df2["HOM_PAC"][int(df2.loc[df2['SNP'] == str(list(wf['ID'][:1])).split("'")[1]].index.values)] + ' ' + k
    for i in set(list(df2['SNP'])):
        wf = vcf.loc[vcf['ID'] == i]
        for k  in wf.columns:
            if sys.argv[4] not in k:
                if ' 0/1' in str(wf[k]) or ' 0|1' in str(wf[k]) or ' 0/2' in str(wf[k]):
                  df2["HET_CON"][int(df2.loc[df2['SNP'] == str(list(wf['ID'][:1])).split("'")[1]].index.values)] = df2["HET_CON"][int(df2.loc[df2['SNP'] == str(list(wf['ID'][:1])).split("'")[1]].index.values)] + ' ' + k
                elif ' 1/1' in str(wf[k]) or ' 1|1' in str(wf[k]):
                  df2["HOM_CON"][int(df2.loc[df2['SNP'] == str(list(wf['ID'][:1])).split("'")[1]].index.values)] = df2["HOM_CON"][int(df2.loc[df2['SNP'] == str(list(wf['ID'][:1])).split("'")[1]].index.values)] + ' ' + k
    return df2

def mergeStatAndEnsembl(df, df2):
    for i in range(len(df2)):
        num = (len(df.loc[df['#Uploaded_variation'] == df2['SNP'][i]])-1)
        df2 = df2.append([df2.loc[df2['SNP'] == df2['SNP'][i]]]*num)
    df2 = df2.sort_values(by=['SNP'], ignore_index=True )
    df = df.sort_values(by=['#Uploaded_variation'], ignore_index=True)
    del df2['CHR']
    df = df.rename(columns={"#Uploaded_variation": "SNP"})
    final_frame = pd.merge(df, df2, on='SNP').drop_duplicates()
    return final_frame

def cleanFinalFrame(final_frame):
    n_col = ['SYMBOL',
    'BIOTYPE',
    'HGVSc',
    'HGVSp',
    'DISTANCE',
    'STRAND',
    'FLAGS',
    'SYMBOL_SOURCE',
    'HGNC_ID',
    'TSL',
    'APPRIS',
    'SIFT',
    'SOMATIC',
    'PHENO',
    'MOTIF_NAME',
    'MOTIF_POS',
    'HIGH_INF_POS',
    'MOTIF_SCORE_CHANGE',
    'Allele',
    'TEST',
    'BP',
    'O(HET)',
    'E(HET)',
    'NCHROBS_A',
    'NCHROBS_U',
    'Gene',
    'Feature_type',
    'TRANSCRIPTION_FACTORS',
    'Location']
    trans_dict = { "ENST00000591111.1": "N2BA",
      "ENST00000460472.2": "N2B",
      "ENST00000342992.6": "N2A",
      "ENST00000359218.5": "Novex-1",
      "ENST00000342175.6": "Novex-2",
      "ENST00000360870.5": "Novex-3",
      "ENST00000589042.1": "Meta"}
    final_frame = final_frame.reset_index(drop=True)
    for key in trans_dict: 
        for i in range(len(final_frame)):
            final_frame['Feature'][i] = final_frame['Feature'][i].replace(key, trans_dict[key])
    for i in n_col:
        del final_frame[i]
    final_frame = final_frame.reset_index(drop=True)
    return final_frame



def addBands(final_frame):
    bands_file = '/home/roma/Chakova/bands_coords.csv'
    bdf = pd.read_csv(bands_file)
    final_frame['Bands'] = ''
    for i in range(len(bdf)):
        for k in range(len(final_frame)):
            if bdf ['0'][i] > int(final_frame['SNP'][k].split(':')[1]) > bdf ['1'][i]:
                final_frame['Bands'][k] = bdf['2'][i]
    return final_frame

def addDomains(final_frame):            
    dom_file = '/home/roma/Chakova/N2BA_domain_coords.csv'
    ddf = pd.read_csv(dom_file)
    final_frame['Domains'] = ''
    for i in range(len(ddf)):
        ddf['ftype'][i] = str(ddf['ftype'][i]).replace('PROSITE-ProRule annotation', '')
    for i in range(len(ddf)):
        for k in range(len(final_frame)):
            if 'N2BA' in final_frame['Feature'][k]:
                if final_frame['Protein_position'][k] != '-':
                    if int(ddf ['id1'][i]) < int(final_frame['Protein_position'][k]) < int(ddf ['id2'][i]):
                        final_frame['Domains'][k] = ddf['ftype'][i]
    for i in range(len(final_frame)):
        if final_frame['Domains'][i] != '':
            tmp_frame = final_frame[final_frame['SNP'].str.contains( final_frame['SNP'][i])]
            for k in list(tmp_frame.index):
                if final_frame['EXON'][k] != '-':
                   final_frame['Domains'][k]=final_frame['Domains'][i]               
    return final_frame

def reorderColumns(final_frame):
    cols = final_frame.columns.tolist()
    order = [0,1,2,3,4,cols.index('Bands'),cols.index('Domains')] + [i for i in range(5,len(cols)-2)] 
    cols = [cols[i] for i in order]
    final_frame=final_frame[cols]
    return final_frame

def processFinalFrame(final_frame):
    final_frame = cleanFinalFrame(final_frame)
    final_frame = addBands(final_frame)
    final_frame = addDomains(final_frame)
    final_frame = reorderColumns(final_frame)
    return final_frame

if __name__ == "__main__":
    df = pd.read_csv(sys.argv[1], sep = '\t')  #ensembl file  
    df2 = pd.read_csv(sys.argv[2], sep = '\t') # csv with snp to annotate
    vcf = readVCFData(sys.argv[3], 2) ## ensembl vcf
    df2 = addPatientsToFrame(df2,vcf)
    final_frame = mergeStatAndEnsembl(df, df2)
    final_frame = processFinalFrame(final_frame)    
    final_frame.to_excel(os.path.dirname(os.path.realpath(sys.argv[1]))+ '/combined_frame_' + sys.argv[4]+'.xlsx', index=False) 

