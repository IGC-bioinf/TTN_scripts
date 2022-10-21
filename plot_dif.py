#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 15 17:13:43 2022

@author: roma
"""

anev_file='/home/roma/Chakova/to_mv/new_all/done/combined_frame_anev.xlsx'
arh_file='/home/roma/Chakova/to_mv/new_all/done/combined_frame_arh.xlsx'
gkmp_file='/home/roma/Chakova/to_mv/new_all/done/combined_frame_gkmp.xlsx'
mio_file='/home/roma/Chakova/to_mv/new_all/done/combined_frame_mio.xlsx'

import matplotlib.ticker as ticker
from itertools import combinations
import pandas as pd
import matplotlib.pyplot as plt

anev=pd.read_excel(anev_file)
arh=pd.read_excel(arh_file)
mio=pd.read_excel(mio_file)
gkmp=pd.read_excel(gkmp_file)
anev=anev[~anev['SNP'].str.contains('179481926')]
arh=arh[~arh['SNP'].str.contains('179481926')]
mio=mio[~mio['SNP'].str.contains('179481926')]
gkmp=gkmp[~gkmp['SNP'].str.contains('179481926')]

frames = [anev, arh, mio,gkmp]
result = pd.concat(frames).drop_duplicates().reset_index(drop=True)
con = result[(result['HOM_CON'].fillna('').str.contains(r'.') | result['HET_CON'].fillna('').str.contains(r'.')) & ~result['SNP'].str.contains('179604264') & ~result['SNP'].str.contains('179481926') ]
bands = [x for x in list(result.sort_values(by=['SNP'])['Bands'].drop_duplicates().dropna()) if str(x) != 'nan']
domains = [x for x in list(result.sort_values(by=['SNP'])['Domains'].drop_duplicates().dropna()) if str(x) != 'nan']
dict_frames = {'Anev':anev, 'Arh':arh, 'Mio':mio, 'Gkmp':gkmp, 'Con':con}
columns = [i for i in dict_frames.keys()]
L_combines = [",".join(map(str, comb)) for comb in combinations(columns, 2)]

def drop_AF(frame):
    for i in range(len(frame)):
        if str(frame['AF'][i]) == '-' :
            frame = frame.drop(index = i)
        elif float(frame['AF'][i]) > 0.01:
            frame = frame.drop(index = i)
            
    frame = frame.reset_index(drop=True)
    return frame

def neg_to_pos(n, position):
    n = int(n)
    if n < 0:
        return str(n * -1)
    else:
        return str(n)

def createPlot(frame1, frame2, Name1, Name2, pattern_list, pattern):
    columns = [Name1+'_missence', Name1 + '_nonsence', Name2 + '_missence', Name2 + '_nonsence', Name1 +'_common',  Name2 +'_common']

    df = pd.DataFrame(columns=columns, index=pattern_list)
    
    if Name1 == 'Con':
        for i in df.index:
            df[Name1 + '_missence'][i] = len(frame1[(frame1[pattern].fillna('') == i)\
            & (frame1['HOM_CON'].fillna('').str.contains(r'.') \
            | frame1['HET_CON'].fillna('').str.contains(r'.')) \
            & (frame1['Consequence'].fillna('').str.contains('missense_variant'))].loc[:, frame1.columns=='SNP'].drop_duplicates())
            
            df[Name1 + '_nonsence'][i] = len(frame1[(frame1[pattern].fillna('') == i)\
            & (frame1['HOM_CON'].fillna('').str.contains(r'.') \
            | frame1['HET_CON'].fillna('').str.contains(r'.')) \
            & (frame1['Consequence'].fillna('').str.contains('stop_gained'))].loc[:, frame1.columns=='SNP'].drop_duplicates())
            
            df[Name2 + '_missence'][i] = len(frame2[(frame2[pattern].fillna('') == i) \
            & (frame2['HOM_PAC'].fillna('').str.contains(r'.') \
            | frame2['HET_PAC'].fillna('').str.contains(r'.')) \
            & (frame2['Consequence'].fillna('').str.contains('missense_variant'))].loc[:, frame2.columns=='SNP'].drop_duplicates())*-1
            
            df[Name2 + '_nonsence'][i] = len(frame2[(frame2[pattern].fillna('') == i) \
            & (frame2['HOM_PAC'].fillna('').str.contains(r'.') \
            | frame2['HET_PAC'].fillna('').str.contains(r'.')) \
            & frame2['Consequence'].fillna('').str.contains('stop_gained')].loc[:, frame2.columns=='SNP'].drop_duplicates())*-1
            
            df[Name1 +'_common'][i] = len(pd.merge(frame1[(frame1[pattern].fillna('') == i)\
            & (frame1['HOM_CON'].fillna('').str.contains(r'.') \
            | frame1['HET_CON'].fillna('').str.contains(r'.')) \
            & (frame1['Consequence'].fillna('').str.contains('missense_variant'))].loc[:, frame1.columns=='SNP'].drop_duplicates() ,\
            frame2[(frame2[pattern].fillna('') == i) \
            & (frame2['HOM_PAC'].fillna('').str.contains(r'.') \
            | frame2['HET_PAC'].fillna('').str.contains(r'.')) \
            & (frame2['Consequence'].fillna('').str.contains('missense_variant'))].loc[:, frame2.columns=='SNP'].drop_duplicates(),\
            how='inner'))
            
            df[Name2 +'_common'][i] = len(pd.merge(frame1[(frame1[pattern].fillna('') == i)\
            & (frame1['HOM_CON'].fillna('').str.contains(r'.') \
            | frame1['HET_CON'].fillna('').str.contains(r'.')) \
            & (frame1['Consequence'].fillna('').str.contains('missense_variant'))].loc[:, frame1.columns=='SNP'].drop_duplicates() ,\
            frame2[(frame2[pattern].fillna('') == i) \
            & (frame2['HOM_PAC'].fillna('').str.contains(r'.') \
            | frame2['HET_PAC'].fillna('').str.contains(r'.')) \
            & (frame2['Consequence'].fillna('').str.contains('missense_variant'))].loc[:, frame2.columns=='SNP'].drop_duplicates(),\
            how='inner'))*-1
            
                
        df.to_excel('/home/roma/Chakova/'+Name1 + '_' + Name2 + '_comparsion.xlsx')        
        for i in df.index:
            n1_all=len(frame1[(frame1['HOM_CON'].fillna('').str.contains(r'.') \
            | frame1['HET_CON'].fillna('').str.contains(r'.'))].loc[:, frame1.columns=='SNP'].drop_duplicates())
            n2_all=len(frame2[(frame2['HOM_PAC'].fillna('').str.contains(r'.') \
            | frame2['HET_PAC'].fillna('').str.contains(r'.'))].loc[:, frame2.columns=='SNP'].drop_duplicates())
            #miss_all= df[Name1 + '_missence'][i] + df[Name2 + '_missence'][i]*-1 - df[Name2 +'_common'][i]
            #non_all = df[Name1 + '_nonsence'][i] + df[Name2 + '_nonsence'][i]*-1 - len(pd.merge(frame1[(frame1[pattern].fillna('') == i)\
            #& (frame1['HOM_PAC'].fillna('').str.contains(r'.') \
            #| frame1['HET_PAC'].fillna('').str.contains(r'.')) \
            #& (frame1['Consequence'].fillna('').str.contains('stop_gained'))].loc[:, frame1.columns=='SNP'].drop_duplicates() ,\
            #frame2[(frame2[pattern].fillna('') == i) \
            #& (frame2['HOM_PAC'].fillna('').str.contains(r'.') \
            #| frame2['HET_PAC'].fillna('').str.contains(r'.')) \
            #& (frame2['Consequence'].fillna('').str.contains('stop_gained'))].loc[:, frame2.columns=='SNP'].drop_duplicates(),\
            #how='inner'))        
            df[Name1 + '_missence'][i] =  df[Name1 + '_missence'][i] / n1_all * 100
            df[Name1 + '_nonsence'][i] = df[Name1 + '_nonsence'][i] / n1_all * 100
            df[Name2 + '_missence'][i] = df[Name2 + '_missence'][i] / n2_all * 100
            df[Name2 + '_nonsence'][i] = df[Name2 + '_nonsence'][i]/ n2_all * 100
            df[Name1 +'_common'][i] = df[Name1 +'_common'][i] / n1_all * 100
            df[Name2 +'_common'][i] =   df[Name2 +'_common'][i] / n2_all * 100
            #if miss_all == 0:
             #   miss_all = 0.01
            #if non_all == 0:
            #    non_all = 0.01
            
    elif Name2 == 'Con':
        for i in df.index:
            df[Name1 + '_missence'][i] = len(frame1[(frame1[pattern].fillna('') == i)\
            & (frame1['HOM_PAC'].fillna('').str.contains(r'.') \
            | frame1['HET_PAC'].fillna('').str.contains(r'.')) \
            & (frame1['Consequence'].fillna('').str.contains('missense_variant'))].loc[:, frame1.columns=='SNP'].drop_duplicates())
            
            df[Name1 + '_nonsence'][i] = len(frame1[(frame1[pattern].fillna('') == i)\
            & (frame1['HOM_PAC'].fillna('').str.contains(r'.') \
            | frame1['HET_PAC'].fillna('').str.contains(r'.')) \
            & (frame1['Consequence'].fillna('').str.contains('stop_gained'))].loc[:, frame1.columns=='SNP'].drop_duplicates())
            
            df[Name2 + '_missence'][i] = len(frame2[(frame2[pattern].fillna('') == i) \
            & (frame2['HOM_CON'].fillna('').str.contains(r'.') \
            | frame2['HET_CON'].fillna('').str.contains(r'.')) \
            & (frame2['Consequence'].fillna('').str.contains('missense_variant'))].loc[:, frame2.columns=='SNP'].drop_duplicates())*-1
            
            df[Name2 + '_nonsence'][i] = len(frame2[(frame2[pattern].fillna('') == i) \
            & (frame2['HOM_CON'].fillna('').str.contains(r'.') \
            | frame2['HET_CON'].fillna('').str.contains(r'.')) \
            & frame2['Consequence'].fillna('').str.contains('stop_gained')].loc[:, frame2.columns=='SNP'].drop_duplicates())*-1
            
            df[Name1 +'_common'][i] = len(pd.merge(frame1[(frame1[pattern].fillna('') == i)\
            & (frame1['HOM_PAC'].fillna('').str.contains(r'.') \
            | frame1['HET_PAC'].fillna('').str.contains(r'.')) \
            & (frame1['Consequence'].fillna('').str.contains('missense_variant'))].loc[:, frame1.columns=='SNP'].drop_duplicates() ,\
            frame2[(frame2[pattern].fillna('') == i) \
            & (frame2['HOM_CON'].fillna('').str.contains(r'.') \
            | frame2['HET_CON'].fillna('').str.contains(r'.')) \
            & (frame2['Consequence'].fillna('').str.contains('missense_variant'))].loc[:, frame2.columns=='SNP'].drop_duplicates(),\
            how='inner'))
            
            df[Name2 +'_common'][i] = len(pd.merge(frame1[(frame1[pattern].fillna('') == i)\
            & (frame1['HOM_PAC'].fillna('').str.contains(r'.') \
            | frame1['HET_PAC'].fillna('').str.contains(r'.')) \
            & (frame1['Consequence'].fillna('').str.contains('missense_variant'))].loc[:, frame1.columns=='SNP'].drop_duplicates() ,\
            frame2[(frame2[pattern].fillna('') == i) \
            & (frame2['HOM_CON'].fillna('').str.contains(r'.') \
            | frame2['HET_CON'].fillna('').str.contains(r'.')) \
            & (frame2['Consequence'].fillna('').str.contains('missense_variant'))].loc[:, frame2.columns=='SNP'].drop_duplicates(),\
            how='inner'))*-1
            #miss_all= df[Name1 + '_missence'][i] + df[Name2 + '_missence'][i]*-1 - df[Name2 +'_common'][i]
            #non_all = df[Name1 + '_nonsence'][i] + df[Name2 + '_nonsence'][i]*-1 - len(pd.merge(frame1[(frame1[pattern].fillna('') == i)\
            #& (frame1['HOM_PAC'].fillna('').str.contains(r'.') \
            #| frame1['HET_PAC'].fillna('').str.contains(r'.')) \
            #& (frame1['Consequence'].fillna('').str.contains('stop_gained'))].loc[:, frame1.columns=='SNP'].drop_duplicates() ,\
            #frame2[(frame2[pattern].fillna('') == i) \
            #& (frame2['HOM_PAC'].fillna('').str.contains(r'.') \
            #| frame2['HET_PAC'].fillna('').str.contains(r'.')) \
            #& (frame2['Consequence'].fillna('').str.contains('stop_gained'))].loc[:, frame2.columns=='SNP'].drop_duplicates(),\
            #how='inner'))
            #if miss_all == 0:
            #    miss_all = 0.01
            #if non_all == 0:
            #    non_all = 0.01
            #df[Name1 + '_missence'][i] =  df[Name1 + '_missence'][i] / miss_all * 100
            #df[Name1 + '_nonsence'][i] = df[Name1 + '_nonsence'][i] / non_all * 100
            #df[Name2 + '_missence'][i] = df[Name2 + '_missence'][i] / miss_all * 100
            #df[Name2 + '_nonsence'][i] = df[Name2 + '_nonsence'][i]/ non_all * 100
            #df[Name1 +'_common'][i] = df[Name1 +'_common'][i] / miss_all * 100
            #df[Name2 +'_common'][i] =   df[Name2 +'_common'][i] / miss_all * 100
        df.to_excel('/home/roma/Chakova/'+Name1 + '_' + Name2 + '_comparsion.xlsx')        
        for i in df.index:
            n1_all= len(frame1[(frame1['HOM_PAC'].fillna('').str.contains(r'.') \
            | frame1['HET_PAC'].fillna('').str.contains(r'.'))].loc[:, frame1.columns=='SNP'].drop_duplicates())
            n2_all=len(frame2[(frame2['HOM_CON'].fillna('').str.contains(r'.') \
            | frame2['HET_CON'].fillna('').str.contains(r'.'))].loc[:, frame2.columns=='SNP'].drop_duplicates())
            df[Name1 + '_missence'][i] =  df[Name1 + '_missence'][i] / n1_all * 100
            df[Name1 + '_nonsence'][i] = df[Name1 + '_nonsence'][i] / n1_all * 100
            df[Name2 + '_missence'][i] = df[Name2 + '_missence'][i] / n2_all * 100
            df[Name2 + '_nonsence'][i] = df[Name2 + '_nonsence'][i]/ n2_all * 100
            df[Name1 +'_common'][i] = df[Name1 +'_common'][i] / n1_all * 100
            df[Name2 +'_common'][i] =   df[Name2 +'_common'][i] / n2_all * 100
            #if miss_all == 0:
    else:
        for i in df.index:
            df[Name1 + '_missence'][i] = len(frame1[(frame1[pattern].fillna('') == i)\
            & (frame1['HOM_PAC'].fillna('').str.contains(r'.') \
            | frame1['HET_PAC'].fillna('').str.contains(r'.')) \
            & (frame1['Consequence'].fillna('').str.contains('missense_variant'))].loc[:, frame1.columns=='SNP'].drop_duplicates())
            
            df[Name1 + '_nonsence'][i] = len(frame1[(frame1[pattern].fillna('') == i)\
            & (frame1['HOM_PAC'].fillna('').str.contains(r'.') \
            | frame1['HET_PAC'].fillna('').str.contains(r'.')) \
            & (frame1['Consequence'].fillna('').str.contains('stop_gained'))].loc[:, frame1.columns=='SNP'].drop_duplicates())
            
            df[Name2 + '_missence'][i] = len(frame2[(frame2[pattern].fillna('') == i) \
            & (frame2['HOM_PAC'].fillna('').str.contains(r'.') \
            | frame2['HET_PAC'].fillna('').str.contains(r'.')) \
            & (frame2['Consequence'].fillna('').str.contains('missense_variant'))].loc[:, frame2.columns=='SNP'].drop_duplicates())*-1
            
            df[Name2 + '_nonsence'][i] = len(frame2[(frame2[pattern].fillna('') == i) \
            & (frame2['HOM_PAC'].fillna('').str.contains(r'.') \
            | frame2['HET_PAC'].fillna('').str.contains(r'.')) \
            & frame2['Consequence'].fillna('').str.contains('stop_gained')].loc[:, frame2.columns=='SNP'].drop_duplicates())*-1
            
            df[Name1 +'_common'][i] = len(pd.merge(frame1[(frame1[pattern].fillna('') == i)\
            & (frame1['HOM_PAC'].fillna('').str.contains(r'.') \
            | frame1['HET_PAC'].fillna('').str.contains(r'.')) \
            & (frame1['Consequence'].fillna('').str.contains('missense_variant'))].loc[:, frame1.columns=='SNP'].drop_duplicates() ,\
            frame2[(frame2[pattern].fillna('') == i) \
            & (frame2['HOM_PAC'].fillna('').str.contains(r'.') \
            | frame2['HET_PAC'].fillna('').str.contains(r'.')) \
            & (frame2['Consequence'].fillna('').str.contains('missense_variant'))].loc[:, frame2.columns=='SNP'].drop_duplicates(),\
            how='inner'))
            
            df[Name2 +'_common'][i] = len(pd.merge(frame1[(frame1[pattern].fillna('') == i)\
            & (frame1['HOM_PAC'].fillna('').str.contains(r'.') \
            | frame1['HET_PAC'].fillna('').str.contains(r'.')) \
            & (frame1['Consequence'].fillna('').str.contains('missense_variant'))].loc[:, frame1.columns=='SNP'].drop_duplicates() ,\
            frame2[(frame2[pattern].fillna('') == i) \
            & (frame2['HOM_PAC'].fillna('').str.contains(r'.') \
            | frame2['HET_PAC'].fillna('').str.contains(r'.')) \
            & (frame2['Consequence'].fillna('').str.contains('missense_variant'))].loc[:, frame2.columns=='SNP'].drop_duplicates(),\
            how='inner'))*-1
            #miss_all= df[Name1 + '_missence'][i] + df[Name2 + '_missence'][i]*-1 - df[Name2 +'_common'][i] 
            #non_all = df[Name1 + '_nonsence'][i] + df[Name2 + '_nonsence'][i]*-1 - len(pd.merge(frame1[(frame1[pattern].fillna('') == i)\
            #& (frame1['HOM_PAC'].fillna('').str.contains(r'.') \
            #| frame1['HET_PAC'].fillna('').str.contains(r'.')) \
            #& (frame1['Consequence'].fillna('').str.contains('stop_gained'))].loc[:, frame1.columns=='SNP'].drop_duplicates() ,\
            #frame2[(frame2[pattern].fillna('') == i) \
            #& (frame2['HOM_PAC'].fillna('').str.contains(r'.') \
            #| frame2['HET_PAC'].fillna('').str.contains(r'.')) \
            #& (frame2['Consequence'].fillna('').str.contains('stop_gained'))].loc[:, frame2.columns=='SNP'].drop_duplicates(),\
            #how='inner'))
            #if miss_all == 0:
            #    miss_all = 0.01
            #if non_all == 0:
            #    non_all = 0.01
        df.to_excel('/home/roma/Chakova/'+Name1 + '_' + Name2 + '_comparsion.xlsx')        
        for i in df.index:                
            n1_all = len(frame1[(frame1['HOM_PAC'].fillna('').str.contains(r'.') \
            | frame1['HET_PAC'].fillna('').str.contains(r'.'))].loc[:, frame1.columns=='SNP'].drop_duplicates())
            n2_all = len(frame2[(frame2['HOM_PAC'].fillna('').str.contains(r'.') \
            | frame2['HET_PAC'].fillna('').str.contains(r'.'))].loc[:, frame2.columns=='SNP'].drop_duplicates())
            df[Name1 + '_missence'][i] =  df[Name1 + '_missence'][i] / n1_all * 100
            df[Name1 + '_nonsence'][i] = df[Name1 + '_nonsence'][i] / n1_all * 100
            df[Name2 + '_missence'][i] = df[Name2 + '_missence'][i] / n2_all * 100
            df[Name2 + '_nonsence'][i] = df[Name2 + '_nonsence'][i]/ n2_all * 100
            df[Name1 +'_common'][i] = df[Name1 +'_common'][i] / n1_all * 100
            df[Name2 +'_common'][i] =   df[Name2 +'_common'][i] / n2_all * 100
            #if miss_all == 0:
            #df[Name1 + '_missence'][i] =  df[Name1 + '_missence'][i] / miss_all * 100
            #df[Name1 + '_nonsence'][i] = df[Name1 + '_nonsence'][i] / non_all * 100
            #df[Name2 + '_missence'][i] = df[Name2 + '_missence'][i] / miss_all * 100
            #df[Name2 + '_nonsence'][i] = df[Name2 + '_nonsence'][i]/ non_all * 100
            #df[Name1 +'_common'][i] = df[Name1 +'_common'][i] / miss_all * 100
            #df[Name2 +'_common'][i] =   df[Name2 +'_common'][i] / miss_all * 100
    
    df = df.loc[~(df==0).all(axis=1)]
    
    
    max_v = df.abs().max().max()
    fig, ax = plt.subplots(1, figsize=(18, len(df)*0.8))
    
    plt.barh(df.index, df[df.columns[0]], color = '#337AE3', height = 1)
    plt.barh(df.index, df[df.columns[1]], left = df[df.columns[0]], color = '#5E96E9', height = 1)
    plt.barh(df.index, df[df.columns[2]], color = '#DB4444', height =1)
    plt.barh(df.index, df[df.columns[3]], left = df[df.columns[2]], color = '#E17979', height = 1)
    plt.barh(df.index, df[df.columns[4]], color = 'green', height = 1)
    plt.barh(df.index, df[df.columns[5]], color = 'green', height = 1)
    
    
    plt.xlim((max_v*-1)-10, max_v+10)
    
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.set_axisbelow(True)
    #px= ax.twiny()
    ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
    ax.xaxis.set_major_formatter(ticker.PercentFormatter())
    #ax.xaxis.set_major_formatter(ticker.FuncFormatter(neg_to_pos))
    
    
    plt.xticks(rotation=90)
    ax.xaxis.grid(color='gray', linestyle='dashed', alpha=0.7)
    ax.tick_params(axis='both', which='major', labelsize=15)
    #plt.xticks(df.index)
    legend_label = list(df.columns)
    for i in range(len(legend_label)):
        legend_label[i] = legend_label[i].replace('_', ' ').replace('Gkmp', 'HCM').replace('Anev', 'TAAD').replace('Mio', 'LVNC').replace('Con', 'Unaffected').replace('Arh', 'LQTS+BrS')
    plt.legend(legend_label, ncol=2, frameon = False, prop={'size': 15}, bbox_to_anchor=(1.04,1), loc="upper left")
    plt.title(Name1.replace('_', ' ').replace('Gkmp', 'HCM').replace('Anev', 'TAAD').replace('Mio', 'LVNC').replace('Con', 'Unaffected').replace('Arh', 'LQTS+BrS') \
    + '  ' + Name2.replace('_', ' ').replace('Gkmp', 'HCM').replace('Anev', 'TAAD').replace('Mio', 'LVNC').replace('Con', 'Unaffected').replace('Arh', 'LQTS+BrS') + ' comparsion', loc = 'left', fontsize=15)
    

    plt.savefig('/home/roma/Chakova/'+Name1 + '_' + Name2 + '_comparsion.png', dpi=100,  bbox_inches = 'tight')
    plt.show()



anev_l = []
arh_l = []
mio_l = [] 
gkmp_l = []

list_phy = [mio_l, gkmp_l]
frames = [mio, gkmp]
for i,v in zip(frames, list_phy):
    tmp = i[(i['HOM_PAC'].fillna('').str.contains(r'.') \
        | i['HET_PAC'].fillna('').str.contains(r'.')) \
        & (i['Consequence'].fillna('').str.contains('missense_variant'))].loc[:, i.columns=='SNP'].drop_duplicates().reset_index(drop=True)
    for r in range(len((tmp))):
        v.append(tmp['SNP'][r].split(':')[1])


con_l = []  
tmp_con = con[(con['HOM_CON'].fillna('').str.contains(r'.') \
            | con['HET_CON'].fillna('').str.contains(r'.')) \
            & (con['Consequence'].fillna('').str.contains('missense_variant'))].loc[:, con.columns=='SNP'].drop_duplicates().reset_index(drop=True)
for i in range(len(tmp_con)):
    con_l.append(tmp_con['SNP'][i].split(':')[1])
from pyvenn import venn
list_phy = [mio_l, gkmp_l, con_l]
labels = venn.get_labels(list_phy)
fig, ax = venn.venn3(labels)
fig.show()
import numpy as np
import seaborn as sns
def createManh(i):
    min_v = 0.05 / len(i[i['P'] < 0.05]['SNP'].drop_duplicates())
    max_v = 0.05/len(i['SNP'].drop_duplicates())
    i['CHR'] = ''
    i['BP'] = ''
    for k in range(len(i)):
        i['CHR'][k] = int(i['SNP'][k].split(':')[0])
        i['BP'][k] = int(i['SNP'][k].split(':')[1])
    i['-logp'] = -np.log10(i.P); i = i.sort_values(['CHR','BP'])
    i.reset_index(inplace=True, drop=True); i['i'] = i.index
    
    # Generate Manhattan plot: (#optional tweaks for relplot: linewidth=0, s=9)
    plot = sns.relplot(data=i, x='i', y='-logp', aspect=3.7, 
                       hue='BP', palette = 'bright', legend=None) 
    chrom_i=i.groupby('CHR')['i'].median()
    plot.ax.set_xlabel('Position'); plot.ax.set_xticks(chrom_i);
    plot.ax.set_xticklabels(chrom_i.index)
    #plot.fig.suptitle(n)
    plt.plot(i['i'], [-np.log10(min_v)] * len(i['i']))
    plt.plot(i['i'], [-np.log10(max_v)] * len(i['i']))
    plt.ylim(0,9)

for i in L_combines:
    createPlot(dict_frames[i.split(',')[0]], dict_frames[i.split(',')[1]], i.split(',')[0], i.split(',')[1], bands, 'Bands')

