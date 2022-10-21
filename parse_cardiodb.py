#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 13:33:01 2022

@author: roma
"""

from bs4 import BeautifulSoup
import requests
import pandas
import os
import sys
import  re
url = 'https://www.cardiodb.org/titin/titin_transcripts.php'
resp =requests.get(url)
soup = BeautifulSoup(resp.text, "lxml")
quotes = soup.find_all('a')
np = []
for i in quotes:
        if (bool(re.search(r'\d\d\d\d\d\d\d\d\d', str(i).split('>')[1]))) or ('png' in str(i)):
            if ('ENS' not in str(i)) and ('NM' not in str(i)) and ('NP' not in str(i)):
                np.append(str(i)[str(i).find('>')+1:str(i).find('</a')])
                
def maybeMakeNumber(s):

    if not s:
        return s
    try:
        f = float(s)
        i = int(f)
        return i if f == i else f
    except ValueError:
        return s

fc = []
sc = []
band = []
np = list(map(maybeMakeNumber, np))
tmp_list = np
result_list = []

str_index = list()

for index, item in enumerate(np): 
     if type(item) == str: str_index.append(index)

str_index.sort(reverse=True) 

for item in str_index: 
    
    sub_list  = np[item-2:item+1]
    
    search_list = re.findall(r'title="(.*?)"', sub_list[2])
    sub_list[2] = ', '.join(search_list)

    result_list.append(sub_list)
    
    del tmp_list[item-2:item+1]
    
for item in range(1, len(tmp_list), 2): 
    result_list.append(tmp_list[item-1:item+1]+[''])

result_frame = pandas.DataFrame(result_list, columns=['id1', 'id2', 'ftype'])
result_frame.to_csv("/home/roma/Chakova/domain_coords.csv", index=False, sep='\t')
