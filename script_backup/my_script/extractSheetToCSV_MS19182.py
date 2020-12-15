#! /usr/bin/env python

import pandas as pd
import numpy as np

df = pd.read_excel('MS19182_human.xlsx', sheet_name=None)

#df = pd.ExcelFile('MS19182_mice.xlsx')

for key, value in df.items(): 
    df[key].to_csv('%s.csv' %key, sep="\t", index=0)


#for i in range(1,73):
#    sheet1 = xls.parse(i)
#    print (sheet1)
