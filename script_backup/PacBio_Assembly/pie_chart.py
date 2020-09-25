#! /usr/bin/env python

import os
import sys
import pandas as pd
#import matplotlib.pyplot as plt
import plotly.express as px

df = pd.read_csv('4_Annotation/drug_class.sig.plotly.txt', header=None, names=['count','drug_class'])
#labels = df['anti']
#count = df['count']
fig = px.pie(df, values='count', names='drug_class', title='CARD: High Relevance AR Drug Class')
fig.update_layout(legend=dict(orientation="v", yanchor="bottom", y=-0.3, xanchor="right", x=1))
fig.write_html('4_Annotation/drug_class.sig.plotly.html')


df = pd.read_csv('4_Annotation/resistance_mechanism.sig.plotly.txt', header=None, names=['count','mechanism'])
#labels = df['anti']
#count = df['count']
fig = px.pie(df, values='count', names='mechanism', title='CARD: High Relevance Resistance Mechanism')
fig.update_layout(legend=dict(orientation="v", yanchor="bottom", y=-0.3, xanchor="right", x=1))
fig.write_html('4_Annotation/resistance_mechanism.sig.plotly.html')


df = pd.read_csv('4_Annotation/AR_family.sig.plotly.txt', header=None, names=['count','AR_family'])
#labels = df['anti']
#count = df['count']
fig = px.pie(df, values='count', names='AR_family', title='CARD: High Relevance AR Family')
fig.update_layout(legend=dict(orientation="v", yanchor="bottom", y=-0.3, xanchor="right", x=1))
fig.write_html('4_Annotation/AR_family.sig.plotly.html')

#plt.pie(count, labels=labels, autopct='%1.1f%%')
#plt.axis('equal')
#plt.show()
