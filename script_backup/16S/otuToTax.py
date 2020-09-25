import os
import sys
import json
import re

inFile = sys.argv[1]
outFile = sys.argv[2]

if not os.path.exists(inFile):
  print('input file not found!')
  exit()

with open(inFile, 'r') as fp:
  dfInfo = json.load(fp)

taxonomyKey = {}
for i in range(len(dfInfo['rows'])):
  taxTmp = []
  for j in dfInfo['rows'][i]['metadata']['taxonomy']:
    if j.find('unclassified') == -1:
      taxTmp.append(j)

  knownTax = ';'.join(taxTmp)
  knownTax = re.sub('[\^<=>\s\-,:?/\.\'\(\)\[\]#\+\|]{1,}','_',knownTax)

  if knownTax not in taxonomyKey:
    taxonomyKey[knownTax] = {}
    taxonomyKey[knownTax]['rowIds'] = []

  taxonomyKey[knownTax]['rowIds'].append(i)

#refIdTax = {}
with open('/export/EC1680U/DataBase/16S/SILVA_132_release/taxonomy/16S_only/99/raw_taxonomy.txt', 'r') as fp:
  myLine = fp.readline().rstrip('\n')
  while myLine:
    myArray = myLine.split('\t')
    for i in taxonomyKey:
      if 'taxId' in taxonomyKey[i]:
        continue

      if myArray[1].startswith(i):
        taxonomyKey[i]['taxId'] = myArray[0]
#        refIdTax[myArray[0]] = myArray[1]

    myLine = fp.readline().rstrip('\n')

for i in taxonomyKey:
  if 'rowIds' not in taxonomyKey[i] or 'taxId' not in taxonomyKey[i]:
    continue

  for rowId in taxonomyKey[i]['rowIds']:
    dfInfo['rows'][rowId]['id'] = taxonomyKey[i]['taxId']

'''

  taxId = '--'
  refTax = '--'
  if 'taxId' in taxonomyKey[i]:
    taxId = taxonomyKey[i]['taxId']
    refTax = refIdTax[taxonomyKey[i]['taxId']]

  rowIds = '--'
  if 'rowIds' in taxonomyKey[i]:
    rowIds = ','.join(taxonomyKey[i]['rowIds'])

  print('\t'.join([i, taxId, refTax, rowIds]))

#print(dfInfo['rows'][0])
'''

with open(outFile, 'w') as fw:
  json.dump(dfInfo, fw, indent=4)
