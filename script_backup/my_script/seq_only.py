#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Hands-on Python v1: create sequencing data symbolic links
# 2020/07/07

import os
import sys
import subprocess
import re


try:
    runID = sys.argv[1]
except:
    print ("### ERROR: please give a parameter, e.g. MS20148!\n")
    sys.exit()

### example 1: 正規表達式寫法 & match寫法
#pwd = subprocess.getoutput('pwd | grep "USERDownLoad" | grep ' + runID + ' | wc -l')
#pwd = subprocess.getoutput('pwd')

#if re.match("(.*)USERDownLoad(.*)", pwd):
#    print("found")
#else:
#    print(pwd)


### example 2: 直接用bash串輸出
pwd = subprocess.getoutput('pwd')
check_pwd = subprocess.getoutput(f'pwd | grep "USERDownLoad" | grep {runID} | wc -l')
if check_pwd == "0":
    print("### ERROR: wrong path, e.g. '/export/EC1680U/USERDownLoad/NCHU_CHLaiLab/MS20148'\n")
    sys.exit()


### 使用python內建語言 create symbolic link
fastq_path = '/mnt/NFS/EC2480U-P/scratch/QCResults/fastq'


raw_path = f'{fastq_path}/{runID}/rawData'
if not os.path.exists(raw_path): # 檢查資料夾是否存在可用 os.path.exists
    print(f"### ERROR: please make sure {runID} has been sequencing complete!\n")
    sys.exit()
elif os.path.exists('rawData'):
    print("### WARNINGS: rawData/ has been linked to this folder!\n")
else:
    subprocess.getoutput(f'ln -s {raw_path} .')


clean_path = f'{fastq_path}/{runID}/cleanData'
if not os.path.exists(clean_path):
    print("""### WARNINGS: 2 possibilities for why there is no cleanData, please check...
1.) iQC failed
2.) Seq-only project\n""") # 多行需使用triple quote
elif os.path.exists('cleanData'):
    print("### WARNINGS: cleanData/ has been linked to this folder!\n")
else:
	subprocess.getoutput(f'ln -s {clean_path} .')


raw_report = f'{fastq_path}/{runID}/reports/multiqc_raw/{runID}.raw.multiqc_report.html'
if not os.path.exists(raw_report): # 檢查資料夾是否存在可用 os.path.exists
    print(f"### ERROR: please make sure {runID} has been sequencing complete!\n")
    sys.exit()
elif os.path.exists(f'{runID}.raw.multiqc_report.html'):
    print("### WARNINGS: rawdata multiqc html has been linked to this folder!\n")
else:
    subprocess.getoutput(f'ln -s {raw_report} .')


clean_report = f'{fastq_path}/{runID}/reports/multiqc_clean/{runID}.clean.multiqc_report.html'
if not os.path.exists(clean_report):
    print("""### WARNINGS: 2 possibilities for why there is no multiqc cleandata report, please check...
1.) iQC failed
2.) Seq-only project\n""") # 多行需使用triple quote
elif os.path.exists(f'{runID}.clean.multiqc_report.html'):
    print("### WARNINGS: cleanData multiqc html has been linked to this folder!\n")
else:
    subprocess.getoutput(f'ln -s {clean_report} .')

print('### Complete <(_ _)>')
