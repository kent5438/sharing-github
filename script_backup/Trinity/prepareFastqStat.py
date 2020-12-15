#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import json
import pandas as pd
import csv

# processing data for downstream report

# settings
stat_dir = sys.argv[1]


def parse_fastqStat(fileName):
    reads = total_bases = len_min = len_mean = len_max = \
                                gc_ratio = q20_ratio = \
                                q30_ratio = qual_mean = phred = ""
    with open(fileName) as f:
        for line in f:
            if line.startswith("reads"):
                reads = line.strip().split("\t")[1]
                continue
            if line.startswith("total bases"):
                total_bases = line.strip().split("\t")[1]
                continue
            if line.startswith("len max"):
                len_max = line.strip().split("\t")[1]
                continue
            if line.startswith("len min"):
                len_min = line.strip().split("\t")[1]
                continue
            if line.startswith("len mean"):
                len_mean = line.strip().split("\t")[1]
                continue
            if line.startswith("%GC"):
                gc_ratio = line.strip().split("\t")[1]
                continue
            if line.startswith("%Q20"):
                q20_ratio = line.strip().split("\t")[1]
                continue
            if line.startswith("%Q30"):
                q30_ratio = line.strip().split("\t")[1]
                continue
            if line.startswith("qual mean"):
                qual_mean = line.strip().split("\t")[1]
                continue
            if line.startswith("phred"):
                phred = line.strip().split("\t")[1]
                continue

    if len_min != len_max:
        length = len_min+"-"+len_max
    else:
        length = len_min
    
    print(fileName,total_bases)
    return {'Total_Bases': int(total_bases),
            'GC': float(gc_ratio),
            'Phred': phred,
            'Qual_Mean': float(qual_mean),
            'Q20_Ratio': float(q20_ratio),
            'Q30_Ratio': float(q30_ratio),
            'Length': length,
            'Length_Mean': float(len_mean),
            'Reads': int(reads),
    }

fastq_stats = {}
for f in os.listdir(stat_dir):
    if f != '.' and f != ".." and f != "summary.csv" and f.endswith("fastq-stats.txt"):
      sample_name = ".".join(f.split(".")[:-2])
      tmp = parse_fastqStat(os.path.join(stat_dir, f))
      #tmp['File_Name'] = sample_name+".fastq.gz"
      fastq_stats[sample_name+".fastq.gz"] = tmp



fastqSTATFile = os.path.join(stat_dir, "summary.csv")
df = pd.DataFrame(data=fastq_stats).T
df.index.name = "File_Name"
df.to_csv(fastqSTATFile)



project_info = dict()
project_info['project'] = dict()
project_info['samples'] = dict()
for fn, v in fastq_stats.items():
    if "_" in fn:
        sample_name = fn.split("_")[0]
    else:
        sample_name = fn.split(".")[0]


    try:
        #print(fn)
        project_info['samples'][sample_name]['Reads'] += v['Reads']
        #print("TotalBases: {}".format(project_info['samples'][sample_name]['Total_Bases'] ))
        #print("TotalQ20Ratio: {}".format(project_info['samples'][sample_name]['Q20_Ratio'] ))
        #print("FileTotalBase: {}".format(v['Total_Bases']))
        #print("FileQ20Ratio: {}".format(v['Q20_Ratio']))
        project_info['samples'][sample_name]['Q30_Ratio'] = (project_info['samples'][sample_name]['Total_Bases'] * project_info['samples'][sample_name]['Q30_Ratio'] + v['Total_Bases'] * v['Q30_Ratio']/100)/( v['Total_Bases'] + project_info['samples'][sample_name]['Total_Bases'])
        project_info['samples'][sample_name]['Q20_Ratio'] = (project_info['samples'][sample_name]['Total_Bases'] * project_info['samples'][sample_name]['Q20_Ratio'] + v['Total_Bases'] * v['Q20_Ratio']/100)/( v['Total_Bases'] + project_info['samples'][sample_name]['Total_Bases'])
        project_info['samples'][sample_name]['Length_Mean'] = (project_info['samples'][sample_name]['Reads'] * project_info['samples'][sample_name]['Length_Mean'] + v['Reads'] * v['Length_Mean']) / (project_info['samples'][sample_name]['Reads'] + v['Reads'])
        project_info['samples'][sample_name]['GC'] = (project_info['samples'][sample_name]['Total_Bases'] * project_info['samples'][sample_name]['GC'] + v['Total_Bases'] * v['GC']/100) / (project_info['samples'][sample_name]['Total_Bases'] + v['Total_Bases'])
        project_info['samples'][sample_name]['Qual_Mean'] = (project_info['samples'][sample_name]['Total_Bases'] * project_info['samples'][sample_name]['Qual_Mean'] + v['Total_Bases'] * v['Qual_Mean']) / (project_info['samples'][sample_name]['Total_Bases'] + v['Total_Bases'])
        project_info['samples'][sample_name]['Total_Bases'] += v['Total_Bases']
        project_info['samples'][sample_name]['details'][fn] = v

        #print(project_info['samples'][sample_name])
    except KeyError:
        #print(e.message)
        #print(fn)
        project_info['samples'][sample_name] = {"Total_Bases": v['Total_Bases'],
                                                "Reads": v['Reads'],
                                                "Q30_Ratio": v['Q30_Ratio']/100,
                                                "Q20_Ratio": v['Q20_Ratio']/100,
                                                "Length_Mean": v['Length_Mean'],
                                                "GC": v['GC']/100,
                                                "Qual_Mean": v['Qual_Mean']}
        project_info['samples'][sample_name]['details'] = dict()
        project_info['samples'][sample_name]['details'][fn] = v
        #print("TotalBases: {}".format(project_info['samples'][sample_name]['Total_Bases'] ))
        #print("TotalQ20Ratio: {}".format(project_info['samples'][sample_name]['Q20_Ratio'] ))
        #print("FileTotalBase: {}".format(v['Total_Bases']))
        #print("FileQ20Ratio: {}".format(v['Q20_Ratio']))

        #print(project_info['samples'][sample_name])

    try:
        project_info['project']['Total_Bases'] += v['Total_Bases']
    except:
        project_info['project']['Total_Bases'] = v['Total_Bases']


print(project_info)
for sn, v in project_info['samples'].items():
    read_type = "SE"

    for fn in v['details'].keys():
        if "R2" in fn:
            read_type = "PE"

    project_info['samples'][sn]['Read_Type'] = read_type
    project_info['samples'][sn]['Total_Bases'] = project_info['samples'][sn]['Total_Bases']
    project_info['samples'][sn]['GC'] = project_info['samples'][sn]['GC'] * 100
    project_info['samples'][sn]['Q20_Ratio'] = project_info['samples'][sn]['Q20_Ratio'] * 100
    project_info['samples'][sn]['Q30_Ratio'] = project_info['samples'][sn]['Q30_Ratio'] * 100

    if read_type == "PE":
        pass
        #project_info['samples'][sn]['Reads'] = project_info['samples'][sn]['Reads'] / 2
        #'reads' means "total reads"" not "read pairs"


fastqSTATJsonFile = os.path.join(stat_dir, "summary.json")
with open(fastqSTATJsonFile, "w") as fw:
    json.dump(project_info, fw, indent=4)


totalBasesFile = os.path.join(stat_dir, "total_bases.txt")
with open(totalBasesFile, "w") as fw:
    fw.write("{}".format(project_info['project']['Total_Bases']))

sampleSummaryFile = os.path.join(stat_dir, "samples.summary.csv")
with open(sampleSummaryFile, "w") as fw:
    f_writer = csv.writer(fw, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    f_writer.writerow(["Sample_Name", "Reads", "Total_Bases", "Q20_Ratio", "Q30_Ratio", "Qual_Mean", "GC", "Length_Mean"])
    #headers = ["Sample_Name", "Reads", "Total_Bases", "Q20_Ratio", "Q30_Ratio", "Qual_Mean", "GC", "Length_Mean"]
    #fw.write(",".join(headers)+"\n")
    for sn, v in  project_info['samples'].items():
        tmps = ["'{}'".format(sn), v["Reads"], v["Total_Bases"], v["Q20_Ratio"], v["Q30_Ratio"], v["Qual_Mean"], v["GC"], v["Length_Mean"]]
        f_writer.writerow(tmps)
        #fw.write(",".join(str(x) for x in tmps)+"\n")

