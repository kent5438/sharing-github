#! /bin/bash

ccs -j 32 --min-rq 0.9 --min-passes 1 --min-length 100 --max-length 2000 --report-file m54126_191025_082136.ccs.report.txt ../rawData/m54126_191024_220212.subreads.bam m54126_191024_220212.ccs.1.bam --chunk 1/10
ccs -j 32 --min-rq 0.9 --min-passes 1 --min-length 100 --max-length 2000 --report-file m54126_191025_082136.ccs.report.txt ../rawData/m54126_191024_220212.subreads.bam m54126_191024_220212.ccs.2.bam --chunk 2/10
ccs -j 32 --min-rq 0.9 --min-passes 1 --min-length 100 --max-length 2000 --report-file m54126_191025_082136.ccs.report.txt ../rawData/m54126_191024_220212.subreads.bam m54126_191024_220212.ccs.3.bam --chunk 3/10
ccs -j 32 --min-rq 0.9 --min-passes 1 --min-length 100 --max-length 2000 --report-file m54126_191025_082136.ccs.report.txt ../rawData/m54126_191024_220212.subreads.bam m54126_191024_220212.ccs.4.bam --chunk 4/10
ccs -j 32 --min-rq 0.9 --min-passes 1 --min-length 100 --max-length 2000 --report-file m54126_191025_082136.ccs.report.txt ../rawData/m54126_191024_220212.subreads.bam m54126_191024_220212.ccs.5.bam --chunk 5/10
ccs -j 32 --min-rq 0.9 --min-passes 1 --min-length 100 --max-length 2000 --report-file m54126_191025_082136.ccs.report.txt ../rawData/m54126_191024_220212.subreads.bam m54126_191024_220212.ccs.6.bam --chunk 6/10
ccs -j 32 --min-rq 0.9 --min-passes 1 --min-length 100 --max-length 2000 --report-file m54126_191025_082136.ccs.report.txt ../rawData/m54126_191024_220212.subreads.bam m54126_191024_220212.ccs.7.bam --chunk 7/10
ccs -j 32 --min-rq 0.9 --min-passes 1 --min-length 100 --max-length 2000 --report-file m54126_191025_082136.ccs.report.txt ../rawData/m54126_191024_220212.subreads.bam m54126_191024_220212.ccs.8.bam --chunk 8/10
ccs -j 32 --min-rq 0.9 --min-passes 1 --min-length 100 --max-length 2000 --report-file m54126_191025_082136.ccs.report.txt ../rawData/m54126_191024_220212.subreads.bam m54126_191024_220212.ccs.9.bam --chunk 9/10
ccs -j 32 --min-rq 0.9 --min-passes 1 --min-length 100 --max-length 2000 --report-file m54126_191025_082136.ccs.report.txt ../rawData/m54126_191024_220212.subreads.bam m54126_191024_220212.ccs.10.bam --chunk 10/10

pbmerge -o m54126_191024_220212.ccs.bam m54126_191024_220212.ccs.*.bam
pbindex m54126_191024_220212.ccs.bam
