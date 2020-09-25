#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Hands-on Python v2: hifi assembly with racon polishing pipeline
# 2020/09/25

import os
import sys
import subprocess
import re
import argparse

import time
def curr_time():
    current_time = time.strftime("%Y-%m-%d %H:%M:%S")
    return current_time

def process_command():
    parser = argparse.ArgumentParser(description="### HIFI Assembly with Racon Polishing Pipeline ###")
    parser.add_argument('-c', '--ccs', help='ccs fastq path', required=True)
    parser.add_argument('-p', '--prefix', help='output prefix name', required=True)
    return parser.parse_args()
if __name__ == '__main__':
    print("\n### HIFI Assembly with Racon Polishing Pipeline ###\n")
    args = process_command()
    print('- input CCS fastq: ', args.ccs)
    print('- output prefix: ', args.prefix)
    pwd = subprocess.getoutput('pwd')
    print('- CWD: ', pwd)
    print(f"\n----- Start Analysis ({curr_time()}) -----\n")


asm_fasta = f'assembly.fasta'
if not os.path.exists(asm_fasta):
    print(f"### ERROR: '{asm_fasta}' not found! Flye assembly may not be completed !!!")
    sys.exit()

print(f"*** [{curr_time()}]: mapped HiFi reads to assembled contigs by pbmm2")
asm_pbmm2_bam = f'assembly.pbmm2.bam'
if not os.path.exists(asm_pbmm2_bam):
    os.system(f"/export/EC1680U/kentchen/miniconda3/envs/pbbioconda/bin/pbmm2 align --preset CCS --sort -j 32 {asm_fasta} {args.ccs} assembly.pbmm2.bam")

print(f"*** [{curr_time()}]: convert bam to sam")
asm_pbmm2_sam = f'assembly.pbmm2.sam'
if not os.path.exists(asm_pbmm2_sam):
    os.system(f"samtools view -F 1796 -q 20 {asm_pbmm2_bam} > {asm_pbmm2_sam}")

print(f"*** [{curr_time()}]: polishing raw assembled contigs by Racon")
polished_ctgs = f"{args.prefix}.polished.fasta"
if not os.path.exists(polished_ctgs):
    os.system(f"/export/EC1680U/kentchen/miniconda3/envs/racon/bin/racon -u -t 32 {args.ccs} {asm_pbmm2_sam} assembly.fasta > {polished_ctgs}")

print(f"\n----- Complete ({curr_time()}) -----\n")
 
