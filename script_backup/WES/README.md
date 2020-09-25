# WES pipeline (customized for 台大醫胡務亮)

> __[Notice]__

> __1. Using SGE cluster, please make sure SGE is keeping alive (Chicken / Sheep / Horse).__

> __2. The mainly shell script executes "GATK_WES.ntuh.pl" automatically by nohup and distributes jobs by SGE__

## [Job running]

### 1. Run "WES_pipe.sh" at anywhere (using NP17027 as example)
`$ sh /export/EC1680U/perl/bin/WES/WES_pipe.sh -i NP17027`

> __If you want to terminate jobs, using "qdel -f" instead of killing nohup.__

### 2. Files deployment

```
NP17027
├── WES-WJY-022A
│   ├── GATK.cfg
│   ├── GATK_WES.ntuh.pl
│   ├── list.txt
│   └── reads
│       ├── WES-WJY-022A.R1.fastq.gz -> /export/EC2480U/WES/NP17027/cleanData/WES-WJY-022A_R1_001.clean.fastq.gz
│       └── WES-WJY-022A.R2.fastq.gz -> /export/EC2480U/WES/NP17027/cleanData/WES-WJY-022A_R2_001.clean.fastq.gz
├── WES-WJY-022B
│   ├── GATK.cfg
│   ├── GATK_WES.ntuh.pl
│   ├── list.txt
│   └── reads
│       ├── WES-WJY-022B.R1.fastq.gz -> /export/EC2480U/WES/NP17027/cleanData/WES-WJY-022B_R1_001.clean.fastq.gz
│       └── WES-WJY-022B.R2.fastq.gz -> /export/EC2480U/WES/NP17027/cleanData/WES-WJY-022B_R2_001.clean.fastq.gz
└── WES-WJY-022C
    ├── GATK.cfg
    ├── GATK_WES.ntuh.pl
    ├── list.txt
    └── reads
        ├── WES-WJY-022C.R1.fastq.gz -> /export/EC2480U/WES/NP17027/cleanData/WES-WJY-022C_R1_001.clean.fastq.gz
        └── WES-WJY-022C.R2.fastq.gz -> /export/EC2480U/WES/NP17027/cleanData/WES-WJY-022C_R2_001.clean.fastq.gz

6 directories, 15 files
```

## 3. Log-in to VNC (Rabbit or Horse) and upload file to NTUH ftp server

### a.) Log-in information
- Host: colo-ftp1.ntuh.gov.tw
- Username: genomicsupload
- Password: onlyupload

### b.) default root path: /genomics
- 右建Create Directory -> WES-xxx -> cd /genomics/WES-xxx
(所有的新創建的資料夾都不會自動顯示，必須要在路徑的地方打上才能進入並出現)

[images.png](https://i.imgur.com/UyHGyQZ.png)

### c.) Create Folder Name & Upload
**- BAM (12 files)**

- *.sorted.dedup.realign.recal.bam|bai (3.bqsr)

- *.bamout.bam|bai (4.hc)

**- fastq (6 files)**

- *.R1|R2.fastq.gz (../reads)

**- vcf (6 files)**

- *.vcf|vcf.idx (4.hc)

**- report (3 files)**

- *_report.zip (8.report)
    
    (need to compress 8.report)