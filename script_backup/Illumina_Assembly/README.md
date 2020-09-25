# Illumina bacterial genome de-novo assembly pipeline

> __[Notice]__

> __1. Run on Snake only (SSD + large local space).__

***

## Run supernova_pipe.sh for Assembly
__This is very naive shell script version. Modify it before you run the pipe__

`$ sh /export/EC1680U/perl/bin/10x_Supernova/supernova_pipe.sh`

***

## Supernova Output
```
HS17026_NTNU-SHLi               <--- project root folder
├── 10x_supernova_report.pl
├── assembly                    <--- 'supernova run' prefix (recommended)
├── Bird_assembly_v2
├── CMD_history.txt
├── HLGCNBBXX_8
├── __HLGCNBBXX_8.mro
├── log
├── mkfastq                     <--- 'supernova mkfastq' prefix (recommended)
├── report
└── supernova_pipe.sh
```

***

## Generate Report by Rmarkdown
`$ perl /export/EC1680U/perl/bin/10x_Supernova/10x_supernova_report.pl [-a assembly] [-d mkfastq]`

-a: assembly (default)

-d: mkfastq (default)

***

## Report Output
```
report
├── css
├── files
├── html
├── images
├── index.html  <--- web-based report
└── site_libs
```

***

## Demo Report

[Demo Report](report_demo.zip)
