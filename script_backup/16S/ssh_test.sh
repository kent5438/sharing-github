#! /bin/bash

ssh -t -t qiime@192.168.0.30<<EOF
cd /home/qiime/Desktop/OTU_analysis/MS17024-2
./ssh_test.sh
exit
EOF

touch aaaaaaaaaaaaaaaaaaaa.txt
