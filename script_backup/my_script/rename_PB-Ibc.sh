#! /bin/bash

while read -r old new
do
    arr=( ${old}.ccs.fastq )
    if (( ${#arr[@]} == 1 ))
    then
        mv "${arr[0]}" "$new.ccs.fastq"
    else
        echo "Error: Multiple files found for $old: ${arr[@]}"
    fi    
done < list.txt
