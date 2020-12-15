#!/bin/sh

# multijoin - join multiple files

# Usages: sh multijoin.sh *.genome.count.txt > all.genome.count.txt

join_rec() {
    if [ $# -eq 1 ]; then
        join - "$1"
    else
        f=$1; shift
        join - "$f" | join_rec "$@"
    fi
}

if [ $# -le 2 ]; then
    join "$@"
else
    f1=$1; f2=$2; shift 2
    join "$f1" "$f2" | join_rec "$@"
fi
