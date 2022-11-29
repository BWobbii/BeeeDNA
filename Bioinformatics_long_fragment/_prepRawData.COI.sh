#!/bin/bash

conda activate vsearch__2.17.1

mkdir logs

echo
echo ====================================
echo Demultiplexing
echo ====================================

mkdir demux

tagF=($(grep "tagF" config.COI.txt | cut -f2 -d"=" | awk '{print $1*2}'))
tagR=($(grep "tagR" config.COI.txt | cut -f2 -d"=" | awk '{print $1*2}'))

head -n $tagF tags_fwd.fasta > tags_fwd.fa
head -n $tagR tags_rev.fasta > tags_rev.fa

cutadapt \
  -e 0.15 --no-indels \
  -O 7 \
  -g file:tags_fwd.fa \
  -G file:tags_rev.fa \
  -o demux/{name1}-{name2}.1.fq \
  -p demux/{name1}-{name2}.2.fq seqs_fwd.fq seqs_rev.fq 2>&1 | tee -a logs/_demux.log



  mkdir unknown
  mv demux/*unknown*.fq unknown/
  mkdir trim
  mkdir untrimmed


  echo
  echo ====================================
  echo Trimming primers
  echo ====================================


prFn=($(grep "prFn" config.COI.txt | cut -f2 -d"=" ))
prRn=($(grep "prRn" config.COI.txt | cut -f2 -d"=" ))



for f in demux/*.1.fq; do

    r=$(sed -e "s/.1.fq/.2.fq/" <<< "$f")
    s=$(sed -e "s/.*\/\(.*\).1.fq/\1/" <<< "$f")


    echo
    echo ====================================
    echo Trimming sample $s
    echo ====================================


    cutadapt \
    -g $prFn \
    -G $prRn \
    -o trim/$s.1.trim.fq \
    -p trim/$s.2.trim.fq \
    -O 23 $f $r --untrimmed-output untrimmed/$s.1.untrim.fq \
    --untrimmed-paired-output untrimmed/$s.2.untrim.fq 2>&1 | tee -a logs/_trim.log


done
