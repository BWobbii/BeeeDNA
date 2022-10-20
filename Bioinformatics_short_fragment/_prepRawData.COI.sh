#!/bin/bash

#script for COI data, preprocessing
#demultiplexing and primer trimming
#short and long possible
conda activate cutadapt__3.7

mkdir logs

echo
echo ====================================
echo Demultiplexing
echo ====================================

mkdir demux

tagF=($(grep "tagF" config.COI.txt | cut -f2 -d"=" | awk '{print $1*2}'))
tagR=($(grep "tagR" config.COI.txt | cut -f2 -d"=" | awk '{print $1*2}'))

head -n $tagF /mnt/data/homes/sickel/Dokumente/bin/tags/tags_fwd.fasta > tags_fwd.fa
head -n $tagR /mnt/data/homes/sickel/Dokumente/bin/tags/tags_rev.fasta > tags_rev.fa

#input file names - make sure they are same in subsetting step!
#or change them?
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
mkdir trim2
mkdir untrimmed2


  echo
  echo ====================================
  echo Trimming primers
  echo ====================================

prFn=($(grep "prFn" config.COI.txt | cut -f2 -d"=" ))
prRn=($(grep "prRn" config.COI.txt | cut -f2 -d"=" ))

prFr=($(grep "prFr" config.COI.txt | cut -f2 -d"=" ))
prRr=($(grep "prRr" config.COI.txt | cut -f2 -d"=" ))

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

    cutadapt \
    -a $prRr  \
    -A $prFr \
    -o trim2/$s.1.trim.fq \
    -p trim2/$s.2.trim.fq \
    -O 23 trim/$s.1.trim.fq trim/$s.2.trim.fq --untrimmed-output untrimmed2/$s.1.untrim.fq \
    --untrimmed-paired-output untrimmed2/$s.2.untrim.fq 2>&1 | tee -a logs/_trim.log

done
