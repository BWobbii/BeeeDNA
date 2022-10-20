#!/bin/bash

#script for COI data, VSEARCH
#read joining, filtering, denoising, taxonomic classification, community table
#path to VSEARCH
conda activate vsearch__2.17.1


echo
echo ====================================
echo Merging reads
echo ====================================

mkdir merged
mkdir notmerged



for f in trim2/*.1.trim.fq; do

    r=$(sed -e "s/.1.trim.fq/.2.trim.fq/" <<< "$f")
    s=$(sed  -e "s/.*\/\(.*\).1.trim.fq/\1/" <<< "$f")


    echo
    echo ====================================
    echo Merging sample $s | tee -a logs/_merging.log
    echo ====================================


    vsearch --fastq_mergepairs $f \
        --reverse $r \
        --fastq_allowmergestagger \
        --fastq_minovlen  20 \
        --fastq_maxdiffs 10 \
        --fastqout merged/$s.merged.fq \
        --relabel R1+2-$s \
        --fastq_eeout \
        --fastqout_notmerged_fwd notmerged/$s.notmerged.fwd.fq \
        --fastqout_notmerged_rev notmerged/$s.notmerged.rev.fq 2>&1 | tee -a logs/_merging.log

done


echo
echo ====================================
echo Quality filtering
echo ====================================


mkdir qual
mkdir qual-disc



for f in merged/*.fq; do

    s=$(sed  -e "s/.*\/\(.*\).merged.fq/\1/" <<< "$f")

    echo
    echo ====================================
    echo Filtering sample $s
    echo ====================================
#relabel?
  vsearch --fastq_filter $f \
    --fastq_maxee 1.0 \
    --fastq_minlen 200 \
    --fastq_maxlen 400 \
    --fastq_maxns 0 \
    --fastaout qual/$s.qual.fa \
    --fasta_width 0 \
    --fastqout_discarded qual-disc/$s.disc.fq \
    --relabel $s. 2>&1 | tee -a logs/_qual.log


done


echo
echo ====================================
echo Dereplication
echo ====================================


mkdir derep

for f in qual/*.fa; do

  s=$(sed  -e "s/.*\/\(.*\).qual.fa/\1/" <<< "$f")

 vsearch --derep_fulllength $f \
  --strand plus \
  --output derep/$s.derep.fa \
  --sizeout \
  --relabel sample=$s. \
  --fasta_width 0 2>&1 | tee -a logs/_derep.log

done

#merge sample derep together
cat derep/*.fa > derep/all.derep.fa


vsearch --derep_fulllength derep/all.derep.fa \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --output all.fa 2>&1 | tee -a logs/_derep.log


    echo
    echo ====================================
    echo Denoising
    echo ====================================

#uses UNOISE version 3 algorithm
#minsize - default 8.0
#discarded seqs can be mapped back using otutab (VSEARCH?)
#unoise_alpha default 2.0
#increasing alpha trades sensitivity to differences
#against increase in number of bad seqs, which are wronlgy predicted to be good
#see UNOISE2 paper
#uchime3 denovo afterwards

vsearch --cluster_unoise all.fa \
        --sizein \
        --sizeout \
        --fasta_width 0 \
        --relabel zOTU. \
        --centroids all.denoise.fa 2>&1 | tee -a logs/_denoise.log


  #chimera check
  echo
  echo ====================================
  echo Chimera check
  echo ====================================

mkdir chimeras

#if troubleshooting needed, add
#--uchimealns filename

vsearch --uchime3_denovo all.denoise.fa \
      --nonchimeras all.chim.fa \
      --sizein \
      --sizeout \
      --chimeras chimeras/chimeras.fa \
      --fasta_width 0 \
      --uchimeout chim.results.txt 2>&1 | tee -a logs/_chimeras.log


#make community table
echo
echo ====================================
echo Make community table
echo ====================================


sed -i 's/\./;/g' derep/all.derep.fa

vsearch --usearch_global derep/all.derep.fa \
    --db all.chim.fa \
    --id 0.97 \
    --strand plus \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --otutabout asvtab.txt \
    --biomout asvtab.biom 2>&1 | tee -a logs/_asvtab.log


#taxonomic assignment, three steps
#1st: direct assignments with NBCI data (prepared with BC-databaser)
#2nd: dierct assignments with BOLD data
#both are Arthropod species in Germany
#3rd: hierarchical classification with SINTAX and MIDORI_UNIQ_GB241_CO1_SINTAX
echo
echo ====================================
echo Taxonomic classification
echo ====================================



refDBs=($(grep "refdb" config.COI.txt | cut -f2 -d"=" | sed 's/\"//g'))
hieDBs=($(grep "hiedb" config.COI.txt | cut -f2 -d"=" | sed 's/\"//g'))

threshold=97

echo ",kingdom,phylum,class,order,family,genus,species" > taxonomy.vsearch

countdb=0
cp  all.chim.fa zotus.direct.$countdb.uc.nohit.fasta
prevDB=$countdb

touch taxonomy.vsearch


for db in "${refDBs[@]}"
  do :
    countdb=$((countdb+1))
    echo "\n\n#### Direct VSEARCH Classification level: $countdb";

    conda activate vsearch__2.17.1

   vsearch --usearch_global zotus.direct.$prevDB.uc.nohit.fasta --db $db --id 0.$threshold --uc zotus.direct.$countdb.uc --fastapairs zotus.direct.$countdb.fasta --strand both 2>  logs/_direct.$countdb.log
    grep "^N[[:space:]]" zotus.direct.$countdb.uc | cut -f 9 > zotus.direct.$countdb.uc.nohit

    conda activate seqkit__2.2.0

    seqkit grep -w 0 -n -f zotus.direct.$countdb.uc.nohit all.chim.fa -o zotus.direct.$countdb.uc.nohit.fasta
    cut -f 9,10 zotus.direct.$countdb.uc  | grep -v "*" | sed "s/[A-Za-z0-9_]*;tax=//" >> taxonomy.vsearch
    prevDB=$countdb
  done



echo "\n\n#### Hierarchical VSEARCH classification";
conda activate vsearch__2.17.1

vsearch --sintax zotus.direct.$countdb.uc.nohit.fasta -db $hieDBs -tabbedout zotus.uc.merge.nohit.sintax -strand plus -sintax_cutoff 0.8 2>  logs/_sintax.log

cut -f1,4 zotus.uc.merge.nohit.sintax | sed -E -e "s/\_[0-9]+//g" -e "s/,s:.*$//"  >> taxonomy.vsearch
sed -i 's/;.*\t/,/' taxonomy.vsearch
sed -i 's/;$//' taxonomy.vsearch
