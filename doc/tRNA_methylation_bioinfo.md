
# tRNA methylation bioinformatics

- [tRNA methylation bioinformatics](#trna-methylation-bioinformatics)
  - [introduction](#introduction)
  - [genome and annotation](#genome-and-annotation)
    - [Gencode](#gencode)
    - [ENSEMBL](#ensembl)
    - [masking tRNAs in the genome](#masking-trnas-in-the-genome)
    - [extracting the tRNA sequences](#extracting-the-trna-sequences)
    - [getting unique sequences only](#getting-unique-sequences-only)
    - [creating an artificial genome](#creating-an-artificial-genome)
  - [mapping](#mapping)
    - [create genomic index](#create-genomic-index)
    - [input fastq files](#input-fastq-files)
    - [clustering and low complexity filter](#clustering-and-low-complexity-filter)
    - [quality check and optical replicates filtering](#quality-check-and-optical-replicates-filtering)
    - [optional step: fqgrep primer masking](#optional-step-fqgrep-primer-masking)
    - [mapping with lastal](#mapping-with-lastal)
  - [parsing MAF format files](#parsing-maf-format-files)
    - [getting tRNA contig matches with counts](#getting-trna-contig-matches-with-counts)
    - [getting two aligment lines (tRNA contig and sequencing read matches)](#getting-two-aligment-lines-trna-contig-and-sequencing-read-matches)
  - [priority tRNA set processing](#priority-trna-set-processing)
    - [extract priority tRNAs](#extract-priority-trnas)
  - [tRNA position numbering](#trna-position-numbering)
    - [html to STOCKHOLM aligment format](#html-to-stockholm-aligment-format)
    - [inferal for database of profiles](#inferal-for-database-of-profiles)
      - [create a ```.cm``` file using ```cmbuild```](#create-a-cm-file-using-cmbuild)
      - [calibrate profiles using ```cmcalibrate```](#calibrate-profiles-using-cmcalibrate)
      - [convert to cmsearch compatible db with ```cmpress```](#convert-to-cmsearch-compatible-db-with-cmpress)
      - [search profiles with filetered read matches](#search-profiles-with-filetered-read-matches)
      - [pre-parse cmscan tabular output](#pre-parse-cmscan-tabular-output)
      - [filtration of the cmscan matches](#filtration-of-the-cmscan-matches)

<small><i><a href='http://ecotrust-canada.github.io/markdown-toc/'>Table of contents generated with markdown-toc</a></i></small>




## introduction

For a successful and timely analysis we used computing cluster managed using SLURM. 
While vast majority of the steps can be performed on a mainstream workstation creation of the genomic database for LAST and fastq mapping should be done on machines with > 32G RAM.


## genome and tRNAs

### T2T genome assembly

```
# download
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_
set/chm13v2.0.fa.gz

# decompress and fix chromosome names

pigz -d --stdout chm13v2.0.fa.gz | sed 's/^>chr/>/g; s/^>M/>MT/' > chm13v2.0.no_chr
.fa
```

### mature tRNAs
```
wget http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-mature-tRNAs.fa

sed '/^>/s/Homo_sapiens_//g; /^>/s/tRNAscan-SE ID: //g ;/^>/!s/U/T/g' hg38-mature-t
RNAs.fa 

```

### tRNAscan-SE predictions


* create region files for chromosomes/MT extraction
``` 
rg -v '^MT' chm13v2.0.no_chr.fa.fai | awk '{print $1":1-"$2}' > chrom_1-22,XY.extract.
list

rg '^MT' chm13v2.0.no_chr.fa.fai | awk '{print $1":1-"$2}' > chrom_MT.extract.list
```

* extract 24 chromosomes and MT to separate fasta files

```
samtools faidx chm13v2.0.no_chr.fa --region-file chrom_MT.extract.list > MT.fa

samtools faidx chm13v2.0.no_chr.fa --region-file chrom_1-22,XY.extract.list > chrom
_1-22,XY.fa

```

* run tRNAscan

```
tRNAscan-SE --threads 16 \
--detail \
--output chm13_2.0.no_chr.trnas.out \
--bed chm13_2.0.no_chr.trnas.bed \
--fasta chm13_2.0.no_chr.trnas.fa \
--log chm13_2.0.no_chr.trnas.log chm13_2.0.no_chr.fa
```

```
tRNAscan-SE --threads 2 \
--detail \
-M mammal \
--output MT_.trnas.out \
--bed MT.trnas.bed \
--fasta MT.trnas.fa \
--log MT.trnas.log chr_MT.fa
```

### Gencode 

Release used: 40 

* Obtained genome fasta and tRNA annotation GTF files:

```sh

wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/GRCh38.primary_assembly.genome.fa.gz
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.tRNAs.gtf.gz

```

* gunzip

```sh
gunzip GRCh38.primary_assembly.genome.fa.gz
gunzip gencode.v40.tRNAs.gtf.gz
```

* fixed chromosome names and removed annotation for MT in GTF

```sh
# ripgrep (rg) can be replaced with grep

rg -v '^chrM' gencode.v40.tRNAs.gtf | sed 's/^chr//g' > gencode.v40.tRNAs.no_chrM.gtf
```

* fixed chr names in genomic fasta

separate Python/pypy3 script
pyfaidx need to be installed, see: https://pypi.org/project/pyfaidx/

```sh
./fix_chr_names.py > GRCh38.primary_assembly.genome.names_fix.fa

```

### ENSEMBL

Release used: 106

* obtain GTF annotation

```sh
wget http://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz
```

* extract mitochondrial tRNA annotations

```sh

rg -z '^MT' Homo_sapiens.GRCh38.106.gtf.gz | rg '\tgene\t' | rg -w Mt_tRNA > ensembl_106.MT_tRNA.gtf

```

### masking tRNAs in the genome

```sh
# create combined GTF for masking

cat gencode.v40.tRNAs.no_chrM.gtf ensembl_106.MT_tRNA.gtf > gencode_ensembl.combined_tRNA.gtf 

# mask the genome

bedtools maskfasta -fi GRCh38.primary_assembly.genome.names_fix.fa \
-bed gencode_ensembl.combined_tRNA.gtf \
-fo GRCh38.primary_assembly.genome.names_fix.masked_tRNAs.fa  

```

### extracting the tRNA sequences

* modify the GTF to get gene names in the second column

separate script Python

```sh
./preprocess_gtf_4_fasta.py gencode_ensembl.combined_tRNA.gtf > gencode_ensembl.combined_tRNA.4_extraction.gtf
```

* extract tRNA genes and pseudogenes from the unmasked genome using bedtools

The '''-s''' switch ensures that the extracted sequence is in the correct orientation

```sh
 bedtools getfasta -s -name \
 -fi GRCh38.primary_assembly.genome.names_fix.fa \
 -bed gencode_ensembl.combined_tRNA.4_extraction.gtf  > gencode_ensembl.combined_tRNA.all.fa

```

* result

670 tRNA seq
556 non Pseudo tRNAs
22 MT tRNAs

Names have following names schemes:

* Gencode

```txt
>Pseudo_tRNA::1:7930278-7930348(-)
>Asn_tRNA::1:16520584-16520658(-)
>Asn_tRNA::1:16532397-16532471(-)
>Glu_tRNA::1:16535278-16535350(-)

```

* ENSEMBL MT tRNAs
  
```txt
>MT-TD::MT:7517-7585(+)
>MT-TK::MT:8294-8364(+)
>MT-TG::MT:9990-10058(+)
>MT-TR::MT:10404-10469(+)
>MT-TH::MT:12137-12206(+)
>MT-TS2::MT:12206-12265(+)
>MT-TL2::MT:12265-12336(+)
>MT-TE::MT:14673-14742(-)
>MT-TT::MT:15887-15953(+)
>MT-TP::MT:15955-16023(-)

```

### getting unique sequences only


Use usearch for clustering identical sequences.



```sh
usearch -cluster_fast gencode_ensembl.combined_tRNA.all.fa  -id 1.00 -centroids gencode_ensembl.combined_tRNA.uniq.fa
```

After clustering all 22 MT tRNA names are present. Pseudo_tRNAs drop from 114 to 109.

**Caveat**
Since there was a possibility that some expressed, true tRNAs will cluster to some Pseudo_tRNA I did a check:

```sh
#extract the Pseudo_tRNAs
rg -A1 '^>Pseudo' gencode_ensembl.combined_tRNA.all.fa | rg -v '^\-' > gencode_ensembl.combined_tRNA.pseudo.fa

# cluster Psuedo_tRNAs only:
usearch -cluster_fast gencode_ensembl.combined_tRNA.pseudo.fa  -id 1.00 -centroids gencode_ensembl.combined_tRNA.pseudo.uniq.fa

# check the number of seq

rg -c '^>'  gencode_ensembl.combined_tRNA.pseudo.uniq.fa
```

After Psuedo_tRNA sequence clustering the number of sequences is also 109  

```sh
# extract non Pseudo tRNAs
./extract_nonPseudo_fasta.py gencode_ensembl.combined_tRNA.all.fa > gencode_ensembl.combined_tRNA.non_pseudo.fa

# cluster 
usearch -cluster_fast gencode_ensembl.combined_tRNA.non_pseudo.fa -id 1.00 -centroids gencode_ensembl.combined_tRNA.non_pseudo.uniq.fa

#check the number of sequences after clustering 

rg -c '^>' gencode_ensembl.combined_tRNA.non_pseudo.uniq.fa

```

393 unique non-Pseudo tRNAs => no true tRNA clusters with Pseudo_tRNA and "vanishes".


### fixing tRNA sequence names

The contig identifiers containing i.e. round brackets are not compatible with IGV viewer. Also the mitochondrial tRNAs names have a different naming convention.
To resolve issues downstream these names were changed:

```
sed -e '/^>MT/ s/_MT:/_tRNA_MT:/g'  gencode_ensembl.combined_tRNA.uniq.fa | sed  '/^>/ s/::/_/g; s/(+)/F/g ; s/(-)/R/g;  s/SeC(e)/SeCe/g' >  
gencode_ensembl.combined_tRNA.uniq.name_fix.fa

#rename file for the compability with the downstream docs 
mv gencode_ensembl.combined_tRNA.uniq.name_fix.fa gencode_ensembl.combined_tRNA.uniq.fa
```


### creating an artificial genome

This is to combina the genomic sequence with masked (using tRNA GTF annotations) with extracted unique tRNA (and Pseudo_tRNA sequences just in case if these are also expresed).

```sh
# concatenate two fasta files
cat GRCh38.primary_assembly.genome.names_fix.masked_tRNAs.fa  gencode_ensembl.combined_tRNA.uniq.fa >  GRCh38_plus_tRNAs.fa


# fix different sequence lenghts for downstream samtools faidx
# used fastareformat from exonerate
### fastareformat GRCh38_plus_tRNAs.ordered_cca.fa > GRCh38_plus_tRNAs.ordered_cca_ref.fa

fastareformat GRCh38_plus_tRNAs.fa > GRCh38_plus_tRNAs.fa.tmp
mv GRCh38_plus_tRNAs.fa.tmp GRCh38_plus_tRNAs.fa

# compress for transfer to HPC
pigz -8 GRCh38_plus_tRNAs.fa
```


## mapping

Used LAST mapper. Installed from source:

```sh
module load gcc/11.2.0
make 
```


### create genomic index

* unpack

```sh
# gunzip fasta
pigz -d GRCh38_plus_tRNAs.fa.gz

```

* slurm script

```sh
#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=24
#SBATCH --partition=mem

module load gcc/11.2.0

export PATH=/scratch/dkedra/soft/bin:$PATH

lastdb -P24 -uNEAR GRCh38_plus_tRNAs.last GRCh38_plus_tRNAs.fa

```

### input fastq files

* sizes/names

```sh
# command
ls -l *fastq.gz 

# output
-r--r--r-- 1 darked darked 421363708 Apr 25 16:25 D_1.fastq.gz
-r--r--r-- 1 darked darked 442198262 Apr 25 16:25 D_BH4_1.fastq.gz

-rw-r--r-- 1 darked darked 458572848 Apr 26 13:16 NoD_1.fastq.gz
-r--r--r-- 1 darked darked 410977493 Apr 25 15:31 NoD_BH4_1.fastq.gz
```

* md5 checksums

```sh
# command
md5sum *fastq.gz 

# result
bb2c0c03a31c47077b4136b599374f51  D_1.fastq.gz
ec37b58216ecc8958df2b4e876133762  D_BH4_1.fastq.gz
68ecdb9735be757c7d53cb556e5258fa  NoD_1.fastq.gz
f8d3140c32e32b7df4888f9697de9923  NoD_BH4_1.fastq.gz
```

* location on HPC cluster

```sh
/scratch/dkedra/proj/trna_20220426/FQ_data/ORIG
```

### clustering and low complexity filter

To shrink size and speed up subsequent processing of the fastq data used clumpify.sh from BBMap 

example SLURM shell script:

```sh
#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=16
#SBATCH --partition=mem

#SBATCH --job-name=job_name

# module load gcc/11.2.0
module load java/11

export PATH=$HOME/soft/bin:$PATH

/scratch/dkedra/soft/progs/bbmap_current/clumpify.sh \
dedupe=t \
optical=t \
reorder=a \
shortname=shrink \
blocksize=2048 \
ziplevel=8 \
pigz=6 \
unpigz=6 \
lowcomplexity=t \
in=./ORIG/D_1.fastq.gz \
out=.//CLUMP_20220426/D_1.clump_opt_dedup.fq.gz
```

Python script for SLURM scripts creation:

```sh
#dir: /scratch/dkedra/proj/trna_20220426/FQ_data

./batch_clumpify.py

```

### quality check and optical replicates filtering

Used fastp.

* example SLURM submission script:

```sh
#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --partition=express

#SBATCH --job-name=job_name

module load conda/current
conda activate fastp


fastp \
--in1=./CLUMP_20220426/D_1.clump_opt_dedup.fq.gz \
--out1=./FASTP_20220426/D_1.clump_opt_dedup.fastp.fq.gz \
--thread=4 \
--json=./FASTP_20220426/D_1.clump_opt_dedup.fastp.json \
--html=./FASTP_20220426/D_1.clump_opt_dedup.fastp.html \
--report_title=D_1.clump_opt_dedup_report \
--overrepresentation_analysis \
--low_complexity_filter \
--disable_length_filtering \
--compression=6

```

Python script for SLURM scripts creation:

```sh
#dir: /scratch/dkedra/proj/trna_20220426/FQ_data

./batch_fastp.py
```

### optional step: fqgrep primer masking

In order to improve mappings/not include primer derived base(s) i.e. at the matches ends we masked the primer sequences in the fastq files.

* program used: https://github.com/indraniel/fqgrep

* example SLURM script to create fqgrep report

```sh
#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=00:59:00
#SBATCH --cpus-per-task=1
#SBATCH --partition=express

#SBATCH --job-name=job_name

module load gcc/11.2.0 

export LD_LIBRARY_PATH=/scratch/dkedra/soft/lib:$LD_LIBRARY_PATH
export PATH=/scratch/dkedra/soft/bin:$PATH


fqgrep -r -a -e -p 'TGGAATTCTCGGGTGCCAAGGC|TGGAATTCTCGGGTGCCAAGG|TGGAATTCTCGGGTGCCAAG|TGGAATTCTCGGGTGCCAA|TGGAATTCTCGGGTGCCA|GTTCAGAGTTCTACAGTCCGACGATC|GTTCAGAGTTCTACAGTCCGACGAT|GTTCAGAGTTCTACAGTCCGACGA|GTTCAGAGTTCTACAGTCCGACG|GTTCAGAGTTCTACAGTCCGA
C|GCCTTGGCACCCGAGAATTCCA|GATCGTCGGACTGTAGAACTCTGAAC' -o ./FQGREP_20220426/D_1.clump_opt_dedup.fastp.fqgrep_report_pat2.out  ./FASTP_20220426/D_1.clump_opt_dedup.fastp.fq.gz
```

* python script to create shell scripts for slurm

```batch_fqgrep.py```

* python/pypy script to re-create fastq files from fqgrep report

```src/parse_fqgrep_report.py```



### mapping with lastal

* example SLURM shell code

```sh
#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=48
#SBATCH --partition=mem

#SBATCH --job-name=lastal_D_1.clump_opt_dedup.fastp.fqgrep_mask

module load gcc/11.2.0 

export PATH=/scratch/dkedra/soft/bin:$PATH


lastal -v -P48 -Qkeep -C2 genome_last/hg38_tRNAs/GRCh38_plus_tRNAs.last ./IN_fqgrep/D_1.clump_opt_dedup.fastp.fqgrep_mask.fq.gz | last-split | gzip > ./D_1.clump_opt_dedup.fastp.fqgrep_mask.lastal.hg38-tRNAs.maf.gz 

```

## parsing MAF format files 

**Caveat** 

MAF format has a 0-based numbering scheme for positions

### getting tRNA contig matches with counts

From a given ```*.maf.gz``` file extract the aligment lines containing names (=> ripgrep patterns) specified in 
the file: ```meta_data/trna_names_short.with_s.txt```

**Caveat**
1. contains data for all tRNA contigs with matches, including 
* mitochondrial tRNAs
* Pseudo_tRNAs 
* SeC(e)_tRNA, SeC_tRNA, Sup_tRNA, Undet_tRNA

```txt
count	seq_name	match_start	match_end	seq_len	match_seq
14	Ala_tRNA::11:50274707-50274779(-)	1	26	72	GGGGGTGTAGCTCAGTGGTAGAGCGG
1	Ala_tRNA::11:50274707-50274779(-)	1	40	72	GGGGGTGTAGCTCAGTGGTAGAGCGGATGCTTTGCATGTA
1	Ala_tRNA::12:124939965-124940037(+)	1	37	72	GGGGATGTAGCTCAGTGGTAGAGCGCATGCTTTGCAC
17	Ala_tRNA::12:124939965-124940037(+)	1	40	72	GGGGATGTAGCTCAGTGGTAGAGCGCATGCTTTGCACGTA
21	Ala_tRNA::12:124939965-124940037(+)	1	41	72	GGGGATGTAGCTCAGTGGTAGAGCGCATGCTTTGCACGTAT

```

* python batch script ```src/batch_extract_tRNA_mappings_maf.py```
* example SLURM script

```sh
#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=00:50:00
#SBATCH --cpus-per-task=4
#SBATCH --partition=express

#SBATCH --job-name=rg_trna_D_1.clump_opt_dedup.fastp.lastal.hg38-tRNAs

module load gcc/11.2.0 

export PATH=/scratch/dkedra/soft/bin:$PATH


cp header.txt.gz .//D_1.clump_opt_dedup.fastp.lastal.hg38-tRNAs.trna.frag_cnt.tsv.gz

rg -z --threads 4 -f /scratch/dkedra/proj/trna_20220426/mapping/trna_names_short.with_s.txt ./D_1.clump_opt_dedup.fastp.lastal.hg38-tRNAs.maf.gz \
| awk '{print $2, $3+1, $3+$4, $6, $7}' \
| sort | uniq -c |  tr --squeeze-repeats " " \
| sed 's/ /\t/g' | sed -e 's/^\t//g' 
| gzip >> .//D_1.clump_opt_dedup.fastp.lastal.hg38-tRNAs.trna.frag_cnt.tsv.gz 

```

### getting two aligment lines (tRNA contig and sequencing read matches) 

**Note**
Using ```pigz``` instead of ripgrep's ```rg -z``` is faster. 
In fish shell
```sh
for fn in *maf.gz 
    pigz --stdout -d $fn | rg -A1 --no-context-separator  -f ../../trna_names_short.with_s.txt |  choose --character-wise 2.. > (basename $fn .maf.gz).trna_matches.2lines
end
```


## priority tRNA set processing

### extract priority tRNAs

Get the aligments for a subset of tRNAs specified in ```meta_data/priority_matches_01.txt```


* preprocessing/selection

fish shell command 


```sh
for fn in *.maf.gz 
    rg -z -A 1 -f priority_matches.txt $fn | gzip > (basename $fn .maf.gz).priority.align.gz
end
```

* output format

```txt
s Val_tRNA::3:169772229-169772302(+)        0 35 + 73 GTTTCCGTAGTGTAGTGGTTATCACGTTCGCCTAA
s SND00105:1227:CCV44ANXX:4:2309:7593:87815 4 35 + 51 GTTTCCGTAGTGTAGTGGTTATCAAGTTCGCCTAA
--
s Val_tRNA::3:169772229-169772302(+)         0 35 + 73 GTTTCCGTAGTGTAGTGGTTATCACGTTCGCCTAA
s SND00105:1227:CCV44ANXX:4:2209:13167:79726 4 35 + 51 GTTTCCGTAGTGTAGTGGTTATCAAGTTCGCCTAA
--
```

* convert subset of maf to fasta

* script: ```src/maf_2_fa_with_counts.py```

* command (fish)

```sh
for fn in *align
    echo $fn
    ./parse_maf_trna_alignments.py $fn > (basename $fn .align).align.10count.fa
end

```


* example output

```txt
>Ala_tRNA::12:124939965-124940037(+)_46_73_count747
CCCCGGGTTCAATCCCCGGCATCTCCA
>Ile_tRNA::6:26554121-26554195(+)_48_75_count1149
TCGCGGGTTCGATCCCCGTACGGGCCA
>Ala_tRNA::6:26553502-26553574(+)_45_73_count470
GTCCCGGGTTCGATCCCCGGCATCTCCA
>Ala_tRNA::6:26553502-26553574(+)_46_73_count3559
TCCCGGGTTCGATCCCCGGCATCTCCA
```


## tRNA position numbering

The canonical numebring scheme is based on the aligments from:
http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/Hsapi38-align.html


### html to STOCKHOLM aligment format

Above html was transformed to text and then to Stockholm aligment format in a following way:

* text processing from html to 2 colum aligments

```sh
wget --output-document=gtrna_align.raw.html  http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/Hsapi38-align.html

sed 's/<span ID=b.>//g' gtrna_align.raw.html |sed 's/<span ID=b..>//g' | sed 's/<\/span>//g' | rg -v ">>" | rg -v '^<' | sed 's/<a name//g' | sed 's/><\/a>//g' | tr -d '=' > gtrna_align.ver1.txt

# create file with desired tRNA names to filter aligments
rg -v 'tRNA covariance' gtrna_align.ver1.txt | choose 1 | rg '^tRNA' | choose -c ..12 | rg -v 'Und|SeC|Sup|iMet' | sort | uniq > gtrna_isotypes_names.51.txt

# extract the lines with relevant tRNA isotypes 
rg -f gtrna_isotypes_names.51.txt gtrna_align.ver1.txt | rg -v filtered | choose ..2 > gtrna_align.isotypes.51.2cols
```

This gives a 2 column file ```gtrna_align.isotypes.51.2cols``` with 422 lines

* alternative more inclusive set

```sh
rg -v '^[\t," "]' gtrna_align.ver1.txt | rg -v 'filtered' | rg 'tRNA-' | choose ..2 >  gtrna_align.isotypes.all.2cols

```

output: ```gtrna_align.isotypes.all.2cols``` has 433 lines.


* from 2 column aligment to stockholm

```sh
# canonical 51 set:
./tab_gtrna_to_stockholm.py gtrna_align.isotypes.51.2cols > gtrna_align.isotypes.51.sto

# alternative all
./tab_gtrna_to_stockholm.py gtrna_align.isotypes.all.2cols > gtrna_align.isotypes.all.sto

```


### inferal for database of profiles

#### create a ```.cm``` file using ```cmbuild```

```sh
# canonical 51 set:
cmbuild --noss -F gtrna_align.isotypes.51.cm  gtrna_align.isotypes.51.sto

# alternative all
cmbuild --noss -F gtrna_align.isotypes.all.cm  gtrna_align.isotypes.all.sto
```

#### calibrate profiles using ```cmcalibrate```

**Caveat** 

This is time/CPU intensive step. Takes about 1.5-2mins per aligment on Ryzen 7 with 16 threads.
On mem machines with 48 threads calibrating one aligment requires 20-30s.  

* example commands

```sh
# canonical 51 set:
cmcalibrate --cpu 48 --gtrna_align.isotypes.51.cm

# alternative all
cmcalibrate --cpu 48 gtrna_align.isotypes.all.cm 
```

* SLURM script for mem partition

```sh
#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=48
#SBATCH --partition=mem

module load conda

conda activate infernal


cmcalibrate --cpu 48 gtrna_align.isotypes.all.cm 

```


#### convert to cmsearch compatible db with ```cmpress```

This works in seconds:

```txt
# canonical 51 set:
cmpress gtrna_align.isotypes.51.sto

# alternative all
cmpress gtrna_align.isotypes.all.cm
```

* output example:

```txt
-rw-r--r-- 1 dkedra ubbiomed  625734 May 16 11:54 gtrna_align.isotypes.all.cm.i1p
-rw-r--r-- 1 dkedra ubbiomed 1753630 May 16 11:54 gtrna_align.isotypes.all.cm.i1m
-rw-r--r-- 1 dkedra ubbiomed    2082 May 16 11:54 gtrna_align.isotypes.all.cm.i1i
-rw-r--r-- 1 dkedra ubbiomed  399759 May 16 11:54 gtrna_align.isotypes.all.cm.i1f


```

#### search profiles with filetered read matches 

One can speed up cmscan using ```--cpu {num_of_threads}``` 


* example command

```sh
cmscan -o D_1.clump_opt_dedup.fastp.lastal.hg38-tRNAs.priority.align.cmscan_tmp  --tblout D_1.clump_opt_dedup.fastp.lastal.hg38-tRNAs.priority.align.10count.cmscan_51_out ../gtrna_align.isotypes.51.cm  D_1.clump_opt_dedup.fastp.lastal.hg38-tRNAs.priority.align.10count.fa
```

* fish command to search all 4 fasta files

```sh
for fn in *fa 
    echo $fn
    cmscan -o (basename $fn .fa).cmscan_tmp --tblout (basename $fn .fa).cmscan_isoall.out ./gtrna_align.isotypes.all.cm $fn 
end
```


#### pre-parse cmscan tabular output

Because the cmscan reports by default not just top hits for a given query the output needs to be pre-parsed to get just the top hit

* script: ```src/parse_top_hits_cmscan.py```

* example command

```sh
./parse_top_hits_cmscan.py D_1.clump_opt_dedup.fastp.lastal.hg38-tRNAs.priority.align.10count.cmscan51_out | sort >  D_1.clump_opt_dedup.fastp.lastal.hg38-tRNAs.priority.align.10count.cmscan51_out.top_hits
```

#### filtration of the cmscan matches

There are two requirements to include fragment mappings into the final stats/ heat plot:

1. each fragment was observed at least 10x in the LAST aligner mappings per fastq file
2. the cmscan E-value <= 1e-04

```txt
min_count = 10
min_score = 1e-04
```

