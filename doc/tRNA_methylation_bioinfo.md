# tRNA methylation

## genome and annotation

### Gencode 

Release used: 40 

* Obtained genome fasta and tRNA annotation GTF files:

```(sh)
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/GRCh38.primary_assembly.genome.fa.gz
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.tRNAs.gtf.gz
```

* gunzip

```
gunzip GRCh38.primary_assembly.genome.fa.gz
gunzip gencode.v40.tRNAs.gtf.gz
```

* fixed chromosome names and removed annotation for MT in GTF

```(sh)
# ripgrep (rg) can be replaced with grep

rg -v '^chrM' gencode.v40.tRNAs.gtf | sed 's/^chr//g' > gencode.v40.tRNAs.no_chrM.gtf
```

* fixed chr names in genomic fasta

separate Python/pypy3 script
pyfaidx need to be installed, see: https://pypi.org/project/pyfaidx/

```(sh)
./fix_chr_names.py > GRCh38.primary_assembly.genome.names_fix.fa

```

### ENSEMBL

Release used: 106

* obtain GTF annotation

```(sh)
wget http://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz
```

* extract mitochondrial tRNA annotations

```(sh)
rg -z '^MT' Homo_sapiens.GRCh38.106.gtf.gz | rg '\tgene\t' | rg -w Mt_tRNA > ensembl_106.MT_tRNA.gtf

```

### masking tRNAs in the genome

```
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

```
./preprocess_gtf_4_fasta.py gencode_ensembl.combined_tRNA.gtf > gencode_ensembl.combined_tRNA.4_extraction.gtf
```

* extract tRNA genes and pseudogenes from the unmasked genome using bedtools

The '''-s''' switch ensures that the extracted sequence is in the correct orientation

```
 bedtools getfasta -s -name -fi GRCh38.primary_assembly.genome.names_fix.fa -bed gencode_ensembl.combined_tRNA.4_extraction.gtf  > gencode_ensembl.combined_tRNA.all.fa

```

* result 

670 tRNA seq 
556 non Pseudo tRNAs
22 MT tRNAs

Names have following names schemes:

* Gencode

```
>Pseudo_tRNA::1:7930278-7930348(-)
>Asn_tRNA::1:16520584-16520658(-)
>Asn_tRNA::1:16532397-16532471(-)
>Glu_tRNA::1:16535278-16535350(-)

```

* ENSEMBL MT tRNAs
```
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



### getting uniqe sequences only


Use usearch for clustering identical sequences.

Version used: usearch11.0.667_i86linux32
Obtained from: https://drive5.com/downloads/usearch11.0.667_i86linux32.gz

```
usearch -cluster_fast gencode_ensembl.combined_tRNA.all.fa  -id 1.00 -centroids gencode_ensembl.combined_tRNA.uniq.fa
```

After clustering all 22 MT tRNA names are present. Pseudo_tRNAs drop from 114 to 109.

**Caveat**
Since there was a possibility that some expressed, true tRNAs will cluster to some Pseudo_tRNA I did a check:

```
#extract the Pseudo_tRNAs
rg -A1 '^>Pseudo' gencode_ensembl.combined_tRNA.all.fa | rg -v '^\-' > gencode_ensembl.combined_tRNA.pseudo.fa

# cluster Psuedo_tRNAs only:
usearch -cluster_fast gencode_ensembl.combined_tRNA.pseudo.fa  -id 1.00 -centroids gencode_ensembl.combined_tRNA.pseudo.uniq.fa

# check the number of seq

rg -c '^>'  gencode_ensembl.combined_tRNA.pseudo.uniq.fa
```

After Psuedo_tRNA sequence clustering the number of sequences is also 109  

```
# extract non Pseudo tRNAs
./extract_nonPseudo_fasta.py gencode_ensembl.combined_tRNA.all.fa > gencode_ensembl.combined_tRNA.non_pseudo.fa

# cluster 
usearch -cluster_fast gencode_ensembl.combined_tRNA.non_pseudo.fa -id 1.00 -centroids gencode_ensembl.combined_tRNA.non_pseudo.uniq.fa

#check the number of sequences after clustering 

rg -c '^>' gencode_ensembl.combined_tRNA.non_pseudo.uniq.fa

```

393 unique non-Pseudo tRNAs => no true tRNA clusters with Pseudo_tRNA tand "vanishes".


### creating an artificial genome

This is to combina the genomic sequence with masked (using tRNA GTF annotations) with extracted unique tRNA (and Pseudo_tRNA sequences just in case if these are also expresed).

```(sh)
# concatenate two fasta files
cat GRCh38.primary_assembly.genome.names_fix.masked_tRNAs.fa gencode_ensembl.combined_tRNA.uniq.fa >  GRCh38_plus_tRNAs.fa

# fingerprint
md5sum GRCh38_plus_tRNAs.fa 
1fdead45afd0e95d59bb0bf8e5e92f84  GRCh38_plus_tRNAs.fa

# compress for transfer to HPC
pigz -8 GRCh38_plus_tRNAs.fa
```


## mapping

### create genomic index

* unpack

```
# gunzip fasta
pigz -d GRCh38_plus_tRNAs.fa.gz

```

* slurm script

```(sh)
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

```(sh)
# command
ls -l *fastq.gz 

# output
-r--r--r-- 1 darked darked 421363708 Apr 25 16:25 D_1.fastq.gz
-r--r--r-- 1 darked darked 442198262 Apr 25 16:25 D_BH4_1.fastq.gz
-r--r--r-- 1 darked darked 318421527 Apr 25 16:24 KO6M_1_1.fastq.gz
-r--r--r-- 1 darked darked 363724791 Apr 25 16:24 KO6M_2_1.fastq.gz
-r--r--r-- 1 darked darked 378515297 Apr 25 16:23 KO6M_3_1.fastq.gz
-rw-r--r-- 1 darked darked 391489072 Apr 26 13:15 KO6M_4_1.fastq.gz
-rw-r--r-- 1 darked darked 458572848 Apr 26 13:16 NoD_1.fastq.gz
-r--r--r-- 1 darked darked 410977493 Apr 25 15:31 NoD_BH4_1.fastq.gz
-r--r--r-- 1 darked darked 391744579 Apr 25 15:30 WT6M_1_1.fastq.gz
-r--r--r-- 1 darked darked 603965167 Apr 25 16:08 WT6M_2_1.fastq.gz
-r--r--r-- 1 darked darked 434950525 Apr 25 16:07 WT6M_3_1.fastq.gz
-r--r--r-- 1 darked darked 481868553 Apr 25 15:29 WT6M_4_1.fastq.gz
```

* md5 checksums

```
# command
md5sum *fastq.gz 

# result
bb2c0c03a31c47077b4136b599374f51  D_1.fastq.gz
ec37b58216ecc8958df2b4e876133762  D_BH4_1.fastq.gz
af5c9a9537fcb9c05187a27afbf56757  KO6M_1_1.fastq.gz
cb51a4e81889fd6f81506eec8608ebdf  KO6M_2_1.fastq.gz
6b9d57d51f9e939056f5cebe56aa4c6e  KO6M_3_1.fastq.gz
1b6bb9fcaed937f0664fc8808987386c  KO6M_4_1.fastq.gz
68ecdb9735be757c7d53cb556e5258fa  NoD_1.fastq.gz
f8d3140c32e32b7df4888f9697de9923  NoD_BH4_1.fastq.gz
aba6ed4ec506bc5bc74d8f751fe5cfc5  WT6M_1_1.fastq.gz
1a7ba3c55167587bd9886e691da7e5a3  WT6M_2_1.fastq.gz
d9d8ec2311039e266210366519893d52  WT6M_3_1.fastq.gz
726ac3f7e73d09195f59ac9bf2804c1d  WT6M_4_1.fastq.gz
```

* location on HPC cluster

```
/scratch/dkedra/proj/trna_20220426/FQ_data/ORIG
```

### clustering and low complexity filter

To shrink size and speed up subsequent processing of the fastq data used clumpify.sh from bbmap ver 38.96

example SLURM shell script:

```
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

```
#dir: /scratch/dkedra/proj/trna_20220426/FQ_data

./batch_clumpify.py

```

### quality check and optical replicates filtering

Used ```fastp``` (installed using conda) for that task.

Example SLURM submission script:

```
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

```
#dir: /scratch/dkedra/proj/trna_20220426/FQ_data

./batch_fastp.py
```

