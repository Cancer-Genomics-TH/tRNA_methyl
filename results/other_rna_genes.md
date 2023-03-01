# Other RNA genes

Verifing that the small RNA NGS sequencing (untreated,  Bo-Seq and external TRAC-Seq)  covers other non-tRNA, non-lncRNA genes.

## method

### annotations

*bigBED* with all annotations
```
wget  https://hgdownload.soe.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/bbi/GCA_009914755.4_T2T-CHM13v2.0.catLiftOffGenesV1/catLiftOffGenesV1.bb
./bigBedToBed catLiftOffGenesV1.bb | sed 's/^chr//g' > catLiftOffGenesV1.no_chr.bed
```
*BED with RNA genes subset*

* RNA classes to exclude (in ```included_rnas.txt```):

```
miRNA
misc_RNA
rRNA
sRNA
scRNA
scaRNA
snRNA
snoRNA
vault_RNA 
```
 * excluded RNAs:

```
lncRNA
Mt_tRNA
Mt_rRNA 
```

```
rg -w -f included_rnas.txt catLiftOffGenesV1.no_chr.bed > non_trna_non_lncrna_genes.bed
```

* convert BED to SAF

FIXME


read counts obtained using:
* SAF file from the previous step
* featureCounts from subread


## read counts 

* at least 1 mapped read per gene.

Table represents number of genes for each RNA gene class with at least one read mapped.  

| rna_class | count _annotation | count_DU145_NT | count_DU145_BoSeq | count_MHCC97H_NT_Chen | count_MHCC97H_TRACSeq_Chen |
|---|---|---|---|---|---|
| miRNA | 2046 | 272 | 292 | 216 | 176 |
| misc_RNA | 2231 | 632 | 685 | 523 | 530 |
| rRNA | 1007 | 71 | 71 | 71 | 72 |
| sRNA | 5 | 1 | 2 | 1 | 1 |
| scRNA | 2 | 0 | 0 | 2 | 2 |
| scaRNA | 48 | 19 | 21 | 21 | 22 |
| snRNA | 1902 | 468 | 434 | 522 | 505 |
| snoRNA | 948 | 364 | 355 | 376 | 373 |
| vault_RNA | 1 | 1 | 1 | 1 | 1 |

* at least 10 reads mapper per gene


| rna_class | annotation_counts | count_DU145_NT | count_DU145_BoSeq | count_MHCC97H_NT_Chen | count_MHCC97H_TRACSeq_Chen |
|---|---|---|---|---|---|
| miRNA | 2046 | 140 | 90 | 43 | 43 |
| misc_RNA | 2231 | 242 | 227 | 137 | 144 |
| rRNA | 1007 | 68 | 52 | 69 | 54 |
| sRNA | 5 | 1 | 1 | 1 | 1 |
| scRNA | 2 | 0 | 0 | 2 | 2 |
| scaRNA | 48 | 17 | 16 | 18 | 17 |
| snRNA | 1902 | 247 | 182 | 202 | 201 |
| snoRNA | 948 | 298 | 277 | 293 | 291 |
| vault_RNA | 1 | 0 | 0 | 1 | 1 |
