# bioinformatics workflow

## overview

To filter out not fully processed tRNA sequences we map FASTQ reads twice: 
1. to a constructed genome with unique tRNA added as separate contigs
2. to a regular human genome

We extract read sequences mapping to tRNA contigs from *1* and filter out reads overlapping the tRNAs genomic flanks as determined by mapping *2*. 

## genome and tRNAs for mapping

To filter out not fully processed tRNA sequences we map FASTQ reads twice: 
* to a constructed genome with unique tRNA added as separate contigs
* to a regular human genome



### constructed genome

We use T2T v2.0 genomic sequence to predict tRNAs using ```tRNAscan-SE```. Predicted tRNAs are masked in the genomic sequence. 
We extract both spliced and non-spliced tRNA sequences predicted by ```tRNAscan-SE``` and combine these with mature tRNAs from gtRNAdb.
Identical (100% identity) sequences are removed using vsearch.
To increase size of mappable tRNA sequences we add CCA to the 3 prime ends of all tRNAs in the set.  
Finally we combine  T2T contigs with masked tRNAs with the above super set of mature and predicted unique tRNA sequences.

```mermaid
graph TD

    genome_fasta{ T2T v2.0 genome} -->|fix names| fasta_names_fix[fasta with chr names: 1-22,X,Y and MT]
    fasta_names_fix --> |extract chromosomes 1-22,X,Y| genome_chromosomes[genome fasta, no MT]
    fasta_names_fix --> |extract MT| genome_mt[mitochondrial sequence MT only]
    genome_chromosomes --> |tRNAscan SE predictions genomic | trna_predicted_genomic_fa[predicted tRNAs genomic fasta]
    genome_chromosomes --> |tRNAscan SE predictions genomic | trna_predicted_genomic_bed[predicted tRNAs genomic BED]
    genome_mt --> |tRNAscan-SE predictions MT | trna_predicted_mt_fa[predicted tRNAs MT fasta]
    genome_mt --> |tRNAscan-SE predictions MT | trna_predicted_mt_bed[predicted tRNAs MT BED]
    trna_predicted_genomic_bed --> |combine predictions BEDs| predicted_trnas_all[predicted tRNAs genomic and MT BED]
    trna_predicted_mt_bed --> |combine predictions BEDs| predicted_trnas_all[predicted tRNAs genomic and MT BED]
    fasta_names_fix --> |bedtools masking|masked_genome[fasta with masked tRNA/Pseudo tRNA genes]
    trna_predicted_genomic_fa -->|merge fasta| trna_predicted_all[T2T predicted tRNAs all]
    trna_predicted_mt_fa --> |merge fasta| trna_predicted_all[T2T predicted tRNAs all]
    predicted_trnas_all --> |bedtools masking|masked_genome[fasta with masked tRNA/Pseudo tRNA genes]
    trna_predicted_all --> |fix tRNA names and sort| trna_predicted_all_name_fixed[T2T predicted tRNAs all with gtRNAdb-like naming scheme ordered by name]
    gtRNAdb --> |download|mature_trnas_gtrnadb_fa[gtRNAdb hg38 mature tRNAs]
    mature_trnas_gtrnadb_fa --> |fix tRNA names| mature_trnas_gtrnadb_name_fixed_fa[gtRNAdb hg38 mature tRNAs no Homo_sapiens prefix]  
    mature_trnas_gtrnadb_name_fixed_fa --> |combine fasta | trna_superset[superset of tRNA sequences]
    trna_predicted_all_name_fixed -->  |combine fasta | trna_superset[superset of tRNA sequences] 
    trna_superset -->  |vsearch clustering|unique_trnas_seq[fasta with unique tRNA/pseudo-tRNA sequences]
    unique_trnas_seq --> | sort by name and add CCA to 3prime | unique_trnas_cca_seq[fasta with unique tRNA/pseudo-tRNA sequences with added 3prime CAA]
    masked_genome --> |merge fasta | masked_genome_trnas[masked T2T genome with separate unique tRNAs]
    unique_trnas_cca_seq --> |merge fasta | masked_genome_trnas[masked T2T genome with separate unique tRNAs]
    masked_genome_trnas --> | index fasta | masked_genome_trnas_fai[T2T genome plus tRNAs  index]
```

### human genome

```mermaid
graph TD

    genome_fasta{ T2T v2.0 genome} -->|fix names| t2t_genome[ T2T fasta with chr names: 1-22,X,Y and MT]
    t2t_genome --> | index fasta | T2T_genome_fai[T2T genome index]

```


## FASTQ mapping

We use the constructed as described above artificial genome to map short RNA NGS sequences. Reads in the input FASTQ files are clustered by sequence, then pre-processed using ```fastp```
removing adaptor sequences. 
Mapping is done using a strict mapper (```LAST```). Only the top hits are retained.

### constructed genome

```mermaid
graph TD
    genome_for_mapping[masked T2T genome with separate unique tRNAs]-->|lastdb|last_database[LAST database]
    input_fastq[FASTQ file] --> |reads clustering| clustered_fastq[clustered FASTQ]
    clustered_fastq --> |QC, adaptor filtering| fastped_fastq[clustered FASTQ, adaptors removed ]
    fastped_fastq -->|lastal mapping|maf_result[mapping result MAF file]
    last_database -->|lastal mapping|maf_result
    maf_result --> |convert to BAM, sort and index| bam_result[mapped results BAM]
    bam_result --> |visual QC in IGV| igv_screen_shots[IGV view mappings]
```

### human genome

```mermaid
graph TD
   t2t_genome[ T2T fasta with chr names: 1-22,X,Y and MT]-->|lastdb|last_database[T2T LAST database]
    input_fastq[FASTQ file] --> |reads clustering| clustered_fastq[clustered FASTQ]
    clustered_fastq --> |QC, adaptor filtering| fastped_fastq[clustered FASTQ, adaptors removed ]
    fastped_fastq -->|lastal mapping|t2t_maf_result[T2T mapping result MAF file]
    last_database -->|lastal mapping|t2t_maf_result
    t2t_maf_result --> |convert to BAM, sort and index| t2t_bam_result[T2T mapped results BAM]
    t2t_bam_result --> |visual QC in IGV| igv_screen_shots[IGV view mappings]
```

# preparing tRNA isoforms aligments for cmscan

The primary source of the alignments is a HTML page from gtRNAdb (hg38 tRNAs). Aligments for the spliced tRNAs (tRNA-Arg-TCT
tRNA-Ile-TAT, tRNA-Leu-CAA, tRNA-Tyr-ATA, tRNA-Tyr-GTA) were created using mlocarna from LocaRNA.

```mermaid
graph TD

html_aligment[gtRNAdb html page] --> |parse_gtrnadb_alignments.ipynb|stockholm_align_A[text file with human tRNAs aligments in Stockholm format]
spliced_trnas[spliced tRNAs fasta files] --> |mlocarna alignment| spliced_sto_files[text files with human spliced tRNAs ]
stockholm_align_A --> |concatenate | stockholm_align_all[superset of spliced and unspliced alignments]
spliced_sto_files --> |concatenate | stockholm_align_all[superset of spliced and unspliced alignments]
stockholm_align_all -->|cmbuild from Infernal |cm_file[cm format alingment]
cm_file --> |cmcalibrate from Infernal |cm_calibrated[files for cmscan searches]
cm_calibrated --> |cmpress |cm_files[files for cmscan searches]
```

## parsing NGS mappings

Objective: extract mappings to mature non-mitochondrial, non-psedogenes tRNAs.

Since short fragment reads can map to often nearly tRNA sequences we filter out reads mapping to pseudogenes/MT tRNAs to obtain reliable read assigments.  
Reads with mapping quality (MAPQ) zero are also dropped.

FIXME
use: bam_2_trna_reads_counts_fasta.py

```mermaid
graph TD
    unique_trnas_seq-->|ripgrep,sed|unique_trnas_seq_names[text file with unique by sequence tRNAs names]
    unique_trnas_seq_names-->|ripgrep, maf_2_fa_with_counts.py|trna_matches_with_counts[fasta with fastq-derived tRNA seq matches with counts] 
    indexed_bam-->|ripgrep, maf_2_fa_with_counts.py|trna_matches_with_counts
    
```

## getting mapping positions using Infernal

FIXME 
```mermaid
graph TD
    cmscan_files-->|cmscan|cmscan_mapping_result[tabular mapping results]
    trna_matches_with_counts-->|cmscan|cmscan_mapping_result
    cmscan_mapping_result-->|top hits selection|cmscan_top_hits[TSV file with top hits and their positions]
    cmscan_top_hits-->|normalize coverage and plot|heat_maps[heat maps]
```