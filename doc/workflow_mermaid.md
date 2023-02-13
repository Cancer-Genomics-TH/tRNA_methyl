# bioinformatics workflow


## genome and tRNAs for mapping

We use T2T v2.0 genomic sequence to predict tRNAs using tRNAscan-SE. Predicted tRNAs are masked in the genomic sequence. 
Since the predictions do not provide us with information about tRNA splicing, we combine our predicted tRNA sequences with mature tRNAs from gtRNAdb.
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
## FASTQ mapping

We use the constructed as described above artificial genome to map short RNA NGS sequences. Reads in the input FASTQ files are clustered by sequence, then pre-processed using ```fastp```
removing adaptor sequences. 
Mapping is done using a strict mapper (```LAST```). Only the top hits are retained.

```mermaid
graph TD
    genome_for_mapping{masked T2T genome with separate unique tRNAs}-->|lastdb|last_database[LAST database]
    input_fastq[FASTQ file] --> |reads clustering| clustered_fastq[clustered FASTQ]
    clustered_fastq --> |QC, adaptor filtereing| fastped_fastq[clustered FASTQ, adaptors removed ]
    fastped_fastq -->|lastal mapping|maf_result[mapping result MAF file]
    last_database -->|lastall mapping|maf_result
```


## parsing mappings

```
    unique_trnas_seq-->|ripgrep,sed|unique_trnas_seq_names[text file with unique by sequence tRNAs names]
    unique_trnas_seq_names-->|ripgrep, maf_2_fa_with_counts.py|trna_matches_with_counts[fasta with fastq-derived tRNA seq matches with counts] 
    maf_result-->|ripgrep, maf_2_fa_with_counts.py|trna_matches_with_counts
    gtRNAdb{gtRNAdb aligment}-->|download|html_aligment[html page with human tRNA aligments]
    html_aligment-->|ripgrep, tab_gtrna_to_stockholm.py|stockholm_align[text file with human tRNAs aligments in Stockholm format]
    stockholm_align-->|Infernal tools|cmscan_files[files for cmscan searches]
    cmscan_files-->|cmscan|cmscan_mapping_result[tabular mapping results]
    trna_matches_with_counts-->|cmscan|cmscan_mapping_result
    cmscan_mapping_result-->|top hits selection|cmscan_top_hits[TSV file with top hits and their positions]
    cmscan_top_hits-->|normalize coverage and plot|heat_maps[heat maps]
    fastq{fastq input} -->|clumpify| clumped_fastq(fastq minus optical replicates, clustered reads)
    clumped_fastq -->|fastp| clumped_fastped_fastq[fastq minus low complexity reads]
    clumped_fastped_fastq -->|fqgrep, batch_fqgrep.py script|clumped_fastped_masked_fastq[fastq with masked primer sequences]
    
```