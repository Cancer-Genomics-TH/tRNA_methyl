/*--------------------------------------------------
    *  t2t_trnas pipeline step_1
    *  Nextflow pipeline for constructing artificial genome using human T2T assembly, trnascan-SE predictions and gtRNAdb mature tRNAs (hg38)
    *  uses awscli; samtools, trnascan-SE  bedtools  LAST aligner 
    *  other utils: pypy, ripgrep (rg), hck 
---------------------------------------------------*/


nextflow.enable.dsl=2
debug=true

params.genome_s3 = "s3://human-pangenomics/T2T/CHM13/assemblies/chm13v2.0.fa"
params.gtrnadb_mature_trnas = "http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-mature-tRNAs.fa"
params.gtrnadb_confident_trnas = "http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-tRNAs.fa"

process download_fasta {

  publishDir './genome', mode: 'copy', overwrite: false

  conda "/scratch/dkedra/.conda/envs/awscli/"

  output:
    path genome_fa

  script:
  genome_fa = "chm13_2.0.fa"
  
  """
  aws s3 --no-sign-request cp ${params.genome_s3} ./$genome_fa

  sed -i 's/^>chr/>/g; s/^>M/>MT/g' chm13_2.0.fa
  """
}


process index_fasta {

  publishDir './genome', mode: 'copy', overwrite: false

  executor 'slurm'
  queue 'express'
  cpus 1
  time '30m'

  conda "/scratch/dkedra/.conda/envs/samtools/"
  
  input:
    path genome_fa
  output:
    path(genome_fai)

  script:
  genome_fai = "chm13_2.0.fa.fai"

  """
  samtools faidx --fai-idx $genome_fai $genome_fa
  """
}

process trna_scan_mt {

  publishDir './results', mode: 'copy', overwrite: false

  executor 'slurm'
  queue 'express'
  cpus 1
  time '30m'

  conda "/scratch/dkedra/.conda/envs/trnascan-se/"

  input:
    path(genome_fa)
    path(genome_fai)
    
  output:
    tuple path(mt_trnas_bed), path(mt_trnas_fa), path(mt_trnas_out)
 
  script:
  my_key = "t2t"


  """
  samtools faidx $genome_fa MT > MT.fa
  
  tRNAscan-SE -E \
    --thread 1 \
    --detail \
    -M mammal \
    --output mt_trnas_out \
    --bed mt_trnas_bed \
    --fasta mt_trnas_fa \
    --log mt.trnas_log MT.fa
  """
}



process trna_scan_chr {
  // uses also MT to avoid extracting all chromosomes from the genome
  // MT predictions are removed later

  publishDir './results', mode: 'copy', overwrite: false

  executor 'slurm'
  queue 'std'
  cpus 16 //no real improvement if using more cores
  time '2h'

  conda "/scratch/dkedra/.conda/envs/trnascan-se/"

  input:
    path genome_fa
    
  output:
    tuple path(chr_trnas_bed), path(chr_trnas_fa), path(chr_trnas_out)
  
  script:

  """
    tRNAscan-SE -E \
    --thread 16 \
    --detail \
    --output chr_trnas_out \
    --bed chr_trnas_bed \
    --fasta chr_trnas_fa \
    --log chr_trnas_log $genome_fa

  """
}

process bed_mask {

  publishDir './results', mode: 'copy', overwrite: false, pattern: 't2t.trna_masked.fa'

  executor 'slurm'
  queue 'std'
  cpus 1
  time '1h'

  conda "/scratch/dkedra/.conda/envs/bedtools/"

  input:
    path genome_fa
    path genome_fai
    tuple path(mt_trnas_bed), path(mt_trnas_fa), path(mt_trnas_out)
    tuple path(chr_trnas_bed), path(chr_trnas_fa), path(chr_trnas_out)
  
  output: 
    path genome_trna_masked_fa
  

  script:
    genome_trna_masked_fa = "t2t.trna_masked.fa"
  """
  rg -v '^MT' $chr_trnas_bed  | hck -f=1,2,3 >  tmp_chr_bed
  
  hck -f=1,2,3 $mt_trnas_bed > tmp_mt_bed

  cat  tmp_chr_bed tmp_mt_bed | bedtools sort -faidx $genome_fai -i stdin  > masking_bed

  bedtools maskfasta  -fi  $genome_fa -bed masking_bed -fo $genome_trna_masked_fa

  """
}


process download_gtrnadb_fastas {
  // download mature tRNAs from UCSC and convert to fasta format
  // add CCA tail
  // script uses pypy to speed up the conversion
  // requires pyfasta library 
  
  publishDir './genome', mode: 'copy', overwrite: false

  executor 'slurm'
  queue 'express'
  cpus 1
  time '10m'
  conda "/scratch/dkedra/.conda/envs/pypy_39/"


  output:
    path gtrnadb_mature_trnas_fa
  
  script:
    gtrnadb_mature_trnas_fa = "gtrnadb_mature_trnas.fa"
  
    """
    wget ${params.gtrnadb_mature_trnas}
    wget ${params.gtrnadb_confident_trnas}
    
    sed -i '/^>/! s/U/T/g' hg38-mature-tRNAs.fa
    
    $baseDir/bin/reformat_combine_gtRNAdb_fa.py > tmp_A_fa
    
    $baseDir/bin/sort_trnas.py  tmp_A_fa >  $gtrnadb_mature_trnas_fa

  """
}


process combine_trnascan_fa {

  executor 'slurm'
  queue 'express'
  cpus 1 
  time '10m'

  conda "/scratch/dkedra/.conda/envs/python_310/"
     
  input:

    tuple path(mt_trnas_bed), path(mt_trnas_fa), path(mt_trnas_out)
    tuple path(chr_trnas_bed), path(chr_trnas_fa), path(chr_trnas_out)

  output:
    path trnascan_combined_fa

  script:
    
  """ 
  $baseDir/bin/trna_tscan_names_fix_sort.py $chr_trnas_fa $mt_trnas_fa  $chr_trnas_out > trnascan_combined_fa
  """
}


process combine_mature_trnascan_fa {

  executor 'slurm'
  queue 'express'
  cpus 1 
  time '10m'
 

  conda "/scratch/dkedra/.conda/envs/vsearch/"

  input:
    path gtrnadb_mature_trnas_fa
    path trnascan_combined_fa
   
  output:
    path combined_uniq_trnas_fa

  script:
   //combined_uniq_trnas_fa = "combined_uniq_trnas.fa"
    """
    cat $gtrnadb_mature_trnas_fa $trnascan_combined_fa > combo_tmp_A
    vsearch -fastx_uniques combo_tmp_A  -fastaout  combo_tmp_B
    $baseDir/bin/sort_trnas.py combo_tmp_B > combined_uniq_trnas_fa

   """

}

process final_genome{

  publishDir './genome', mode: 'copy', overwrite: false, pattern: 't2t_trnas.fa*'
  
  executor 'slurm'
  queue 'express'
  cpus 1
  time '30m'

  conda "/scratch/dkedra/.conda/envs/samtools/"


  input:
    path genome_trna_masked_fa
    path combined_uniq_trnas_fa

  output:
    tuple path(genome_for_mapping), path(genome_for_mapping_fai)

  script:
    genome_for_mapping =  "t2t_trnas.fa"
    genome_for_mapping_fai = "t2t_trnas.fa.fai"

    """
    cat $genome_trna_masked_fa $combined_uniq_trnas_fa > $genome_for_mapping
    samtools faidx --fai-idx $genome_for_mapping_fai $genome_for_mapping
    """
}

process lastdb {

  publishDir './genome/lastdb', mode: 'copy', overwrite: false

  executor 'slurm'
  queue 'mem'
  cpus 32
  time '3h'
  
  conda "/scratch/dkedra/.conda/envs/last_aligner/"

  input:
    path genome_for_mapping

  output:
    tuple val(lastdb_prefix), path("lastdb")

  script:
    lastdb_prefix = "t2t_trnas.last_soft"
    
    """
    lastdb   -P32 \
    -uNEAR -v -c \
    ./lastdb/$lastdb_prefix $genome_for_mapping

    """

}

workflow {
  download_fasta()
  
  index_fasta(download_fasta.out)

  download_mature_trnas()
  trna_scan_mt(download_fasta.out, index_fasta.out)
  trna_scan_chr(download_fasta.out)
  combine_trnascan_fa(trna_scan_mt.out, trna_scan_chr.out)

  bed_mask(download_fasta.out, index_fasta.out, trna_scan_mt.out, trna_scan_chr.out)

  combine_mature_trnascan_fa( download_mature_trnas.out,  combine_trnascan_fa.out)

  final_genome(bed_mask.out, combine_mature_trnascan_fa.out)
  lastdb(final_genome.out)

}


