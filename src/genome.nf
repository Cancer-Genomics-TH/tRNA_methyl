nextflow.enable.dsl=2
debug=true

params.genome_s3 = "s3://human-pangenomics/T2T/CHM13/assemblies/chm13v2.0.fa"
//"https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz"

process download_fasta {
   
  publishDir './genome', mode: 'copy', overwrite: false

  output:
    path genome_fa

  script:
  genome_fa = "chm13_2.0.fa"
  
  """
  /dades/ubbiomed/soft/bin/aws s3 --no-sign-request cp ${params.genome_s3} ./$genome_fa
  
  sed -i 's/^>chr/>/g; s/^>M/>MT/g' chm13_2.0.fa*
  """
}


process index_fasta {

  publishDir './genome', mode: 'copy', overwrite: false

  executor 'slurm'
  queue 'express'
  cpus 1
  time '30m'

  conda '/scratch/dkedra/.conda/envs/samtools/'
  
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
    tuple val(my_key), path(mt_trnas_bed), path(mt_trnas_fa), path(mt_trnas_out)
 
  script:
  my_key = "t2t"


  """
  #samtools faidx $genome_fa
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

  publishDir './results', mode: 'copy', overwrite: false

  executor 'slurm'
  queue 'std'
  cpus 16 //no real improvement if using more cores
  time '2h'

  conda "/scratch/dkedra/.conda/envs/trnascan-se/"

  input:
    path genome_fa
    
  output:
    tuple val(my_key), path(chr_trnas_bed), path(chr_trnas_fa), path(chr_trnas_out)
  script:
  my_key = "t2t"

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


process download_mature_trnas {

  publishDir './genome', mode: 'copy', overwrite: false
  
   executor 'slurm'
  queue 'express'
  cpus 1 
  time '10m'
  conda "/scratch/dkedra/.conda/envs/seqtk/"


  output:
    path gtrnadb_mature_trnas_fa
  
  script:
  gtrnadb_mature_trnas_fa = "gtrnadb_mature_trnas.fa"
  //println "debug 1"
  """
  wget -O tmp_A_fa  http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-mature-tRNAs.fa
  
  seqtk seq  -C -l 80  tmp_A_fa | sed '/^>/ s/Homo_sapiens_//g' | sed '/^>/! s/U/T/g' > tmp_B_fa
  
  $baseDir/bin/trna_add_CCA.py  tmp_B_fa >  $gtrnadb_mature_trnas_fa

  """
}

process combine_trnascan_fa {
   
   input:

    tuple path(mt_trnas_bed), path(mt_trnas_fa), path(mt_trnas_out)
    tuple path(chr_trnas_bed), path(chr_trnas_fa), path(chr_trnas_out)

output:
  path trnascan_combined_fa

  script:
    //trnascan_combined_fa = "trnascan_combined.fa"

  """ 
  $baseDir/bin/trna_tscan_names_fix_sort.py $chr_trnas_fa $mt_trnas_fa > trnascan_combined_fa
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
  $baseDir/bin//sort_trnas.py combo_tmp_B > combined_uniq_trnas_fa

   """

}

process final_genome{

 publishDir './genome', mode: 'copy', overwrite: false

input:
path genome_trna_masked_fa
path combined_uniq_trnas_fa

output:
  path genome_for_mapping
script:
 genome_for_mapping =  "t2t_trnas.fa"
"""
cat $genome_trna_masked_fa $combined_uniq_trnas_fa > $genome_for_mapping
"""
}

workflow {

  genome_fa_ch = Channel.fromPath("./genome/chm13_2.0.fa") .ifEmpty('download_fasta()').view() 
  genome_fai_ch = Channel.fromPath("./genome/chm13_2.0.fa.fai") .ifEmpty('index_fasta()').view()
  rna_scan_mt_ch = Channel.fromPath("./results/mt_trnas_*") .ifEmpty('trna_scan_mt(genome_fa_ch, genome_fai_ch)').toSortedList().view()
  rna_scan_chr_ch = Channel.fromPath("./results/chr_trnas_*") .ifEmpty('trna_scan_chr(genome_fa_ch)').toSortedList().view()
  download_mature_trnas()
  
  bed_mask(genome_fa_ch, genome_fai_ch, rna_scan_mt_ch, rna_scan_chr_ch)
  //mature_trnas_ch = Channel.fromPath("./genome/gtrnadb_mature_trnas.fa") .ifEmpty('download_mature_trnas()').view()


 combine_trnascan_fa(rna_scan_mt_ch, rna_scan_chr_ch)
 combine_mature_trnascan_fa( download_mature_trnas.out,  combine_trnascan_fa.out)
 final_genome(bed_mask.out, combine_mature_trnascan_fa.out)


}

