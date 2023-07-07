
/*--------------------------------------------------
    *  t2t_trnas pipeline step_2
    *  Nextflow pipeline for mapping tRNAs to t2t genome with masked tRNAs and clustered tRNAS as contigs
    *  uses fastp, clumpify, lastal, samtools
---------------------------------------------------*/


nextflow.enable.dsl=2
debug=true

params.fq_primers = "$baseDir/meta/primers_Nextflex9.fa"
lastdb = "$baseDir/genome/lastdb/t2t_trnas.last_soft"
    
process fastp_filter {
   
  // basic QC and adapter trimming
  
  executor 'slurm'
  queue 'std'
  cpus 8
  time '60m'

  conda "/scratch/dkedra/.conda/envs/fastp/"
   

  input:
    path fq_fn
  output:
    tuple path(fastp_fq), path(fastp_fail_fq), path(fastp_json), path(fastp_html)

  script:
    
    fn_prefix = fq_fn.getSimpleName()
    fastp_fq      = fn_prefix + ".clump.fq.gz"
    fastp_fail_fq = fn_prefix + ".clump.fail.fq.gz"
    fastp_json    = fn_prefix + ".fastp.json"
    fastp_html    = fn_prefix + ".fastp.html"

  """
  seqtk seq -C -q 0 $fq_fn  | fastp \
   --thread 8 \
   --stdin \
   --overrepresentation_analysis \
   --adapter_fasta ${params.fq_primers} \
   --out1 $fastp_fq \
   --report_title $fq_fn \
   --json $fastp_json \
   --html $fastp_html \
   --failed_out $fastp_fail_fq 

  """
}

process clump {

  executor 'slurm'
  queue 'std'
  cpus 16
  time '2h'

  conda "/scratch/dkedra/.conda/envs/bbmap/"

  publishDir './results_fq',  mode: "copy", pattern: "*.fastp.clump.fq.gz"

    input:
      path fq_fn

    output:
      path clumped_fq

    script:
      clumped_fq = fq_fn.getSimpleName() + ".fastp.clump.fq.gz"
      num_pigz_threads = 6
      
      """
      clumpify.sh  \
      reorder=a \
      dedupe=t \
      optical=t \
      changequality=t \
      blocksize=2048 \
      ziplevel=8 pigz=$num_pigz_threads unpigz=$num_pigz_threads \
      in=$fq_fn \
      out=$clumped_fq
      """
}

process lastal_map {

  executor 'slurm'
  queue 'std'
  cpus 16 // depending on the IO on the cluster(?) one can use 16 threads instead
  memory '48 GB'
  time '120m'

  publishDir './results_fq/',  mode: "copy", pattern: "*.maf.gz"

  //module load conda
  conda  "/scratch/dkedra/.conda/envs/last_aligner"

  input:
    path fq_fn

  output:
      tuple val(fn_prefix), path(result_maf_gz)

  script:
    fn_prefix = fq_fn.getSimpleName()
    result_maf_gz = fn_prefix + ".t2t_trnas.maf.gz"

    """
    lastal -v -P12  -Qkeep  -C2  ${lastdb}  $fq_fn \
    | last-split \
    | pigz -4 --processes 4 > $result_maf_gz
    """
}

process maf_to_bam {

executor 'slurm'
  queue 'std'
  cpus 16 // depending on the IO on the cluster(?) one can use 16 threads instead
  memory '48 GB'
  time '120m'

  publishDir './results_fq/',  mode: "copy", pattern: "*.bam*"
  conda  "/scratch/dkedra/.conda/envs/last_aligner"
  
  input:
    tuple val(fn_prefix), path(maf_gz)
  
  output:
    tuple val(fn_prefix), path(result_bam)
  
  script:
    maf_conv_py = "/scratch/dkedra/.conda/envs/last_aligner/bin/maf-convert"
    result_bam = maf_gz.getSimpleName() + ".t2t_trnas.bam"
    """
    $baseDir/bin/pypy3  $maf_conv_py sam $maf_gz  \
    | samtools view -b --fai-reference $baseDir/genome/t2t_trnas.fa.fai - \
    | samtools sort -@ 12  - > $result_bam
    """

}

process bam_index {

  executor 'slurm'
  queue 'express'
  cpus 4 // depending on the IO on the cluster(?) one can use 16 threads instead
  time '15m'

  publishDir './results_fq/',  mode: "copy", pattern: "*.bam.bai"
  conda  "/scratch/dkedra/.conda/envs/samtools"
  
  input:
    tuple val(fn_prefix), path(bam_fn)
  
  output:
    tuple val(fn_prefix), path(bam_fn), path(result_bai)
  
  script:
    result_bai = bam_fn.getSimpleName() + ".t2t_trnas.bam.bai"
    """
    samtools index -@ 4 -o $result_bai $bam_fn
    """
}

workflow {

  fq_ch = Channel.fromPath("./data/*.fastq.gz")
  fq_ch.view()
  fastp_filter(fq_ch)
  clump(fastp_filter.out.map { it[0] })
  lastal_map(clump.out)
  maf_to_bam(lastal_map.out)
  bam_index(maf_to_bam.out)
  
}

