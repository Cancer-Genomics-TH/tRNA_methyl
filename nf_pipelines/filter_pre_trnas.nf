nextflow.enable.dsl=2
debug=true

params.t2t_trnas_superset_noMT_bed = "$baseDir/genome/trnas_superset.noMT.bed"
params.t2t_trnas_flanks_bed = "$baseDir/genome/trnascan_t2t_chr_trnas.noMT.flanks_2up_6down.bed"
params.selected_trnas_bed = "$baseDir/genome/superset_trnas_selected.bed"
params.trna_hmm_profiles = "$baseDir/hmm/infernal/calibrated_hg38_gtrnadb_combined.cm"


process trna_only_bam {

  executor 'slurm'
  queue 'express'
  cpus 4 
  time '60m'

  //publishDir '$baseDir/results_filtering',  mode: "copy", pattern: "*.list"

  conda  "/scratch/dkedra/.conda/envs/samtools"

  input:
    tuple val(sample_id), path(bam_bai_pairs_list)
  
  output:
      tuple val(sample_prefix), path(trna_only_bam), path(trna_only_bam_bai) 


  script:
    genome_trnas_bam = bam_bai_pairs_list[0]
    sample_prefix = genome_trnas_bam.getSimpleName()
    trna_only_bam = sample_prefix + ".t2t_trnas_superset.train.trnas_only.bam"
    trna_only_bam_bai = sample_prefix + ".t2t_trnas_superset.train.trnas_only.bam.bai"

    """
    samtools view -@ 4 -b -L ${params.t2t_trnas_superset_noMT_bed} $genome_trnas_bam > $trna_only_bam
    samtools index -o $trna_only_bam_bai $trna_only_bam
    """
}


process reads_in_flanks {

  executor 'slurm'
  queue 'express'
  cpus 4 
  time '60m'

  //publishDir '$baseDir/results_filtering',  mode: "copy", pattern: "*.list"

  conda  "/scratch/dkedra/.conda/envs/samtools"


input:
    tuple val(sample_id), path(bam_bai_pairs_list)

output:
      tuple val(sample_prefix), path(reads_in_flanks)


script:
   genomic_bam = bam_bai_pairs_list[0]
   sample_prefix = genomic_bam.getSimpleName()
   reads_in_flanks = sample_prefix + ".reads_in_flanks.list"

""" 
 samtools view -@ 4 -L ${params.t2t_trnas_flanks_bed} $genomic_bam \
  | hck -f=1 | sort > $reads_in_flanks

"""
}


process filter_pre_trna {

executor 'slurm'
  queue 'express'
  cpus 4 
  time '60m'

  publishDir './results/',  mode: "copy", pattern: "*.bam*"
  conda  "/scratch/dkedra/.conda/envs/samtools"
  
  input:
    tuple val(sample_prefix), path(reads_in_flanks), path(trna_only_bam), path(trna_only_bam_bai)
    
  
  output:
    tuple val(sample_prefix), path(filtered_trna_bam), path(filtered_trna_bam_bai)
  
  script:
    filtered_trna_bam = sample_prefix + ".t2t_trnas_superset.train.trnas_only.no_flank_reads.bam"
    filtered_trna_bam_bai = sample_prefix + ".t2t_trnas_superset.train.trnas_only.no_flank_reads.bam.bai"
  
  """

  samtools view -@ 4 --with-header $trna_only_bam \
  | rg -v --fixed-strings -f $reads_in_flanks \
  | samtools view -b - \
  | samtools sort - >  $filtered_trna_bam
  samtools index -@ 4 -o $filtered_trna_bam_bai $filtered_trna_bam   
  """
}


process extract_trna_reads_from_bam {

 
// use: params.t2t_trnas_superset_noMT_bed  params.selected_trnas_bed

publishDir './results/',  mode: "copy", pattern: "*.fa"

executor 'slurm'
  queue 'express'
  cpus 4 
  time '60m'

  conda  "/scratch/dkedra/.conda/envs/samtools"

  input:
    tuple val(sample_prefix), path(filtered_trna_bam), path(filtered_trna_bam_bai)
    
  output:
    tuple val(sample_prefix), path(matches_fasta)

  script:
    matches_fasta = sample_prefix + ".last.t2t_trnas.no_flank_reads.matches.fa"
   """
   $baseDir/bin/bam_2_trna_reads_counts_fasta_nf.py ${params.selected_trnas_bed} $filtered_trna_bam $matches_fasta
   """

}

process cmscan_trna_profiles {
  
  publishDir  './results/',  mode: "copy", pattern: "*.cmscan"
  executor 'slurm'
  queue 'std'
  cpus 12
  time '60m'
  
  conda  "/scratch/dkedra/.conda/envs/infernal"

  input:
    tuple val(sample_prefix), path(matches_fasta)
  
  output:
    tuple val(sample_prefix), path(cmscan_output_tab), path(cmscan_output)
  
  script:

    cmscan_output = sample_prefix + ".last.t2t_trnas.no_flank_reads.matches.cmscan"
    cmscan_output_tab = sample_prefix + ".last.t2t_trnas.no_flank_reads.matches.cmscan.tab"
  
    """
    cmscan --cpu 12 --toponly --tblout $cmscan_output_tab --tblout $cmscan_output_tab $baseDir/genome/tRNAs.cm $matches_fasta
    """

}

process cmscan_top_hits {
  
  publishDir  './results/',  mode: "copy", pattern: "*.cmscan"
  executor 'slurm'
  queue 'std'
  cpus 12
  time '60m'
  
  conda  "/scratch/dkedra/.conda/envs/pypy_39"

  input:
    tuple val(sample_prefix), path(cmscan_output_tab),  path(cmscan_output)
  output:
    tuple val(sample_prefix), path(cmscan_top_hits)
  script:

    cmscan_top_hits = sample_prefix + ".last.t2t_trnas.no_flank_reads.matches.cmscan.top_hits.tsv"
  
    """
    $baseDir/bin/parse_top_hits_cmscan.py $cmscan_output_tab > $cmscan_top_hits 
    
    """

}


workflow {

  genome_bam_ch =  Channel.fromFilePairs("$baseDir/genome_bam_data/*.{bam,bai}", checkIfExists:true) { file -> file.name.replaceAll(/.bam|.bai$/,'') }
  genome_bam_ch.view()
  
  trna_bam_ch = Channel.fromFilePairs("$baseDir/trna_bam_data/*.{bam,bai}", checkIfExists:true)  { file -> file.name.replaceAll(/.bam|.bai$/,'') }
  trna_bam_ch.view()
  trna_only_bam(trna_bam_ch)
  
  reads_in_flanks(genome_bam_ch)

  trna_only_ch = trna_only_bam.out
  reads_in_flanks_ch = reads_in_flanks.out

  combined_ch = reads_in_flanks_ch.combine(trna_only_ch, by: 0)
  combined_ch.view()
  filter_pre_trna(combined_ch)
  extract_trna_reads_from_bam(filter_pre_trna.out)
}

