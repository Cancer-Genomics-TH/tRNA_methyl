nextflow.enable.dsl=2

params.last_db = "/scratch/dkedra/proj/trna_20220426_rsynced_20220920/last_index/t2t_trnas/t2t_sorted_uniq_trnas.ref.last"

params.outdir = "results"


process lastal_map {

  executor 'slurm'
  queue 'std'
  cpus 16
  memory '48 GB'  
  time '120m'
  
  publishDir params.outdir,  mode: "copy", pattern: "*.maf.gz"

  //module load conda
  conda  "/scratch/dkedra/.conda/envs/last_latest"
  
  input:
    path fq_fn

  output:
      tuple val(fn_prefix), path(result_maf_gz)
      
  script:

    
    fn_prefix = fq_fn.getSimpleName()
    result_maf_gz = fn_prefix + ".last.t2t_trnas.maf.gz"
   
    """
    lastal -v -P16  -Qkeep  -C2  ${params.last_db}  $fq_fn  | last-postmask | last-split | pigz  > $result_maf_gz
    """
}


workflow {
    data = channel.fromPath("./data/*.fq.gz") 
    data.view()
    lastal_map(data)
    
}