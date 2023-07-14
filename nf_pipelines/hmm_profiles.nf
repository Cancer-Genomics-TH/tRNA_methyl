/*--------------------------------------------------
    *  Nextflow pipeline for coonstracting HMM profiles using gtRNAdb aligments / mature tRNAs file
    *  uses mlocarna from LOCARNA, and Infernal

    *  FIXME: to fix issues with the same python libraries required for custom scripts some scripts use pypy3 instead of python3  
    
---------------------------------------------------*/


nextflow.enable.dsl=2

debug=true


params.gtrnadb_mature_trnas = "http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-mature-tRNAs.fa"

//unused, hardcoded in the script:"
//params.gtrnadb_html = "http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-tRNAs.html"


process gtrnadb_html_to_stockholm {
  //
  // 1. download html file from gtrnadb
  // 2. parse html file to stockholm format

  publishDir "./hmm", mode: 'copy'

  conda "/scratch/dkedra/.conda/envs/python_310/"

  output:
    path hs38_gtrnadb_fromhtml_sto
  
  script:
    hs38_gtrnadb_fromhtml_sto = "trnas_hs38_gtrnadb_fromhtml.sto"

    """
    $baseDir/bin/parse_gtrnadb_alignments.py  $hs38_gtrnadb_fromhtml_sto

    """

}

process gtrnadb_spliced_to_fasta {
  // 1. download mature tRNAs fasta file from gtrnadb
  // 2. extract spliced tRNAs from fasta file to separate files

  publishDir "./hmm", mode: 'copy'

  conda "/scratch/dkedra/.conda/envs/python_310/"

  output:
    tuple path(spliced_fasta_dir), path(gtrnadb_mature_trnas) 

  script:
    spliced_fasta_dir = "gtrnadb_trnas_spliced"
    gtrnadb_mature_trnas = "hg38-mature-tRNAs.fa"
    
    """
    mkdir $spliced_fasta_dir

    wget -O $gtrnadb_mature_trnas ${params.gtrnadb_mature_trnas}
    
    $baseDir/bin/extract_spliced_trnas_gtrnadb.py  $spliced_fasta_dir
    
    """

}

process locarna_spliced_trnas{
  
  // 1. run mlocarna on spliced tRNAs fasta files
  // 2. merge Stokholm outputs fixing banes 
  // 3. special case for tRNA-Tyr-ATA (one sequence, no mlocarna alignment)
  // 3a. add tRNA-Tyr-ATA to the stockholm file 

  publishDir "./hmm", mode: 'copy'

  conda "/scratch/dkedra/.conda/envs/locarna/"

  input:
    tuple path(spliced_fasta_dir), path(gtrnadb_mature_trnas)

  output:
    path spliced_trna_stockholm

  script:
    """
    $baseDir/bin/spliced_locarna_stockholm.py  $spliced_fasta_dir > spliced_trna_stockholm

    """
}

process combine_sto {
  // 1. combine stockholm files from gtrnadb_html and mlocarna

  publishDir "./hmm", mode: 'copy'
  conda "/scratch/dkedra/.conda/envs/locarna/"

  input:
    path hs38_gtrnadb_fromhtml_sto
    path spliced_trna_stockholm


  output:
    path combined_trnas_stockholm

  script:
  combined_trnas_stockholm = "hg38_gtrnadb_combined.sto"
  """
  cat $hs38_gtrnadb_fromhtml_sto $spliced_trna_stockholm > $combined_trnas_stockholm
  """
}

process cmbuild {

  // 1. build Infernall hmm from stockholm file

  publishDir "./hmm", mode: 'copy'
  conda "/scratch/dkedra/.conda/envs/infernal/"

  input:
    path combined_trnas_stockholm
  
  output:
    path combined_trnas_cm
  
  script:
  combined_trnas_cm = "hg38_gtrnadb_combined.cm"
  """
  cmbuild --noss -F $combined_trnas_cm $combined_trnas_stockholm

  """

}


process cmcalibrate {
  // 1. calibrate Infernal hmm "hg38_gtrnadb_combined.cm"
  // 2. save calibrated hmm to "calibrated_hg38_gtrnadb_combined.cm"
  // 3. CPU intense step, use slurm
  // TODO: resolve issue with saving output file with the same name as input file

  publishDir "./hmm", mode: saveAs: { filename -> "calibrated_$filename" }
 
  executor 'slurm'
  queue 'std'
  cpus 24
  time '8h'


  conda "/scratch/dkedra/.conda/envs/infernal/"

  input:
    path combined_trnas_cm
  
  output:
    path combined_trnas_cm
  
  script:
  """
  cmcalibrate --cpu 24 $combined_trnas_cm
  
  """
}⏎


process cmpress {
  // 1. compress Infernal hmm "hg38_gtrnadb_combined.cm"
  // 2. save compressed hmm to "calibrated_hg38_gtrnadb_combined.cm"
  // files created:
  // calibrated_hg38_gtrnadb_combined.cm.i1f
  // calibrated_hg38_gtrnadb_combined.cm.i1i
  // calibrated_hg38_gtrnadb_combined.cm.i1m
  // calibrated_hg38_gtrnadb_combined.cm.i1p
 
  publishDir "./hmm/infernal/", mode: 'copy'

  
  conda "/scratch/dkedra/.conda/envs/infernal/"
  
  input:
    path combined_trnas_cm
  output:
    path "*.cm*"

  script:

    """
    cmpress $combined_trnas_cm
    #hg38_gtrnadb_combined.cm
    """
}


workflow {
    gtrnadb_html_to_stockholm()
    gtrnadb_spliced_to_fasta()
    locarna_spliced_trnas(gtrnadb_spliced_to_fasta.out)
    combine_sto(gtrnadb_html_to_stockholm.out, locarna_spliced_trnas.out)
    cmbuild(combine_sto.out)
    //commented out to speed up testing
    //cmcalibrate file exists already, hence the fromPath 
    //cmcalibrate(cmbuild.out)
    myFileChannel = Channel.fromPath('./data/calibrated_hg38_gtrnadb_combined.cm') //, checkIfExists=true)
    myFileChannel.view()
    cmpress(myFileChannel)
    
}
