/*--------------------------------------------------
    *  Nextflow pipeline for coonstracting HMM profiles using gtRNAdb aligments / mature tRNAs file

    *  uses LOCARNA, Infernal
    *  
---------------------------------------------------*/


nextflow.enable.dsl=2

debug=true


params.gtrnadb_mature_trnas = "http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-mature-tRNAs.fa"


process gtrnadb_html_to_stockholm {

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
  
  //pattern =  spliced_fasta_dir + "/tRNA*.fa"
    //spliced_trna_stockholm = "foo.bar.sto"
    //println pattern
    
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




process cmpress {
   
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
    //cmcalibrate(cmbuild.out)
    myFileChannel = Channel.fromPath('./data/calibrated_hg38_gtrnadb_combined.cm') //, checkIfExists=true)
    myFileChannel.view()
    cmpress(myFileChannel)
        //cmcalibrate.out)
}
