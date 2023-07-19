## Nextflow pipelines

## prerequisites
For the time being the pipeline requires installation of several conda environments. 

### conda environments

Pipelines use following conda envs

* awscli
* bbmap
* bedtools
* fastp
* infernal
* last_aligner
* locarna
* pypy_39
* python_310
* samtools
* trnascan-se
* vsearch




Apart from that it uses few command line utils for speed/ease of use:
* pypy3
* ripgrep (rg), grep substitute
* hck (cut/awk substitute)

And Python libraries:

* BeautifulSoup
* executor
* polars
* pyfaidx



## order of the commands

```
nextflow run genome.nf   -with-conda true

nextflow run hmm_profiles.nf -with-conda true

nextflow run map_fastq.nf  -with-conda true

nextflow filter_pre_trnas.nf -with-conda true
```

Tested with Nextflow version 23.04.2.5871