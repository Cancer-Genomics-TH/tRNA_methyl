# Software


- [aligners/mappers](#aligners-mappers)
  * [LAST mapper](#last-mapper)
  * [Infernal](#infernal)
- [general tools](#general-tools)
  * [bedtools](#bedtools)
  * [clumpify from BBMap](#clumpify-from-bbmap)
  * [fastp](#fastp)
  * [```usearch```](#---usearch---)
  * [```fqgrep```](#---fqgrep---)
- [system utils](#system-utils)
  * [```ripgrep```](#---ripgrep---)
  * [```choose```](#---choose---)
  * [pigz](#pigz)
  * [fish shell](#fish-shell)

<small><i><a href='http://ecotrust-canada.github.io/markdown-toc/'>Table of contents generated with markdown-toc</a></i></small>



## aligners/mappers

### LAST mapper

Used to map fastq files to the genome

* www:  https://gitlab.com/mcfrith/last
* docs: https://gitlab.com/mcfrith/last/-/tree/main/doc
* version: 1281 
* obtained from: https://gitlab.com/mcfrith/last/-/archive/1281/last-1281.tar.gz
* installation: compiled from source


### Infernal

Used to map tRNA fragments to human tRNAs profiles 

* www: http://eddylab.org/infernal/
* doc (PDF): http://eddylab.org/infernal/Userguide.pdf
* version:  1.1.4 
* obtained from: https://anaconda.org/bioconda/infernal
* installation: ```conda install -c bioconda infernal```


## general tools

### bedtools 

* www/docs: https://bedtools.readthedocs.io/en/latest/
* version: 2.30.0
* obtained from: https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary
* installation: static binary

### clumpify from BBMap
* www:  https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmap-guide/
* doc: https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/clumpify-guide/
* version: 38.96
* obtained from: https://sourceforge.net/projects/bbmap/files/BBMap_38.96.tar.gz/download
* installation: download from the above www and unpack
* depencencies: java

### fastp

Fastq quality checking, low complecity reads filtering

* www: 
* version: 0.23.2
* obtained from: https://anaconda.org/bioconda/fastp
* installation: 

```
conda create --name fastp
conda activate fastp
conda install -c bioconda fastp
```



### ```usearch```

* version: usearch11.0.667_i86linux32
* obtained from: https://drive5.com/downloads/usearch11.0.667_i86linux32.gz
* installation: static binary

### ```fqgrep```

* www: https://github.com/indraniel/fqgrep
* doc: https://github.com/indraniel/fqgrep/blob/master/README
* obtained from: https://github.com/indraniel/fqgrep/archive/refs/tags/v0.4.4.tar.gz
* installation: compiled from source
* depencencies: 
    * libtre https://github.com/laurikari/tre

## system utils

This is a set of command line tools used for speed/convinience. One can  


### ```ripgrep```

Faster and with extra features ```grep``` replacement.

* www: https://github.com/BurntSushi/ripgrep
* doc: https://github.com/BurntSushi/ripgrep/blob/master/README.md
* version: 13.0.0
* obtained from: https://github.com/BurntSushi/ripgrep/releases/download/13.0.0/ripgrep-13.0.0-x86_64-unknown-linux-musl.tar.gz
* installation: static binary linked with musl

### ```choose```

Easy to use ```cut``` and partialy ```awk``` alternative.

* www: https://github.com/theryangeary/choose
* doc: https://github.com/theryangeary/choose/blob/master/readme.md
* version: 1.3.4
* obtained from: https://github.com/theryangeary/choose/releases/download/v1.3.4/choose-x86_64-unknown-linux-musl
* installation: static binary linked with musl

### pigz

Multithreaded ```gzip``` replacement.

* www: https://github.com/madler/pigz
* doc: see ```pigz --help```
* version: 2.7
* obtained from: https://github.com/madler/pigz/archive/refs/tags/v2.7.tar.gz
* installation: compiled from source

### fish shell

Alternative to ```bash`` shell.


### ```pypy3```

Faster Python implementation

* www: https://www.pypy.org/
* version:  Python 3.9.12 / PyPy 7.3.9 
* https://downloads.python.org/pypy/pypy3.9-v7.3.9-linux64.tar.bz2
* installation: executable binary after unpacking

### python libraries

Installed using ```pip```

#### main scripts
* pyfaidx
* numpy
* pandas
* seaborn

#### extras
* executor


