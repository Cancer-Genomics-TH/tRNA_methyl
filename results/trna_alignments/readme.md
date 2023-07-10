# tRNA alignmens

Files in this directory were constructed using gtRNAdb data:

## trnas_hs38_gtrnadb_fromhtml.sto

* parsed aligments from HTML page: http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/Hsapi38-align.html
* script: `parse_gtrnadb_alignments.py`

## trnas_hs38_gtrnadb_algn_spliced.sto

* source: http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-mature-tRNAs.fa
* spliced trnas:
  * spld-Arg-TCT
  * spld-Ile-TAT
  * spld-Leu-CAA
  * spld-Tyr-GTA
  * spld-Tyr-ATA

Sequences corresponding to the above isoforms were extracted as individual files and aligned using mlocarna from LOCARNA package.
Sequence IDs were added to Stockholm alignment files using `fix_stockholm_id.py` script and concatenated afterwards. 