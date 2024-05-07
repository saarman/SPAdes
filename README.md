# SPAdes
Culex species identification and blood meal analysis with de novo assembly of illumina reads from amplicon resequencing (COI, Ace2, cqm1)

# A place to keep track of Culex blood meal bioinformatics steps

## Logging onto CHPC with Terminal on a mac


## Outline of steps:
1. seqtk to subsample reads (tens to hundreds of 1000's?)
2. fastqc to trim/clean/quality control
3. de novo assembly with SPAdes https://biomedicalhub.github.io/genomics/04-part4-denovo-assembly.html
4. map reads back to assembly and confirm which have good support
5. blast all consensus sequences (.fa) - relabel as COi if match COi
6. blast/bold consensus of COi matches only

## Manuals for seqtk, fastqc, SPades 
 - seqtk https://docs.csc.fi/apps/seqtk/#usage   
 - fastqc https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
 - SPAdes https://ablab.github.io/spades/ 

# Task 1: seqtk 










