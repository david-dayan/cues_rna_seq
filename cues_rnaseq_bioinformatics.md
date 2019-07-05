Cues RNAseq Bioinformatics
================

Summary: A log of bioinformatics work to convert raw sequence reads to
read counts for the plasticity project and the predator cues project.

A few notes about directories, permissions and data storage:

1.  Only a few of the code chunks are meant to be run locally. Most of
    them are slurm scripts or short sets of commands meant to run on the
    cluster, so paths and directories refer to cluster, not local
    directories. Be sure to tell R not to run the code chunks (e.g. set
    " eval = FALSE " in the code chunk header)

2.  Any data to be used locally (on our own computers) MUST be stored in
    the github repository (and therefore on your local copy of the
    repository.

3.  Similarly, any local paths must refer to the github repository (./)

4.  Commit your work often. Always pull before pushing any changes.

5.  Knit any changes to Rmarkdown and Rnotebook documents before pushing
    your changes so that others can see the rendered work.

6.  If any local code chunks are computationally intensive, cache the
    results, so that others won’t have to run them again, (e.g. add
    “cache = TRUE” to the code chunk header)

### Outline

The approach has 5 steps:

  - Data Storage
  - Data QC  
  - Demultiplexing
  - Alignment  
  - Gene Read Counting

We roughly follow the approach put forward by Lexogen for steps 4 and 5

#### Detailed approach outline

  - Move all data to appropriate directories and change permissions  
  - QC
      - run fastqc on all raw datasets  
      - trim and clean up reads using bbduk.sh: remove polyA, adapter
        contamination and low quality reads  
      - run fastqc again to look at success of (b)  
  - Demultiplex in one of two strategies
      - using “fixed” fastq files (i.e. with index replaced in header) -
        defq, demuxbyname (bbmap)  
      - using separate fastq files for read 1 and i7 read - deML or the
        q2-demux plugin from QIIME  
      - at this step the files are ready to quality controlled and
        demulitplexed and ready to input into bluebee pipeline if
        desired (or possible)  
  - Alignment
      - align the cleaned, demultiplexed reads to the stickleback
        transcriptome/genome using STAR aligner
      - assess alignment quality using RSEQC, if looks it looks good
        proceed
  - Gene read counting
      - take the individual aligned bam files and move forward with
        making gene counts using HTSeq-count
      - after this we’re ready to move on with the actual analysis from
        the read counts :) using DESeq2 or similar

### Data Storage

Because a lot of data will be shared and room is limited in the /ddayan/
directory, will store raw sequence data in Melissa’s or Dale’s cluster
directory.

Also need to change permissions so RWX for all
    users

    # in directory containing the directory containing the reads (e.g. /home/dalstevens/)
    chmod -R 777 ./rnaseq_data/30-196056229

### Quality Control

#### First look

A first look at the data using fastqc and (zcat | head) on a single lane
(QuantSeq-1) revealed that the data met expectations except no index
read was made\! Contacted genewiz and they were able to send index reads
as separate files.

Other fastqc takeaways: \* Read quality looks consistently high  
\* Begins to hit polyA tail after around read 60-80 for majority of
reads  
\* Adapter sequence picks up after polyA  
\* First 8 or so look nonrandom (per base sequence content looks
biased), but this makes sense given the random primer

A first look of the adapter sequences show they have the same number of
reads as the corresponding forward read files and that the first six
base pairs correspond to i7 sequences from the quantseq user guide
(7001-7080)

#### Index Sequences

Adding index sequences back to reads

May not need to be done if we demultiplex using a utility that can
handle a separate index read file such as deML or QIIME (q2-demux) -
will wait until we plan further

##### Summary

Genewiz failed to make an index read for the RNAseq data, when they
fixed this they made an index only read and sent the resulting fastq
file. We need to add the index back to the right place in fastq header.

##### Approach outline

This will take 3 steps for each file:

1)  truncate the index read seqeunce to 6bp - genewiz gave us an 8bp
    index read, our indexes are only 6bp, so why not clean them up now
2)  remove the current “index” from the fastq files (currently it
    corresponds to the lane number and will always be a number between 1
    and 8
3)  concatenate the index to the end of the header for each read

Here’s the example code for each
    step:

    where "R1.fq.gz" is the forward read and "index1.fq.gz" is the index read

1)  index trimming - using cutadapt

<!-- end list -->

``` 
    /path/to/cutadapt --cut -2 index1.fq.gz index1_trimmed.fq.gz
```

2)  delete sample id / index from
reads

<!-- end list -->

``` 
    awk 'NR%4==1{sub(/[12345678]$/,"")}1' <(gunzip -c R1.fq.gz ) | gzip > tmp_reads
```

3)  add index sequence to
header

<!-- end list -->

``` 
    paste -d '~' <(gunzip -c tmp_reads) <(gunzip -c index1_trimmed.fq.gz) | perl -F'~' -lane 'push(@buffer, $F[0]); if($line == 1){@buffer[0] .= "$F[1]"}; if(($line == 3) && @buffer){print join("\n",@buffer); @buffer = ()}; $line = ($line+1) % 4;' | gzip - > R1_windex.fq.gz
```

##### Scaling

each of these three steps need to be done in sequence for all of the
files and the scripts need to be adapted to run on the server in a job

tasks to complete for slurm script: - figure out path to cutadapt maybe:
/home/ddayan/.local/bin/cutadapt if installed with PIP to instal with
pip: `pip install --user --upgrade cutadapt` - set up variable system to
cover all the files in a directory and loop through - or just write the
command 32 times and run in sequence

slurm script:

``` bash

#!/bin/bash

# set max wall-clock time (D-HH:MM:SS)
#SBATCH --time=0-23:59:00

#SBATCH --cpus-per-task=1

# set partition/queue to use
#SBATCH --partition=day-long-cpu

# set name of job
#SBATCH --job-name=index1

# set name of output file
#SBATCH --output=index1.out

# mail alert at start, end and abortion of execution
#SBATCH --mail-type=ALL

# send mail to this address
#SBATCH --mail-user=ddayan@clarku.edu

cutadapt --cut -2 Quantseq-1_S1_I1_001.fastq.gz index_trimmed.fq.gz
awk 'NR%4==1{sub(/[12345678]$/,"")}1' <(zcat Quantseq-1_R1_001.fastq ) | gzip > tmp_reads
paste -d '~' <(zcat tmp_reads) <(zcat index_trimmed.fq.gz) | perl -F'~' -lane 'push(@buffer, $F[0]); if($line == 1){@buffer[0] .= "$F[1]"}; if(($line == 3) && @buffer){print join("\n",@buffer); @buffer = ()}; $line = ($line+1) % 4;' | gzip - > Quantseq-1_w_index.fq.gz
rm index_trimmed.fq.gz
rm tmp_reads
```
