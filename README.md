# rnaseq_as_analysis
Alternative Splicing analysis using docker and singularity

## Usage:
``` 
perl run_as.pl <options>

-ref_fasta <path to reference .fa/fasta>
-gtf_file <gtf file to be used in rMATS>
-gtf_index_file <gtf index to be used in kallisto>
-annotation_file <annotation file to be used in SpliceTools>
-output_dir <main output directory, default is directory script is current working directory>
-mapping_file <mapping file (see format below)>
-grouping_file <file with sample - grouping information (see format below)>
-comparison_file <file with group comparisions (see format below)>
-alignments_dir <directory with all bams. Bam filenames must end with associated sample id found in mapping file + ".bam">
-strand <options: reverse, forward, unstranded default: unstranded>
```

## mapping_file
Mapping file must be tab delimited, and will contain this format:
```
_1  SampleName  SequencingMachine   PathToFQ    PE/SE
```
**_1** - Unused column for this pipeline. Put anything here  
**SampleName** - Sample name - Will be used to link samples to grouping file, and bams (bams must end in `<sample name>` + ".bam")  
**SequencingMachine** - Unused column for this pipeline. Put anything here  
**PathToFQ** - Path to fastq files for this sample. Inside this directory there should be one or more `fastq.gz` files. If paired end, filenames must contain a `_R1_` and `_R2_` in them.  
**PE or SE** - This column must contain either the letters `PE` or `SE`. `PE` means that the reads are Paired End, while `SE` means they are Single ended. 

## grouping_file
Grouping file must be tab delimited, and will contain two columns. First column is the sample name (as found in mapping file). Second column will contain a group name. 
Example:
```
sample1ctl	exp1ctrl
sample2ctl	exp1ctrl
sample3ctl	exp1ctrl
samp1test	exp1test
samp2test	exp1test
samp3test	exp1test
```

## comarison_file
Comparisons file must be tab delimited. Each line contains two group names that we will be comparing. 
Example:
```
exp1ctrl	exp1test
exp2ctrl	exp2test
exp1ctrl	exp2ctrl
```

## annotation_file
Splicetools annotation file requires the 3rd column to be setup as `geneID_geneSymbol`. I basically just altered `gtf_index_file` in order to get this correct. 

## How tools are run
Currently the tools are run by singularity using a docker container that is hardcoded into the code. Obviously if this code were to be developed more we could fix that, but until it's actually going to be used in the future, we shouldn't waste the time.
