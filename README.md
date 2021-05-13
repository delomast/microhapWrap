# microhapWrap v1
A wrapper for using microTyper (https://github.com/delomast/microTyper) and 
calculating a set of useful summary statistics. 

microhapWrap doesn't do much on its own, it delegates most tasks and 
then collects and outputs the results. It therefore relies on having 
several other programs in the system path. It is written 
in Python (tested in 3.6.10, requires >= 3.5) for Linux systems. 
The programs it calls are

- wc
- samtools
- cut
- mtype2 (part of microTyper)
- bowtie2 (but only if you use it to align reads, if you produce bam files yourself, this is not needed)

microhapWrap is able to provide all or a subset of

- genotypes for microhaplotypes and SNPs
- genotypes for presence/absence markers (typically sex markers)
- proportion of genotypes called per individual (genotyping success)
- proportion of reads that had valid alignments
- total number of reads
- total number of reads with valid alignments
- an attempt at a potential contamination score
- counts of reads aligning to each locus within each individual
- counts of reads beginning with each forward primer within each individual

Call as:

```
  python microhapWrap.py -f ...
```
or if it is executable directly on your system, 
```
  microhapWrap.py -f ...
```

## Inputs
microhapWrap assumes one fastq file for each sample and assumes that each read is represented by four 
lines (so no multi-line sequences). microhapWrap can align reads through bowtie2, or you can align reads 
separately and provide it with bam files. If you provide bam files, ensure that the names are 
identical to the fastq files but with the extension ".bam" (e.g., sampleOne.fastq and sampleOne.bam). 
If you provide just the bam files (no fastqs) microhapWrap will genotype the samples and calculate what it can,
but some summary statistics will not be calculated.

microhapWrap takes the inputs needed for mtype2 and some other optional inputs (see below) depending on exactly what you 
want it to produce.

### Arguments

- `-f` a space separated list of fastq files. If not specified, then all files in the 
  current directory ending in ".fastq" are used. Ex: `-f sample*.fastq`, `-f sample1.fastq sample2.fastq`
- `-bam` For use if you have aligned the reads yourself. A space separated list of bam files to use. If you 
  specify `-bam` but do not provide any files, then all files in the current directory ending in 
  ".bam" are used.
- `-t` number of threads to use. Default is to use as many threads as are available to the process.
- `-p` the -p input to mtype2 giving the position file
- `-r` the -r input to mtype2 giving the reference fasta file
- `-b` the -b input to mtype2 giving the batch size (default 200)
- `-c` the -c input to mtype2 giving the threshold posterior probability (default 0.99)
- `-eps` the -eps input to mtype2 giving the read error rate (default 0.01)
- `-ploidy` the ploidy of the individuals being genotyped (default 2)
- `-d` the -d input to mtype2 giving the minimum depth (default 5 \* ploidy)
- `-pa` if you want to genotype presence/absence markers, this is the input file (see below)
- `-bt2ref` if you want microhapWrap to align your reads, this is the path/prefix to the bowtie2 reference
- `-fwd` if you want microhapWrap to count reads beginning with each forward primer within each individual, this is the input file (see below)
- `--alignByLocus` pass this argument to have microhapWrap write out a file with counts of reads aligning to each locus within each individual

## Method

### Presence/Absence Markers
Presence/absence markers are heuristically genotyped using the ratio of the number of reads that align to the 
microhaplotype/SNP loci (e.g. non-presence/absence loci) compared to the number that align to the locus. microhapWrap 
calculates this ratio and then compares it to user specified thresholds as well as allowing a user to specify a minimum number 
of reads aligning to microhaplotype/SNP loci for a presence absence marker to be called. 

microhapWrap uses a MaxRatioPres and MinRatioAbs threshold with MinRatioAbs > MaxRatioPres. If 
ratio > MinRatioAbs, the absent genotype is called. If ratio < MaxRatioPres, the presence genotype is called. 
If MinRatioAbs > ratio > MaxRatioPres, the genotype is missing. The input file is tab-separated, 
_has a header line_, has one row for each presence/absence marker, and has six columns. The columns are, in order,

1. Name of the locus (should match name in reference fasta)
1. Character string representing a present genotype
1. Character string representing an absence genotype
1. Minimum number of reads aligning to microhaplotype/SNP loci for the genotype to be called
1. MaxRatioPres
1. MinRatioAbs

See the "example_pres_abs_input.txt" file for an example. I suggest you determine the minimum 
number of reads, MaxRatioPres, and MinRatioAbs empirically for a given locus.

### Count forward primer matches
Forward primer sequences are input as a tab-separated file _with_ a header line. The first 
column is the locus name and the second is the forward primer sequence (can be given as a 
regular expression as interpreted by the python re library). microhapWrap will search through the fastq files 
for reach individual and count the number of reads that _start_ with exact matches to the input primers (or regex). 
I find this useful when developing/optimizing a panel. It's a little slow (python isn't the fastest at text operations) 
so I don't usually run it for routine genotyping with established panels. 

### Contamination Score
The contamination score is a heuristic calcualted for each individual using microhaplotype/SNP 
loci with called genotypes. The score for a single locus is calculated 
as |(number of reads of Allele1 / total number of reads aligning to the locus) - (dosage of Allele1 / ploidy)|. 
This is averaged across all loci with called genotypes (i.e., loci with missing genotypes are skipped) to give something that might be 
described as the mean deviation of read ratios from their expectation. The higher the value, the more likely there 
is something wrong with the sample. Samples with low genotyping success (and few reads) will likely also have higher values of 
this metric, but poor performing samples are typically regenotyped or omitted from analyses anyway. We have found 
that utilizing this metric is best performed by comparing the values across samples sequenced with a given panel. It 
is frequently informative to view a scatterplot of contamination score vs number of genotyped loci. Any 
outliers are typically obvious. 

This is _not_ a perfect metric for assessing potential contamination, but it seems to work ok in our limited 
experience. Better metrics could certainly address this, and if you have developed one, please let me 
know - I would be interested in implementing it.


## Output Files

- alignment_output.txt: whatever output bowtie2 prints to stdErr during alignment
- microhap_genotypes.txt: genotypes, output by mtype2
- pres_abs_genotypes.txt: genotypes of presence/absence loci
- fwd_primer_counts.txt: counts of reads beginning with each forward primer sequence. All loci x individual combinations 
  should be present, even if the count is 0.
- summary_stats.txt: For each sample, the proportion of markers genotyped (includes any presence/absence markers), 
  proportion of reads that had valid alignments (number of mapped reads in the bam file / total number of reads), 
  contamination score, number of mapped reads in the bam file, number of reads in the fastq file
- align_by_ind_locus.txt: number of reads aligning to each locus within each individual. If a locus has 0 reads 
  aligning within an individual, _it is not listed in the file_.
