# pathogenseq-scripts

## Installation

#### Prerequisites:
* java<br>
* perl<br>
* samtools<br>
* bwa<br>
* vcftools<br>
* R<br>

#### Perl modules:
* Statistics::Lite<br>
* Term::ProgresBar

These programs and perl modules must be installed for full functionality of the pipeline. Additionally the following programs and scripts much also be added to the path, however we provide a script to download the latest versions:<br>
* trimmomatic<br>
* sambamba<br>
* vcfutils.pl<br>
* GEM-library<br>
* bedops<br>

#### To install please follow these instructions:
	git clone https://github.com/jodyphelan/pathogenseq-scripts.git
	cd pathogenseq-scripts
	bash downloadPrograms.sh
	
##### Before using the pipeline add the scripts and programs to the path each time you start a new session:
	source add_to_path.source 


## Example usage:
To map fastq files to a reference and perform variant calling we can use the `map_call_snps_0.2.pl` script. This can be done either step by step using the different filtering modules individually or can be performed in one run. A simple directory structure and naming scheme must be followed in order for the pipeline to run correctly. Make a directory called `fastq` and place the reads using the naming scheme `<sample-id>_1.fastq.gz` and `<sample-id>_2.fastq.gz` for paired end reads. Navigate back to the directory containing the fastq folder and you are ready to run the pipeline:<br>
```
map_call_snps.pl all <sample> <ref> <threads> <working_dir> <storage_dir>
```

| Option | Description |
| ------ | ----------- |
| sample | The sample ID used in naming the fastq files |
| ref	 | The reference fasta |
| threads | The number of threads to be used |
| working_dir | The directory used to process data and store temporary files (most cases the current directory) |
| storage_dir | The directory used to store data (most cases the current directory) |

The final structure will look like the following:
```
base_dir/
├──fastq/
│  ├──sample1_1.fastq
│  └──sample1_2.fastq
├──coverage/
│  └──sample1.coverage.gz
├──vcf/
│  └──sample1.filt.vcf.gz
└──bam/
   └──sample1.bam
```

To create a SNP matrix from the VCFs and coverage files we will use the `filter_SNPs_MT_0.2.pl` script. This can be done either step by step using the different filtering modules individually or can be performed in one run. Before we can run the main pipeline we should create a mappability file listing the unique regions in the genome:<br>
```
filter_SNPs_MT.pl mappability <ref> <kmer> <threads>
```

| Option | Description |
| ------ | ----------- |
| ref   | The reference fasta to use |
| kmer | The size of the kmer to look for uniqueness (should be roughly the same as the read length) |
| threads | Number of threads to use |

This will create a file named mappability.bed in the current directory which will feed into the pipeline. We will also be filtering based on the total coverage over all samples. Because non-nuclear DNA will have higher coverage we do not want to filter these based on the same cut-offs. We must give a file listing all the chromosomes (one per line) which we would like to apply this filtering to.
The main command requires a number of different parameters which will be used during different stages in the filtering process: <br>
```
filter_SNPs_MT.pl all <sample_file> <base_dir> <%cov> <min_cov> <threads> <samtools|gatk|both> <mappability_file> <nuclear_chromosomes>
```

| Option | Description |
| ------ | ----------- |
| sample_file | File containing the names of all samples (one per line). The naming of vcf and coverage files must be consistent with these names  |
| base_dir | Directory containing two directories names ‘coverage’ and ‘vcf’  |
| %cov | The minimum allelic frequency of the alternate allele for it to be called |
| min_cov | Minimum coverage for a position to be considered. If lower it will be marked as missing |
| threads | Number of threads to use during processing |
| samtools\|gatk\|both | Which source of VCFs to use. Samtools should have the extension .filt.vcf.gz GATK should have the extension .gatk.raw.vcf.gz. Using both will use the union of the results.| 
| mappability_file | The mappability file to use. If you created this using mappability function of the script is will be names mappability.bed. |
| nuclear_chromosomes | The chromosomes which should use the total coverage filter |


This will perform several filters on the dataset interactively. For each of these steps a graph plot will appear and ask you to set cut-off values. Please click (anywhere on the plot) at the corresponding value you would like to set the cut-off at using the indicated axis:<br>
+ Filtering samples with large number of missing or mixed calls and those with low median depth
	- A plot of the number of missing calls as a proportion of the total number of SNP positions. (Y-axis)
	- A plot of the number of mixed calls as a proportion of the total number of SNP positions (Y-axis)
	- A plot of the number of the median depth per sample (Y-axis)
+ Filtering out SNPs with large number of missing or mixed calls
	- A plot of the number of missing calls as a proportion of call per SNP. (Y-axis)
	- A plot of the number of mixed calls as a proportion of call per SNP. (Y-axis)
+ Filtering out SNPs with low coverage over all samples
	- A plot containing the density function of the total coverage per SNP (X-axis)



