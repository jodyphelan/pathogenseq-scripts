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


### Example usage:
To create a SNP matrix from the VCFs and coverage files we will use the `filter_SNPs_MT_0.2.pl` script. This can be done either step by step using the different filtering modules individually or can be performed in one run. Before we can run the main pipeline we should create a mappability file listing the unique regions in the genome.<br><br>
     `filter_SNPs_MT.pl mappability <ref> <kmer> <threads>`<br>

| Option | Description |
| ------ | ----------- |
| ref   | The reference fasta to use |
| kmer | The size of the kmer to look for uniqueness (should be roughly the same as the read length) |
| threads | Number of threads to use |
