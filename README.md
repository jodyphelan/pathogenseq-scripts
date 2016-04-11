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
