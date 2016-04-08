#! /bin/bash
### VCFTOOLS.PL
wget https://raw.githubusercontent.com/lh3/samtools/master/bcftools/vcfutils.pl
chmod 755 vcfutils.pl

### SAMBAMBA
wget https://github.com/lomereiter/sambamba/releases/download/v0.6.0/sambamba_v0.6.0_linux.tar.bz2
tar -xvf sambamba_v0.6.0_linux.tar.bz2
mv sambamba_v0.6.0 sambamba
chmod 755 sambamba

### TRIMMOMATIC
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip
unzip Trimmomatic-0.36.zip 
ln -s Trimmomatic-0.36/trimmomatic-0.36.jar trimmomatic.jar

