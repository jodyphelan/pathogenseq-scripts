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

### GEM
wget http://master.dl.sourceforge.net/project/gemlibrary/gem-library/Binary%20pre-release%202/GEM-binaries-Linux-x86_64-core_2-20121106-022124.tbz2
tar -xvf GEM-binaries-Linux-x86_64-core_2-20121106-022124.tbz2 
ln -s GEM-binaries-Linux-x86_64-core_2-20121106-022124/gem-* .

### BEDOPS
wget https://github.com/bedops/bedops/releases/download/v2.4.16/bedops_linux_x86_64-v2.4.16.tar.bz2
tar -xvf bedops_linux_x86_64-v2.4.16.tar.bz2 
mv bin/ bedops-bin
ln -s bedops-bin/wig2bed .

