#! /usr/bin/perl
#===============================================================================
##
##		 FILE:  map_call_snps.pl
##
##		USAGE:  ./filter_SNPs.pl  
##
##  DESCRIPTION:  map with bwa and call snps
##
##	  OPTIONS:  ---
## REQUIREMENTS:  ---
##		 BUGS:  ---
##		NOTES:  ---
##	   AUTHOR:  Jody Phelan (mn), jody.phelan@lshtm.ac.uk
##	  COMPANY:  LSHTM
##	  VERSION:  1.0
##	  CREATED:  02/05/2015 11:51:02 AM
##	 REVISION:  29/03/2016 
##===============================================================================
#

use strict;
use warnings;
use List::Util 'max';
use Cwd;
use Cwd 'abs_path';

if (scalar @ARGV == 0){ print "\n################# map_call_snps.pl ################\n\n\ttrim - trim reads\n\tmap - map a sample with reads\n\tsamtools - call variants using samtools\n\tcoverage - create coverage file\n\tsamtools_coverage - create coverage file using samtools\n\tall - perform whole mapping pipeline\n\n################################################### \n\n"; exit;}


if ($ARGV[0] eq "trim"){
	if (scalar @ARGV != 3){ print "\nmap_call_snps.pl trim <sample> <threads>\n\n"; exit;}
	my $sample = $ARGV[1];
	my $threads = $ARGV[2];
	trim($sample,$threads);	
} elsif ($ARGV[0] eq "map"){
	if (scalar @ARGV != 4){ print "\nmap_call_snps.pl map <sample> <ref> <threads>\n\n"; exit;}
	my $sample = $ARGV[1];
	my $ref = $ARGV[2];
	my $threads = $ARGV[3];
	bwa_mapping($sample,$ref,$threads);
} elsif ($ARGV[0] eq "samtools"){
	if (scalar @ARGV != 3){ print "\nmap_call_snps.pl samtools <sample> <ref>\n\n"; exit;}
	my $sample = $ARGV[1];
	my $ref = $ARGV[2];
	samtools($sample,$ref);
} elsif ($ARGV[0] eq "coverage"){
	if (scalar @ARGV != 4){ print "\nmap_call_snps.pl coverage <sample> <ref> <threads>\n\n"; exit;}
	my $sample = $ARGV[1];
	my $ref = $ARGV[2];
	my $threads = $ARGV[3];
	extract_cov($sample,$ref,$threads);
} elsif ($ARGV[0] eq "samtools_coverage"){
	if (scalar @ARGV != 3){ print "\nmap_call_snps.pl samtools_coverage <sample> <ref>\n\n"; exit;}
	my $sample = $ARGV[1];
	my $ref = $ARGV[2];
	samtools_cov($sample,$ref);
} elsif ($ARGV[0] eq "all"){
	if (scalar @ARGV != 6){ print "\nmap_call_snps.pl all <sample> <ref> <threads> <working_dir> <storage_dir>\n\n"; exit;}
	my $sample = $ARGV[1];
	my $ref = $ARGV[2];
	my $threads = $ARGV[3];
	my $wd = $ARGV[4];
	my $sd = $ARGV[5];
	pipeline($sample,$ref,$threads,$wd,$sd);
} else { print "\nmap_call_snps.pl\n\n\ttrim - trim reads\n\tmap - map a sample with reads\n\tsamtools - call variants using samtools\n\tcoverage - create coverage file\n\n"; exit;}



#---------------------------------------------------------------------------
#  Trimming reads
#---------------------------------------------------------------------------
sub trim{
	my $script_dir = abs_path($0);
	$script_dir =~ s/map_call_snps_0.2.pl//;
	
	my $sample = $_[0];	
	my $threads = $_[1];
	if (!-e "./fastq/${sample}_1.fastq.gz" or !-e "./fastq/${sample}_2.fastq.gz"){ print "Can't find ./fastq/${sample}_1.fastq.gz\n"; exit;}
	if (-e  "${sample}_1_trimmed_paired.txt"){ print "Found trimmed reads\n"; exit;}
	print "Running trimmomatic on $sample\n";
	`java -jar $script_dir/trimmomatic.jar PE -threads $threads -phred33 ./fastq/${sample}_1.fastq.gz ./fastq/${sample}_2.fastq.gz ${sample}_1_trimmed_paired.txt ${sample}_1_trimmed_unpaired.txt ${sample}_2_trimmed_paired.txt ${sample}_2_trimmed_unpaired.txt LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 2>> $sample.log`;

}


#---------------------------------------------------------------------------
#  Running BWA
#---------------------------------------------------------------------------
sub bwa_mapping{
	
	my $sample = $_[0];
	my $ref = $_[1];
	my $threads = $_[2];
	print "Running BWA for $sample\n";
	my $bwamem = "bwa mem -t $threads -c 100 -R '\@RG\\tID:$sample\\tSM:$sample\\tPL:Illumina' -M -T 50 ";
	print "Mapping pairs...";
	`$bwamem $ref ${sample}_1_trimmed_paired.txt ${sample}_2_trimmed_paired.txt 2>>$sample.log | sambamba view -t $threads -S -f bam /dev/stdin 2>>$sample.log | sambamba sort -o ${sample}_pair.sorted.bam -t $threads /dev/stdin 2> $sample.log`;
	print "Mapping unpaired\n";
	`$bwamem $ref ${sample}_1_trimmed_unpaired.txt 2>> $sample.log | sambamba view -t $threads -S -f bam /dev/stdin 2>>$sample.log | sambamba sort -o ${sample}_single1.sorted.bam -t $threads /dev/stdin 2> $sample.log`;
	`$bwamem $ref ${sample}_2_trimmed_unpaired.txt 2>> $sample.log | sambamba view -t $threads -S -f bam /dev/stdin 2>>$sample.log | sambamba sort -o ${sample}_single2.sorted.bam -t $threads /dev/stdin 2> $sample.log`;

	`sambamba merge -t $threads ${sample}.bam ${sample}_pair.sorted.bam ${sample}_single1.sorted.bam ${sample}_single2.sorted.bam 2>> $sample.log`;

	`sambamba sort -t $threads -o $sample.sorted.bam $sample.bam 2>> $sample.log`;

	`mv $sample.sorted.bam $sample.bam`;
	`sambamba index -t $threads $sample.bam 2>> $sample.log`;

	`rm ${sample}_1_trimmed_paired.txt ${sample}_1_trimmed_unpaired.txt ${sample}_2_trimmed_paired.txt ${sample}_2_trimmed_unpaired.txt ${sample}_pair.sorted.bam ${sample}_single1.sorted.bam ${sample}_single2.sorted.bam `;
	`samtools flagstat $sample.bam > $sample.stats.txt 2>> $sample.log`;
}

#---------------------------------------------------------------------------
#  Samtools SNPs
#---------------------------------------------------------------------------
sub samtools{
	my $sample = $_[0];
	my $ref = $_[1];
	if( !-e "$sample.bam" ){ print "Can't find bam file\n"; exit;}
	if (-e "$sample.filt.vcf.gz"){ print "Found VCF\n"; exit;}
	
	print "Calling samtools variants on $sample\n";
	`samtools mpileup -B -Q 23 -d 2000 -C 50 -ugf $ref $sample.bam 2>> $sample.log | bcftools view -bvcg - > $sample.raw.bcf 2>> $sample.log`;
	`bcftools view $sample.raw.bcf | vcfutils.pl varFilter -d 10 -D 2000 > $sample.filt.vcf 2>> $sample.log`;
	`bgzip $sample.filt.vcf`;
	`tabix -p vcf $sample.filt.vcf.gz`;

}


#---------------------------------------------------------------------------
#  Extract Coverage
#---------------------------------------------------------------------------

sub extract_cov{

	my $sample = $_[0];
	my $reference = $_[1];
	my $threads = $_[2];

	`sambamba depth base -t $threads -q 20 -z -F  'mapping_quality > 30 and not unmapped' $sample.bam > $sample.temp.coverage`;
	
	my ($r1,$r2) = hash_contigs($reference); 
	my %hasContigs = %{$r1}; 
	my @arrContigs = @{$r2}; 
 
	my $x = 0; 
	my %cov; 
	open F, "$sample.temp.coverage" or die; 
	open OUT, ">$sample.coverage" or die;
	while(<F>){ 
	    chomp; 
	    if ($x<1){$x++; print OUT "$_\n"; next;} 
	    my ($chr,$pos) = (split /\s+/,$_)[0,1]; 
	    $pos = $pos; 
	    $cov{$chr}{$pos} = $_; 
	} 
	close(F); 
 
 
	foreach my $chr ( @arrContigs ) { 
	     
	    for (my $pos=0; $pos<=(length($hasContigs{$chr})-1); $pos++){ 
						
    	    if(exists($cov{$chr}{$pos})){ 
	            print OUT "$cov{$chr}{$pos}\n"; 
	        } else { 
	            print OUT "$chr\t$pos\t0\t0\t0\t0\t0\t0\t0\t$sample\n"; 

	        } 
	    } 
	 
	}
	close(OUT);

	`rm $sample.temp.coverage`;



}

#---------------------------------------------------------------------------
#  Hash Contigs
#---------------------------------------------------------------------------
sub hash_contigs {
  if (@_ != 1) {
	print "Usage: hash_contigs contigs_file";
	exit;
  }
  my $contigs_file = shift;
  if( $contigs_file =~ /\.gz$/ ){
	open(CONTIGS, $contigs_file, "gunzip -c $contigs_file |" ) or die "Cant open contigs file: $!\n";
  } else {
	open(CONTIGS, $contigs_file) or die "Cant open contigs file: $!\n";
  }
  my %contigs_hash; # hash to store contig names and sequences
  my $contigName;
  my @ids;
  while (<CONTIGS>) {
	if (/^>(\S+)/) {
	  $contigName=$1;
	  push @ids, $contigName;
	} else {
	  chomp;
	  $contigs_hash{$contigName} .= $_;
	}
  }
  close(CONTIGS);
  return(\%contigs_hash,\@ids);
}



#---------------------------------------------------------------------------
#  Whole pipeline
#---------------------------------------------------------------------------
sub pipeline{
	
	my $sample = $_[0];
	my $threads = $_[2];
	my $ref = $_[1];
	my $wd = $_[3];
	my $sd = $_[4];
	if( !-d "$wd"){print "$wd not found";exit;}
	$wd  = abs_path($wd);
	$sd  = abs_path($sd);
	$ref = abs_path($ref);	

	print "Working directory: $wd\n";
	print "Storage directory: $sd\n";


	if( !-d "$sd/bam"){`mkdir $sd/bam`;}
	if( !-d "$sd/vcf"){`mkdir $sd/vcf`;}
	if( !-d "$sd/coverage"){`mkdir $sd/coverage`;}
	if( !-d "$sd/logs"){`mkdir $sd/logs`;}
	
	if( !-d "$wd/fastq"){`mkdir $wd/fastq`;}
	my $read1 = abs_path("./fastq/${sample}_1.fastq.gz");
	my $read2 = abs_path("./fastq/${sample}_2.fastq.gz");
	chdir("$wd");
	unless ($wd eq $sd){
		`ln -s $read1 fastq/`;
		`ln -s $read2 fastq/`;
	}

	trim($sample,$threads); 
	bwa_mapping($sample,$ref,$threads);
	samtools($sample,$ref);
	extract_cov($sample,$ref,$threads);
	
	

	#---------------------------------------------------------------------------
	#  clean up directory
	#---------------------------------------------------------------------------
	`bgzip $sample.coverage`;
	`rm ${sample}.raw.bcf`;
	`mv ${sample}.bam ${sample}.bam.bai ${sample}.stats.txt $sd/bam/`;
	`mv ${sample}.filt.vcf.gz  ${sample}.filt.vcf.gz.tbi $sd/vcf/`;
	`mv ${sample}.coverage.gz $sd/coverage/`;
	`mv $sample.log $sd/logs/`;
	

}


sub samtools_cov{
	my $sample = $_[0];
	my $ref_seq = $_[1];
	my ($r1, $r2) = hash_contigs($ref_seq);
	my %lhChromosomes  = %{$r1};
	my @laChromosomes = @{$r2};
	
	`samtools mpileup -B -Q 23 -d 2000 -C 50 -f $ref_seq $sample.bam > $sample.pileup`;
	
	if( !-e "$sample.pileup"){ print "Can't find pileup\n"; exit;}
	if( -e "$sample.coverage"){print "Found coveraqge\n"; exit;}

	open HC, ">$sample.coverage";
	print HC "REF\tPOS\tA\tC\tG\tT\tNtDEL\tN\tSAMPLE\n";
	foreach my $chr (@laChromosomes) {
		my %lhCoverage;
		my $com = "awk \' { if (\$1 == \"$chr\") { print } } \' $sample.pileup";
		my @laPileup = `$com`;
		foreach (@laPileup) {
			chomp;
	  		my @laLine = split /\s+/,$_;
			my $liLine = scalar(@laLine);
			my ($chr, $pos, $ref, $total) = @laLine[0..3];
			my $id = $chr."%".$pos;
			$ref = uc($ref);
			if ($liLine >= 5 && $total > 0) {
				my ($lsPile, $ins, $del) = ($laLine[4], 0, 0);
				$lsPile = uc($lsPile);
				$lsPile =~ s/\^.//g;
				$lsPile =~ s/\$//g;
				$lsPile =~ s/\*/N/g;
				$lsPile =~ s/[.,]/$ref/g;
				$ins = () = $lsPile =~ /\+/g;
				$del = () = $lsPile =~ /\-/g;
				my @laIndels = split /[\+\-]/, $lsPile;
				$lsPile = "";
				foreach (@laIndels){
					my @laTemp = /(\d+)(.*)/;
					if (@laTemp){
						$lsPile .= substr($laTemp[1],$laTemp[0]);
					} else {
						$lsPile .= $_;
					}
				}
				my ($a, $c, $g, $t, $n) = (0, 0, 0, 0, 0);
				$a = $lsPile =~ tr/A//;
				$c = $lsPile =~ tr/C//;
				$g = $lsPile =~ tr/G//;
				$t = $lsPile =~ tr/T//;
				$n = $lsPile =~ tr/N//;
			
				if ($total != ($a + $c + $g + $t + $n)) {
					open H, ">>coverage/$sample.error";
					die "$chr\t$pos\t$total\t$a\t$c\t$g\t$t\n";
					close H;
 					$total = $a + $c + $g + $t + $n;
				}

				my ($abs, $rel) = ("0", "0");
				if ($total > 0) {
					$abs = max(($a,$c,$g,$t));
					$rel = $abs/$total;
					$rel = sprintf("%0.3f",$rel);
					$lhCoverage{$id} = "$chr\t$pos\t$total\t$a\t$c\t$g\t$t\t$del\t$n\t$sample";
				}
			}
		}
		my $len = length($lhChromosomes{$chr});
#		my @seq = split //, $lhChromosomes{$chr};
		for (my $i = 1; $i <= $len; $i++) {
			my $id = $chr."%".$i;
			if (defined $lhCoverage{$id}) {
				print HC "$lhCoverage{$id}\n";
			} else {
				print HC "$chr\t$i\t0\t0\t0\t0\t0\t0\t0\t$sample\n";
#				if (defined $seq[$i-1]) {
#					my $ref = uc($seq[$i-1]);
#					print HC "$chr\t$i\t$ref\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n";
#					
#				}
			}
		}
	}
	close HC;
}

