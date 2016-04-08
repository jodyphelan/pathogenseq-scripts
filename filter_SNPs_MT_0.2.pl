#!/usr/bin/perl 
#===============================================================================
#
#		 FILE:  multi_threaded.pl
#
#		USAGE:  ./multi_threaded.pl  
#
#  DESCRIPTION:  
#
#	  OPTIONS:  ---
# REQUIREMENTS:  ---
#		 BUGS:  ---
#		NOTES:  ---
#	   AUTHOR:  Jody Phelan (JP), jody.phelan@lshtm.ac.uk
#	  COMPANY:  LSHTM
#	  VERSION:  1.0
#	  CREATED:  09/02/15 20:11:36
#	 REVISION:  ---
#===============================================================================

use threads;
use strict;
use warnings;
use Statistics::Lite qw(:all);
use POSIX qw(ceil);
use Term::ProgressBar;
use Cwd 'abs_path';

my $usage = 'filter_SNPs_MT.pl

	all - perform the whole pipeline

###### Modules ######

	raw - create a raw matrix and quality metrics
	sampleFilt - Filter using a missingness filter
	snpQual - Get quality of SNPs
	mappability - generate mappability file for genome
	snpFilt - Filter SNPs based on missingness and mixednessi
	covFilt - Filter based on the total coverage over all samples
	
###### Experimental ######

	indelMat - Create matrix with indel calls

';

#---------------------------------------------------------------------------
#  Check inputs
#---------------------------------------------------------------------------
if (scalar @ARGV == 0){ print $usage; exit;}

if ($ARGV[0] eq "raw"){
	if (scalar @ARGV != 7){ print "\nfilter_SNPs_MT.pl raw <sample_file> <base_dir> <%cov> <min_cov> <threads> <samtools|gatk|both>\n\n"; exit;}
	raw($ARGV[1],$ARGV[2],$ARGV[3],$ARGV[4],$ARGV[5],$ARGV[6]);
} elsif ($ARGV[0] eq "sampleFilt"){
	if (scalar @ARGV != 4){ print "\nfilter_SNPs_MT.pl sampleFilt <Missing cutoff> <Mixed cutoff> <Cov cutoff>\n\n"; exit;}
	sampleFilter($ARGV[1],$ARGV[2],$ARGV[3]);
} elsif ($ARGV[0] eq "snpQual"){
	if (scalar @ARGV != 1){ print "\nfilter_SNPs_MT.pl snpQual\n\n"; exit;}
	snpQual();
} elsif ($ARGV[0] eq "snpFilt"){
    if (scalar @ARGV != 4){ print "\nfilter_SNPs_MT.pl <Missing cutoff> <Mixed cutoff> <mappability file>\n\n"; exit;}
    filterSNPs($ARGV[1],$ARGV[2],$ARGV[3]);
} elsif ($ARGV[0] eq "covFilt"){
	if (scalar @ARGV !=4){ print "\nfilter_SNPs_MT.pl <lower-bound> <upper-bound> <nuclear_chromosomes>\n\n"; exit;}
	covFilter($ARGV[1],$ARGV[2],$ARGV[3]);
} elsif ($ARGV[0] eq "indelMat"){
	if (scalar @ARGV != 4){ print "\nfilter_SNPs_MT.pl <sample_file> <base_dir> <threads>\n\n"; exit;}
	createIndelMat($ARGV[1], $ARGV[2], $ARGV[3]);
} elsif ($ARGV[0] eq "mappability"){
	if (scalar @ARGV != 4){ print "\nfilter_SNPs_MT.pl <ref> <kmer> <threads>\n\n"; exit;}
	generateMapFile($ARGV[1],$ARGV[2],$ARGV[3]);
} elsif ($ARGV[0] eq "all"){
	if (scalar @ARGV != 9){ print "\nfilter_SNPs_MT.pl all <sample_file> <base_dir> <%cov> <min_cov> <threads> <samtools|gatk|both> <mappability_file> <nuclear_chromosomes>\n\n"; exit;}
	all($ARGV[1],$ARGV[2],$ARGV[3],$ARGV[4],$ARGV[5],$ARGV[6],$ARGV[7],$ARGV[8]);
} else {
	print $usage; exit;
} 

sub all{

my $sampleFile = $_[0];
my $base_dir = $_[1];
my $percentCov = $_[2];
my $minCov = $_[3];
my $threads = $_[4];
my $vcf = $_[5];
my $mapFile = $_[6];
my $nuclear = $_[7];

open FLOG, ">filtering.log" or die;
writeRfiles();
raw($sampleFile,$base_dir,$percentCov,$minCov,$threads,$vcf);
print "Filtering samples\n\n";
my ($l1,$l2,$l3) = samplePlot();
print FLOG "Sample Missingness = $l1\nSample Mixedness = $l2\nSample Median Coverage = $l3\n";
sampleFilter($l1,$l2,$l3);
print "Calculating snp qualities\n";
snpQual();
my ($r1,$r2) = snpPlot();
print FLOG "SNP Missingness = $r1\nSNP Mixedness = $r2\n";
filterSNPs($r1,$r2,$mapFile);
print "\nFiltering SNPs\n";
my ($q1,$q2) = covPlot();
print FLOG "SNP Total Coverage Lower = $q1\nSNP Total Coverage Higher = $q2\n";
covFilter($q1,$q2,$nuclear);

print "\n\n**********************\nFinished filtering\n**********************\n";

}



sub raw{
### SUB START ###

#---------------------------------------------------------------------------
#  Open VCFs and store all SNP locations
#---------------------------------------------------------------------------

my @samples;
my $base_dir = $_[1];
my $percentTotalcov = $_[2];
my $min_cov = $_[3];
my $threads = $_[4];
my $vcfType = $_[5];



if (($vcfType ne "both") and ($vcfType ne "samtools") and ($vcfType ne "gatk")){
	print "State which VCF to use using both|samools|gatk\n";
	exit;
}
my %lhAll_pos;
my %final_calls;
my %ref_calls;
my %unionINDEL_pos;
my %lhAllINDEL_pos;


open F, $_[0] or die;
while(<F>){
	chomp;
	push @samples, $_;
}
close(F);


open SNP, ">snps.pos.txt" or die;
open INDEL, ">indel.pos.txt" or die;

foreach my $sample (@samples){
	print "Processing $sample VCF\n";	
	my ($r1,$r2,$r3) = parse_vcf($sample,$base_dir,$vcfType);
	my %union_pos = %{$r1};
	my %unionINDEL_pos = %{$r2};
	my %tempRef = %{$r3};

	foreach my $chr (sort keys %union_pos){
		foreach my $pos ( sort {$a<=>$b} keys %{$union_pos{$chr}} ){
			$lhAll_pos{$chr}{$pos} = 1;
			$ref_calls{$chr}{$pos} = $tempRef{$chr}{$pos};
		}
	}
	foreach my $chr (sort keys %unionINDEL_pos){
		foreach my $pos ( sort {$a<=>$b} keys %{$unionINDEL_pos{$chr}}){
			$lhAllINDEL_pos{$chr}{$pos} = 1;
		}
	}
}	


print INDEL "chr\tpos\tref\n";
foreach my $chr (sort keys %lhAllINDEL_pos){
	foreach my $pos (sort {$a<=>$b} keys %{$lhAllINDEL_pos{$chr}}  ) {
		print INDEL "$chr\t$pos\tNA\n";
	}
}

open REF, ">chr.pos.ref.txt" or die;
print REF "chr\tpos\tref\n";
my %snpPos;
foreach my $chr (sort keys %lhAll_pos){
	foreach my $pos (sort {$a<=>$b} keys %{$lhAll_pos{$chr}}){
		print SNP "$chr\t$pos\n";
		print REF "$chr\t$pos\t$ref_calls{$chr}{$pos}\n";
	}
}

my $script_dir = abs_path($0);
$script_dir =~ s/filter_SNPs_MT_0.2.pl//;

#  Extracting coverage
`cat $_[0] | awk '{print \$1 " $base_dir $percentTotalcov $min_cov "}' | tr '\n' ' ' > xargs.txt`;
`cat $ARGV[1] | awk '{print \$1 " $base_dir $percentTotalcov $min_cov "}' | tr '\n' ' ' | xargs -n4 -P $threads $script_dir/extract_cov_0.2.pl`;
`cat *.sample.qual > sample.quals`;

my %popCov;
my $counter=0;
foreach my $sample (@samples){
	$counter = 0;
	open COV, "$sample.COV" or die;
	while(<COV>){
		chomp;
		$popCov{$counter} += $_;

		$counter++;
	}
	close(COV);
}

open COV,">temp.COV" or die;
print COV "Cov\n";
foreach my $counter ( sort {$a<=>$b} keys %popCov ) {

	print COV "$popCov{$counter}\n";

}
close(COV);

`paste chr.pos.ref.txt temp.COV > pop.COV`;

my @tempSamples = @samples;
my $tempNum = (scalar @samples)/1000;
$tempNum = ceil($tempNum);

for (my $i=1; $i<=$tempNum; $i++){
	open OUT,">samples.$i" or die;
	for (my $j=0; $j<=999; $j++){
		if (!@tempSamples){next;}
		my $tempSample = shift @tempSamples;
		print OUT "$tempSample.calls ";
	}
}

open TMP, ">tempMat.txt" or die;
for (my $i=1; $i<=$tempNum; $i++){
	print TMP "tempMat.$i ";
	`paste \`cat samples.$i\` > tempMat.$i`;
}
close TMP;
`paste \`cat tempMat.txt\` > tempMat.final`;
`paste chr.pos.ref.txt tempMat.final > unfiltered.mat`;

if(!-d"temp"){
	`mkdir temp`;
}
`mv xargs.txt *.COV *.qual *.calls temp`;
`mv temp/pop.COV .`;
`rm tempMat.* samples.*`;
### SUB END ###
}

sub sampleFilter{
### SUB START ###

my $missCutoff = $_[0];
my $mixCutoff = $_[1];
my $covCutoff = $_[2];
my %quals;
my @samples;
my %samples;
open F, "sample.quals" or die;
while(<F>){
	chomp;
	my ($sample,$miss,$mix,$cov) = split /\s+/,$_;	
	if ($miss<$missCutoff and $mix < $mixCutoff and $cov > $covCutoff){
		push @samples,$sample;
		$samples{$sample} = 1;
	}
}
close(F);

#progress bar
my $totSNPs = `wc -l < unfiltered.mat`;
chomp $totSNPs;
my $progress = Term::ProgressBar->new({name => 'Filtering', count => $totSNPs, remove => 1});
$progress->minor(0);
my $next_update = 0;

open BIN, ">sample.filt.mat.bin" or die;
open OUT, ">sample.filt.mat" or die;
my $x = 0;
my %keepSamples;
open F, "unfiltered.mat" or die;
while(<F>){
	chomp;
	if ($x<1){
		$x++;
		my $line = "chr\tpos\tref";
		my @a = split /\s+/,$_;
		@a = @a[3..$#a];
		  for (my $i=0; $i<=$#a; $i++){
			if (exists($samples{$a[$i]})){
				$keepSamples{$i} = $a[$i];
				$line .= "\t$a[$i]";
			}
		}
		print OUT "$line\n";
		print BIN "$line\n";
		next;
	}
	if (scalar keys %keepSamples ==0 ){ print "No samples left after filtering\n"; exit;}

	if ($_ =~ /\w+\t(\d+)\t\w+\t(.+)/){
		my $pos = $1;
		my $line = $2;
		my $a = () = $line =~ /A/g;
		my $c = () = $line =~ /C/g;
		my $g = () = $line =~ /G/g;
		my $t = () = $line =~ /T/g;
		my $count = () = "%$a%%$c%%$g%%$t%" =~ /%0%/g;
		if ($count == 3){
			$x++;
			next;            # this removes polmorhic snps
		}
	}

	

	my @a = split /\s+/,$_;
	my ($chr,$pos,$ref) = @a[0..2];
	@a = @a[3..$#a];
	my %alt_present;
	my $line = "$chr\t$pos\t$ref";
	my $binLine = "$chr\t$pos\t$ref";
	for (my $i=0; $i<=$#a; $i++){
		if (exists($keepSamples{$i})){
			$line .= "\t$a[$i]";
			
			if ($a[$i] eq "-" || $a[$i] eq "N"){
				$binLine .= "\tNA";
			} elsif ($a[$i] eq $ref){
				$binLine .= "\t0";
			} elsif ($a[$i] ne $ref){
				if (length($a[$i])==1){
					$binLine.= "\t1";
				} elsif (length($a[$i])>1){
					$binLine.= "\tMX";
				} else {
					die;
				}
			} else {
				die;
			}

			if ($a[$i] eq "N" or $a[$i] eq "-"){next;}
			$alt_present{$a[$i]} = 1;
		}
	}
	if (scalar keys %alt_present > 1){
		print OUT "$line\n";
		print BIN "$binLine\n";
	}

	$x++;
	if ($x >= $next_update){
		$next_update = $progress->update($x);
	}
}
if ($x >= $next_update){
	$next_update = $progress->update($x);
}
close(F);
close(OUT);
close(BIN);

### SUB END ###
}

sub snpQual{
	
my $snpMat = "sample.filt.mat.bin";
my @samples;

my $totSNPs = `wc -l < sample.filt.mat.bin`;
chomp $totSNPs;
my $progress = Term::ProgressBar->new({name => 'Filtering', count => $totSNPs, remove => 1});
$progress->minor(0);
my $next_update = 0;


my $x = 0;
open OUT, ">snps.qual" or die;
open F, $snpMat or die;
while(<F>){
	chomp;
	if ($x<1){
		$x++;
		my $line = "chr\tpos\tref";
		my @a = split /\s+/,$_;
		@a = @a[3..$#a];
		for (my $i=0; $i<=$#a; $i++){
			push @samples,$a[$i];
		}
		next;
	}
	my $line;
	my $numSamples = $#samples+1;
	if ($_ =~ /([^\s]+)\t(\d+)\t\w+\t(.+)/){
		my $chr = $1;
		my $pos = $2;
		$line = $3;
		my $zeros = () = $line =~ /0/g;
		my $ones = () = $line =~ /1/g;
		my $mix = () = $line =~ /MX/g;
		my $NA = () = $line =~ /NA/g;
		print OUT "$chr\t$pos\t",$zeros/$numSamples,"\t",$ones/$numSamples,"\t",$NA/$numSamples,"\t",$mix/$numSamples,"\n";
	}
	$x++;
	if ($x >= $next_update){
		$next_update = $progress->update($x);
	}
}
close(F);
close(OUT);

}

sub filterSNPs{


my $cutMiss = $_[0];
my $cutMix = $_[1];
my $mapFile = $_[2];

my $totSNPs = `wc -l < $mapFile`;
chomp $totSNPs;
my $progress = Term::ProgressBar->new({name => 'Loading Mappability', count => $totSNPs, remove => 1});
$progress->minor(0);
my $next_update = 0;

my $x = 0;
my %hasMappability;
open F, "$mapFile" or die;
while(<F>){
	chomp;
	my ($chr,$start,$end,$map) = (split /\s+/,$_)[0,1,2,4];
	for (my $i = $start+1; $i<=$end+1; $i++){
		$hasMappability{$chr}{$i} = $map;
	}
	$x++;
	if ($x >= $next_update){
		$next_update = $progress->update($x);
	}
}
close(F);

my %snpQual;
open F, "snps.qual" or die;
while(<F>){
	chomp;
	my ($chr,$pos,$missing,$mix) = (split /\s+/,$_)[0,1,4,5];
	$snpQual{$chr}{$pos}{"miss"} = $missing;
	$snpQual{$chr}{$pos}{"mix"} = $mix;
}
close(F);

$totSNPs = `wc -l < sample.filt.mat`;
chomp $totSNPs;
$progress = Term::ProgressBar->new({name => 'Filtering', count => $totSNPs, remove => 1});
$progress->minor(0);
$next_update = 0;


open OUT, ">snp.map.sample.filt.mat" or die;
open LOG, ">snps.log" or die;
$x = 0;
open F, "sample.filt.mat" or die;
while(<F>){
    if ($x<1){$x++; print OUT $_; next;}
	my ($chr,$pos) = (split /\s+/,$_)[0,1];
	if (!exists($hasMappability{$chr}{$pos})){print "$chr\t$pos\tNAAAAAAA\n";exit;}
	if ($hasMappability{$chr}{$pos} < 1 or $hasMappability{$chr}{$pos} == 0){
		print LOG "$chr\t$pos\t$hasMappability{$chr}{$pos}\t$snpQual{$chr}{$pos}{'miss'}\t$snpQual{$chr}{$pos}{'mix'}\tNON_UNIQ\n";
		$x++;
		next;
	}
	if (!defined($snpQual{$chr}{$pos}{"miss"} ) or !defined($snpQual{$chr}{$pos}{"miss"})){die"$chr $pos\n"};
	if ($snpQual{$chr}{$pos}{"miss"} > $cutMiss or $snpQual{$chr}{$pos}{"mix"} > $cutMix){
		print LOG "$chr\t$pos\t$hasMappability{$chr}{$pos}\t$snpQual{$chr}{$pos}{'miss'}\t$snpQual{$chr}{$pos}{'mix'}\tLOW_QUAL\n";
		$x++;
		next;
	} 
	print LOG "$chr\t$pos\t$hasMappability{$chr}{$pos}\t$snpQual{$chr}{$pos}{'miss'}\t$snpQual{$chr}{$pos}{'mix'}\tOK\n";
	print OUT  $_;
	$x++;
	if ($x >= $next_update){
		$next_update = $progress->update($x);
	}
}
 if ($x >= $next_update){
	$next_update = $progress->update($x);
}		


}








sub covFilter{


my $low = $_[0];
my $high = $_[1];
my $covFile = $_[2];

my %snpCov;
open F, "pop.COV" or die;
while(<F>){
	chomp;
	my ($chr,$pos,$cov) = (split /\s+/,$_)[0,1,3];
	$snpCov{$chr}{$pos} = $cov;
}
close(F);

my %nuclear;
open F, $covFile or die;
while(<F>){
	chomp;
	$nuclear{$_} = 1;
}
close(F);

my $totSNPs = `wc -l < snp.map.sample.filt.mat`;
chomp $totSNPs;
my $progress = Term::ProgressBar->new({name => 'Filtering', count => $totSNPs, remove => 1});
$progress->minor(0);
my $next_update = 0;


open OUT, ">cov.snp.map.sample.filt.mat" or die;
open LOG, ">cov.snps.log" or die;
print LOG "chr\tpos\tcov\n";
my $x = 0;
open F, "snp.map.sample.filt.mat" or die;
while(<F>){
    if ($x<1){$x++; print OUT $_; next;}
	my ($chr,$pos) = (split /\s+/,$_)[0,1];
	if (!exists($nuclear{$chr})){
		print LOG "$chr\t$pos\t$snpCov{$chr}{$pos}\tOK\n";
		print OUT  $_;
		$x++;
		next;
	}	
	if ($snpCov{$chr}{$pos} < $low or $snpCov{$chr}{$pos} > $high){
		print LOG "$chr\t$pos\t$snpCov{$chr}{$pos}\tBAD_COV\n";
		$x++;
		next;
	} 
	print LOG "$chr\t$pos\t$snpCov{$chr}{$pos}\tOK\n";
	print OUT  $_;
	$x++;
	if ($x >= $next_update){
		$next_update = $progress->update($x);
	}
}
if ($x >= $next_update){
                $next_update = $progress->update($x);
        }
	

}



################## SUBROUTINES #######################

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

sub parse_vcf{
	my $sample = $_[0];
	my $base_dir = $_[1];
	my $vcfType = $_[2];
	
	my $gatk_vcf = "$base_dir/vcf/$sample.gatk.raw.vcf.gz";
	my $samtools_vcf = "$base_dir/vcf/$sample.filt.vcf.gz";

	my %lhSamtools;
	my %lhSamtoolsINDEL;
	my %lhGatk;
	my %lhGatkINDEL;
	my %ref;


	if ($vcfType eq "both"){
		if (!-e $gatk_vcf || !-e $samtools_vcf){
			print "Cant't find VCF files for $sample\n";
			exit;
		}
		open V, "gunzip -c $samtools_vcf |";
		while (<V>){
			if ($_ =~ /#/){next;}
			my ($chr,$pos,$ref,$alt) = (split/\s+/, $_)[0,1,3,4];
			if (length($ref) == length($alt)){
				if (exists($lhSamtools{$chr}{$pos})){ 
					$lhSamtools{$chr}{$pos} = $lhSamtools{$chr}{$pos}."%".$alt;
				} else {
					$lhSamtools{$chr}{$pos} = $alt;
				}
				$ref{$chr}{$pos} = $ref;
			} else {
				if (exists($lhSamtoolsINDEL{$chr}{$pos})){
					$lhSamtoolsINDEL{$chr}{$pos} = $lhSamtoolsINDEL{$chr}{$pos}."%".$alt;
				} else {
					$lhSamtoolsINDEL{$chr}{$pos} = $alt;
				}
			}
		}
		close V;
	
		open V, "gunzip -c $gatk_vcf |";
		while (<V>){
			if ($_ =~ /#/){next;}
			my ($chr,$pos,$ref,$alt) = (split/\s+/, $_)[0,1,3,4];
			if (length($ref) == length($alt)){
				if (exists($lhGatk{$chr}{$pos})){ 
					$lhGatk{$chr}{$pos} = $lhGatk{$chr}{$pos}."%".$alt;
				} else {
					$lhGatk{$chr}{$pos} = $alt;
				}
				if (exists($ref{$chr}{$pos})){
					if ($ref{$chr}{$pos} ne $ref){ 
						print "Samtools and GATK disagree on reference allele at $pos in $sample\n";
						exit;
					}
				} else {
					$ref{$chr}{$pos} = $ref;
				}
			} else {
				if (exists($lhGatkINDEL{$chr}{$pos})){
					$lhGatkINDEL{$chr}{$pos} = $lhGatkINDEL{$chr}{$pos}."%".$alt;
				} else {
					$lhGatkINDEL{$chr}{$pos} = $alt;
				}
			}
		}
		close V;

		my %unionINDEL_pos;
		my %union_pos;
	
		
		foreach my $chr (sort keys %lhSamtools) {
			foreach my $pos ( sort {$a<=>$b} keys %{$lhSamtools{$chr}} ) {
				$union_pos{$chr}{$pos} = $lhSamtools{$chr}{$pos};
			}
		}
		foreach my $chr (sort keys %lhGatk){
			foreach my $pos ( sort {$a<=>$b} keys %{$lhGatk{$chr}} ) {
				$union_pos{$chr}{$pos} = $lhGatk{$chr}{$pos};
			}
		}
	
		foreach my $chr (sort keys %lhSamtoolsINDEL ) {
			foreach my $pos ( sort {$a<=>$b} keys %{$lhSamtoolsINDEL{$chr}} ) {
				$unionINDEL_pos{$chr}{$pos} = $lhSamtoolsINDEL{$chr}{$pos};
			}
		}
		foreach my $chr (sort keys %lhGatkINDEL ) {
			foreach my $pos ( sort {$a<=>$b} keys %{$lhGatkINDEL{$chr}} ) {
				$unionINDEL_pos{$chr}{$pos} = $lhGatkINDEL{$chr}{$pos};
			}
		}	
	
		return(\%union_pos,\%unionINDEL_pos,\%ref);	
		

	} elsif ($vcfType eq "samtools"){
		if (!-e $samtools_vcf){
			print "Cant't find $samtools_vcf\n";
			exit;
		}
		open V, "gunzip -c $samtools_vcf |";
		while (<V>){
			if ($_ =~ /#/){next;}
			my ($chr,$pos,$ref,$alt) = (split/\s+/, $_)[0,1,3,4];
			if (length($ref) == length($alt)){
				if (exists($lhSamtools{$chr}{$pos})){ 
					$lhSamtools{$chr}{$pos} = $lhSamtools{$chr}{$pos}."%".$alt;
				} else {
					$lhSamtools{$chr}{$pos} = $alt;
				}
				$ref{$chr}{$pos} = $ref;
			} else {
				if (exists($lhSamtoolsINDEL{$chr}{$pos})){
					$lhSamtoolsINDEL{$chr}{$pos} = $lhSamtoolsINDEL{$chr}{$pos}."%".$alt;
				} else {
					$lhSamtoolsINDEL{$chr}{$pos} = $alt;
				}
			}
		}
		close V;
			return(\%lhSamtools,\%lhSamtoolsINDEL,\%ref);
	} elsif ($vcfType eq "gatk"){
		if (!-e $gatk_vcf){
			print "Cant't find VCF files for $sample\n";
			exit;
		}
		open V, "gunzip -c $gatk_vcf |";
		while (<V>){
			if ($_ =~ /#/){next;}
			my ($chr,$pos,$ref,$alt) = (split/\s+/, $_)[0,1,3,4];
			if (length($ref) == length($alt)){
				if (exists($lhGatk{$chr}{$pos})){ 
					$lhGatk{$chr}{$pos} = $lhGatk{$chr}{$pos}."%".$alt;
				} else {
					$lhGatk{$chr}{$pos} = $alt;
				}
				if (exists($ref{$chr}{$pos})){
					if ($ref{$chr}{$pos} ne $ref){ 
						print "Samtools and GATK disagree on reference allele at $pos in $sample\n";
						exit;
					}
				} else {
					$ref{$chr}{$pos} = $ref;
				}
			} else {
				if (exists($lhGatkINDEL{$chr}{$pos})){
					$lhGatkINDEL{$chr}{$pos} = $lhGatkINDEL{$chr}{$pos}."%".$alt;
				} else {
					$lhGatkINDEL{$chr}{$pos} = $alt;
				}
			}
		}
		close V;
		return(\%lhGatk,\%lhGatkINDEL,\%ref);

	}


				
}

sub createIndelMat{
	my $scriptDir = abs_path($0);
	$scriptDir =~ s/filter_SNPs_MT_0.2.pl//;

	my $sampleFile = $_[0];
	my $base_dir = $_[1];
	my $threads = $_[2];

	my @samples;
	open OUT, ">xargs.indel.txt" or die; 
	open F, $sampleFile or die;
	while(<F>){
		chomp;
		push @samples,$_;
		print OUT "$_ $base_dir ";
	}
	close(F);
	close(OUT);
	
	`cat xargs.indel.txt | xargs -n2 -P $threads $scriptDir/makeIndel_calls.pl`;
	
	my @tempSamples = @samples;
	my $tempNum = (scalar @samples)/1000;
	$tempNum = ceil($tempNum);

	for (my $i=1; $i<=$tempNum; $i++){
	    open OUT,">samples.$i" or die;
	    for (my $j=0; $j<=999; $j++){
	        if (!@tempSamples){next;}
	        my $tempSample = shift @tempSamples;
	        print OUT "$tempSample.indel.calls ";
	    }
		close(OUT);
	}
	
	open TMP, ">tempMat.txt" or die;
	for (my $i=1; $i<=$tempNum; $i++){
	    print TMP "tempMat.$i ";
	    `paste \`cat samples.$i\` > tempMat.$i`;
	}
	close TMP;
	`paste \`cat tempMat.txt\` > tempMat.final.indel`;
	`paste indel.pos.txt tempMat.final.indel > unfiltered.indel.mat`;
	
	my $i = 0;
	open OUT, ">filtered.indel.mat" or die;
	open F, "unfiltered.indel.mat" or die;
	while(<F>){
		chomp;
		if ($i<1){print OUT "$_\n"; $i++; next;}
		my $line = $_;
		$_ =~ s/\w+\t\d+\tNA\t//;
		my $alts = () = $_ =~ /1/g;
		if ($alts > 0 ){
			print OUT "$line\n";
		}
	}
	close(F);
	close(OUT);
}



sub generateMapFile{

my $ref = $_[0];
my $kmer = $_[1];
my $threads = $_[2];

my ($r1,$r2) = hash_contigs($ref);
my %hasContigs = %{$r1};
my @arrContigs = @{$r2};

open OUT, ">processed.ref.fa" or die;
for (sort @arrContigs){
	print OUT ">$_\n$hasContigs{$_}\n";
}
close(OUT);

`gem-indexer -i processed.ref.fa -o index -T $threads`;
`gem-mappability -T $threads -I index.gem -l $kmer -m 0.00 -o genome`;
`gem-2-wig -I index.gem -i genome.mappability -o genome`;
`wig2bed < genome.wig > mappability.temp.bed`;

open F, "mappability.temp.bed" or die;
open OUT, ">mappability.bed" or die;
my @lines = <F>;
foreach my $chr (@arrContigs){
	for (my $i = 0; $i<=$#lines;$i++){
		if ($lines[$i] =~ /$chr/){
			if (defined($lines[$i+1])){
				if ($lines[$i+1] =~ /$chr/){
					print OUT $lines[$i];
				} else {
					my $end = (split /\s+/,$lines[$i])[2];
					print OUT $lines[$i];
					my $chrEnd = length($hasContigs{$chr});
					print OUT "$chr\t$end\t$chrEnd\t.\t1.000000\n";
				}
			} else {
				my $end = (split /\s+/,$lines[$i])[2];
				print OUT $lines[$i];
				my $chrEnd = length($hasContigs{$chr});
				print OUT "$chr\t$end\t$chrEnd\t.\t1.000000\n";
			}
		}
	}
}
close(OUT);
}




sub samplePlot{


print "Please select cutoffs for mixed missing and median coverage per sample:\n";
my $x = `Rscript samplePlot.R`;
chomp $x;
my ($l1,$l2,$l3) = (split /\s+/,$x)[1,2,3];
print "Mixed call fraction:$l1\n";
print "Missing call fraction:$l2\n";
print "Median coverage:$l3\n\n";
print "Do you want to use these values? [Y/n]";
my $ans1 = <STDIN>;
chomp $ans1;
$ans1 = uc($ans1);

while ($ans1 ne "" && $ans1 ne "Y" && $ans1 ne "N"){
	print "Do you want to use these values? [Y/n]";
	$ans1 = <STDIN>;
	chomp $ans1;
	$ans1 = uc($ans1);
}
if ($ans1 eq "" or $ans1 eq "Y"){
	print "Using previously defined values\n";
} else {
	print "Please provide new missiness cutoff:\n";
	$l1 = <STDIN>;
	chomp $l1;
	print "Please provide new mixedness cutoff:\n";
	$l2 = <STDIN>;
	chomp $l2;
	print "Please provide new median coverage cutoff:\n";
        $l3 = <STDIN>;
        chomp $l3;
}

print "\n-----------------------------------------------------------\n";
print "Filtering snps using values $l1, $l2 and $l3\n";
print "-----------------------------------------------------------\n";
return($l1,$l2,$l3);
	
}


sub snpPlot{


print "Please select cutoffs for missing and mexed fraction per sample:\n";
my $x = `Rscript snpPlot.R`;
chomp $x;
my ($l1,$l2) = (split /\s+/,$x)[1,2];
print "Mixed call fraction:$l1\n";
print "Missing call fraction:$l2\n";
print "Do you want to use these values? [Y/n]";
my $ans1 = <STDIN>;
chomp $ans1;
$ans1 = uc($ans1);

while ($ans1 ne "" && $ans1 ne "Y" && $ans1 ne "N"){
	print "Do you want to use these values? [Y/n]";
	$ans1 = <STDIN>;
	chomp $ans1;
	$ans1 = uc($ans1);
}
if ($ans1 eq "" or $ans1 eq "Y"){
	print "Using previously defined values\n";
} else {
	print "Please provide new missiness cutoff:\n";
	$l1 = <STDIN>;
	chomp $l1;
	print "Please provide new mixedness cutoff:\n";
	$l2 = <STDIN>;
	chomp $l2;

}

print "\n-----------------------------------------------------------\n";
print "Filtering snps using values $l1 and $l2\n";
print "-----------------------------------------------------------\n";
return($l1,$l2);
	
}

sub covPlot{

my %okSNPs;
open F, "snps.log" or die;
while(<F>){
	chomp;
	if ($_ =~ /OK/){
		my ($chr,$pos) = (split /\s+/,$_)[0,1];
		$okSNPs{$chr}{$pos} = 1;
	}
}
close(F);


open F, "pop.COV" or die;
open OUT, ">filtered.cov" or die;
while(<F>){
	chomp;
	my ($chr,$pos) = (split /\s+/,$_)[0,1];
	if (exists($okSNPs{$chr}{$pos})){
		print OUT "$_\n";
	}
}
close(F);

print "Please select cutoffs based on total coverage per SNP:\n";

my $x = `Rscript covPlot.R`;
chomp $x;
my ($l1,$l2) = (split /\s+/,$x)[1,2];
print "Lower-bound:$l1\n";
print "Higher-bound:$l2\n";

print "Do you want to use these values? [Y/n]";
my $ans1 = <STDIN>;
chomp $ans1;
$ans1 = uc($ans1);

while ($ans1 ne "" && $ans1 ne "Y" && $ans1 ne "N"){
    print "|$ans1|\n";
    print "Do you want to use these values? [Y/n]";
    $ans1 = <STDIN>;
    chomp $ans1;
    $ans1 = uc($ans1);
}

if ($ans1 eq "" or $ans1 eq "Y"){
    print "Using previously defined values\n";
} else {
    print "Please provide new lower bound:\n";
    $l1 = <STDIN>;
    chomp $l1;
    print "Please provide new higher bound:\n";
    $l2 = <STDIN>;
    chomp $l2;
}


print "\n-----------------------------------------------------------\n";
print "Filtering snps based on $l1 lower bound and $l2 upper bound\n";
print "-----------------------------------------------------------\n";

return($l1,$l2);

}

sub writeRfiles{

open OUT, ">covPlot.R" or die;
print OUT 'x<-read.table("filtered.cov")
qs<-quantile(x$V4,probs=c(0,0.995))
x11()
plot(density(x$V4),xlim=qs,xlab="Coverage",pch=21,main="Density of total coverage over all SNPs")
y<-locator(2)
pdf("covPlot.pdf")
plot(density(x$V4),xlim=qs,xlab="Coverage",pch=21,main="Density of total coverage over all SNPs")
abline(v=y$x,col="red")
err<-dev.off()
print(sort(y$x))';
close(OUT);

open OUT, ">samplePlot.R" or die;
print OUT 'x<-read.table("sample.quals")
x11()
plot(sort(x$V2),main="Proportion of missing calls over total SNP positions",ylab="Missing proportion")
y<-locator(1)
plot(sort(x$V3),main="Proportion of mixed calls over total SNP positions",ylab="Mixed proportion")
z<-locator(1)
plot(sort(x$V4),main="Median coverage of samples",ylab = "Median coverage")
a<-locator(1)
pdf("sample-missing.pdf")
plot(sort(x$V2),main="Proportion of missing calls over total SNP positions",ylab="Missing proportion")
abline(h=y$y,col="red",lty=2)
abline(v=y$x,col="blue",lty=2)
err<-dev.off()
pdf("sample-mixed.pdf")
plot(sort(x$V3),main="Proportion of mixed calls over total SNP positions",ylab="Mixed proportion")
abline(h=z$y,col="red",lty=2)
abline(v=z$x,col="blue",lty=2)
err<-dev.off()
pdf("sample-cov.pdf")
plot(sort(x$V4),main="Median coverage of samples",ylab = "Median coverage")
abline(h=a$y,col="red",lty=2)
abline(v=a$x,col="blue",lty=2)
err<-dev.off()
print(c(y$y,z$y,a$y))';
close(OUT);


open OUT, ">snpPlot.R" or die;
print OUT 'x<-read.table("snps.qual")
x11()
plot(sort(x$V5),main="Proportion of missing calls over each SNP",ylab="Missing proportion")
y<-locator(1)
plot(sort(x$V6),main="Proportion of mixed calls over each SNP",ylab="Mixed proportion")
z<-locator(1)

pdf("snp-missing.pdf")
plot(sort(x$V5),main="Proportion of missing calls over each SNP",ylab="Missing proportion")
abline(h=y$y,col="red",lty=2)
abline(v=y$x,col="blue",lty=2)
err<-dev.off()
pdf("snp-mixed.pdf")
plot(sort(x$V6),main="Proportion of mixed calls over each SNP",ylab="Mixed proportion")
abline(h=z$y,col="red",lty=2)
abline(v=z$x,col="blue",lty=2)
err<-dev.off()

print(c(y$y,z$y))';
close(OUT);
}

sub fileSummary{

open OUT, ">summary.txt" or die;

my $unfilt =  `wc -l unfiltered.mat | awk '{print \$1}'`;
chomp $unfilt;
my $sampleFilt = `wc -l sample.filt.mat | awk '{print \$1}'`;
chomp $sampleFilt;
my $snpFilt = `wc -l snp.map.sample.filt.mat | awk '{print \$1}'`;
chomp $snpFilt;
my $covFilt = `wc -l cov.snp.map.sample.filt.mat | awk '{print \$1}'`;
chomp $covFilt;

print OUT "Matrix filtering:\n";
print OUT "Untiltered\t$unfilt\nSample_filtered\t$sampleFilt\nSNP_filtered\t$snpFilt\nCoverage_filtered\t$covFilt\n";

my $unfiltNoSamp = `head -1 unfiltered.mat | awk '{print NF-3}'`;
chomp $unfiltNoSamp;
my $filtNoSamp = `head -1 sample.filt.mat | awk '{print NF-3}'`;
chomp $filtNoSamp;
print OUT "\nSample Filtering:\nPre\t$unfiltNoSamp\nPost\t$filtNoSamp\n";

my $lowQual = `grep LOW_QUAL snps.log | wc -l`;
chomp $lowQual;
my $nonUniq = `grep NON_UNIQ snps.log | wc -l`;
chomp $nonUniq;
my $ok = `grep OK snps.log | wc -l`;
chomp $ok;
print OUT "\nSNP Qual Filtering:\nLow-Qual\t$lowQual\nNon-Uniq\t$nonUniq\nOK\t$ok\n";

close(OUT);
}
