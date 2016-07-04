#!/usr/bin/perl 
#===============================================================================
#
#		 FILE:  verifyIndels.pl
#
#		USAGE:  ./verifyIndels.pl  
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
#	  CREATED:  06/23/2016 03:14:31 PM
#	 REVISION:  ---
#===============================================================================

use strict;
use warnings;
use Cwd 'abs_path';
use Cwd;
use Statistics::Lite qw(:all);
use POSIX;

if ($#ARGV+1 < 4){print "\nverifyIndels.pl <sample> <base_dir> <ref> <careful|fast>\n\n"; exit;}

my $minKmer = 19;
my $maxKmer = 79;
my @kmers;
my @expCov;
my @covCut;
my $minIndelSize  = 100;
my $minFlankSize = 20;
my $minID = 95;
my $ref = $ARGV[2];
$ref = abs_path($ref);
my $sample = $ARGV[0];
my $base_dir = $ARGV[1];
my $alg = $ARGV[3];

if ($alg ne "careful" and $alg ne "fast"){
	print "Algorithm must be careful or fast\n";
	exit;
}

$base_dir = abs_path($base_dir);
my $script_dir = abs_path($0);
$script_dir =~ s/verifyIndels.pl//;
my $samtools = "$script_dir/samtools-1.3.1/samtools";
my $velvetOpt = "/usr/local/src/VelvetOptimiser-2.2.5/VelvetOptimiser.pl";

my $r1 = parseDelly($sample,$base_dir); 
my %delly = %{$r1};
my @dellyIndelNo = sort {$a<=>$b} keys %delly;
`mkdir $sample`; 
chdir("$sample");

our ($expCov,$covCut,$kmer) = calibrateAssembly($base_dir,$sample,"Chromosome",20000,30000);
print "Expected coverage:$expCov\tCoverage cutoff:$covCut\tKmer:$kmer\n";

my $pwd = `pwd`;
print $pwd;

open KMER,">kmer.txt" or die;
open EXPCOV, ">expCov.txt" or die;
open COVCUT, ">covCut.txt" or die;
open RESULTS, ">sv.results.txt" or die;
open POS, "../indelGenes.positions.txt" or die;
while(<POS>){
	
	print "In kmers:";
	for (@kmers){
		print "\t$_";
	}
	print "\n";

	print KMER "$kmer\n";
	print EXPCOV "$expCov\n";
	print COVCUT "$covCut\n";
	chomp;
	my ($name,$chr,$start,$end,$minStart,$maxEnd) = (split /\s+/,$_);
	my $regionStart = $start - 2000;
	my $regionEnd = $end + 2000;

	my $dellyPresent = 0;
	for (@dellyIndelNo){
#		print "$delly{$_}{'start'} < $start && $delly{$_}{'end'} > $end\n";
#		if ($delly{$_}{'start'} < $start && $delly{$_}{'end'} > $end){
#			print "Reassigning coordinates based on sample specific delly call:\n\t$start => $delly{$_}{'start'}\n\t$end => $delly{$_}{'end'}\n";
#			$start = $delly{$_}{'start'};
#			$end = $delly{$_}{'end'};
#			$dellyPresent ++;
#			last;
#		} 
	}	
	if ($dellyPresent == 0){
		print "Reassigning coordinates based on maximum possible delly call:\n\t$start => $minStart\n\t$end => $maxEnd\n";
		$start = $minStart;
		$end = $maxEnd;
	}

	print RESULTS "$name\t$chr\t$start\t$end";
	print "$name\n";
	`mkdir \'$name\'`;
	chdir("$name");
	my $LA = localAssembly($base_dir,$sample,$chr,$regionStart,$regionEnd);
	if ($LA eq "FAIL"){
		print RESULTS "\tNA\n";
		chdir("../");
		next;
	}
	my $res = bwaAlign($ref,$start,$end);
	if ($res eq "NA"){
		if ($alg eq "careful"){
			print "Not possible to determine...Trying full assembly\n";
			fullLocalAssembly($base_dir,$sample,$chr,$regionStart,$regionEnd);
			$res = bwaAlign($ref,$start,$end);
		} else {
			print RESULTS "\tABSENT\n";
			print "Indel absent\n";
			`echo "NO INDELS" > results.txt`;
			sam2graph($chr,$regionStart,$regionEnd,$start,$end);
			chdir("../");
			next;
		}
	}
	if ($res eq "ABSENT"){
		print RESULTS "\tABSENT\n";
		print "Indel absent\n";
		`echo "NO INDELS" > results.txt`
	} elsif ($res eq "NA"){
		print "Not possible to determine\n";
		print RESULTS "\tNA\n";
		`echo "NO CONTIGS COVERING" > results.txt`
	} else {	
		if ($alg eq "careful"){
			print "Possible indel...runing AGE\n";
			my $ageResult = ageAlign($ref,$chr,$regionStart,$regionEnd,$minFlankSize,$minID);
			if ($ageResult eq "NA"){
				print RESULTS "\tMAYBE\n";
			} elsif ($ageResult eq "PRESENT"){
				print RESULTS "\tPRESENT\n";
			} else {
				die;
			}
		} else {
			print RESULTS "\tPRESENT\n";
		}
	}

	sam2graph($chr,$regionStart,$regionEnd,$start,$end);	
	chdir("../");
}
close(POS);
close RESULTS;
close KMER;
close EXPCOV;
close COVCUT;

#---------------------------------------------------------------------------
#  Sub Routines
#---------------------------------------------------------------------------

sub calibrateAssembly{


	my ($base_dir,$sample,$chr,$start,$end) = @_;
	`mkdir calibration`;
	chdir("calibration");

	my $numReads = `sambamba view -F "not (unmapped or mate_is_unmapped) and mapping_quality >=30" $base_dir/bam/$sample.bam $chr:$start-$end | wc -l`;
	chomp $numReads;
	if ($numReads < 1000){
	    die"Calibration failed";
	}

	`sambamba view -F "not (unmapped or mate_is_unmapped) and mapping_quality >=30" -o filt.bam $base_dir/bam/$sample.bam $chr:$start-$end -f bam`;
	
	`velveth test $minKmer,$maxKmer,2 -shortPaired -bam filt.bam`;
	
	
	my %assembly;
	my %expCov;
	my %covCut;
	for (my $i=$minKmer; $i<$maxKmer; $i=$i+2){
		my $res = `velvetg test_$i -cov_cutoff auto -exp_cov auto -clean yes| tail -3 | tr '\n' ' '`;
		if ($res =~ /EMPTY/){
			next;
		}
		$res =~ m/Estimated Coverage = ([\d\.]+).+Estimated Coverage cutoff = ([\d\.]+).+n50 of (\d+)/;
		my ($exp_cov,$cov_cut,$n50) = ($1,$2,$3);
		$assembly{$n50} = $i;
		$expCov{$i} = $exp_cov;
		$covCut{$i} = $cov_cut;
		print "$i\t$n50\t$exp_cov\t$cov_cut\n";
	} 
	
	
	my $bestn50 = (sort {$b<=>$a} keys %assembly)[0];
	my $best = $assembly{$bestn50};
	`mv test_$best k$best`;
	#`rm -r test*`;
	`ln -s k$best/contigs.fa .`;

	my $exp_cov = $expCov{$best};
	my $cov_cut = $covCut{$best};
	
	
	push @kmers,$best;
	push @expCov,$exp_cov;
	push @covCut,$cov_cut;	
	
	my $newCovCut = mean @covCut;
	my $newExpCov = mean @expCov;
	my $newKmer = mean @kmers;
	
	$kmer = ceil $newKmer;
	$expCov = ceil $newExpCov;
	$covCut = $newCovCut;
	chdir("../");
	return($exp_cov,$cov_cut,$best);	

}

sub localAssembly{

#### start #####
if (-e "contigs.fa"){
	print "Contigs already exists...Exiting\n"; exit;
}


my ($base_dir,$sample,$chr,$start,$end) = @_;


my $numReads = `sambamba view -F "not (unmapped or mate_is_unmapped) and mapping_quality >=30" $base_dir/bam/$sample.bam $chr:$start-$end | wc -l`;
chomp $numReads;
if ($numReads < 1000){
	return("FAIL");
}
`sambamba view -F "not (unmapped or mate_is_unmapped) and mapping_quality >=30" -o filt.bam $base_dir/bam/$sample.bam $chr:$start-$end -f bam`;
`velveth k$kmer $kmer -shortPaired -bam filt.bam 2>> err; velvetg k$kmer -exp_cov $expCov -cov_cutoff $covCut -very_clean yes 2>> err`;
my $folder = `ls | grep ^k`;
chomp $folder;
if (!-e "$folder/contigs.fa"){
	print "Local assembly failed...Exiting\n"; exit;
}
`ln -s $folder/contigs.fa contigs.fa`;
return("TRUE");
#### end ####
}

sub fullLocalAssembly{

if (-e "contigs.fa"){
    `rm -r *`;
    print "Revmoing old run\n";
}

my ($base_dir,$sample,$chr,$start,$end) = @_;
`sambamba view -F "not (unmapped or mate_is_unmapped) and mapping_quality >=30" -o filt.bam $base_dir/bam/$sample.bam $chr:$start-$end -f bam`;


`velveth test $minKmer,$maxKmer,2 -shortPaired -bam filt.bam`;
	
	
	my %assembly;
	my %expCov;
	my %covCut;
	for (my $i=$minKmer; $i<$maxKmer; $i=$i+2){
		my $res = `velvetg test_$i -cov_cutoff auto -exp_cov auto -clean yes| tail -3 | tr '\n' ' '`;
		$res =~ m/Estimated Coverage = ([\d\.]+).+Estimated Coverage cutoff = ([\d\.]+).+n50 of (\d+)/;
		my ($exp_cov,$cov_cut,$n50) = ($1,$2,$3);
		$assembly{$n50} = $i;
		$expCov{$i} = $exp_cov;
		$covCut{$i} = $cov_cut;
		print "$i\t$n50\t$exp_cov\t$cov_cut\n";
	} 
	
	
	my $bestn50 = (sort {$b<=>$a} keys %assembly)[0];
	my $best = $assembly{$bestn50};
	`mv test_$best k$best`;
	#`rm -r test*`;
	`ln -s k$best/contigs.fa .`;

	my $exp_cov = $expCov{$best};
	my $cov_cut = $covCut{$best};
	
	
	push @kmers,$best;
	push @expCov,$exp_cov;
	push @covCut,$cov_cut;	
	
	my $newCovCut = mean @covCut;
	my $newExpCov = mean @expCov;
	my $newKmer = mean @kmers;
	print "Reassigning kmer $kmer => ";
	$kmer = ceil $newKmer;
	$expCov = ceil $newExpCov;
	$covCut = $newCovCut;
	print "$kmer\n";




}

sub OPTIMlocalAssembly{

#### start #####
if (-e "contigs.fa"){
	`rm -r *`;
	print "Revmoing old run\n";
}


my $velvetOpt = "/usr/local/src/VelvetOptimiser-2.2.5/VelvetOptimiser.pl";
my ($base_dir,$sample,$chr,$start,$end) = @_;
`sambamba view -F "not (unmapped or mate_is_unmapped) and mapping_quality >=30" -o filt.bam $base_dir/bam/$sample.bam $chr:$start-$end -f bam`;
`$velvetOpt --s $minKmer --e $maxKmer --x 2 -f '-shortPaired -bam filt.bam' 2>>err`;
my $folder = `ls | grep auto`;
chomp $folder;
if (!-e "$folder/contigs.fa"){
	print "Local assembly failed...Exiting\n"; exit;
}
`ln -s $folder/contigs.fa contigs.fa`;

my $exp_cov = `tail -n20 *Log* | grep Velvetg | awk '{print \$8}'`;
my $cov_cut = `tail -n20 *Log* | grep Velvetg |awk  '{print \$10}'`;
my $tempkmer = `ls -d auto_data*`;
chomp $tempkmer;
$tempkmer =~ s/auto_data_//;
chomp $exp_cov;
chomp $cov_cut;
push @kmers,$tempkmer;
push @expCov,$exp_cov;
push @covCut,$cov_cut;	

my $newCovCut = mean @covCut;
my $newExpCov = mean @expCov;
my $newKmer = mean @kmers;

print "Reassigning kmer $kmer => ";
$kmer = ceil $newKmer;
$expCov = ceil $newExpCov;
$covCut = $newCovCut;
print "$kmer\n";

#### end ####
}


sub bwaAlign{
##### start #####
my ($ref,$bp1,$bp2) = @_;
`bwa bwasw $ref contigs.fa 2>> err | $samtools view -S - > bwaContigAln.sam 2>> err`;


my ($c1,$c2) = hash_contigs("contigs.fa");
my %hasContigs = %{$c1};
my @arrContigs = @{$c2};

my %hits;
open F, "bwaContigAln.sam" or die;
my $i = 0;
while (<F>){
	chomp;
	my ($contig) = (split /\s+/,$_)[0];
	$hits{$contig} ++;
#   print "Segment found:\t$contig\n";
}
close F;

my @hits;
for( keys %hits){
	if ($hits{$_} >1){
		push @hits,$_;
#		print "Analysing $_ further\n";
	}
}

my $result = checkContigCov($bp1,$bp2);
if ($result eq "ABSENT"){
	return("ABSENT");
} elsif ($result eq "NA"){
	return("NA");
} else {

	if (scalar @hits < 1){
		die"No split contig\n";
	} else {
		open OUT, ">hits.fasta" or die;
		for (@hits){
			print OUT ">$_\n$hasContigs{$_}\n";
		}
	}
	close OUT;
	return("POSSIBLE");
}



######### end ###########
}













sub ageAlign{
##### start #####

my ($ref,$chr,$start,$end,$minFlankSizem,$minID) = @_;
my ($r1,$r2) = hash_contigs("$ref");
my %hasRef = %{$r1};
my @arrRef = @{$r2};

open REF, ">refChr.fa" or die;
print REF ">$chr\n";
my $refSeq = $hasRef{$chr} ;
$refSeq =~ s/(.{1,60})/$1\n/gs;
print REF "$refSeq\n";
close(REF);

my $res;
my $hits = `age_align -indel -both refChr.fa hits.fasta -coor1=$start-$end`;
`age_align -indel -both refChr.fa hits.fasta -coor1=$start-$end>report.txt`;
my @agehits = split /Alignment time is/,$hits;
for (my $i=0; $i<$#agehits; $i++){
	my $hit = $agehits[$i];
	my ($hitName,$identity,$indelsize,$bp1,$bp2,$flank1start,$flank1end,$flank2start,$flank2end) = parseAge($hit);
	if ($hitName eq "INSERTION"){
		print "FAIL assembly quality\n";
		next;
	}
#	print "$hitName\t$identity\t$indelsize\t$bp1-$bp2\t$flank1start\t$flank1end\t$flank2start\t$flank2end\n";
	my $flank1Len = abs($flank1end-$flank1start);
	my $flank2Len = abs($flank2end-$flank2start);
	if ($flank1Len < $minFlankSize or $flank2Len <$minFlankSize or $identity <$minID){
		print "FAIL length or identity\n"
	} else {
		print "PASS\n";
		$res .= "$hitName\t$indelsize\t$flank1Len\t$flank2Len\t$bp1\t$bp2\n";
	}
	
}

open RES, ">results.txt" or die;
if (!$res){
	print RES "NO INDELS\n";
	return("NA");
} else {
	print RES $res;
	return("PRESENT");

}
close RES;
##### end #####
}












sub parseAge {
### start###

my $hit = $_[0];


if ($hit !~ /Second.+\'(.+)\'/){
#	print "Hit name not found: Alignment malformed\n";
	exit;
}
$hit =~ m/Second.+\'(.+)\'/;
my $hitName = $1;


if ($hit !~ /Identic:\s+\d+\s+\(\s*(\d+)\%/){
	print "Indentity not found: Alignment malformed\n";
	exit;
}
$hit =~ m/Identic:\s+\d+\s+\(\s*(\d+)\%/;
my $identity = $1;


if ($hit !~ /EXCISED REGION\(S\):\n first  seq =>\s+(\d+)\s+nucs\s+\[(\d+),(\d+)\]/){
#	print "Indel size and breakpoints not found: Alignment malformed\n";
#	exit;
	return("INSERTION");
}
$hit =~ m/EXCISED REGION\(S\):\n first  seq =>\s+(\d+)\s+nucs\s+\[(\d+),(\d+)\]/;
my ($indelsize,$bp1,$bp2) = ($1,$2,$3);


if ($hit !~ / second seq =>\s+\[\s*(\d+),\s*(\d+)\]\s+EXCISED\s+REGION\s+\[\s*(\d+),\s*(\d+)\]/){
	print "Flanking regions not found: Alignment malformed\n";
	exit;
}
$hit =~ m/ second seq =>\s+\[\s+(\d+),\s+(\d+)\]\s+EXCISED\s+REGION\s+\[\s+(\d+),\s+(\d+)\]/;
my ($flank1start,$flank1end,$flank2start,$flank2end) = ($1,$2,$3,$4);

return ($hitName,$identity,$indelsize,$bp1,$bp2,$flank1start,$flank1end,$flank2start,$flank2end);
###end###
}







sub checkContigCov{
##### start ######

my ($bp1,$bp2) = @_;
my %contigs;
open F, "bwaContigAln.sam" or die;
while(<F>){
	chomp;
	my ($contig,$chr,$start,$cigar) = (split /\s+/,$_)[0,2,3,5];
	if ($cigar eq "*"){next;}
	if ($cigar =~ /S/){
		$cigar =~ s/\d+S//g;
	}
	if ($cigar =~ /I/){
		$cigar =~ s/I/M/g;
	}
	if ($cigar =~ /D/){
		$cigar =~ s/D/M/g;
	}

	my @m = split /M/, $cigar;
	my $m;
	for (@m){
		$m += $_;
	}
	my $end = $start+$m;
#	print "$contig\t$chr\t$start\t$m\n";
	for (my $i = $start; $i<=$end; $i++){
		$contigs{$contig}{$i}++;
	}
}
close(F);


foreach my $contig (keys %contigs){
	my $delLen = $bp2-$bp1;
	my $delCover=0;
	for (my $i=$bp1; $i<$bp2; $i++){
		if (exists($contigs{$contig}{$i})){
			$delCover ++;
		}
	}
#	print "\n\n$contig\t$delLen\t$delCover\n";
	if ($delLen eq $delCover){return("ABSENT");}
	
	my $flank1Start = $bp1 - 50;
	my $flank1End = $bp1;
	my $flank2Start = $bp2;
	my $flank2End = $bp2 + 50;
	
	my $flank1Cover=0;
	my $flank2Cover=0;
	
	for (my $i=$flank1Start; $i<$flank1End; $i++){
		if (exists($contigs{$contig}{$i})){
			$flank1Cover ++;
		}
	}
	for (my $i=$flank2Start; $i<$flank2End; $i++){
		if (exists($contigs{$contig}{$i})){
			$flank2Cover ++;
		}
	}
	my $flank1CoverPCT = $flank1Cover/50;
	my $flank2CoverPCT = $flank2Cover/50;
#	print "$contig\t$flank1CoverPCT\t$flank2CoverPCT\n";
	if ($flank1CoverPCT>0.50 && $flank2CoverPCT>0.50){return("POSSIBLE")}
}
return("NA");



########## end ############
}







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


















sub sam2graph{

my ($regionChr,$regionStart,$regionEnd,$bp1,$bp2) = @_;
`echo "$bp1 $bp2" > delly.coords.txt`;
open F, "bwaContigAln.sam" or die;
open OUT, ">samGraph.txt" or die;
$regionStart -= 1000;
$regionEnd += 1000;
while(<F>){
	chomp;
	my ($contig,$chr,$start,$cigar) = (split /\s+/,$_)[0,2,3,5];
	if ($regionChr ne $chr){next;}
	if ($cigar eq "*"){next;}
	if ($cigar =~ /S/){
		$cigar =~ s/\d+S//g;
	}
	if ($cigar =~ /I/){
		$cigar =~ s/I/M/g;
	}
	if ($cigar =~ /D/){
		$cigar =~ s/D/M/g;
	}

	my @m = split /M/, $cigar;
	my $m;
	for (@m){
		$m += $_;
	}
	my $end = $start+$m;
	if ($start<$regionStart or $end>$regionEnd){next;}
	print OUT "$contig\t$start\t$end\n";
}
close(F);

plotResults();
#### close ####
}





sub plotResults{
#### start #####

`$samtools view filt.bam| cut -f4,9 > insertSizes.txt`;


open R, ">sam2graph.r" or die;
print R '
contigs<-read.table("samGraph.txt",stringsAsFactors=F)
xmax<-max(c(contigs$V2,contigs$V3))
xmin<-min(c(contigs$V2,contigs$V3))

ymax<-dim(contigs)[1]
ymin<-0-ymax;
png("results.png", 800,600)
plot (xmax,ymax,xlim=c(xmin,xmax),ylim=c(ymin,ymax),pch="",yaxt="n")

uniq.contigs<-unique(contigs$V1)

for (i in 1:dim(contigs)[1]){
		height<-match(contigs$V1[i],uniq.contigs)
		rect(contigs$V2[i],height-1,contigs$V3[i],height,col="light grey")

}


yaxis.points<-seq(1:length(uniq.contigs))-0.5
yaxis.names<-substr(uniq.contigs,1,7)
axis(2,yaxis.points,yaxis.names)

segments(xmin,0,xmax,0,lwd=5)


ins<-read.table("insertSizes.txt")
ins$V2<-abs(ins$V2)
if(length(which(ins$V2>10000))>0){
	ins<-ins[-(which(ins$V2>10000)),]
}

par(new=T)
ymax<-max(ins$V2)*2.5

plot (xmax,ymax,xlim=c(xmin,xmax),ylim=c(ymin,ymax),pch="",yaxt="n")
points(ins$V1,ins$V2,pch=20)
axis(2,seq(0,max(ins$V2),500),seq(0,max(ins$V2),500))

del<-read.table("delly.coords.txt")
rect(del$V1,ymin,del$V2,ymax,col=paste(substr("#00FF00FF",1,7),40,sep=""),border=F)

res<-read.table("results.txt")
if (res$V1[1] != "NO"){
		for (i in 1:dim(res)[1]){
				rect(res$V5[i],ymin,res$V6[i],ymax,col=paste(substr("#FF00FFFF",1,7),40,sep=""),border=F)
		}
}


dev.off()
';
close R;
`Rscript sam2graph.r 2>> err`;

#### end ####
}










sub parseDelly{
#### start ####

print "Parsing $_[0] indels\n";
	my $sample = $_[0];
	my $base_dir = $_[1];

	my %delly;
	my $i=0;
	open F, "$base_dir/vcf/$sample.vcf" or die "Can't find $sample.vcf\n";
	while(<F>){
		chomp;
		if ($_ =~ /^#/){next;}
		my ($chr,$start,$qual,$temp) = (split /\s+/,$_)[0,1,6,7];
		if ($qual ne "PASS"){next;}
		$temp =~ m/END=(\d+)/;
		my $end = $1;
		$delly{$i}{'chr'} = $chr;
		$delly{$i}{'start'} = $start;
		$delly{$i}{'end'} = $end;

		$i++;
	}
	close(F);
	return(\%delly);
}


