#!/usr/bin/perl 
#===============================================================================
#
#		 FILE:  filterIndels.pl
#
#		USAGE:  ./filterIndels.pl  
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
#	  CREATED:  06/23/2016 12:06:52 PM
#	 REVISION:  ---
#===============================================================================

use strict;
use warnings;
use List::Util qw(max min);
use Cwd 'abs_path';



if ($#ARGV+1 < 4){print "\nfilterIndels.pl <sample file> <base_dir> <reference> <outfile>\n\n"; exit;}

my $sampleFile = $ARGV[0];
my $baseDir = $ARGV[1];
my $refFile = $ARGV[2];
my $outfile = $ARGV[3];
$baseDir = abs_path($baseDir);
#parseVCF($sampleFile,$baseDir);

my $err = `cat $sampleFile | xargs -i -P20 sh -c "verifyIndels.pl {} $baseDir $refFile"`;
print  $err;
parseAssemblyResults($sampleFile);











#---------------------------------------------------------------------------
#  Sub Routines
#---------------------------------------------------------------------------

sub parseAssemblyResults{
### start ###
my $samplefile = $_[0];
my @samples;

open F, $samplefile or die;
while(<F>){
	chomp;
	push @samples,$_;
}
close(F);

my %meta;
my @indels;
open F, "indelGenes.positions.txt" or die;
while (<F>){
	chomp;
	my ($name) = (split /\s+/,$_);
	$meta{$name} = $_;
	push @indels,$name;
}

my %results;

foreach my $sample ( @samples ) {
	open F, "$sample/sv.results.txt" or die;
	while(<F>){
		chomp;
		my ($name,$result) = (split /\s+/,$_)[0,4];
		if ($result eq "ABSENT"){
			$results{$name}{$sample} = 0;
		} elsif ($result eq "PRESENT"){
			$results{$name}{$sample} = 1;
		} elsif ($result eq "NA"){
			$results{$name}{$sample} = "NA";
		} elsif ($result eq "MAYBE"){
			$results{$name}{$sample} = "NA";
		} else {
			die"Sort it out $sample $name\n";
		}
	}
	close(F);
}


open OUT, ">$outfile" or die;
print OUT "chr\tpos\tgenes";
foreach my $sample ( @samples ){
	print OUT "\t$sample";
}
print OUT "\n";

foreach my $indelName ( sort keys %results ) {
	my ($chr,$start) = (split /\s+/,$meta{$indelName})[1,4];
	print OUT "$chr\t$start\t$indelName";
	foreach my $sample ( @samples ) {
		print OUT "\t$results{$indelName}{$sample}";
	}
	print OUT "\n";
}

### end ###
}

sub parseVCF{

######### start ###########



my $samplefile = $_[0];
my $base_dir = $_[1];
my @samples;
open F, $samplefile or die;
while(<F>){
	chomp;
	push @samples,$_;
}
close(F);

my %genes;
open GFF, "$base_dir/Mtb.gff" or die;
while(<GFF>){
	chomp;
	my ($chr,$type,$start,$end,$temp) = (split /\t/,$_)[0,2,3,4,8];
	if ($type ne "gene"){next;}
	my $gene;
	if ($temp =~ /locus_tag/){
		$temp =~ m/locus_tag (\S+)/;
		$gene = $1;
	}
#	print "$chr\t$type\t$start\t$end";
	if (!$gene){
		$gene  = "NA";		
	}
	$genes{$gene}{'chr'} = $chr;
	$genes{$gene}{'start'} = $start;	
	$genes{$gene}{'end'} = $end;
}
close(GFF);


my %indels;
foreach my $sample ( @samples ) {
	my ($chr,$qual,$start,$end) = parseDelly($sample,$base_dir);
	my @chr = @{$chr};
	my @qual = @{$qual};
	my @start = @{$start};
	my @end = @{$end};
	
	for (my $i=0; $i<=$#chr; $i++){
		my $indelLen = $end[$i] - $start[$i];
		if ($indelLen > 10000){next;}
#		print "$sample\t$chr[$i]\t$qual[$i]\t$start[$i]\t$end[$i]\n";
		$indels{$chr[$i]}{$start[$i]} = $end[$i];
	}
}

#print "\n\n****************************************************************\n\n";
my %maxSize;
my %indelsUniq;
foreach my $chr ( sort keys %indels ) {
	foreach my $start ( sort {$a<=>$b} keys %{$indels{$chr}} ) {
		my $end = $indels{$chr}{$start};
		my $indelgenes = findGenes($chr,$start,$end,\%genes);
		print "$chr\t$start\t$end\t$indelgenes\n";
		$indelsUniq{$indelgenes} ++;
		if (!exists($maxSize{$indelgenes})){
			$maxSize{$indelgenes}{'start'} = $start;
			$maxSize{$indelgenes}{'end'} = $end;
		} else {
			if ($maxSize{$indelgenes}{'start'}<$start){
				$start = $maxSize{$indelgenes}{'start'};
			}
			if ($maxSize{$indelgenes}{'end'}=$end){
				$end = $maxSize{$indelgenes}{'end'};
			}
		}
	}
}

#print "\n\n****************************************************************\n\n";

open OUT, ">indelGenes.positions.txt" or die;
foreach my $indelgenes ( sort keys %indelsUniq ) {
	my @indelgenes = split /%/,$indelgenes;
	my @bp;
	my $tempChr;
	for (@indelgenes){
		if (!exists($genes{$_})){
			my ($chr,$start,$end) = (split /-/,$_);
			push @bp, $start;
			push @bp, $end;					
			$tempChr = $chr;
			next;
		}
		push @bp,$genes{$_}{'start'};
		push @bp,$genes{$_}{'end'};
		$tempChr = $genes{$_}{'chr'};
	}
	my $max = max @bp;
	my $min = min @bp;
	if (!$min){
		die"$indelgenes\n";
	}
	print OUT "$indelgenes\t$tempChr\t$min\t$max\t$maxSize{$indelgenes}{'start'}\t$maxSize{$indelgenes}{'end'}\n";
}
close OUT;


############# end #############
}


sub parseDelly{
#	print "Parsing $_[0] indels\n";
	my $sample = $_[0];
	my $base_dir = $_[1];
	my @chr;
	my @qual;
	my @start;
	my @end;
	open F, "$base_dir/vcf/$sample.vcf" or die "Can't find $sample.vcf\n";
	while(<F>){
		chomp;
		if ($_ =~ /^#/){next;}
		my ($chr,$start,$qual,$temp) = (split /\s+/,$_)[0,1,6,7];
		if ($qual ne "PASS"){next;}
		$temp =~ m/END=(\d+)/;
		my $end = $1;
		push @chr,$chr;
		push @qual,$qual;
		push @start,$start;
		push @end,$end;
	}
	close(F);
	return(\@chr,\@qual,\@start,\@end);
}


sub findGenes{

	my ($chr,$start,$end,$pointGenes) = @_;
	my %genes = %{$pointGenes};
	my $tempGenes;
	foreach my $gene ( keys %genes ){
		if (($start > $genes{$gene}{'start'} && $start < $genes{$gene}{'end'}) or ($end > $genes{$gene}{'start'} && $start < $genes{$gene}{'end'})){
			$tempGenes .= "%$gene";
#			print "$chr\t$start\t$end\t$gene\n";
		}
	}
	if (!$tempGenes){
		$tempGenes = "%$chr-$start-$end";
	} 
	my $res = substr($tempGenes,1);
	return($res);

}
