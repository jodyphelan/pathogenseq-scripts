#!/usr/bin/perl 
#===============================================================================
#
#			FILE:  extract_cov.pl
#
#		  USAGE:  ./extract_cov.pl  
#
#  DESCRIPTION:  
#
#		OPTIONS:  ---
# REQUIREMENTS:  ---
#			BUGS:  ---
#		  NOTES:  ---
#		 AUTHOR:  Dr. Fritz Mehner (mn), mehner@fh-swf.de
#		COMPANY:  FH Südwestfalen, Iserlohn
#		VERSION:  1.0
#		CREATED:  11/30/2015 12:25:03 PM
#	  REVISION:  ---
#===============================================================================

use strict;
use warnings;

if (scalar @ARGV != 4){ print "\nextract_cov.pl <sample> <base_dir> <%cov> <min_cov>  \n\n"; exit;}

my %lhAll_pos;
my $percentTotalcov = $ARGV[2];
my $min_cov = $ARGV[3];

open F, "snps.pos.txt" or die;
while(<F>){
	chomp;
	my ($chr,$pos) = (split /\s+/,$_);
	$lhAll_pos{$chr}{$pos} = 1;
}
close(F);

my $sample = $ARGV[0];
my $base_dir = $ARGV[1];
my $cov_file = "$base_dir/coverage/$sample.coverage.gz";
my @laCov_lines;

print "Extracting coverage for $sample\n";
if (!-e $cov_file){
	 print "Cant't find coverage file  for $sample\n";
	 exit;
}

my $x=0;
open F, "zcat $cov_file |";
open OUT, ">$sample.calls" or die;
open COV, ">$sample.COV" or die;
print OUT "$sample\n";

my %calls;
my %cov;
open QUAL, ">$sample.sample.qual" or die;
my $NA = 0;
my $count = 0;
my $mix = 0;
my @totalCov;
while(<F>){
	$count++;
	
	if ($x<1){$x++; next;}
	my ($test1,$test2) = (split /\s+/,$_)[0,1];
	$test2=$test2+1;
	if (!exists($lhAll_pos{$test1}{$test2})){next;}
	my ($chr, $pos, $tot, $cov_A, $cov_C, $cov_G, $cov_T) = (split /\s+/,$_)[0,1,2,3,4,5,6,8];
	$pos = $pos+1;
	push @totalCov,$tot;
	my $p = $percentTotalcov*$tot;
	my $cutoff;
	if ($p>=$min_cov){$cutoff = $p;} else {$cutoff = $min_cov;}
	my $c = "";
	if ($tot == 0){
		$c.="-";
	} else {
		if ($cov_A >= $cutoff){ $c.="A";}
		if ($cov_C >= $cutoff){ $c.="C";}
		if ($cov_G >= $cutoff){ $c.="G";}
		if ($cov_T >= $cutoff){ $c.="T";}
	}
#				print "Cutoff: $cutoff | A: $cov_A | C: $cov_C | G: $cov_G | T: $cov_T ......... ";
	if ($c eq ""){ $c.="N";}
	if ($c eq "N" or $c eq "-"){$NA++;}
	if (length($c) > 1){$mix++;}
#				print "$c\n";
	$calls{$chr}{$pos} = $c;
	$cov{$chr}{$pos} = $tot;
}

foreach my $chr (sort keys %calls){
	foreach my $pos (sort {$a<=>$b} keys %{$calls{$chr}}){
		print OUT "$calls{$chr}{$pos}\n";
		print COV "$cov{$chr}{$pos}\n";
	}	
}
my $medCov = median(@totalCov);
my $fNA = $NA/$count;
my $fMix = $mix/$count;
print QUAL "$sample\t$fNA\t$fMix\t$medCov\n";
close(F);
close(OUT);


sub median
{
        return unless @_;
        return $_[0] unless @_ > 1;
        @_= sort{$a<=>$b}@_;
        return $_[$#_/2] if @_&1;
        my $mid= @_/2;
        return ($_[$mid-1]+$_[$mid])/2;
}

