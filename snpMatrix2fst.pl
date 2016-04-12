#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  snpMatrix2fst.pl
#
#        USAGE:  ./snpMatrix2fst.pl  
#
#  DESCRIPTION:  
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  Jody Phelan (mn), jody.phelan@lshtm.ac.uk
#      COMPANY:  LSHTM
#      VERSION:  1.0
#      CREATED:  04/12/2016 12:31:27 PM
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;


if (scalar @ARGV != 3){ print "\nsnpMatrix2fst.pl <snp_mat.bin> <categories> <out>\n\n"; exit;}

my %category;
my %categories;
open F, $ARGV[1] or die;
open OUT, ">$ARGV[2]" or die;

while(<F>){
	chomp;
	my ($sample,$category) = split /\s+/,$_;
	$category{$sample} = $category;
	$categories{$category} = 1;
}
close(F);

my @categories = sort keys %categories;
print OUT "chr\tpos";
for (@categories) {
	print OUT "\t$_";
}
print OUT "\n";

my %fsts;
my %keepSamples;
open F, $ARGV[0] or die;
my $x=0;
while(<F>){
	chomp;
	if ($x<1){
		$x++;
		my @a = split /\s+/,$_;
		my $chr = shift @a;
		my $pos = shift @a;
		my $ref = shift @a;
		for (my $i=0; $i<=$#a; $i++){
			if ($a[$i] eq "NA"){next;}
			if (exists($category{$a[$i]})){
				$keepSamples{$i} = $category{$a[$i]}; #keepSamples{1} = CAT1
			}
		}
		next;
	}
	my @a = split /\s+/,$_;
	my $chr = shift @a;
	my $pos = shift @a;
	my $ref = shift @a;
	print "$chr\t$pos\n";

	my %gt;
	my %N;
	my $totAlt =0;
	for (my $i=0; $i<=$#a; $i++){
		if (exists($keepSamples{$i})){
			if ($a[$i] eq "NA"){next;}
			$gt{$keepSamples{$i}} += $a[$i];
			$N{$keepSamples{$i}} ++;
			$totAlt += $a[$i];
		}
	}
	if ($totAlt == 0){print "No alts\n"; next;}
	my %af;
	foreach my $cat (keys %gt){
		$af{$cat} = $gt{$cat}/$N{$cat};
	} 
	my @cats = keys %gt;
	
	my @gt;
	my @N;
	print OUT "$chr\t$pos";
	foreach my $masterCat (@categories) {
		print "Analysing $masterCat\n";
		foreach my $cat (keys %gt){
			if($cat eq $masterCat){
				$gt[0] = $gt{$cat};
				$N[0] = $N{$cat};
			} else {
				$gt[1] += $gt{$cat};
				$N[1] += $N{$cat};
			}
		}
		my @af;
		$af[0] = $gt[0]/$N[0];
		$af[1] = $gt[1]/$N[1];
	
		my $totSamps = $N[0] + $N[1];
		my $nBar = ($N[0]/2) + ($N[1]/2);
		my $pBar = ($N[0]*$af[0]/$totSamps) + ($N[1]*$af[1]/$totSamps);
		
		my $fst = (($N[0]*(($af[0]-$pBar)**2)/$nBar) + ($N[1]*(($af[1]-$pBar)**2)/$nBar))/$pBar*(1-$pBar);
		printf("%f\t%f", $af[0], $af[1]);
		print "\t$fst\n";
		print OUT "\t$fst";
	}
	print OUT "\n";


}
close(F);






