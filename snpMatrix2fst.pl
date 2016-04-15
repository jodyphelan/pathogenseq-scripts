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
	if ($category eq "NA"){next;}
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
#			print "Adding $a[$i] to $keepSamples{$i}\n";
			$N{$keepSamples{$i}} ++;
			$totAlt += $a[$i];
		}
	}
	if ($totAlt == 0){print "No alts\n"; next;}
	my %af;
	foreach my $cat (@categories){
		$af{$cat} = $gt{$cat}/$N{$cat};
#		print "Allele freq of $cat is $af{$cat}\n"
	} 
	my @cats = keys %gt;
	

	print OUT "$chr\t$pos";
	foreach my $masterCat (@categories) {
#		print "Analysing $masterCat\n";
		my @gt;
		my @N;
		foreach my $cat (@categories){
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

		my @H;
		$H[0] = 2*($af[0]*(1-$af[0]));
		$H[1] = 2*($af[1]*(1-$af[1]));
		my $Hs = ($H[0] + $H[1])/2;
		my $af1 = ($af[0] + $af[1])/2;
		my $af2 = (1-$af[0] + 1-$af[1])/2;
		my $Ht = 2*$af1*$af2;
		my $Fst = ($Ht-$Hs)/$Ht;
#		printf("AF: %f\t%f\n", $af[0], $af[1]); 		
#		print "H1: $H[0]\nH2: $H[1]\n";
#		print "Hs = $Hs\nHt = $Ht\nFst = $Fst\n";
		print OUT "\t$Fst";
	}
	print OUT "\n";


}
close(F);






