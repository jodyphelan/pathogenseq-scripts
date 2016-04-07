#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  bigSnpMatrix_subset.pl
#
#        USAGE:  ./bigSnpMatrix_subset.pl  
#
#  DESCRIPTION:  
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  Dr. Fritz Mehner (mn), mehner@fh-swf.de
#      COMPANY:  FH Südwestfalen, Iserlohn
#      VERSION:  1.0
#      CREATED:  12/02/2015 01:19:42 PM
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;


if (scalar @ARGV != 3){ print "\nbigSnpMatrix_subset.pl <samples> <snp_mat> <outfile>  \n\n"; exit;}

my %samples;
open F, $ARGV[0] or die;
while(<F>){
	chomp;
	$samples{$_} = 1
}
close(F);

open OUT, ">$ARGV[2]" or die;
my $x = 0;
my %keepSamples;
open F, $ARGV[1] or die;
while(<F>){
	chomp;
	if ($x<1){
		$x++;
		my @a = split /\s+/,$_;
		my $chr = shift @a;
		my $pos = shift @a;
		my $ref = shift @a;
		print OUT "$chr\t$pos\t$ref";
	
		for (my $i=0; $i<=$#a; $i++){
			if (exists($samples{$a[$i]})){
				$keepSamples{$i} = $a[$i];
				print OUT "\t$a[$i]";
			}
		}
		print OUT "\n";
		next;

	}
	my @a = split /\s+/,$_;
	my $chr = shift @a;
	my $pos = shift @a;
	my $ref = shift @a;
	print "$chr\t$pos\n";
	my $line = "$chr\t$pos\t$ref";
	my %alt_present;
	for (my $i=0; $i<=$#a; $i++){
		if (exists($keepSamples{$i})){
			$line .= "\t$a[$i]";
			if ($a[$i] eq "N" or $a[$i] eq "-"){next;}
			$alt_present{$a[$i]} = 1;
		}
	}
	if (scalar keys %alt_present > 1){print OUT "$line\n";}
}
close(F);
