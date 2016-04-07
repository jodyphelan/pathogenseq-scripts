#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  snpMatrix_removeMulti.pl
#
#        USAGE:  ./snpMatrix_removeMulti.pl  
#
#  DESCRIPTION:  remove multi calls from a snp matrix
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  Jody Phelan (mn), jody.phelan@lshtm.ac.uk
#      COMPANY:  LSHTM
#      VERSION:  1.0
#      CREATED:  03/31/2015 12:13:23 PM
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;


if ($#ARGV+1 < 2 ){print "\n./snpMatrix_removeMulti.pl <snp_matrix> <outfile>\n\n"; exit;}

my %final_calls;
my @pos;
my @samples;
my $i = 0;
my %ref_calls;
my $outfile = $ARGV[1];


open OUT, ">$outfile" or die;
print OUT "chr\tpos\tref";

print "Loading snp matrix\n";
open F, $ARGV[0] or die"Can't open $ARGV[0]";
while (<F>){	
	if ($i < 1){
		$i ++; 
		my @a = split(/\s+/, $_);
		for ( my $x = 3; $x<=$#a; $x++){
			print OUT "\t$a[$x]";
			push @samples, $a[$x];
		}
		print OUT "\n";
		next;
	}	
	chomp;
	my @a = split(/\s+/, $_);
	my $chr = shift @a;
	my $pos = shift @a;
	push @pos, $pos;
	my $ref = shift @a;
	print "$pos\n";
	$ref_calls{$pos} = $ref;
	
	my $line .= "$chr\t$pos\t$ref_calls{$pos}";
	my %alt_present;

	for (my $j = 0; $j<=$#a; $j++){
		if (length($a[$j])>1){ 
			$a[$j] = "N";
		}

		if ($a[$j] ne "N" && $a[$j] ne "-"){
			$alt_present{$a[$j]} = 1
		}
		$line .= "\t$a[$j]";
	}
	if (scalar keys %alt_present > 1){print OUT "$line\n";}
}
close(F);



