#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  snpMatrix2fasta.pl
#
#        USAGE:  ./snpMatrix2fasta.pl  
#
#  DESCRIPTION:  convert snp_matrix to fasta file
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  Jody Phelan (JP), jody.phelan@lshtm.ac.uk
#      COMPANY:  LSHTM
#      VERSION:  1.0
#      CREATED:  02/17/2015 04:56:12 PM
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;

if (scalar @ARGV != 2){print "\nsnpMatrix2fasta.pl <snp_matrix> <out_file>\n\n"; exit;}
my @samples;
my %seqs;
my @pos;

my $i = 0;
open F, $ARGV[0] or die"Can't open $ARGV[0]";
print "Loading SNP matrix\n";
while (<F>){
	chomp;	
	if ($i < 1){
		$i ++; 
		my @a = split(/\s+/, $_);
		shift @a; shift @a; shift @a;
		for ( my $x = 0; $x<=$#a; $x++){
			$seqs{$x}{"seq"} = "";
			$seqs{$x}{"sample"} = $a[$x];
		}
		next;
	}	
	my @a = split(/\s+/, $_);
	shift @a;
	my $pos = shift @a;
	push @pos, $pos;
	print "Loading position $pos\n";
	my $ref = shift @a;
	for (my $j = 0; $j<=$#a; $j++){
		if ($a[$j] eq "-"){$a[$j] = "N";}
		$seqs{$j}{"seq"} .= $a[$j];
	}
}
close(F);


open OUT, ">$ARGV[1]" or die;
foreach my $key ( sort {$a<=>$b} keys %seqs){
	my $sample = $seqs{$key}{'sample'};
	print "Writing sequence for $sample\n";
	print OUT ">$sample\n";
	print OUT "$seqs{$key}{'seq'}\n";
}
