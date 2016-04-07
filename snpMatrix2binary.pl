#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  bigSnpMatrix2binary.pl
#
#        USAGE:  ./bigSnpMatrix2binary.pl  
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
#      CREATED:  01/12/15 17:40:41
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;


if (scalar @ARGV != 2){ print "\nbigSnpMatrix2binary.pl <snp_mat> <outfile>  \n\n"; exit;}
my $i = 0;
open OUT, ">$ARGV[1]" or die;
open F, $ARGV[0] or die;
while(<F>){
	chomp;
	if ($i<1){$i++; print OUT "$_\n"; next;}
	my @a = split /\s+/,$_; # hello archie! what are you doing?
	my $chr = shift @a;
	my $pos = shift @a;
	my $ref = shift @a;
	print "$pos\n";
	print OUT "$chr\t$pos\t$ref";
	for (my $i=0; $i<=$#a; $i++){
		if ($a[$i] eq "-" || $a[$i] eq "N"){
			print OUT "\tNA";
		} elsif ($a[$i] eq $ref){
			print OUT "\t0";
		} elsif ($a[$i] ne $ref){
			print OUT "\t1";
		} else {
			die;
		}
	}
	print OUT "\n";
}
close(F);
