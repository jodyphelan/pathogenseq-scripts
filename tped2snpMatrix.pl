#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  recode.pl
#
#        USAGE:  ./recode.pl  
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
#      CREATED:  04/07/2016 10:25:02 AM
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;


if (scalar @ARGV != 2){ print "\nrecode.pl <prefix> <outfile>\n\n"; exit;}

open F, "$ARGV[0].tfam" or die;
open OUT, ">$ARGV[1]" or die;
print OUT "chr\tpos\tid";
while(<F>){
	my $id = (split /\s+/,$_)[1];
	print OUT "\t$id";
}
close(F);
print OUT "\n";

open F, "$ARGV[0].tped" or die;
while(<F>){
	chomp;
	$_ =~ s/(\S+\s\S+\s\d+\s\d+)\s//;
	my $info = $1;
	$info =~ s/ /\t/g;
	my ($chr,$id,$pos) = (split /\s+/, $info)[0,1,3];
	print OUT "$chr\t$pos\t$id";
	print "$chr $pos\n";
	$_ =~ s/0/NA/g;
	$_ =~ s/1/0/g;
	$_ =~ s/2/1/g;
	my @a = split /\s+/,$_;
	for (my $i = 0; $i<=$#a; $i=$i+2){
		if ($a[$i] eq "NA" and $a[$i+1] eq "NA"){print OUT "\tNA"; next;}
		
		my $gt = ($a[$i]+$a[$i+1])/2;
		print OUT "\t$gt";
	}
	print OUT "\n";
}
close(F);
