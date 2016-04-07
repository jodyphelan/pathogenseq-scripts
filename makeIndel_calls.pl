#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  makeIndel_calls.pl
#
#        USAGE:  ./makeIndel_calls.pl  
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
#      CREATED:  01/14/2016 03:27:21 PM
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;


my $base_dir = $ARGV[1];
my $sample = $ARGV[0];

my $i=0;
my %indelPos;
open F, "indel.pos.txt" or die;
while(<F>){
	if ($i<1){$i++;next;}
	chomp;
	my $pos = (split /\s+/,$_)[1];
	
	$indelPos{$pos} = 0;
}
close(F);

open F, "zcat $base_dir/vcf/$sample.vcf.gz |" or die;  ######## THIS SHOULD BE GATK ###########
while(<F>){
	chomp;
	if ($_ =~ /#/){next;}
	my ($chromosome,$pos,$ref,$alt) = (split /\s+/,$_)[0,1,3,4];
	if (exists($indelPos{$pos})){
		if (length($ref) != length($alt)){
			$indelPos{$pos} = 1;
		}
	}		
}
close(F);

open OUT, ">$sample.indel.calls" or die;
print OUT "$sample\n";
foreach my $pos ( sort {$a<=>$b} keys %indelPos ) {
	print OUT "$indelPos{$pos}\n";
}
