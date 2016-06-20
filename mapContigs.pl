#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  coords2overlap.pl
#
#        USAGE:  ./coords2overlap.pl  
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
#      CREATED:  06/15/2016 04:21:08 PM
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;

my $sample = $ARGV[0];
my $ref = $ARGV[1];
`nucmer --prefix=$sample $ref contigs/$sample.contigs 2> $sample.log`;
`show-snps -ClrT $sample.delta > $sample.snps`;
`show-coords -rT $sample.delta > $sample.coords`;

if (!-d "rawFiles"){`mkdir rawFiles`;}
`mv *.coords *.snps *.log *.delta rawFiles`;

