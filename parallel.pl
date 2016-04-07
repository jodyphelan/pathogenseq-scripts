#! /usr/bin/perl

use strict;
use warnings;

if (scalar @ARGV != 6){print "parallel.pl <sample_list> <ref> <threads_per_run>  <temp_dir> <storage_dir> <total_threads>\n";exit;}

my @samples;
open F, $ARGV[0] or die;
while(<F>){
	chomp;
	push @samples,$_;
}
close(F);

my $line = "";
my $hostname = `hostname`;
chomp $hostname;
open OUT, ">$hostname.xargs" or die;
foreach my $sample (@samples){
	$line .= " all $sample $ARGV[1] $ARGV[2] $ARGV[3] $ARGV[4]";
}
print OUT $line;
close(OUT);
`cat $hostname.xargs | xargs -n6 -P $ARGV[5] map_call_snps_0.2.pl`;
