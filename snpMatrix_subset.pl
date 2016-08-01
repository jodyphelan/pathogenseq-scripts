#!/usr/bin/perl 

use strict;
use warnings;
use Term::ProgressBar;


if (scalar @ARGV != 3){ print "\nbigSnpMatrix_subset.pl <samples> <snp_mat> <outfile>  \n\n"; exit;}

my %samples;
open F, $ARGV[0] or die;
while(<F>){
	chomp;
	$samples{$_} = 1
}
close(F);

my $totSNPs = `wc -l < $ARGV[1]`;
chomp $totSNPs;
my $progress = Term::ProgressBar->new({name => 'Filtering', count => $totSNPs, remove => 1});
$progress->minor(0);
my $next_update = 0;
my $x = 0;


open OUT, ">$ARGV[2]" or die;
my $q = 0;
my %keepSamples;
open F, $ARGV[1] or die;
while(<F>){
	chomp;
	if ($q<1){
		$q++;
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
#	print "$chr\t$pos\n";
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
	
	$x++;
	if ($x >= $next_update){
		$next_update = $progress->update($x);
	}
}
close(F);

if ($x >= $next_update){
        $next_update = $progress->update($x);
}
