#!/usr/bin/perl 

use strict;
use warnings;
use Term::ProgressBar;

if (scalar @ARGV != 2){ print "\nbigSnpMatrix2binary.pl <snp_mat> <outfile>  \n\n"; exit;}

my $totSNPs = `wc -l < $ARGV[0]`;
chomp $totSNPs;
my $progress = Term::ProgressBar->new({name => 'Filtering', count => $totSNPs, remove => 1});
$progress->minor(0);
my $next_update = 0;
my $x = 0;

my $i = 0;
open OUT, ">$ARGV[1]" or die;
open F, $ARGV[0] or die;
while(<F>){
	chomp;
	if ($i<1){$x++;$i++; print OUT "$_\n"; next;}
	my @a = split /\s+/,$_; # hello archie! what are you doing?
	my $chr = shift @a;
	my $pos = shift @a;
	my $ref = shift @a;
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
	$x++;
	if ($x >= $next_update){
		$next_update = $progress->update($x);
	}
}
close(F);

if ($x >= $next_update){
	$next_update = $progress->update($x);
}
