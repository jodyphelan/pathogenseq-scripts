#!/usr/bin/perl 
use strict;
use warnings;
use Term::ProgressBar;
if (scalar @ARGV != 4){ print "\nsnpMatrix_filter.pl <snpMatrix> <bed_file> <out> <exclude/extract>\n\n"; exit;}

my $alg = $ARGV[3];
if ($alg ne "exclude" and $alg ne "extract"){
	print "Filtering should be \"include\" or \"exclude\"\n";
	exit;
}

my %bed;
open F, $ARGV[1] or die;
while(<F>){
	chomp;
	my ($chr,$start,$end) = (split /\s+/,$_)[0,1,2];
	for (my $i = $start; $i<=$end; $i++){
		$bed{$chr}{$i} = 1;
	}
}
close(F);

my $totSNPs = `wc -l < $ARGV[0]`;
chomp $totSNPs;
my $progress = Term::ProgressBar->new({name => 'Filtering', count => $totSNPs, remove => 1});
$progress->minor(0);
my $next_update = 0;
my $x = 0;

my $q = 0;
open F, $ARGV[0] or die;
open OUT, ">$ARGV[2]" or die;
while(<F>){
	chomp;
	if ($q<1){$x++;$q++;print OUT "$_\n";next;}
	my ($chr,$pos) = (split /\s+/,$_)[0,1];
	if ($alg eq "exclude"){
		if (!exists($bed{$chr}{$pos})){
			print OUT "$_\n";
		}
	} elsif ($alg eq "extract"){
		if (exists($bed{$chr}{$pos})){
			print OUT "$_\n";
		}
	}
	$x++;
	if ($x >= $next_update){
		$next_update = $progress->update($x);
	}
}
close(F);

if ($x >= $next_update){
        $next_update = $progress->update($x);
}


