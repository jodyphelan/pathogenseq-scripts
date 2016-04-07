#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  snpMatrix2genome.pl
#
#        USAGE:  ./snpMatrix2genome.pl  
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
#      CREATED:  08/12/2015 12:04:18 PM
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;

if (scalar @ARGV != 3){print "\nsnpMatrix2fasta.pl <snp_matrix> <ref.fa> <out_file>\n\n"; exit;}
my @samples;
my %seqs;
my @pos;

my $i = 0;
open F, $ARGV[0] or die"Can't open $ARGV[0]";
print "Loading SNP matrix\n";
while (<F>){
    if ($i < 1){
        $i ++;
        my @a = split(/\s+/, $_);
        for ( my $x = 3; $x<=$#a; $x++){
            push @samples, $a[$x];
        }
        next;
    }
    chomp;
    my @a = split(/\s+/, $_);
    shift @a;
    my $pos = shift @a;
    push @pos, $pos;
    my $ref = shift @a;
    for (my $j = 0; $j<=$#a; $j++){
        $seqs{$samples[$j]}{$pos} = $a[$j];
    }
}
close(F);

my ($r1,$r2) = hash_contigs("$ARGV[1]");

my %hasContigs = %{$r1};
my @arrContigs = @{$r2};
my $ref = $hasContigs{$arrContigs[0]};

my %hasGenomes;

open OUT,">$ARGV[2]" or die;
foreach my $sample (@samples){
	print "Analysing $sample\n";
	my @arrGenome = split //,$ref;

	foreach my $pos (sort {$a<=>$b} keys %{$seqs{$sample}}){
		$pos = $pos-1;
		print "$arrGenome[$pos] $pos $seqs{$sample}{($pos+1)}\n";
		$arrGenome[$pos] = $seqs{$sample}{($pos+1)};
	}

	my $strGenome = "";
	for (@arrGenome){
		$strGenome = $strGenome.$_;
	}
	print OUT ">$sample\n$strGenome\n";
}



sub hash_contigs {
	if (@_ != 1) {
		print "Usage: hash_contigs contigs_file";
		exit;
	}
	my $contigs_file = shift;
	if( $contigs_file =~ /\.gz$/ ){
		open(CONTIGS, $contigs_file, "gunzip -c $contigs_file |" ) or die "Cant open contigs file: $!\n";
	} else {
		open(CONTIGS, $contigs_file) or die "Cant open contigs file: $!\n";
	}
	my %contigs_hash; # hash to store contig names and sequences
		my $contigName;
	my @ids;
	while (<CONTIGS>) {
		if (/^>(\S+)/) {
			$contigName=$1;
			push @ids, $contigName;
		} else {
			chomp;
			$contigs_hash{$contigName} .= $_;
		}
	}
	close(CONTIGS);
	return(\%contigs_hash,\@ids);
}
