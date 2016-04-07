#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  permuteRef2vcf.pl
#
#        USAGE:  ./permuteRef2vcf.pl  
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
#      CREATED:  03/29/2016 01:38:31 PM
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;


if (scalar @ARGV != 2){ print "\npermuteRef2vcf.pl <ref.fa> <output>\n\n"; exit;}

my ($r1,$r2) = hash_contigs($ARGV[0]);
my %hasContigs = %{$r1};
my @arrContigs = @{$r2};

open OUT, ">$ARGV[1]" or die;
print OUT "##fileformat=VCFv4.1\n";
print OUT '##FORMAT=<ID=GT,Number=1,Type=Integer,Description="Genotype">)';
print OUT "\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n";

foreach my $contig (@arrContigs){
	for (my $i = 0; $i<=(length($hasContigs{$contig})-1); $i++){
		my $nuc = substr($hasContigs{$contig},$i,1);
		$nuc = uc($nuc);
		my $pos = $i+1;
		if ($nuc eq "A"){
		print OUT "$contig\t$pos\t.\t$nuc\tC\t30\tPASS\t.\tGT\t1/1\n";
		print OUT "$contig\t$pos\t.\t$nuc\tG\t30\tPASS\t.\tGT\t1/1\n";
		print OUT "$contig\t$pos\t.\t$nuc\tT\t30\tPASS\t.\tGT\t1/1\n";
		} elsif ($nuc eq "C"){
		print OUT "$contig\t$pos\t.\t$nuc\tA\t30\tPASS\t.\tGT\t1/1\n";
		print OUT "$contig\t$pos\t.\t$nuc\tG\t30\tPASS\t.\tGT\t1/1\n";
		print OUT "$contig\t$pos\t.\t$nuc\tT\t30\tPASS\t.\tGT\t1/1\n";
		} elsif ($nuc eq "G"){
		print OUT "$contig\t$pos\t.\t$nuc\tA\t30\tPASS\t.\tGT\t1/1\n";
		print OUT "$contig\t$pos\t.\t$nuc\tC\t30\tPASS\t.\tGT\t1/1\n";
		print OUT "$contig\t$pos\t.\t$nuc\tT\t30\tPASS\t.\tGT\t1/1\n";
		} elsif ($nuc eq "T"){
		print OUT "$contig\t$pos\t.\t$nuc\tA\t30\tPASS\t.\tGT\t1/1\n";
		print OUT "$contig\t$pos\t.\t$nuc\tC\t30\tPASS\t.\tGT\t1/1\n";
		print OUT "$contig\t$pos\t.\t$nuc\tG\t30\tPASS\t.\tGT\t1/1\n";
		}
	}
}

close(OUT);

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
