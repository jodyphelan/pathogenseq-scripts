#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  snpMatrix2summaryVCF.pl
#
#        USAGE:  ./snpMatrix2summaryVCF.pl  
#
#  DESCRIPTION: convert a snp matrix to summary vcf one line per snp 
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  Jody Phelan (mn), jody.phelan@lshtm.ac.uk
#      COMPANY:  LSHTM
#      VERSION:  1.0
#      CREATED:  03/27/2015 10:58:48 AM
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;

if (scalar @ARGV != 2){print "\nsnpMatrix2fasta.pl <snp_matrix> <out_file>\n\n"; exit;}

my %final_calls;
my @pos;
my @samples;
my $i = 0;
my %ref_calls;
my $outfile = $ARGV[1];

print "Loading snp matrix\n";
open F, $ARGV[0] or die"Can't open $ARGV[1]";
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
	my $chr = shift @a;
	my $pos = shift @a;
	push @pos, $pos;
	my $ref = shift @a;
	$ref_calls{$chr}{$pos} = $ref;
	for (my $j = 0; $j<=$#a; $j++){
		$final_calls{$chr}{$pos}{$samples[$j]} = $a[$j];
	}
}
close(F);

print "Printing VCF file\n";
open OUT, ">$outfile" or die;
my $line1 = "##fileformat=VCFv4.1";
my $line2 = '##FORMAT=<ID=GT,Number=1,Type=Integer,Description="Genotype">)';
my $line3 = "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample";

print OUT $line1."\n".$line2."\n".$line3."\n";

foreach my $chr (sort keys %final_calls){
	foreach my $pos (sort {$a<=>$b} keys %{$final_calls{$chr}}){
		print "Loading position $chr $pos\n";
		my %genotypes;
		foreach my $sample (@samples){
			if (length($final_calls{$chr}{$pos}{$sample})>1){
				for (my $i = 0; $i<length($final_calls{$chr}{$pos}{$sample}); $i++){
					my $nuc = substr($final_calls{$chr}{$pos}{$sample}, $i, 1);
					if ($nuc eq "-" or $nuc eq "N"){next;}
					if (!exists $genotypes{$nuc} && $nuc ne $ref_calls{$chr}{$pos}){
						$genotypes{$nuc} = 1;
					} elsif (exists $genotypes{$nuc} && $nuc ne $ref_calls{$chr}{$pos}) {
						$genotypes{$nuc} ++;
					}
				}
			} else {	
				my $nuc = $final_calls{$chr}{$pos}{$sample};
				if ($nuc eq "-" or $nuc eq "N"){next;}
				if (!exists $genotypes{$nuc} && $nuc ne $ref_calls{$chr}{$pos}){
					$genotypes{$nuc} = 1;
				} elsif(exists $genotypes{$nuc} && $genotypes{$nuc} ne $ref_calls{$chr}{$pos}) {
					$genotypes{$nuc} ++;
				}
			}
		}
		print OUT "$chr\t$pos\t.\t$ref_calls{$chr}{$pos}\t";
		my $alts = "";
		foreach my $alt ( sort keys %genotypes){
			$alts .= ",$alt";
		}	
		print OUT substr($alts,1);
		print OUT "\t.\tPASS\t.\tGT\t1/1\n";
	}
}

close OUT;
