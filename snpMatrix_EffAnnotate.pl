#!/usr/bin/perl 
use strict;
use warnings;

if ($#ARGV+1 < 3 ){print "\nsubset_snpMatrix.pl <snp_matrix> <snpEFF VCF> <outfile> \n\n"; exit;}

my %final_calls;
my @pos;
my @samples;
my $i = 0;
my %ref_calls;
my $outfile = $ARGV[2];

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

my %hasVCF;

print "Loading VCF\n";
open F, $ARGV[1] or die"Can't open $ARGV[1]";
while(<F>){
	chomp;
	if ($_ =~/#/){next;}
	my ($chr,$pos,$alt,$info) = (split /\s+/,$_)[0,1,4,7];
	if (!exists($ref_calls{$chr}{$pos})){next;}
	my ($type,$gene,$aa_change) = (split /\|/, $info)[1,3,10];
	$hasVCF{$chr}{$pos}{$alt}{type} = $type;	
	$hasVCF{$chr}{$pos}{$alt}{gene} = $gene;
#	$hasVCF{$chr}{$pos}{$alt}{aa_change} = $aa_change;
}
close(F);


#---------------------------------------------------------------------------
#  Print final matrix
#---------------------------------------------------------------------------

print "Writing snp matrix\n";
open OUT, ">$outfile" or die;
print OUT "chr\tpos\tref";

foreach my $sample ( @samples ) {
	print OUT "\t$sample";
}
print OUT "\n";

foreach my $chr ( sort keys %final_calls){
	foreach my $pos ( sort {$a<=>$b} keys %{$final_calls{$chr}}){
		my $line;
		$line .= "$chr\t$pos\t$ref_calls{$chr}{$pos}";
		my %genotypes;
		foreach my $sample ( @samples ) {
			$line .= "\t$final_calls{$chr}{$pos}{$sample}";
			$genotypes{$final_calls{$chr}{$pos}{$sample}} = 1;
		}
		my $type = "";
		my $gene = "";
		my $aa_change = "";
		foreach my $alt (sort keys %genotypes){
			if ($alt eq $ref_calls{$chr}{$pos}){next;}
			if ($alt eq "N"){next;}
			if ($alt eq "-"){next;}
			if (!exists($hasVCF{$chr}{$pos}{$alt})){ print "Can't find $chr $pos $alt in hash\n"; exit; }
			$type .= "%".$hasVCF{$chr}{$pos}{$alt}{type};
#			if (!exists($hasVCF{$chr})){ die"$chr $pos";}
			$gene .= "%".$hasVCF{$chr}{$pos}{$alt}{gene};
#			if ($hasVCF{$chr}{$pos}{$alt}{aa_change} ne ""){
#				$aa_change .= "%".$hasVCF{$chr}{$pos}{$alt}{aa_change};
#			}
#			print "Chr: $chr Pos: $pos Alt: $alt Type: $type\n";
		}
#		if ($aa_change eq ""){$aa_change = "-";} else {$aa_change = substr($aa_change,3);}
		$type = substr($type,1);

		$gene = substr($gene,1);

		$line .= "\t$type\t$gene";		
		print OUT "$line\n";
	}
}



close OUT;

