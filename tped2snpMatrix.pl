#!/usr/bin/perl 

use strict;
use warnings;
use Term::ProgressBar;

if (scalar @ARGV != 2){ print "\nrecode.pl <prefix> <outfile>\n\n"; exit;}

open F, "$ARGV[0].tfam" or die;
open OUT, ">$ARGV[1]" or die;
print OUT "chr\tpos\tid";
while(<F>){
	my $id = (split /\s+/,$_)[1];
	print OUT "\t$id";
}
close(F);
print OUT "\n";


my $totSNPs = `wc -l < $ARGV[0].tped`;
chomp $totSNPs;
my $progress = Term::ProgressBar->new({name => 'Converting', count => $totSNPs, remove => 1});
$progress->minor(0);
my $next_update = 0;
my $x = 0;




open F, "$ARGV[0].tped" or die;
while(<F>){
	chomp;
	$_ =~ s/(\S+\s\S+\s\d+\s\d+)\s//;
	my $info = $1;
	$info =~ s/ /\t/g;
	my ($chr,$id,$pos) = (split /\s+/, $info)[0,1,3];
	print OUT "$chr\t$pos\t$id";
#	print "$chr $pos\n";
	$_ =~ s/0/NA/g;
#	$_ =~ s/1/0/g;
	$_ =~ s/2/0/g;
	my @a = split /\s+/,$_;
	for (my $i = 0; $i<=$#a; $i=$i+2){
		if ($a[$i] eq "NA" and $a[$i+1] eq "NA"){print OUT "\tNA"; next;}
		
		my $gt = ($a[$i]+$a[$i+1])/2;
		print OUT "\t$gt";
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
