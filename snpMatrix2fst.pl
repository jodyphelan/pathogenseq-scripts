#!/usr/bin/perl 

use strict;
use warnings;
use Term::ProgressBar;
use Getopt::Long qw(GetOptions);

my $matFile;
my $annFile;
my $outFile;
my $type;

my $usage = '
snpMatrix2fst.pl
	--mat|-m	Matrix
	--ann|-a	Annotation
	--out|-o	Output
	--type|-t	Method [wright|nei]

';

GetOptions(
	'mat|m=s' => \$matFile,
	'ann|a=s' => \$annFile,
	'out|o=s' => \$outFile,
	'type|t=s' => \$type,
) or die $usage;

if (!$matFile or !$annFile or !$outFile or !$type){
	print $usage; exit;
}

if (($type ne "wright") and ($type ne "nei")){
	print "Method should be wright or nei\n";exit;
}

if ($type eq "wright"){
	wright($matFile,$annFile,$outFile);
} elsif ($type eq "nei"){
	nei($matFile,$annFile,$outFile);
}


sub nei{

##### start

my $totSNPs = `wc -l < $_[0]`;
chomp $totSNPs;
my $progress = Term::ProgressBar->new({name => 'Computing', count => $totSNPs, remove => 1});
$progress->minor(0);
my $next_update = 0;
my $x = 0;


my %category;
my %categories;
open F, $_[1] or die;
open OUT, ">$_[2]" or die;

while(<F>){
	chomp;
	my ($sample,$category) = split /\s+/,$_;
	$category{$sample} = $category;
	if ($category eq "NA"){next;}
	$categories{$category} = 1;
}
close(F);

my @categories = sort keys %categories;
print OUT "chr\tpos\tFst";
for (@categories){
	print OUT "\t$_";
}
print OUT "\n";
my %fsts;
my %keepSamples;
open F, $_[0] or die;
my $q=0;
while(<F>){
	chomp;
	if ($q<1){
		$x++;
		$q++;
		my @a = split /\s+/,$_;
		my $chr = shift @a;
		my $pos = shift @a;
		my $ref = shift @a;
		for (my $i=0; $i<=$#a; $i++){
			if ($a[$i] eq "NA"){next;}
			if (exists($category{$a[$i]})){
				$keepSamples{$i} = $category{$a[$i]}; #keepSamples{1} = CAT1
			}
		}
		next;
	}
	my @a = split /\s+/,$_;
	my $chr = shift @a;
	my $pos = shift @a;
	my $ref = shift @a;
#	print "$chr\t$pos\n";

	my %gt;
	my %N;
	my $totAlt =0;
	my %values; ### stores number of values per binary key (0,1)
	for (my $i=0; $i<=$#a; $i++){
		if (exists($keepSamples{$i})){
			if ($a[$i] eq "NA"){next;}
			$gt{$keepSamples{$i}} += $a[$i];
#			print "Adding $a[$i] to $keepSamples{$i}\n";
			$N{$keepSamples{$i}} ++;
			$totAlt += $a[$i];
			$values{$a[$i]} ++;
		}
	}
	if ($totAlt == 0){
#		print "No alts\n";
		next;
	}
	if (scalar keys %values ==1){ 
#		print "All alts\n"; 
		next;
	}
	my $Af1;
	my $Af2;
	my $Hs;
	my %af;
	foreach my $cat (@categories){
		my $af1 = $gt{$cat}/$N{$cat};
		my $af2 = 1-($gt{$cat}/$N{$cat});
		$Af1 += $af1;
		$Af2 += $af2;
		my $H = 2*$af1*$af2;
#		print "$cat\tAf1=$af1\tAf2=$af2\tH=$H\n";
		$Hs += $H;
#		print "Allele freq of $cat is $af{$cat}\n"
		$af{$cat} = $af1;
	} 
	my $numCats = scalar @categories;
	$Hs = $Hs/$numCats;
	$Af1 = $Af1/$numCats;
	$Af2 = $Af2/$numCats;
	my $Ht = 2*$Af1*$Af2;
	my $Fst = ($Ht-$Hs)/$Ht;
#	print "Hs=$Hs\nHt=$Ht\nFst=$Fst\n";
	print OUT "$chr\t$pos\t$Fst";
	foreach my $cat (@categories){
		print OUT "\t$af{$cat}";
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


##### end 

}


sub wright{

#### Start ####

my %category;
my %categories;
open F, $_[1] or die;
open OUT, ">$_[2]" or die;

while(<F>){
	chomp;
	my ($sample,$category) = split /\s+/,$_;
	$category{$sample} = $category;
	if ($category eq "NA"){next;}
	$categories{$category} = 1;
}
close(F);

my @categories = sort keys %categories;
print OUT "chr\tpos\tFst";
for (@categories){print OUT "\t$_";}
print OUT "\n";


my $totSNPs = `wc -l < $_[0]`;
chomp $totSNPs;
my $progress = Term::ProgressBar->new({name => 'Computing', count => $totSNPs, remove => 1});
$progress->minor(0);
my $next_update = 0;
my $x = 0;




my %fsts;
my %keepSamples;
open F, $_[0] or die;
my $q=0;
while(<F>){
	chomp;
	if ($q<1){
		$q++;
		$x++;
		my @a = split /\s+/,$_;
		my $chr = shift @a;
		my $pos = shift @a;
		my $ref = shift @a;
		for (my $i=0; $i<=$#a; $i++){
			if ($a[$i] eq "NA"){next;}
			if (exists($category{$a[$i]})){
				$keepSamples{$i} = $category{$a[$i]}; #keepSamples{1} = CAT1
			}
		}
		next;
	}
	my @a = split /\s+/,$_;
	my $chr = shift @a;
	my $pos = shift @a;
	my $ref = shift @a;
#	print "$chr\t$pos\n";

	my %gt;
	my %N;
	my $totSamps =0;
	for (my $i=0; $i<=$#a; $i++){
		if (exists($keepSamples{$i})){
			if ($a[$i] eq "NA"){next;}
			$gt{$keepSamples{$i}} += $a[$i];
#			print "Adding $a[$i] to $keepSamples{$i}\n";
			$N{$keepSamples{$i}} ++;
			$totSamps ++ ;
		}
	}

	my $tempNumCats = scalar keys %N;
	my $tempNumCats2 = scalar @categories;
	if ($tempNumCats != $tempNumCats2){
		print OUT "$chr\t$pos\tNA";
		for (@categories){
			print OUT "\tNA";
		}
		print OUT "\n";
	
		next;
	}


	if ($totSamps == 0){
#		print "No alts\n"; 
		print OUT "$chr\t$pos\tNA"; for (@categories){print OUT "\tNA";} print OUT "\n"; $x++; next;}
	my $sumGT;
	for (@categories){
			$sumGT+=$gt{$_};
	}
	if ($sumGT == 0){
#		print "No alts\n"; 
		print OUT "$chr\t$pos\tNA"; for (@categories){print OUT "\tNA";} print OUT "\n"; $x++;next;}

	my $pbar;
	foreach my $cat (@categories){
		$pbar += $gt{$cat};
#		print "Allele freq of $cat is $gt{$cat}\n"
	} 
	my $numCats = scalar @categories;
#	print "temp pba  = $pbar\n";
#	print "NumCats: $numCats\n";

	$pbar = $pbar/$totSamps;
	if ($pbar == 1){print OUT "$chr\t$pos\tNA"; for (@categories){print OUT "\tNA";} print OUT "\n"; $x++; next;}
#	print "pbar = $pbar\n";
	my $nbar = $totSamps/$numCats;
	my $term1;
#	print "nbar = $nbar\n";
	foreach my $cat (@categories){
		my $af = $gt{$cat}/$N{$cat};
		$term1 += $N{$cat}*(($af-$pbar)**2);
	}
	my $topTerm = $term1/(($numCats*$nbar));
########	my $topTerm = $term1/(($numCats-1)*$nbar);
	my $Fst = $topTerm/($pbar*(1-$pbar)); 
#	print "Pbar = $pbar\nnbar = $nbar\nterm1 = $term1\ntopTerm = $topTerm\n";
#	print "Fst = $topTerm/($pbar*(1-$pbar)) = $Fst\n";

	print OUT "$chr\t$pos\t$Fst";
	for (@categories){
		my $af = $gt{$_}/$N{$_};
		print OUT "\t$af";
	}
	print OUT "\n";

	$x++;
	if ($x >= $next_update){
        $next_update = $progress->update($x);
   }



}
close(F);

$x++;
if ($x >= $next_update){
       $next_update = $progress->update($x);
}




#### End
}

