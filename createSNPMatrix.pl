#!/usr/bin/perl 
#===============================================================================
#
#		 FILE:  createSNPMatrix.pl
#
#		USAGE:  ./createSNPMatrix.pl  
#
#  DESCRIPTION:  
#
#	  OPTIONS:  ---
# REQUIREMENTS:  ---
#		 BUGS:  ---
#		NOTES:  ---
#	   AUTHOR:  Jody Phelan (mn), jody.phelan@lshtm.ac.uk
#	  COMPANY:  LSHTM
#	  VERSION:  1.0
#	  CREATED:  06/15/2016 06:05:58 PM
#	 REVISION:  ---
#===============================================================================

use strict;
use warnings;
use POSIX qw(ceil);
use Term::ProgressBar;
use Statistics::Lite qw(:all);
 

if (scalar @ARGV != 3){ print "\n$0 <samples> <base_dir> <threads>\n\n"; exit;}

my $samplesFile = $ARGV[0];
my $base_dir = $ARGV[1];
my $threads = $ARGV[2];

raw($samplesFile,$base_dir,$threads);
print "Filtering samples\n\n";
writeRfiles();
my ($l1) = samplePlot();
sampleFilter($l1);
snpQual();
my ($r1) = snpPlot();
filterSNPs($r1);



sub raw {
########### START ###########
my ($samplesFile,$base_dir,$threads)  = @_;
my @samples;
open F, $ARGV[0] or die;
while(<F>){
	chomp;
	push @samples,$_;
}
close(F);

my %snpPos;
foreach my $sample (@samples){
	print "Extracting SNPs for $sample\n";
	my $x=0;
	open F, "$base_dir/rawFiles/$sample.snps" or die;
	while(<F>){
		chomp;
		if ($x<4){$x++;next;}
		my ($pos,$ref,$alt,$chr) = (split /\s+/,$_)[0,1,2,10];
		if ($alt eq "." or $ref eq "."){next;} 		######### this skips the indels!!!! #####
#		print "$sample\t$chr\t$pos\t$ref\t$alt\n";
		$snpPos{$chr}{$pos} = $ref;
	}
	close(F);
}


open OUT, ">chr.pos.ref.txt" or die;
print OUT "chr\tpos\tref\n";
foreach my $chr (sort keys %snpPos) {	
	foreach my $pos ( sort {$a<=>$b} keys %{$snpPos{$chr}} ) {
		print OUT "$chr\t$pos\t$snpPos{$chr}{$pos}\n";
	}
}
close OUT;


writeCallsScript();

print "Generating calls\n";
`cat $ARGV[0] | xargs -i -P $threads sh -c "perl generateAssemblyCalls.pl {} $base_dir"`;


my @tempSamples = @samples;
my $tempNum = (scalar @samples)/1000;
$tempNum = ceil($tempNum);

for (my $i=1; $i<=$tempNum; $i++){
	open OUT,">samples.$i" or die;
	for (my $j=0; $j<=999; $j++){
		if (!@tempSamples){next;}
		my $tempSample = shift @tempSamples;
		print OUT "$tempSample.calls ";
	}
}

print "Generating final matrix\n";
open TMP, ">tempMat.txt" or die;
for (my $i=1; $i<=$tempNum; $i++){
	print TMP "tempMat.$i ";
	`paste \`cat samples.$i\` > tempMat.$i`;
}
close TMP;
`paste \`cat tempMat.txt\` > tempMat.final`;
`paste chr.pos.ref.txt tempMat.final > unfiltered.mat`;

`rm samples.* tempMat*`;
`cat *.sample.qual > sample.quals`;
`mv *.sample.qual *.calls rawFiles`;
########## END #########
}



sub sampleFilter{
### SUB START ###

my $missCutoff = $_[0];
my %quals;
my @samples;
my %samples;
open F, "sample.quals" or die;
while(<F>){
        chomp;
        my ($sample,$miss) = split /\s+/,$_;
        if ($miss<$missCutoff){
                push @samples,$sample;
                $samples{$sample} = 1;
        }
}
close(F);

#progress bar
my $totSNPs = `wc -l < unfiltered.mat`;
chomp $totSNPs;
my $progress = Term::ProgressBar->new({name => 'Filtering', count => $totSNPs, remove => 1});
$progress->minor(0);
my $next_update = 0;

open BIN, ">sample.filt.mat.bin" or die;
open OUT, ">sample.filt.mat" or die;
my $x = 0;
my %keepSamples;
open F, "unfiltered.mat" or die;
while(<F>){
        chomp;
        if ($x<1){
                $x++;
                my $line = "chr\tpos\tref";
                my @a = split /\s+/,$_;
                @a = @a[3..$#a];
                  for (my $i=0; $i<=$#a; $i++){
                        if (exists($samples{$a[$i]})){
                                $keepSamples{$i} = $a[$i];
                                $line .= "\t$a[$i]";
                        }
                }
                print OUT "$line\n";
                print BIN "$line\n";
                next;
        }
        if (scalar keys %keepSamples ==0 ){ print "No samples left after filtering\n"; exit;}

        if ($_ =~ /\w+\t(\d+)\t\w+\t(.+)/){
                my $pos = $1;
                my $line = $2;
                my $a = () = $line =~ /A/g;
                my $c = () = $line =~ /C/g;
                my $g = () = $line =~ /G/g;
                my $t = () = $line =~ /T/g;
                my $count = () = "%$a%%$c%%$g%%$t%" =~ /%0%/g;
                if ($count == 3){
                        $x++;
                        next;            # this removes polmorhic snps
                }
        }


        my @a = split /\s+/,$_;
        my ($chr,$pos,$ref) = @a[0..2];
        @a = @a[3..$#a];
        my %alt_present;
        my $line = "$chr\t$pos\t$ref";
        my $binLine = "$chr\t$pos\t$ref";
        for (my $i=0; $i<=$#a; $i++){
                if (exists($keepSamples{$i})){
                        $line .= "\t$a[$i]";

                        if ($a[$i] eq "-" || $a[$i] eq "N"){
                                $binLine .= "\tNA";
                        } elsif ($a[$i] eq $ref){
                                $binLine .= "\t0";
                        } elsif ($a[$i] ne $ref){
                                if (length($a[$i])==1){
                                        $binLine.= "\t1";
                                } elsif (length($a[$i])>1){
                                        $binLine.= "\tMX";
                                } else {
                                        die;
                                }
                        } else {
                                die;
                        }

                        if ($a[$i] eq "N" or $a[$i] eq "-"){next;}
                        $alt_present{$a[$i]} = 1;
                }
        }
        if (scalar keys %alt_present > 1){
                print OUT "$line\n";
                print BIN "$binLine\n";
        }

        $x++;
        if ($x >= $next_update){
                $next_update = $progress->update($x);
        }
}
if ($x >= $next_update){
        $next_update = $progress->update($x);
}
close(F);
close(OUT);
close(BIN);

### SUB END ###
}



sub snpQual{

my $snpMat = "sample.filt.mat.bin";
my @samples;

my $totSNPs = `wc -l < sample.filt.mat.bin`;
chomp $totSNPs;
my $progress = Term::ProgressBar->new({name => 'Filtering', count => $totSNPs, remove => 1});
$progress->minor(0);
my $next_update = 0;


my $x = 0;
open OUT, ">snps.qual" or die;
open F, $snpMat or die;
while(<F>){
        chomp;
        if ($x<1){
                $x++;
                my $line = "chr\tpos\tref";
                my @a = split /\s+/,$_;
                @a = @a[3..$#a];
                for (my $i=0; $i<=$#a; $i++){
                        push @samples,$a[$i];
                }
                next;
        }
        my $line;
        my $numSamples = $#samples+1;
        if ($_ =~ /([^\s]+)\t(\d+)\t\w+\t(.+)/){
                my $chr = $1;
                my $pos = $2;
                $line = $3;
                my $zeros = () = $line =~ /0/g;
                my $ones = () = $line =~ /1/g;
                my $NA = () = $line =~ /NA/g;
                print OUT "$chr\t$pos\t",$zeros/$numSamples,"\t",$ones/$numSamples,"\t",$NA/$numSamples,"\n";
        }
        $x++;
        if ($x >= $next_update){
                $next_update = $progress->update($x);
        }
}
close(F);
close(OUT);

}



sub filterSNPs{


my $cutMiss = $_[0];

my %snpQual;
open F, "snps.qual" or die;
while(<F>){
        chomp;
        my ($chr,$pos,$missing,$mix) = (split /\s+/,$_)[0,1,4,5];
        $snpQual{$chr}{$pos}{"miss"} = $missing;
}
close(F);

my $totSNPs = `wc -l < sample.filt.mat`;
chomp $totSNPs;
my $progress = Term::ProgressBar->new({name => 'Filtering', count => $totSNPs, remove => 1});
$progress->minor(0);
my $next_update = 0;
open OUT, ">snp.map.sample.filt.mat" or die;
open LOG, ">snps.log" or die;
my $x = 0;
open F, "sample.filt.mat" or die;
while(<F>){
    if ($x<1){$x++; print OUT $_; next;}
        my ($chr,$pos) = (split /\s+/,$_)[0,1];
        if (!defined($snpQual{$chr}{$pos}{"miss"} ) or !defined($snpQual{$chr}{$pos}{"miss"})){die"$chr $pos\n"};
        if ($snpQual{$chr}{$pos}{"miss"} > $cutMiss){
                print LOG "$chr\t$pos\t$snpQual{$chr}{$pos}{'miss'}\tLOW_QUAL\n";
                $x++;
                next;
        }
        print LOG "$chr\t$pos\t$snpQual{$chr}{$pos}{'miss'}\tOK\n";
        print OUT  $_;
        $x++;
        if ($x >= $next_update){
                $next_update = $progress->update($x);
        }
}
 if ($x >= $next_update){
        $next_update = $progress->update($x);
}


}






sub writeRfiles{
open OUT, ">samplePlot.R" or die;
print OUT 'x<-read.table("sample.quals")
x11()
plot(sort(x$V2),main="Proportion of missing calls over total SNP positions",ylab="Missing proportion")
y<-locator(1)
pdf("sample-missing.pdf")
plot(sort(x$V2),main="Proportion of missing calls over total SNP positions",ylab="Missing proportion")
abline(h=y$y,col="red",lty=2)
abline(v=y$x,col="blue",lty=2)
err<-dev.off()
print(c(y$y))';
close(OUT);


open OUT, ">snpPlot.R" or die;
print OUT 'x<-read.table("snps.qual")
x11()
plot(sort(x$V5),main="Proportion of missing calls over each SNP",ylab="Missing proportion")
y<-locator(1)

pdf("snp-missing.pdf")
plot(sort(x$V5),main="Proportion of missing calls over each SNP",ylab="Missing proportion")
abline(h=y$y,col="red",lty=2)
abline(v=y$x,col="blue",lty=2)
err<-dev.off()

print(c(y$y))';
close(OUT);
}








sub writeCallsScript{

open SC, ">generateAssemblyCalls.pl" or die;
print SC '
#!/usr/bin/perl 
use strict;
use warnings;


my $sample = $ARGV[0];
my $base_dir = $ARGV[1];


open F, "rawFiles/$sample.snps" or die;
my %snpPos;
my $x=0;
while(<F>){
	chomp;
	if ($x<4){$x++;next;}
	my ($pos,$ref,$alt,$chr) = (split /\s+/,$_)[0,1,2,10];
    if ($alt eq "." or $ref eq "."){$snpPos{$chr}{$pos} = "N"; next;}      ######### this skips the indels!!!! #####

	$snpPos{$chr}{$pos} = $alt;
}
close(F);


print "Extracting contig overlap\n";
my %overlap;
$x = 0;
open F, "rawFiles/$sample.coords" or die;
while(<F>){
	chomp;
	if ($x<4){$x++;next;}
	my ($start,$stop,$chr) = (split /\s+/,$_)[0,1,7];
	for (my $i=$start; $i<=$stop; $i++){
		$overlap{$chr}{$i} ++;
	}
}
close(F);

my $missing = 0;
my $positions = 0;
open F, "chr.pos.ref.txt" or die;
open OUT, ">$sample.calls" or die;
print OUT "$sample\n";
my $q=0;
while(<F>){
	chomp;
    if ($q<1){$q++;next;}
	my ($chr,$pos,$ref) = split /\s+/,$_;
    if (!exists($overlap{$chr}{$pos})){$overlap{$chr}{$pos} = 0;$missing++;$positions++;}
	if ($overlap{$chr}{$pos} != 1){
		print OUT "N\n";$missing++;$positions++;next;
	}
	$positions ++;
	if (exists($snpPos{$chr}{$pos})){
			print OUT "$snpPos{$chr}{$pos}\n";
	} else {
			print OUT "$ref\n";
	}
}
close(F);

open QV,">$sample.sample.qual" or die;
my $pct = $missing/($positions);
print QV "$sample\t$pct\n";
close QV;

';
}






sub samplePlot{


print "Please select cutoffs for mixed missing and median coverage per sample:\n";
my $x = `Rscript samplePlot.R`;
chomp $x;
my ($l1) = (split /\s+/,$x)[1];
print "Mixed call fraction:$l1\n";
print "Do you want to use these values? [Y/n]";
my $ans1 = <STDIN>;
chomp $ans1;
$ans1 = uc($ans1);

while ($ans1 ne "" && $ans1 ne "Y" && $ans1 ne "N"){
    print "Do you want to use these values? [Y/n]";
    $ans1 = <STDIN>;
    chomp $ans1;
    $ans1 = uc($ans1);
}
if ($ans1 eq "" or $ans1 eq "Y"){
    print "Using previously defined values\n";
} else {
    print "Please provide new missiness cutoff:\n";
    $l1 = <STDIN>;
    chomp $l1;
}

print "\n-----------------------------------------------------------\n";
print "Filtering snps using value: $l1\n";
print "-----------------------------------------------------------\n";
return($l1);

}


sub snpPlot{


print "Please select cutoffs for missing and mexed fraction per sample:\n";
my $x = `Rscript snpPlot.R`;
chomp $x;
my ($l1) = (split /\s+/,$x)[1];
print "Mixed call fraction:$l1\n";
print "Do you want to use these values? [Y/n]";
my $ans1 = <STDIN>;
chomp $ans1;
$ans1 = uc($ans1);

while ($ans1 ne "" && $ans1 ne "Y" && $ans1 ne "N"){
    print "Do you want to use these values? [Y/n]";
    $ans1 = <STDIN>;
    chomp $ans1;
    $ans1 = uc($ans1);
}
if ($ans1 eq "" or $ans1 eq "Y"){
    print "Using previously defined values\n";
} else {
    print "Please provide new missiness cutoff:\n";
    $l1 = <STDIN>;
    chomp $l1;
}

print "\n-----------------------------------------------------------\n";
print "Filtering snps using value: $l1\n";
print "-----------------------------------------------------------\n";
return($l1);

}

