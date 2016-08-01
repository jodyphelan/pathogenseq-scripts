#!/usr/bin/perl 

use strict;
use warnings;
use Getopt::Long qw(GetOptions);


my $samples;
my $sampleFile;
my $baseDir;
my $pos;

my $usage = '
	generateCovPlot.pl
		--samples|-s	Samples seperated by commas
		--file|-f		File containing sample names
		--base|-b		Base directory
		--pos|p			Position

';

GetOptions(
	'samples|s=s' => \$samples,
	'file|f=s' => \$sampleFile,
	'base|b=s' => \$baseDir,
	'pos|p=s' => \$pos,
) or die $usage;

if (!$pos){	print $usage;exit;}

if (!$samples and !$sampleFile and $pos){print "Provide sample or sampleFile\n"; exit;}
if ($samples and $sampleFile and $pos){print "Provide sample or sampleFile\n"; exit;}

my @samples;
if ($samples){
	 @samples = split /,/,$samples;
} else {
	open F, $sampleFile or die;
	while(<F>){
		chomp;
		push @samples,$_;
	}
	close(F);
}



open OUT, ">tempCovPlot/files.txt" or die;
if (!-d "tempCovPlot"){`mkdir tempCovPlot`;}
foreach my $sample ( @samples ) {
	if (!-e "$baseDir/coverage/$sample.coverage.gz.tbi"){
		die "$baseDir/coverage/$sample.coverage.gz.tbi";
	}
	`tabix $baseDir/coverage/$sample.coverage.gz $pos > tempCovPlot/$sample.cov`;
	print OUT "tempCovPlot/$sample.cov\n";
}

writeRfiles();
`Rscript covplot.r`;

sub writeRfiles{

open R, ">covplot.r" or die;
print R "x11()\n";
if (scalar @samples > 1){	print R "par(mfrow=c(",scalar @samples,",1))\n";}
print R "
library(zoo)
samples<-read.table(\"tempCovPlot/files.txt\")[,1]
par(mar=c(3,4,1,1)+0.1)

for (s in samples){
	sname<-gsub(\"tempCovPlot/\",\"\",s)
	sname<-gsub(\".cov\",\"\",sname)
	dat<-read.table(s)
#	y<-rollmean(dat[,3],k=100)
	y<-dat[,3]
	x<-dat[1:length(y),2]
	xx<-c(x,rev(x))
	yy<-c(y,rep(0,length(y)))
	plot(y ~ x, type=\"n\",ylim=c(0,max(y)),ylab=\"Total coverage\",main=sname)	
	polygon(xx,yy,col=\"light blue\",ylim=c(0,60),border=\"NA\")

}
locator(1)
";
}

