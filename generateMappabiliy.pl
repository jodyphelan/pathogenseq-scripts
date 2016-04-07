#! /usr/bin/perl

if (scalar @ARGV != 4){ print "\ngenerateMappability.pl <base_dir> <prefix> <kmer> <threads>\n\n";exit;}

my $bd = $ARGV[0];
my $prefix = $ARGV[1];
my $kmer = $ARGV[2];
my $threads = $ARGV[3];
my $ref = "$bd/$prefix.fasta";
`mkdir $bd/mappability`;
chdir("$bd/mappability");
if (!-e $ref){print "Cannot find $ref\n";exit;}
`gem-indexer -i $ref -o index -T $threads`;
if (!-e "index.gem"){print "Gem indexer failed\n";exit;}
`gem-mappability -T $threads -I index.gem -l $kmer -o $prefix`;
if (!-e "$prefix.mappability"){print "Gem mappability failed\n";exit;}
`gem-2-wig -I index.gem -i $prefix.mappability -o $prefix`;

`wig2bed < $prefix.wig > $prefix.mappability.bed`;

