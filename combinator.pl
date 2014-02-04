#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

if (@ARGV != 2 ) { die "./combinator.pl <GFF3> <FPKM>\n";}

my %contig = ();

#open mannotatored.gff file

open(MAN,$ARGV[0]) || die "Couldn't open file $ARGV[0]\n";
while (my $mannotator = <MAN>)
{
        next if $mannotator =~ /#/;
        chomp $mannotator;
#       print $mannotator;
        my @columns = split (/\t/, $mannotator);
        my $gene = "$columns[0],$columns[3],$columns[4]";
#       print $gene, "\t", $columns[8], "\n";
        $contig{$gene} = $columns[8];
}
close MAN;
#print Dumper(%contig);

#open fpkm file
#print "Opening $ARGV[1]\n";

open my $FPKM, "<", $ARGV[1] or die "Couldn't open file $ARGV[1]\n";
while (<$FPKM>)
{
        chomp $_;
        my @names = split (/\t/, $_);
#       print "**$names[0]**\n";
        if (exists $contig{$names[0]})
        {
        #       print $contig{$names[0]};
                print "$names[0]\t$contig{$names[0]}\t$names[1]\n";
        }
}

exit;

