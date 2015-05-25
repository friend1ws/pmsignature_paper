#! /usr/local/bin/perl

use strict;


print STDERR "Start reading refGene.txt\n";
my $n = 0;
open(IN, "refGene.txt") || die "cannot open $!";
while(<IN>) {
    s/[\r\n\"]//g;
    my @F = split("\t", $_);

    my $chr = $F[2]; 
    my @starts = split(",", $F[9]);
    my @ends = split(",", $F[10]);
    my $strand = $F[3];
    my $exonNum = $F[8];
    my $gene = $F[1];
    my $symbol = $F[12];

    for (my $i = 0; $i <= $#starts; $i++) {

        my $key = $chr . "\t" . $starts[$i] . "\t" . $ends[$i];
        if ($strand eq "+") {
            print $key . "\t" . $symbol . "(" . $gene . ")" . "." . $i . "\t" . "0" . "\t" . "+" . "\n";
        } else {
            print $key . "\t" . $symbol . "(" . $gene . ")" . "." . ($exonNum - $i - 1) . "\t" . "0" . "\t" . "-" . "\n";
        }

    }


    $n = $n + 1;
    if ($n % 1000 == 0) {
        print STDERR "$n genes completed.\n";
    }

}
close(IN);
print STDERR "Reading refGene.txt completed.\n";


