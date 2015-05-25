#! /usr/local/bin/perl

use strict;

my $input = $ARGV[0];

my $mid = 0;
open(IN, $input) || die "cannot open $!";
while(<IN>) {
    s/[\r\n]//g;
    my @F = split("\t", $_);
    next if ($F[8] ne "SBS");

    my $strand = $F[4];
    (my $ntref, my $ntalt) = split(">", $F[11]);

    my $fseq = substr($F[12], 8, 5);
    if ($strand eq "-") {
        $ntref = &complementSeq($ntref);
        $ntalt = &complementSeq($ntalt); 
        $fseq = &complementSeq($fseq);
    }

    # if (substr($F[12], 10, 1) ne $ntref) {
    #     print $ntref . "\t" . $ntalt . "\n";
    # }
    
    $F[9] =~ s/^g\.//g;
    $F[9] =~ s/[ACGT]>[ACGT]$//g;
    (my $chr, my $pos) = split(":", $F[9]);

    print "mutation_" . $mid . "\t" . $F[2] . "\t" . $fseq . "\t" . $ntalt . "\t" . $chr . "\t" . $pos . "\t" . $strand . "\n";
    $mid = $mid + 1;

}
close(IN);



sub complementSeq {

    my $tseq = reverse($_[0]);

    $tseq =~ s/A/S/g;
    $tseq =~ s/T/A/g;
    $tseq =~ s/S/T/g;

    $tseq =~ s/C/S/g;
    $tseq =~ s/G/C/g;
    $tseq =~ s/S/G/g;

    return $tseq;
}
    
