#! /usr/local/bin/perl

use strict;

my $input = $ARGV[0];

open(IN, $input) || die "cannot open $!";
while(<IN>) {
    s/[\r\n]//g;
    my @F = split("\t", $_);

    $F[2] =~ s/23/X/;
    $F[2] =~ s/24/Y/;
    next if ($F[2] eq "25");
    next if ($F[2] eq "MT");

    next if ($F[1] ne "subs");

    print $F[0] . "\t" . "chr" . $F[2] . "\t" . $F[3] . "\t" . $F[5] . "\t" . $F[6] . "\n";
}
close(IN);


