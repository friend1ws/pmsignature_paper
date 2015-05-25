#! /usr/local/bin/perl

use strict;

my $input = $ARGV[0];
my $numBase = $ARGV[1];
my $dirOrNot = $ARGV[2];

if (not ($numBase == 3 or $numBase == 5)) {
    die "the number of bases should be 3 or 5\n";
}

if (not ($dirOrNot == 0 or $dirOrNot == 1)) {
    die "the third parameter (direction ? or not) should be 0 or 1\n";
}

my %sample2index = ();
my $index = 1;

my %base2ind = ();
$base2ind{"A"} = 1;
$base2ind{"C"} = 2;
$base2ind{"G"} = 3;
$base2ind{"T"} = 4;

my %mut2ind = ();
$mut2ind{"C>A"} = 1;
$mut2ind{"C>G"} = 2;
$mut2ind{"C>T"} = 3;
$mut2ind{"T>A"} = 4;
$mut2ind{"T>C"} = 5;
$mut2ind{"T>G"} = 6;


# my %sample2exists = ();
# my $sampleNum = 0;
# open(IN, $input) || die "cannot open $!";
# while(<IN>) {
#     s/[\r\n]//g;
#     my @F = split("\t", $_);

#     next if (length($F[3]) >= 2);
#     next if ($F[2] =~ /N/);
#     next if ($F[3] =~ /N/);

#     if ($dirOrNot == 1) {
#         next if ($F[6] ne "+" and $F[6] ne "-");
#     }

#     if (not exists $sample2exists{$F[1]}) {
#         $sample2exists{$F[1]} = 1;
#         $sampleNum = $sampleNum + 1;
#     }
# }


 
open(IN, $input) || die "cannot open $!";

while(<IN>) {
    s/[\r\n]//g;
    my @F = split("\t", $_);

    next if (length($F[3]) >= 2);
    next if ($F[2] =~ /N/);
    next if ($F[3] =~ /N/);

    if ($dirOrNot == 1) {
        next if ($F[6] ne "+" and $F[6] ne "-");
    }


    # if (not exists $sample2index{$F[1]}) {
    #     $sample2index{$F[1]} = $index;
    #     $index = $index + 1;
    # }        

    my $fSeq = "";
    my $altBase = "";
    my $fStrand = "";
    if (substr($F[2], 2, 1) eq "A" or substr($F[2], 2, 1) eq "G") {
        $fSeq = complementSeq($F[2]);
        $altBase = complementSeq($F[3]);

        if ($F[6] eq "+") {
            $fStrand = 2;
        } elsif ($F[6] eq "-") {
            $fStrand = 1;
        } else {
            $fStrand = 3;
        }

    } else {
        $fSeq = $F[2];
        $altBase = $F[3];

        if ($F[6] eq "+") {
            $fStrand = 1;
        } elsif ($F[6] eq "-") {
            $fStrand = 2;
        } else {
            $fStrand = 3;
        }


    }

    next if ($dirOrNot == 1 and $fStrand == 3);

    print $F[1];
    print "\t" . $mut2ind{substr($fSeq, 2, 1) . ">" . $altBase};

    if ($numBase == 3) {
        for(my $i = 1; $i < length($fSeq) - 1; $i++) {
            next if ($i == 2);
            print "\t" . $base2ind{substr($fSeq, $i, 1)};
        }
    } elsif ($numBase == 5) {
        for(my $i = 0; $i < length($fSeq); $i++) {
            next if ($i == 2);
            print "\t" . $base2ind{substr($fSeq, $i, 1)};
        }
    } 

    if ($dirOrNot == 0) {
        print "\n";
    } elsif ($dirOrNot == 1) {
        print "\t" . $fStrand . "\n";
    } 

}


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


