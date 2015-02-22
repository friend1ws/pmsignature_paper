#! /usr/local/bin/perl

use strict;

my $input = $ARGV[0];
my $bed = $ARGV[1];

my $numBase = $ARGV[2];
my $dirOrNot = $ARGV[3];

if (not ($numBase == 3 or $numBase == 5)) {
    die "the number of bases should be 3 or 5\n";
}

if (not ($dirOrNot == 0 or $dirOrNot == 1)) {
    die "the third parameter (direction ? or not) should be 0 or 1\n";
}


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

my %dir2ind = ();
$dir2ind{"+"} = 1;
$dir2ind{"-"} = 2;

my %mutID2count = ();
my $mutSum = 0;

my %pos2dir = ();
open(IN, $bed) || die "cannot open $!";
while(<IN>) {
    s/[\r\n]//g;
    my @F = split("\t", $_);
    $pos2dir{$F[0] . ":" . $F[1] . "-" . $F[2]} = $F[5];
}
close(IN);

$DB::single = 1;
open(IN, $input) || die "cannot open $!";

while(<IN>) {

    s/[\r\n]//g;
    $_ =~ s/^>//;
    my $pos = $_;
    my $dir = $pos2dir{$pos};

    if ($dir ne "+" and $dir ne "-") {
        print $pos . "\t" . $dir . "\n";
        die "something is wrong.\n";
    }
    $_ = <IN>;
    s/[\r\n]//g;

    if ($pos =~ /79919006/) {
        $DB::single = 1;
    }

    my $key = uc($_);
    next if ($key =~ /N/);

    my $fSeq = "";
    if (substr($key, 2, 1) eq "A" or substr($key, 2, 1) eq "G") {
        $fSeq = complementSeq($key);
        if ($dir eq "-") {
            $dir = "+";
        } else {
            $dir = "-";
        }
    } else {
        $fSeq = $key;
    }

    my @mutFeatures = ();
   
    if ($numBase == 3) {
        push @mutFeatures, $base2ind{substr($fSeq, 1, 1)};
        push @mutFeatures, $base2ind{substr($fSeq, 3, 1)};
    } elsif ($numBase == 5) {
        push @mutFeatures, $base2ind{substr($fSeq, 0, 1)};
        push @mutFeatures, $base2ind{substr($fSeq, 1, 1)};
        push @mutFeatures, $base2ind{substr($fSeq, 3, 1)};
        push @mutFeatures, $base2ind{substr($fSeq, 4, 1)};
    }

    if ($dirOrNot == 1) {
        push @mutFeatures, $dir2ind{$dir};
    }


    if (substr($fSeq, 2, 1) eq "C") {
        $mutID2count{join(",", (1, @mutFeatures))} = $mutID2count{join(",", (1, @mutFeatures))} + 1;
        $mutID2count{join(",", (2, @mutFeatures))} = $mutID2count{join(",", (2, @mutFeatures))} + 1;
        $mutID2count{join(",", (3, @mutFeatures))} = $mutID2count{join(",", (3, @mutFeatures))} + 1;
    } elsif (substr($fSeq, 2, 1) eq "T") { 
        $mutID2count{join(",", (4, @mutFeatures))} = $mutID2count{join(",", (4, @mutFeatures))} + 1;        
        $mutID2count{join(",", (5, @mutFeatures))} = $mutID2count{join(",", (5, @mutFeatures))} + 1;
        $mutID2count{join(",", (6, @mutFeatures))} = $mutID2count{join(",", (6, @mutFeatures))} + 1;
    } else {
        print $pos . "\t" . $key . "\t" . $fSeq . "\t" . $dir . "\n";
        die "center base is neither C nor T\n";
    } 

    $mutSum = $mutSum + 3;
}

foreach my $mutID (sort keys %mutID2count) {
    print $mutID . "\t" . sprintf("%.10f", $mutID2count{$mutID} / $mutSum) . "\n";
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


