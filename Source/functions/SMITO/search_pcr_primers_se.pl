###########################################################################
## Author: Youtao Lu <luyoutao@sas.upenn.edu>
## 
## Copyright (c) 2021-2023, Youtao Lu and Junhyong Kim, Department of Biology, University of Pennsylvania
## Copyright (c) 2021-2023, James Eberwine, Perelman School of Medicine, University of Pennsylvania
## Copyright (c) 2021-2023, David Issadore, Penn Engineering, University of Pennsylvania
## 
## This program is free software; you can redistribute it and/or modify it under the terms of the Artistic License (2.0). You may obtain a copy of the full license at: 
## https://opensource.org/license/artistic-2-0/
###########################################################################
## v0.5

## Search the PCR forward and reverse primers for SNV1-12 successively for each mito
## for single-end data only

use strict;
use warnings;
use IO::Zlib;
use File::Path qw( make_path );
use File::Copy;
use Cwd qw( abs_path );
use String::Approx qw( amatch aindex );
use POSIX qw ( strftime );

my $fh;
my %loci;

my $inDir = $ARGV[0]; #' eg. Data/E.smito/Sample_L18R29P2/cutdemux/BC_M1
my $inFilePrefix = $ARGV[1]; #' eg. unaligned 

print STDERR "[" . strftime("%Y-%m-%d %H:%M:%S %z", localtime) . "] Starting...\n";

open $fh, "<", "Data/snv_loci_v2.csv" or die "Cannot open SNV loci file!";
while (<$fh>) {
    chomp;
    next if $. == 1;
    my @f = split(',', $_);
    $loci{$f[0]}{fw} = $f[11];
    $loci{$f[0]}{rv} = $f[12];
}
close $fh;

my $inFile1 = "$inDir/${inFilePrefix}_1.fq.gz";

my $input1 = "$inDir/input_1.fq.gz"; #' the file to search
unlink $input1;
symlink abs_path($inFile1), $input1 or die "Cannot symlink $inFile1 $input1: $!";

sub nextReadPair {
    my ($inFh1, $inFh2) = @_;
    my $rname1  = readline($inFh1);
    if ( !defined($rname1) ) { return () }
    my $seq1   = readline($inFh1);
    my $sep1   = readline($inFh1);
    my $phred1 = readline($inFh1);
    my $rname2 = undef;
    my $seq2   = undef;
    my $sep2   = undef;
    my $phred2 = undef;
    return ( { rname => $rname1, seq => $seq1, sep => $sep1, phred => $phred1 }, 
             { rname => $rname2, seq => $seq2, sep => $sep2, phred => $phred2 } );
}

sub editRName0 {
    my ($rname, $patt, $ind) = @_;
    my @f = split(' ', $rname); 
    return "$f[0]::$patt~$ind $f[1]\n";
}

for my $i (1..12) {
    my $snv = "SNV$i";
    for my $end ("fw", "rv") {
        my $seq = $loci{$snv}{$end};
        print STDERR "[" . strftime("%Y-%m-%d %H:%M:%S %z", localtime) . "] " . "$snv\t$end\t$seq\n";
        my ( $inputFh1, $inputFh2 );
        $inputFh1 = IO::Zlib->new( $input1, "r" ) or die "Cannot read $input1!";
        my $outDir = "$inDir/${snv}_${end}_1";
        make_path $outDir or warn "Failed creating $outDir: $!";
        my $outfound1 = "$outDir/unaligned_1.fq.gz";
        my $outnotfound1 = "$inDir/notfound_1.fq.gz";
        my $outfoundFh1 = IO::Zlib->new( $outfound1, "w" ) or die "Cannot open $outfound1: $!";
        my $outnotfoundFh1 = IO::Zlib->new( $outnotfound1, "w" ) or die "Cannot open $outnotfound1: $!";
        while (1) {
            my @rp = &nextReadPair($inputFh1, undef);
            last if $#rp == -1;
            my ($r1, $r2) = @rp;
            my $ind = -1;
            $ind = index($r1->{seq}, $seq); # if exact match works, skip approximate match
            $ind = aindex($seq, [ "I1", "D1", "S1" ], $r1->{seq}) if $ind == -1;
            if ($ind == -1) { # not found
                print $outnotfoundFh1 $r1->{rname}, $r1->{seq}, $r1->{sep}, $r1->{phred};
            } else {
                $r1->{rname} = &editRName0($r1->{rname}, "R1_${snv}_${end}", $ind);
                print $outfoundFh1 $r1->{rname}, $r1->{seq}, $r1->{sep}, $r1->{phred};
            }
        }
        close $outfoundFh1 or warn "Cannot close $outfoundFh1: $!";
        close $outnotfoundFh1 or warn "Cannot close $outnotfoundFh1: $!";
        close $inputFh1 or warn "Cannot close $inputFh1: $!";
#' keep the last round output as "notfound_[12].fq.gz"
        unless ( $i == 12 && $end eq "rv" ) {
            move $outnotfound1, $input1 or die "Cannot mv $outnotfound1 $input1: $!";
        } else {
            unlink $input1 or warn "Cannot remove $input1: $!"; 
        }
    }
}

print STDERR "[" . strftime("%Y-%m-%d %H:%M:%S %z", localtime) . "] Done.\n";
