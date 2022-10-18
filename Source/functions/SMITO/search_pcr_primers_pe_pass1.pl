###########################################################################
## Youtao Lu <luyoutao@sas.upenn.edu>
## Copyright (c) 2021, Youtao Lu and Junhyong Kim, Department of Biology, University of Pennsylvania
## Copyright (c) 2021, James Eberwine, Perelman School of Medicine, University of Pennsylvania
## Copyright (c) 2021, David Issadore, Penn Engineering, University of Pennsylvania
## All Rights Reserved.

## You may not use this file except in compliance with the Kim Lab License
## located at
##
##     http://kim.bio.upenn.edu/software/LICENSE
##
## Unless required by applicable law or agreed to in writing, this
## software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
## CONDITIONS OF ANY KIND, either express or implied.  See the License
## for the specific language governing permissions and limitations
## under the License.
###########################################################################
## v0.7

## Search the PCR forward and reverse primers for SNV1-12 successively for each mito
## for paired-end seqeuencing, pass one (searching R1)

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
my $inFile2 = "$inDir/${inFilePrefix}_2.fq.gz";

my $input1 = "$inDir/input_1.fq.gz"; #' the file to search
my $input2 = "$inDir/input_2.fq.gz";
unlink $input1;
unlink $input2;
symlink abs_path($inFile1), $input1 or die "Cannot symlink $inFile1 $input1: $!";
symlink abs_path($inFile2), $input2 or die "Cannot symlink $inFile2 $input2: $!";

sub nextReadPair {
    my ($inFh1, $inFh2) = @_;
    my $rname1  = readline($inFh1);
    if ( !defined($rname1) ) { return () }
    my $seq1   = readline($inFh1);
    my $sep1   = readline($inFh1);
    my $phred1 = readline($inFh1);
    my $rname2 = readline($inFh2);
    my $seq2   = readline($inFh2);
    my $sep2   = readline($inFh2);
    my $phred2 = readline($inFh2);
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
        $inputFh2 = IO::Zlib->new( $input2, "r" ) or die "Cannot read $input2!";
        my $outDir = "$inDir/${snv}_${end}_1";
        make_path $outDir or warn "Failed creating $outDir: $!";
        my $outfound1 = "$outDir/unaligned_1.fq.gz";
        my $outfound2 = "$outDir/unaligned_2.fq.gz";
        my $outnotfound1 = "$inDir/notfound_1.fq.gz";
        my $outnotfound2 = "$inDir/notfound_2.fq.gz";
        my $outfoundFh1 = IO::Zlib->new( $outfound1, "w" ) or die "Cannot open $outfound1: $!";
        my $outfoundFh2 = IO::Zlib->new( $outfound2, "w" ) or die "Cannot open $outfound2: $!";
        my $outnotfoundFh1 = IO::Zlib->new( $outnotfound1, "w" ) or die "Cannot open $outnotfound1: $!";
        my $outnotfoundFh2 = IO::Zlib->new( $outnotfound2, "w" ) or die "Cannot open $outnotfound2: $!";
        while (1) {
            my @rp = &nextReadPair($inputFh1, $inputFh2);
            last if $#rp == -1;
            my ($r1, $r2) = @rp;
            my $ind = -1;
            $ind = index($r1->{seq}, $seq); # if exact match works, skip approximate match
            $ind = aindex($seq, [ "I1", "D1", "S1" ], $r1->{seq}) if $ind == -1;
            if ($ind == -1) { # not found
                print $outnotfoundFh1 $r1->{rname}, $r1->{seq}, $r1->{sep}, $r1->{phred};
                print $outnotfoundFh2 $r2->{rname}, $r2->{seq}, $r2->{sep}, $r2->{phred};
            } else {
                $r1->{rname} = &editRName0($r1->{rname}, "R1_${snv}_${end}", $ind);
                $r2->{rname} = &editRName0($r2->{rname}, "R1_${snv}_${end}", $ind);
                print $outfoundFh1 $r1->{rname}, $r1->{seq}, $r1->{sep}, $r1->{phred};
                print $outfoundFh2 $r2->{rname}, $r2->{seq}, $r2->{sep}, $r2->{phred};
            }
        }
        close $outfoundFh1 or warn "Cannot close $outfoundFh1: $!";
        close $outfoundFh2 or warn "Cannot close $outfoundFh2: $!";
        close $outnotfoundFh1 or warn "Cannot close $outnotfoundFh1: $!";
        close $outnotfoundFh2 or warn "Cannot close $outnotfoundFh2: $!";
        close $inputFh1 or warn "Cannot close $inputFh1: $!";
        close $inputFh2 or warn "Cannot close $inputFh2: $!";
        # keep the last round output as "notfound_[12].fq.gz"
        unless ( $i == 12 && $end eq "rv" ) {
            move $outnotfound1, $input1 or die "Cannot mv $outnotfound1 $input1: $!";
            move $outnotfound2, $input2 or die "Cannot mv $outnotfound2 $input2: $!";
        } else {
            unlink $input1 or warn "Cannot remove $input1: $!"; 
            unlink $input2 or warn "Cannot remove $input2: $!"; 
        }
    }
}

print STDERR "[" . strftime("%Y-%m-%d %H:%M:%S %z", localtime) . "] Done.\n";
