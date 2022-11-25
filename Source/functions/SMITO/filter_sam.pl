###########################################################################
## Author: Youtao Lu <luyoutao@sas.upenn.edu>
## 
## Copyright (c) 2021, Laboratory of Junhyong Kim, Department of Biology, University of Pennsylvania
## Copyright (c) 2021, Laboratory of James Eberwine, Perelman School of Medicine, University of Pennsylvania
## Copyright (c) 2021, Laboratory of David Issadore, Penn Engineering, University of Pennsylvania
## All Rights Reserved.
## 
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
## Filter the BAM file reads for reads (1) that are mapped to the positive strand,
## (2) that start at approximately the same position as the 5' end of the PCR forward
## primer, (3) with length above certain threshold.

use strict;
use warnings;
use POSIX qw ( strftime );
use Getopt::Long;
use Inline Config => DIRECTORY => $ENV{HOME};
use Inline 'C'    => <<CODE;
/********************************************************************************* 
*** How to pass Perl reference to C routine is very tricky, see:
*** https://stackoverflow.com/questions/18291827/perl-inline-c-passing-arrayref-to-c-function 
*** Most baffling is that if the C routine is written wrong, the error will be like
*** Undefined subroutine &main::is_pos_ok called at filter_sam.pl line 138, <infh> line 10.
***********************************************************************************/
int is_pos_ok(int p, AV *s, int k, int u, int d) {
    int i = 0, j = k;
    if (p <= SvNV(*av_fetch(s, 0, 0))) {
        return SvNV(*av_fetch(s, 0, 0)) - p <= u ? 1 : 0;
    }
    if (p >= SvNV(*av_fetch(s, k, 0))) {
        return p - SvNV(*av_fetch(s, k, 0)) <= d ? 1 : 0;
    } 

    while (i < j - 1) {
        int m = (int) (i + j) / 2; 
        if (p < SvNV(*av_fetch(s, m, 0))) { 
            j = m;
            continue;
        } else if (p > SvNV(*av_fetch(s, m, 0))) {
            i = m;
            continue;
        }
        return 1;
    }
    return p - SvNV(*av_fetch(s, i, 0)) <= d || SvNV(*av_fetch(s, j, 0)) - p <= u ? 1 : 0;
}
CODE

my ( $infile,  $infh );
my ( $outfile, $outfh );
my ( $snvfile, $snvfh ) = ( "Data/snv_loci_v2.csv", undef );
my ( $offup, $offdn )   = ( 0, 0 );
my $minlen = 135;
my $increv = 0;
my @start  = ();
my ( $n, $i ) = ( 0, 0 );
my $help    = 0;
my $version = 0;
our $VERSION = "0.3";

sub usage {
    print <<DOC;
Summary:
    Filter for reads (1) that are mapped to the positive strand, (2) that start at
    approximately the same position as the 5' end of the PCR forward primer, (3) 
    with mapped length above certain threshold. 

Usage:
    perl $0 --infile unfiltered.bam --outfile filtered.bam [--snvfile Data/snv_loci_v2.csv] [--offup 0] [--offdn 0] [minlen 135] [--increv]

Options:
    --infile    input BAM/SAM
    --outfile   output BAM/SAM
    --snvfile   SNV list file (default: "Data/snv_loci_v2.csv")
    --offup     max. bases (>=0) upstream the expected start to allow (default: 0) 
    --offdn     max. bases (>=0) downstream the expected start to allow (default: 0) 
    --minlen    min. bases of the mapped length (default: 135) 
    --increv    whether to include reads mapped to the reverse strand
DOC
}

GetOptions(
    "infile|i=s"  => \$infile,
    "outfile|o=s" => \$outfile,
    "snvfile=s"   => \$snvfile,
    "offup=i"     => \$offup,
    "offdn=i"     => \$offdn,
    "minlen=i"    => \$minlen,
    "increv"      => \$increv,
    "help"        => \$help,
    "version"     => \$version,
) or &usage() && exit(-1);

&usage                     && exit(0) if $help;
( print "$0 v$VERSION\n" ) && exit(0) if $version;

die "--infile missing!\n"         if !defined($infile);
die "--outfile missing!\n"        if !defined($outfile);
die "--snvfile missing!\n"        if !defined($snvfile);
die "--offup cannot be negative!" if $offup < 0;
die "--offdn cannot be negative!" if $offdn < 0;

print STDERR "["
  . strftime( "%Y-%m-%d %H:%M:%S %z", localtime )
  . "] { infile = $infile, outfile = $outfile, snvfile = $snvfile, offup = $offup, $offdn = $offdn, minlen = $minlen, increv = $increv, VERSION = $VERSION }\n";

open $snvfh, "<", $snvfile or die "Cannot open $snvfile for read!\n";
while (<$snvfh>) {
    next if $. == 1;
    my @F = split( ',', $_ );
    push @start, int( $F[6] );
}
@start = sort { $a <=> $b } @start;
close $snvfh;

print STDERR "["
  . strftime( "%Y-%m-%d %H:%M:%S %z", localtime )
  . "] Loci found: @start.\n";

print STDERR "["
  . strftime( "%Y-%m-%d %H:%M:%S %z", localtime )
  . "] Filter begins.\n";

open $infh, "samtools view -h $infile | "
  or die "Cannot open $infile for read!\n";
open $outfh, "| samtools view -b > $outfile"
  or die "Cannot open $outfile for write!\n";

while (<$infh>) {
    ( print $outfh $_ ) && next if /^@/;
    $n++;
    my @F = split( "\t", $_ );
    my ( $flag, $pos, $cigar ) = @F[ ( 1, 3, 5 ) ];
    next if !$increv && $flag & 0x10;
    next if !&is_pos_ok( $pos, \@start, $#start, $offup, $offdn );
    my $m = 0;
    while ( $cigar =~ s/(\d+)M// ) {
        $m += $1;
    }
    next if $m < $minlen;
    print $outfh $_;
    $i++;
}
close $outfh;
close $infh;

print STDERR "[" . strftime( "%Y-%m-%d %H:%M:%S %z", localtime ) . "] Done.\n";
print STDERR "["
  . strftime( "%Y-%m-%d %H:%M:%S %z", localtime )
  . "] $n reads checked, $i reads passed\n";

