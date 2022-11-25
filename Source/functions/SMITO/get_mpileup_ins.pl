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
## v0.2

## Extract insertion from samtools mpileup READ BASE column
## Input:
## chrM	12798	G	15	.,.-1A,,,,....,,..	EEEEEEEEEE/AEEA
##                       EEE   EEEEEEE/AEEA
## chrM	12798	G	15	.,.,,,,....,,..	EEEEEEEEEE/AEEA 
## chrM	9452	G	5	.,+2TGA.a+2TG..+3TTT	EEEEEAA
##                      EE    EEE    AA
## chrM	9452	G	5	.,A.a	EEEEE
##                      EEEEE
## chrM	12799	A	15	.,*,,,,....,,..	/E!EEEEEEEEEEEE
##                      /E!EEEEEEEEEEEE
## Output: (column 6 is the respective frequency of 5, separated by comma; column 7 total insertions)
## chrM 9452    G   GTG,GTTT    2,1    3
## 
## Example:
## perl Source/functions/SMITO/get_mpileup_ins.pl <(gzcat < Data/E.smito/analyzed/Sample_L2R53P1_b/Sample_L2R53P1_b_M1/star/SNV1/Sample_L2R53P1_b_M1.star.unique.mpileup.sub500k.tsv.gz)
## chrM    9360    C       12051   CC      4       4
## chrM    9362    A       12054   AA      1       1
## chrM    9366    G       12061   GG      3       3
## chrM    9372    A       12065   AA,AG   1,1     2
## chrM    9376    A       12066   AG      3       3


## Note, since `samtools mpileup` use lowercase and uppercase to distinguish 
## the strandedness (ACGT means the forward strand, acgt means the reverse
## strand), here we should first make all inserted bases to uppercase before
## we count the insertion frequency. 
###########################################################################
use strict;
use warnings;
use List::Util qw( reduce );

sub parse_ins {
    my $x = shift;
    my %h;
    while ( $x =~ /\+(\d+)/) {
        my $n = $1;
        $x =~ s/\+$n([ACGTNacgtn]{$n})//;
        my $k = $1;
        ++$h{uc($k)};
    }
    return \%h;
}

my $infile = $ARGV[0];
my $fh;
if ( defined($infile) && -e $infile ) {
    open $fh, "<", $infile or die("Cannot open $infile!$!");
} else {
    $fh = *STDIN;
}

while ( <$fh> ) {
    my @F = split("\t", $_);
    next if $F[3] == 0;
    
    if ( $F[4] =~ /\+/ ) {
        my $h = &parse_ins($F[4]);
        my @k = sort keys %{$h};
        print $F[0], "\t", $F[1], "\t", uc($F[2]), "\t", $F[3], "\t", 
              join(',', map { uc($F[2]) . $_ } @k), "\t", 
              join(',', map { ${$h}{$_} } @k), "\t", 
              (reduce { $a + $b } (0, values %{$h})),             
              "\n";
    }
}

if ( defined($infile) && -e $infile ) {
    close $fh;
}
