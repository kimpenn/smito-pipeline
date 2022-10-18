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
## v0.3

## Remove indel markups from samtools mpileup READ BASE column so that we can count
## the frequency of each base variant (substitution or deletion). Insertion will be 
## treated separately by `get_mpileup_ins.pl`.
##
## chrM	9452	G	5	.,+2TGA.a	EEEEE                         becomes
##                       EE    EEE
## chrM	9452	G	5	.,A.a	    EEEEE
##
## chrM	12798	G	15	.,.-1A,,,,....,,..	EEEEEEEEEE/AEEA       becomes
##                      EEE   EEEEEEE/AEEA
## chrM	12798	G	15	.,.,,,,....,,..	EEEEEEEEEE/AEEA 
##
## chrM	12799	A	15	.,*,,,,....,,..	/E!EEEEEEEEEEEE
##                      /E!EEEEEEEEEEEE
##                        ^
## Note, the removal of the deletion marker (e.g. "-1A") won't cause data loss as it 
## can be also found in the next base position (i.e. `*`), therefore we still have 
## access to its frequency. 

## Note, the deletion's phred score is just a placeholder. It is confirmed by 
## samtools developers:
## https://github.com/samtools/samtools/issues/1468
## 
## $ samtools view test.bam 
## 1       0       chr1    10      255     2S4M1S  *       0       0       AGTGGTA EDBA@?>
## 2       0       chr1    50      60      1M2I3M  *       0       0       CGTGCC  IHGFEC
## 3       0       chr1    90      30      3M2D1M  *       0       0       GACT    DCBA
## 
## $ samtools mpileup test.bam 
## [mpileup] 1 samples in 1 input files
## chr1    10      N       1       ^~T     B
## chr1    11      N       1       G       A
## chr1    12      N       1       G       @
## chr1    13      N       1       T$      ?
## chr1    50      N       1       ^]C+2GT I
## chr1    51      N       1       G       F
## chr1    52      N       1       C       E
## chr1    53      N       1       C$      C
## chr1    90      N       1       ^?G     D
## chr1    91      N       1       A       C
## chr1    92      N       1       C-2NN   B
## chr1    93      N       1       *       A
## chr1    94      N       1       *       A
## chr1    95      N       1       T$      A
 
## My puzzle is why the two deletions (position 93-94) got score A but not B, if we 
## consider the padded lines as associated with position 92 (the origin of the deletions) 
## according to https://www.biostars.org/p/61965.
###########################################################################
use strict;
use warnings;

sub rm_indel_markup {
    my $x = shift;
    while ( $x =~ /[+-](\d+)/) {
        my $n = $1;
        $x =~ s/[+-]$n[ACGTNacgtn]{$n}//;
    }
    $x;
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
    $F[4] =~ s/\^.//g; # the ASCII of the character following `^' minus 33 gives the mapping quality
    $F[4] =~ s/\$//g;
    next if $F[3] == 0;
    
    if ( $F[4] =~ /\d/) {
        $F[4] = &rm_indel_markup($F[4]);
    }
    print join("\t", @F);
}

if ( defined($infile) && -e $infile ) {
    close $fh;
}
