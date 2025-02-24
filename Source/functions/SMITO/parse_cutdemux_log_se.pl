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
## v0.1

## Parse CUTDEMUX LOG and generate a table tallying the full and partial matches
## For single-end data

print "Total\tSeq1\tSeq2\tR1\tR2\tR1_N9\tR1_N10\tR1_N11\tR1_N12\tR1_N13\tR2_N9\tR2_N10\tR2_N11\tR2_N12\tR2_N13\n";
$done = 0;
$begin = 0;
$end = 1;
while (<>) {
    chomp;
    if (!$done) {
        if (/-g X([CGAT]+)/) { ($seq1, $seq2) = ($1, $2); }
        if (/Total reads processed: \s+ ([0-9,]+)$/) { 
            $tot = $1; 
            $tot =~ s/,//g;
            print "$tot\t$seq1\t$seq2\t";
        }
        if (/Reads with adapters:\s+([0-9,]+) \(/) { 
            $R1 = $1; 
            $R1 =~ s/,//g; 
            $R2 = "";
            print "$R1\t$R2";
            $done = 1;
        }
    } else {
        if (/^=== /) {
            $begin = 1;
            $end = 0;
            ($n9, $n10, $n11, $n12, $n13) = (0, 0, 0, 0, 0);
        }
        if ($begin && /^9\s+(\d+)\s+/) { $n9 = $1; $end = 1; }
        if ($begin && /^10\s+(\d+)\s+/) { $n10 = $1; $end = 1; }
        if ($begin && /^11\s+(\d+)\s+/) { $n11 = $1; $end = 1; }
        if ($begin && /^12\s+(\d+)\s+/) { $n12 = $1; $end = 1; }
        if ($begin && /^13\s+(\d+)\s+/) { $n13 = $1; $end = 1; }
        if ($begin && $end && /^$/) {
            $begin = 0;
            print "\t$n9\t$n10\t$n11\t$n12\t$n13";
        }
    }
}
if ($begin) {
    $begin = 0;
    print "\t$n9\t$n10\t$n11\t$n12\t$n13";
}
print "\n";
