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
## v0.2

if (!exists("SMITO") || is.environment(SMITO)) SMITO <- new.env(parent = emptyenv())
local({


    ## For each line of the `samtools mpileup` output, parse the base composition and 
    ## the Phred score in pairs. The output is a long table. 
    ## Input: 
    ## chrM	9534	C	6	,.....	B/>A6>	~~~~~~	136,122,120,134,133,67	NB501328:332:HGL3VBGXJ:3:23404:1861:7397,NB501328:332:HGL3VBGXJ:3:22504:3980:13205,NB501328:332:HGL3VBGXJ:3:12410:9299:1395,NB501328:332:HGL3VBGXJ:3:11411:6750:13356,NB501328:332:HGL3VBGXJ:3:23606:12446:9361,NB501328:332:HGL3VBGXJ:1:13112:10451:13828	16,0,0,0,0,0	9399,9415,9415,9433,9442,9491

    ## Output: 
    ##   chrom  pos ref depth obs phred
    ## 1  chrM 9534   C     6   ,    33
    ## 2  chrM 9534   C     6   .    14
    ## 3  chrM 9534   C     6   .    29
    ## 4  chrM 9534   C     6   .    32
    ## 5  chrM 9534   C     6   .    21
    ## 6  chrM 9534   C     6   .    29
    parse_sam_mpileup_perbase <- function(X) {
        bases <- X[,5]
        bases <- strsplit(bases, "")[[1]]
        phreds <- X[,6]
        phreds <- as.integer(charToRaw(phreds)) - 33
        data.frame(
            chrom = X[,1], 
            pos = X[,2], 
            ref = X[,3], 
            depth = X[,4], 
            obs = bases, 
            phred = phreds,
            stringsAsFactors = FALSE
        )
    }

    ## Given `samtools mpileup` output, parse the base composition and the Phred score
    ## in pairs. The input is the output of rm_mpileup_indel.pl and the output is a 
    ## long, VCF-like table.
    parse_sam_mpileup <- function(X) {
        X <- subset(X, V4 > 0)
        if (nrow(X) == 0) {
            return(data.frame(
                chrom = character(0), 
                pos = integer(0), 
                ref = character(0), 
                depth = integer(0), 
                obs = character(0), 
                phred = integer(0),
                stringsAsFactors = FALSE
            ))
        }
        ret <- lapply(1:nrow(X), function(i) {
            parse_sam_mpileup_perbase(X[i, ])
        })
        do.call(rbind, ret)
    }

    for (obj in ls()) { assign(obj, get(obj), envir = SMITO) }
})
