###########################################################################
## Author: Youtao Lu <luyoutao@sas.upenn.edu>
## 
## Copyright (c) 2021-2023, Youtao Lu and Junhyong Kim, Department of Biology, University of Pennsylvania
## Copyright (c) 2021-2023, James Eberwine, Perelman School of Medicine, University of Pennsylvania
## Copyright (c) 2021-2023, David Issadore, Penn Engineering, University of Pennsylvania
## 
## This program is free software; you can redistribute it and/or modify it under the terms of the the Artistic License (2.0). You may obtain a copy of the full license at: 
## https://opensource.org/license/artistic-2-0/
###########################################################################
## v0.6
source("Source/functions/SMITO.R")
library("parallel")

mito_barcodes <- read.csv("Data/mito_barcodes.csv", as.is = TRUE)
mitoIDs <- mito_barcodes[, "ID"] 

snv_info <- read.csv("Data/snv_loci_v2.csv", as.is = TRUE)
snvIDs <- snv_info[, "SNVID"]

LibraryInfo <- read.csv("Data/LibraryInfo.csv", as.is = TRUE, check.names = FALSE)
rownames(LibraryInfo) <- LibraryInfo[, "LibraryID"]
exptIDs <- sort(unique(LibraryInfo[, "ExptID"]))

ncores <- 10

for (exptID in exptIDs) {
    libraryIDs <- subset(LibraryInfo, ExptID == exptID)[, "LibraryID"]

    ## 1. Process substitutions and deletions. 
    ## 1.1 Load the markup fixed data after running `rm_mpileup_indel.pl`.
    mpileups_cutdemux <- lapply(libraryIDs, function(x) {
        x1 <- paste0("Sample_", x)
        X <- mclapply(mitoIDs, function(y) {
            sapply(snvIDs, function(z) {
                f <- sprintf("Data/E.smito/analyzed/%s/%s_%s/star/%s/%s_%s.star.unique.mpileup.sub500k.noindel.tsv.gz", x1, x1, y, z, x1, y)
                message(f)
                tryCatch(read.csv(f, as.is = TRUE, head = FALSE, sep = "\t", quote = ""), error = function(e) {
                    warning("empty file: ", f, "!")
                    data.frame(V1=character(0), V2=integer(0), V3=character(0), V4=character(0), V5=character(0), V6=character(0), V7=character(0), V8=character(0), V9=character(0), V10=character(0), V11=character(0), stringsAsFactors = FALSE) 
                })
            }, simplify = FALSE)
        }, mc.cores = ncores)
        names(X) <- mitoIDs
        X
    })
    names(mpileups_cutdemux) <- libraryIDs

    dir.create("Data/SNVs", FALSE, TRUE)
    saveRDS(mpileups_cutdemux, file = sprintf("Data/SNVs/mpileups_cutdemux_E%s_sub500k.RDS", exptID))
    mpileups_cutdemux <- readRDS(file = sprintf("Data/SNVs/mpileups_cutdemux_E%_sub500k.RDS", exptID))
        
    ## 1.2 Calculate each base's frequency for each position, strand aware.
    basedifflevels <- c(".", ",", "A", "a", "C", "c", "G", "g", "T", "t", "N", "n", "*", "#")
    basedifffreq_cutdemux_q30 <- sapply(libraryIDs, function(x) {
        X <- mclapply(mitoIDs, function(y) {
            Y <- mclapply(snvIDs, function(z) {
                message(x, " ", y, " ", z)
                mpileup <- mpileups_cutdemux[[x]][[y]][[z]]
                basediff <- SMITO$parse_sam_mpileup(mpileup)
                if (nrow(basediff) > 0) { ## deletions don't have phred quality
                    basediff[basediff[, "obs"] == "*" | basediff[, "obs"] == "#", "phred"] <- NA
                    basediff <- subset(basediff, phred >= 30 | is.na(phred)) # no quality filtering on deletion
                    if (nrow(basediff) > 0) {
                        basediff$obs <- factor(basediff$obs, levels = basedifflevels)
                        basedifffreq <- with(basediff, tapply(obs, pos, table, simplify = FALSE))
                        refbase <- with(basediff, tapply(ref, pos, '[', 1))
                        basedifffreq <- do.call(rbind, basedifffreq)
                        basedifffreq <- data.frame(pos = as.integer(rownames(basedifffreq)), ref = refbase, basedifffreq, stringsAsFactors = FALSE, row.names = NULL, check.names = FALSE)
                    } else {
                        basedifffreq <- data.frame(pos = integer(0), ref = character(0), sapply(basedifflevels, function(i) integer(0)), stringsAsFactors = FALSE, row.names = NULL, check.names = FALSE)
                    }
                } else {
                    basedifffreq <- data.frame(pos = integer(0), ref = character(0), sapply(basedifflevels, function(i) integer(0)), stringsAsFactors = FALSE, row.names = NULL, check.names = FALSE)
                }
                basedifffreq
            }, mc.cores = 10)
            names(Y) <- snvIDs
            Y
        }, mc.cores = 1)
        names(X) <- mitoIDs
        X
    }, simplify = FALSE)

    dir.create("Data/SNVs", FALSE, TRUE)
    saveRDS(basedifffreq_cutdemux_q30, file = sprintf("Data/SNVs/basedifffreq_cutdemux_E%s_sub500k_q30.RDS", exptID))
    basedifffreq_cutdemux_q30 <- readRDS(file = sprintf("Data/SNVs/basedifffreq_cutdemux_E%s_sub500k_q30.RDS", exptID))

    basedifffreq_cutdemux_q30_df <- do.call(rbind, lapply(libraryIDs, function(x) {
        Z <- do.call(rbind, lapply(mitoIDs, function(y) {
            Y <- do.call(rbind, lapply(snvIDs, function(z) {
                message(x, " ", y, " ", z)
                X <- basedifffreq_cutdemux_q30[[x]][[y]][[z]]
                if (nrow(X) > 0) {
                    data.frame(SNVID = z, X, stringsAsFactors = FALSE, check.names = FALSE)
                } else {
                    data.frame(SNVID = character(0), pos = integer(0), ref = character(0), sapply(basedifflevels, function(i) integer(0)), stringsAsFactors = FALSE, row.names = NULL, check.names = FALSE)
                }
            }))
            if (nrow(Y) > 0) {
                data.frame(MitoID = y, Y, stringsAsFactors = FALSE, check.names = FALSE)
            } else {
                data.frame(MitoID = character(0), SNVID = character(0), pos = integer(0), ref = character(0), sapply(basedifflevels, function(i) integer(0)), stringsAsFactors = FALSE, row.names = NULL, check.names = FALSE)
            }
        }))
        if (nrow(Z) > 0) {
            data.frame(ExptID = LibraryInfo[x, "ExptID"], LibraryID = x, Z, stringsAsFactors = FALSE, check.names = FALSE)
        } else {
            data.frame(ExptID = character(0), LibraryID = character(0), MitoID = character(0), SNVID = character(0), pos = integer(0), ref = character(0), sapply(basedifflevels, function(i) integer(0)), stringsAsFactors = FALSE, row.names = NULL, check.names = FALSE)
        }
    }))

    dir.create("Report/SNVs", FALSE, TRUE)
    write.csv(basedifffreq_cutdemux_q30_df, file = gzfile(sprintf("Report/SNVs/basedifffreq_cutdemux_E%s_sub500k_q30.csv.gz", exptID)), row.names = FALSE)
    basedifffreq_cutdemux_q30_df <- read.csv(sprintf("Report/SNVs/basedifffreq_cutdemux_E%s_sub500k_q30.csv.gz", exptID), as.is = TRUE, check.names = FALSE)

    ## 1.3 Calculate each base's frequency for each position, strand merged
    basedifffreq_cutdemux_q30_unstranded <- data.frame(
        basedifffreq_cutdemux_q30_df[, 1:6],
        "=" = basedifffreq_cutdemux_q30_df[, "."] + basedifffreq_cutdemux_q30_df[, ","],
        "A" = basedifffreq_cutdemux_q30_df[, "A"] + basedifffreq_cutdemux_q30_df[, "a"],
        "C" = basedifffreq_cutdemux_q30_df[, "C"] + basedifffreq_cutdemux_q30_df[, "c"],
        "G" = basedifffreq_cutdemux_q30_df[, "G"] + basedifffreq_cutdemux_q30_df[, "g"],
        "T" = basedifffreq_cutdemux_q30_df[, "T"] + basedifffreq_cutdemux_q30_df[, "t"],
        "D" = basedifffreq_cutdemux_q30_df[, "*"] + basedifffreq_cutdemux_q30_df[, "#"],
        stringsAsFactors = FALSE, check.names = FALSE
    )
    basedifffreq_cutdemux_q30_unstranded$depth <- rowSums(basedifffreq_cutdemux_q30_unstranded[, 7:12])
    basedifffreq_cutdemux_q30_unstranded <- basedifffreq_cutdemux_q30_unstranded[, c(1:6, 13, 7:12)]
    write.csv(basedifffreq_cutdemux_q30_unstranded, file = gzfile(sprintf("Report/SNVs/basedifffreq_cutdemux_E%s_sub500k_q30_unstranded.csv.gz", exptID)), row.names = FALSE)
    basedifffreq_cutdemux_q30_unstranded <- read.csv(sprintf("Report/SNVs/basedifffreq_cutdemux_E%s_sub500k_q30_unstranded.csv.gz", exptID), as.is = TRUE, check.names = FALSE)

    ## 2. Process insertions.
    mpileups_ins_cutdemux <- lapply(libraryIDs, function(x) {
        x1 <- paste0("Sample_", x)
        X <- lapply(mitoIDs, function(y) {
            sapply(snvIDs, function(z) {
                f <- sprintf("Data/E.smito/analyzed/%s/%s_%s/star/%s/%s_%s.star.unique.mpileup.sub500k.ins.tsv.gz", x1, x1, y, z, x1, y)
                message(f)
                tryCatch(read.csv(f, as.is = TRUE, head = FALSE, sep = "\t", quote = "", colClasses = c("character", "integer", "character", "integer", "character", "character", "integer")), error = function(e) 
                    data.frame(V1=character(0), V2=integer(0), V3=character(0), V4=integer(0), V5=character(0), V6=character(0), V7=integer(0), stringsAsFactors = FALSE))
            }, simplify = FALSE)
        })
        names(X) <- mitoIDs
        X
    })
    names(mpileups_ins_cutdemux) <- libraryIDs
    saveRDS(mpileups_ins_cutdemux, file = sprintf("Data/SNVs/mpileups_ins_cutdemux_E%s_sub500k.RDS", exptID))
    mpileups_ins_cutdemux <- readRDS(file = sprintf("Data/SNVs/mpileups_ins_cutdemux_E%s_sub500k.RDS", exptID))

    mpileups_ins_cutdemux <- lapply(libraryIDs, function(x) {
        X <- lapply(mitoIDs, function(y) {
            Y <- lapply(snvIDs, function(z) {
                message(x, " ", y, " ", z)
                mpileup <- mpileups_ins_cutdemux[[x]][[y]][[z]]
                if (nrow(mpileup) > 0) { ## insersions don't have phred quality
                    return(with(mpileup, data.frame(ExptID = exptID, LibraryID = x, MitoID = y, SNVID = z, pos = V2, ref = V3, depth = V4, alt = V5, altfreq = V6, altfreqtot = V7, stringsAsFactors = FALSE, row.names = NULL, check.names = FALSE)))
                } else {
                    return(data.frame(ExptID = character(0), LibraryID = character(0), MitoID = character(0), SNVID = character(0), pos = integer(0), ref = character(0), depth = integer(0), alt = character(0), altfreq = character(0), altfreqtot = integer(0), stringsAsFactors = FALSE, row.names = NULL, check.names = FALSE))
                }
            })
            do.call(rbind, Y)
        })
        do.call(rbind, X)
    })
    mpileups_ins_cutdemux <- do.call(rbind, mpileups_ins_cutdemux)

    dir.create("Report/SNVs", FALSE, TRUE)
    write.csv(mpileups_ins_cutdemux, file = gzfile(sprintf("Report/SNVs/ins/mpileups_ins_cutdemux_E%s_sub500k.csv.gz", exptID)), row.names = FALSE)
}
