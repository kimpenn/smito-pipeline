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
source("Source/functions/SMITO.R")

if (!require("optparse")) { 
    install.packages("optparse")
    library("optparse") 
}

parser <- OptionParser()
parser <- add_option(parser, "--exptID", action = "store", type = "character", default = NULL, help = "Experiment run ID (e.g. '749'")
parser <- add_option(parser, "--libraryID", action = "store", type = "character", default = NULL, help = "Library ID (e.g. 'L1R45P1'")
parser <- add_option(parser, "--mitoID", action = "store", type = "character", default = NULL, help = "Mito ID (e.g. M1, M2, ..., M10)")
parser <- add_option(parser, "--qth", action = "store", type = "integer", default = 30, help = "Phred score cutoff (default: 30)")
parser <- add_option(parser, "--starbase", action = "store", type = "character", default = "unique", help = "Part of the STAR filename that specifies which filter was applied (default: unique)")
parser <- add_option(parser, "--cache_dir", action = "store", type = "character", default = "~/Workspace.cache", help = "Directory for temporary results on cluster (default: ~/Workspace.cache)")
parser <- add_option(parser, "--data_dir", action = "store", type = "character", default = "Data", help = "Data folder (default: 'Data')")
parser <- add_option(parser, "--expt_dir", action = "store", type = "character", default = "E.smito", help = "Pooled experiment folder (default: 'E.smito')")
parser <- add_option(parser, "--result_dir", action = "store", type = "character", default = "Report", help = "Result folder (default: 'Report')")
parser <- add_option(parser, "--version_dir", action = "store", type = "character", default = "20210510", help = "Version folder (default: '20210510')")


opts <- parse_args(parser, args = commandArgs(trailingOnly = TRUE))

is_online <- function(host = "\\.upenn\\.edu$") {
    grepl(host, system("hostname", intern = TRUE), ignore.case = TRUE)
}

parse_basediff_bymito <- function(exptID, libraryID, mitoID, snvIDs = paste0("SNV", 1:12), qth = 30, starbase = "unique", cache_dir = "~/Workspace.cache/SMITO", data_dir = "Data", expt_dir = "E.smito", result_dir = "Report", version_dir = "20210510") {

## 1.1 Load the markup fixed data after running `rm_mpileup_indel.pl`.
    libraryName <- paste0("Sample_", libraryID)
    libraryMitoID <- paste0(libraryID, "_", mitoID)
    libraryMitoName <- paste0("Sample_", libraryMitoID)
    mpileups_cutdemux <- sapply(snvIDs, function(snvID) {
        f <- sprintf("%s/%s/analyzed/%s/%s/star/%s/%s.star.%s.mpileup.sub500k.noindel.tsv.gz", data_dir, expt_dir, libraryName, libraryMitoName, snvID, libraryMitoName, starbase)
            message("Reading ", f)
            tryCatch(read.csv(f, as.is = TRUE, head = FALSE, sep = "\t", quote = ""), error = function(e) {
                warning("empty file: ", f, "!")
                data.frame(V1=character(0), V2=integer(0), V3=character(0), V4=character(0), V5=character(0), V6=character(0), V7=character(0), V8=character(0), V9=character(0), V10=character(0), V11=character(0), stringsAsFactors = FALSE) 
            })
    }, simplify = FALSE)

    data_dir_cache <- sprintf("%s/%s", cache_dir, data_dir)
    result_dir_cache <- sprintf("%s/%s", cache_dir, result_dir)

    datadir <- sprintf("%s/%s/SNVs", ifelse(is_online(), data_dir_cache, data_dir), version_dir)
    dir.create(datadir, FALSE, TRUE)
    saveRDS(mpileups_cutdemux, file = sprintf("%s/mpileups_cutdemux_%s_sub500k.RDS", datadir, libraryMitoID))
    mpileups_cutdemux <- readRDS(file = sprintf("%s/mpileups_cutdemux_%s_sub500k.RDS", datadir, libraryMitoID))
        
    ## 1.2 Calculate each base's frequency for each position, strand aware.
    basedifflevels <- c(".", ",", "A", "a", "C", "c", "G", "g", "T", "t", "N", "n", "*", "#")
    basedifffreq_cutdemux <- sapply(snvIDs, function(snvID) {
        message("Processing substitute/del: ", exptID, " ", libraryID, " ", mitoID, " ", snvID)
        mpileup <- mpileups_cutdemux[[snvID]]
        basediff <- SMITO$parse_sam_mpileup(mpileup)
        if (nrow(basediff) > 0) { ## deletions don't have phred quality
            basediff[basediff[, "obs"] == "*" | basediff[, "obs"] == "#", "phred"] <- NA
            basediff <- subset(basediff, phred >= qth | is.na(phred)) # no quality filtering on deletion
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
    }, simplify = FALSE)

    saveRDS(basedifffreq_cutdemux, file = sprintf("%s/basedifffreq_cutdemux_%s_sub500k_q%s.RDS", datadir, libraryMitoID, qth))
    basedifffreq_cutdemux <- readRDS(file = sprintf("%s/basedifffreq_cutdemux_%s_sub500k_q%s.RDS", datadir, libraryMitoID, qth))

    basedifffreq_cutdemux_df <- do.call(rbind, lapply(snvIDs, function(snvID) {
        message(exptID, " ", libraryID, " ", mitoID, " ", snvID)
        X <- basedifffreq_cutdemux[[snvID]]
        if (nrow(X) > 0) {
            data.frame(SNVID = snvID, X, stringsAsFactors = FALSE, check.names = FALSE)
        } else {
            data.frame(SNVID = character(0), pos = integer(0), ref = character(0), sapply(basedifflevels, function(i) integer(0)), stringsAsFactors = FALSE, row.names = NULL, check.names = FALSE)
        }
    }))
    if (nrow(basedifffreq_cutdemux_df) > 0) {
        basedifffreq_cutdemux_df <- data.frame(ExptID = exptID, LibraryID = libraryID, MitoID = mitoID, basedifffreq_cutdemux_df, stringsAsFactors = FALSE, check.names = FALSE)
    } else {
        basedifffreq_cutdemux_df <- data.frame(ExptID = character(0), LibraryID = character(0), MitoID = character(0), SNVID = character(0), pos = integer(0), ref = character(0), sapply(basedifflevels, function(i) integer(0)), stringsAsFactors = FALSE, row.names = NULL, check.names = FALSE)
    }

    resultdir <- sprintf("%s/%s/SNVs", ifelse(is_online(), result_dir_cache, result_dir), version_dir)
    dir.create(resultdir, FALSE, TRUE)
    write.csv(basedifffreq_cutdemux_df, file = gzfile(sprintf("%s/basedifffreq_cutdemux_%s_sub500k_q%s.csv.gz", resultdir, libraryMitoID, qth)), row.names = FALSE)
    basedifffreq_cutdemux_df <- read.csv(sprintf("%s/basedifffreq_cutdemux_%s_sub500k_q%s.csv.gz", resultdir, libraryMitoID, qth), as.is = TRUE, check.names = FALSE)

    ## 1.3 Calculate each base's frequency for each position, strand merged
    basedifffreq_cutdemux_unstranded <- data.frame(
        basedifffreq_cutdemux_df[, 1:6],
        "=" = basedifffreq_cutdemux_df[, "."] + basedifffreq_cutdemux_df[, ","],
        "A" = basedifffreq_cutdemux_df[, "A"] + basedifffreq_cutdemux_df[, "a"],
        "C" = basedifffreq_cutdemux_df[, "C"] + basedifffreq_cutdemux_df[, "c"],
        "G" = basedifffreq_cutdemux_df[, "G"] + basedifffreq_cutdemux_df[, "g"],
        "T" = basedifffreq_cutdemux_df[, "T"] + basedifffreq_cutdemux_df[, "t"],
        "del" = basedifffreq_cutdemux_df[, "*"] + basedifffreq_cutdemux_df[, "#"],
        stringsAsFactors = FALSE, check.names = FALSE
    )
    basedifffreq_cutdemux_unstranded$depth <- rowSums(basedifffreq_cutdemux_unstranded[, 7:12])
    basedifffreq_cutdemux_unstranded <- basedifffreq_cutdemux_unstranded[, c(1:6, 13, 7:12)]
    write.csv(basedifffreq_cutdemux_unstranded, file = gzfile(sprintf("%s/basedifffreq_cutdemux_%s_sub500k_q%s_unstranded.csv.gz", resultdir, libraryMitoID, qth)), row.names = FALSE)
    basedifffreq_cutdemux_unstranded <- read.csv(sprintf("%s/basedifffreq_cutdemux_%s_sub500k_q%s_unstranded.csv.gz", resultdir, libraryMitoID, qth), as.is = TRUE, check.names = FALSE)

    ## 2. Process insertions.
    mpileups_ins_cutdemux <- sapply(snvIDs, function(snvID) {
        f <- sprintf("%s/%s/analyzed/%s/%s/star/%s/%s.star.%s.mpileup.sub500k.ins.tsv.gz", data_dir, expt_dir, libraryName, libraryMitoName, snvID, libraryMitoName, starbase)
        message("Reading ", f)
        tryCatch(read.csv(f, as.is = TRUE, head = FALSE, sep = "\t", quote = "", colClasses = c("character", "integer", "character", "integer", "character", "character", "integer")), error = function(e) 
            data.frame(V1=character(0), V2=integer(0), V3=character(0), V4=integer(0), V5=character(0), V6=character(0), V7=integer(0), stringsAsFactors = FALSE))
    }, simplify = FALSE)

    datadir <- sprintf("%s/%s/SNVs/ins", ifelse(is_online(), data_dir_cache, data_dir), version_dir)
    dir.create(datadir, FALSE, TRUE)
    saveRDS(mpileups_ins_cutdemux, file = sprintf("%s/mpileups_ins_cutdemux_%s_sub500k.RDS", datadir, libraryMitoID))
    mpileups_ins_cutdemux <- readRDS(file = sprintf("%s/mpileups_ins_cutdemux_%s_sub500k.RDS", datadir, libraryMitoID))

    mpileups_ins_cutdemux <- do.call(rbind, lapply(snvIDs, function(snvID) {
        message("Processing insertion: ", exptID, " ", libraryID, " ", mitoID, " ", snvID)
        mpileup <- mpileups_ins_cutdemux[[snvID]]
        if (nrow(mpileup) > 0) { ## insersions don't have phred quality
            return(with(mpileup, data.frame(ExptID = exptID, LibraryID = libraryID, MitoID = mitoID, SNVID = snvID, pos = V2, ref = V3, depth = V4, alt = V5, altfreq = V6, altfreqtot = V7, stringsAsFactors = FALSE, row.names = NULL, check.names = FALSE)))
        } else {
            return(data.frame(ExptID = character(0), LibraryID = character(0), MitoID = character(0), SNVID = character(0), pos = integer(0), ref = character(0), depth = integer(0), alt = character(0), altfreq = character(0), altfreqtot = integer(0), stringsAsFactors = FALSE, row.names = NULL, check.names = FALSE))
        }
    }))

    resultdir <- sprintf("%s/%s/SNVs/ins", ifelse(is_online(), result_dir_cache, result_dir), version_dir)
    dir.create(resultdir, FALSE, TRUE)
    write.csv(mpileups_ins_cutdemux, file = gzfile(sprintf("%s/mpileups_ins_cutdemux_%s_sub500k.csv.gz", resultdir, libraryMitoID)), row.names = FALSE)
    mpileups_ins_cutdemux <- read.csv(sprintf("%s/mpileups_ins_cutdemux_%s_sub500k.csv.gz", resultdir, libraryMitoID), as.is = TRUE, check.names = FALSE)
}
parse_basediff_bymito(exptID = opts$exptID, libraryID = opts$libraryID, mitoID = opts$mitoID, snvIDs = paste0("SNV", 1:12), qth = opts$qth, starbase = opts$starbase, version_dir = opts$version_dir)
