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
source Source/functions/SMITO.bash

ln -sf ../../../Datasets Data/
mkdir -p Data/Datasets/lab/repo/E.smito
ln -sf ../../../Datasets/lab/repo/E.smito Data/
## limit to 64 cores
awk -vFS=',' '$6=="SE" { gsub("M", "", $4); gsub(";", ",", $4); print }' Data/LibraryInfo.csv | parallel -j 9 --colsep ' ' "run_pipe --exptID {2} --sampleID {1} --endType {6} --subsam_n 500000 --subsam_l 500k --genome Data/mm10.mito/chrM.fa --bcidx {4} --ncores 7 --verbose" 
