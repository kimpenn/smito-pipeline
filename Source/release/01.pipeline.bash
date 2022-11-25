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
source Source/functions/SMITO.bash

ln -sf ../../../Datasets Data/
mkdir -p Data/Datasets/lab/repo/E.smito
ln -sf ../../../Datasets/lab/repo/E.smito Data/
## limit to 64 cores
awk -vFS=',' '$6=="SE" { gsub("M", "", $4); gsub(";", ",", $4); print }' Data/LibraryInfo.csv | parallel -j 9 --colsep ' ' "run_pipe --exptID {2} --sampleID {1} --endType {6} --subsam_n 500000 --subsam_l 500k --genome /lab/repo/resources/src/mm10.mito/chrM.fa --bcidx {4} --ncores 7 --verbose" 
