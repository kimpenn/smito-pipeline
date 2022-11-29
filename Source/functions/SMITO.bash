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
## v1.3

## Prepare folders to transfer data from the repo to user workspace
## repoDir: the repo folder (default: /lab/repo)
## baseDir: the top-level user folder (default: ./Data)
## exptDir: the user folder that follows the PennSCAP-T structure (for the specification and more details see https://github.com/kimpenn/ngs-pipeline) 
## exptID: the library's run ID in /lab/repo
## sampleID: the library name
## bcidx: the mito-barcode indices, separated by space (default: '1 2 3 4 5 6 7 8 9 10')
function init_repo {
    local repoDir="/lab/repo"
    local baseDir="Data"
    local exptDir="E.smito"
    local verbose=false
    local bcidx='1 2 3 4 5 6 7 8 9 10'
    while [[ $# -gt 0 ]]; do
        case $1 in
            --repoDir)
                repoDir=$2
                shift; shift
                ;;
            --baseDir)
                baseDir=$2
                shift; shift
                ;;
            --exptDir)
                exptDir=$2
                shift; shift
                ;;
            --exptID)
                exptID=$2
                shift; shift
                ;;
            --sampleID)
                sampleID=$2
                shift; shift
                ;;
            --bcidx)
                bcidx="$2"
                shift; shift
                ;;
            --verbose)
                verbose=true
                shift
                ;;
            *)
                echo "init_repo [--repoDir /lab/repo] [--baseDir Data] [--exptDir E.smito] --exptID 730 --sampleID L19R23P3 [--bcidx '1 2 3 4 5 6 7 8 9 10'] [--verbose]"
                return 1
                ;;
        esac
    done

    sampleName=Sample_$sampleID
    exptName=E.$exptID
    for d in blast fastqc log; do 
        if [[ $verbose == "true" ]]; then
            echo mkdir -p $baseDir/$exptDir/analyzed/$sampleName/$d
            echo ln -sf $repoDir/$exptName/analyzed/$sampleName/$d/* $baseDir/$exptDir/analyzed/$sampleName/$d/
        fi
        mkdir -p $baseDir/$exptDir/analyzed/$sampleName/$d
        ln -sf $repoDir/$exptName/analyzed/$sampleName/$d/* $baseDir/$exptDir/analyzed/$sampleName/$d/
    done

    for m in $bcidx; do
        sampleNameM=${sampleName}_M$m
        bcM=BC_M$m
        for d in fastqc fastqc.trim log star trim verse verse-SMITO; do 
            if [[ $verbose == "true" ]]; then
                echo mkdir -p $baseDir/$exptDir/analyzed/$sampleName/$sampleNameM/$d
                echo ln -sf $repoDir/$exptName/analyzed/$sampleName/$sampleNameM/$d/* $baseDir/$exptDir/analyzed/$sampleName/$sampleNameM/$d/
            fi
            mkdir -p $baseDir/$exptDir/analyzed/$sampleName/$sampleNameM/$d
            ln -sf $repoDir/$exptName/analyzed/$sampleName/$sampleNameM/$d/* $baseDir/$exptDir/analyzed/$sampleName/$sampleNameM/$d/
        done
        if [[ $verbose == "true" ]]; then
            echo mkdir -p $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$bcM
            echo ln -sf $repoDir/$exptName/analyzed/$sampleName/cutdemux/$bcM/* $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$bcM/
        fi
        mkdir -p $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$bcM
        ln -sf $repoDir/$exptName/analyzed/$sampleName/cutdemux/$bcM/* $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$bcM/
    done

    if [[ $verbose == "true" ]]; then
        echo mkdir -p $baseDir/$exptDir/analyzed/$sampleName/cutdemux/notfound
        echo ln -sf $repoDir/$exptName/analyzed/$sampleName/cutdemux/notfound/* $baseDir/$exptDir/analyzed/$sampleName/cutdemux/notfound/
    fi
    mkdir -p $baseDir/$exptDir/analyzed/$sampleName/cutdemux/notfound
    ln -sf $repoDir/$exptName/analyzed/$sampleName/cutdemux/notfound/* $baseDir/$exptDir/analyzed/$sampleName/cutdemux/notfound/

    if [[ $verbose == "true" ]]; then
        echo mkdir -p $baseDir/$exptDir/raw/$sampleName
        echo ln -sf $repoDir/$exptName/raw/$sampleName/${sampleID}_* $baseDir/$exptDir/raw/$sampleName/
    fi
    mkdir -p $baseDir/$exptDir/raw/$sampleName
    ln -sf $repoDir/$exptName/raw/$sampleName/${sampleID}_* $baseDir/$exptDir/raw/$sampleName/
}

## Search and count PCR forward and reverse primers for SNV1-12 successively in single-end reads
## baseDir: the top-level user folder (default: ./Data)
## exptDir: the user folder that follows the PennSCAP-T structure (for the specification and more details see https://github.com/kimpenn/ngs-pipeline) 
## sampleName: the library name with "Sample_" prefix
## mitoID: the mito-barcode name, choose from "M1", "M2", ..., "M10" 
function primer_search_se {
    local baseDir="Data"
    local exptDir="E.smito"
    local verbose=false
    while [[ $# -gt 0 ]]; do
        case $1 in
            --baseDir)
                baseDir=$2
                shift; shift
                ;;
            --exptDir)
                exptDir=$2
                shift; shift
                ;;
            --sampleName)
                sampleName=$2
                shift; shift
                ;;
            --mitoID)
                m=$2
                shift; shift
                ;;
            --verbose)
                verbose=true
                shift
                ;;
            *)
                echo "primer_search_se [--baseDir Data] [--exptDir E.smito] --sampleName Sample_L19R23P3 --mitoID M1 [--verbose]"
                return 1
                ;;
        esac
    done
    if [[ $verbose == "true" ]]; then
        echo perl Source/functions/SMITO/search_pcr_primers_se.pl $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m unaligned 2\> $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/search_pcr_primers_pass_se.log
    fi
    perl Source/functions/SMITO/search_pcr_primers_se.pl $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m unaligned 2> $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/search_pcr_primers_pass_se.log
    if [[ $verbose == "true" ]]; then
        echo mkdir -p $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/notfound_1
    fi
    mkdir -p $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/notfound_1
    if [[ $verbose == "true" ]]; then
        echo mv $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/notfound_1.fq.gz $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/notfound_1/unaligned_1.fq.gz
    fi
    mv $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/notfound_1.fq.gz $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/notfound_1/unaligned_1.fq.gz

    for vi in SNV{1..12}; do 
        for si in fw_1 rv_1; do 
            snvi=${vi}_${si}
            if [[ $verbose == "true" ]]; then
                echo gzip -cd $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/$snvi/unaligned_1.fq.gz \| wc -l \| awk \'{ print \$1 / 4 }\' \> $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/$snvi/unaligned_1.nreads.txt
            fi
            gzip -cd $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/$snvi/unaligned_1.fq.gz | wc -l | awk '{ print $1 / 4 }' > $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/$snvi/unaligned_1.nreads.txt
        done
    done
    
    if [[ $verbose == "true" ]]; then
        echo gzip -cd $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/notfound_1/unaligned_1.fq.gz \| wc -l \| awk \'{ print \$1 / 4 }\' \> $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/notfound_1/unaligned_1.nreads.txt
    fi
    gzip -cd $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/notfound_1/unaligned_1.fq.gz | wc -l | awk '{ print $1 / 4 }' > $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/notfound_1/unaligned_1.nreads.txt
}

## Search and count PCR forward and reverse primers for SNV1-12 successively in paired-end reads
## baseDir: the top-level user folder (default: ./Data)
## exptDir: the user folder that follows the PennSCAP-T structure (for the specification and more details see https://github.com/kimpenn/ngs-pipeline) 
## sampleName: the library name with "Sample_" prefix
## mitoID: the mito-barcode name, choose from "M1", "M2", ..., "M10" 
function primer_search_pe {
    local baseDir="Data"
    local exptDir="E.smito"
    local verbose=false
    while [[ $# -gt 0 ]]; do
        case $1 in
            --baseDir)
                baseDir=$2
                shift; shift
                ;;
            --exptDir)
                exptDir=$2
                shift; shift
                ;;
            --sampleName)
                sampleName=$2
                shift; shift
                ;;
            --mitoID)
                m=$2
                shift; shift
                ;;
            --verbose)
                verbose=true
                shift
                ;;
            *)
                echo "primer_search_pe [--baseDir Data] [--exptDir E.smito] --sampleName Sample_L19R23P3 --mitoID M1 [--verbose]"
                return 1
                ;;
        esac
    done
    if [[ $verbose == "true" ]]; then
        echo perl Source/functions/SMITO/search_pcr_primers_pe_pass1.pl $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m unaligned 2\> $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/search_pcr_primers_pe_pass1.log
    fi
    perl Source/functions/SMITO/search_pcr_primers_pe_pass1.pl $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m unaligned 2> $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/search_pcr_primers_pe_pass1.log

    if [[ $verbose == "true" ]]; then
        echo mkdir -p $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/notfound_1
    fi
    mkdir -p $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/notfound_1

    if [[ $verbose == "true" ]]; then
        echo mv $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/notfound_1.fq.gz $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/notfound_1/unaligned_1.fq.gz
        echo mv $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/notfound_2.fq.gz $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/notfound_1/unaligned_2.fq.gz
    fi
    mv $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/notfound_1.fq.gz $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/notfound_1/unaligned_1.fq.gz
    mv $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/notfound_2.fq.gz $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/notfound_1/unaligned_2.fq.gz

    if [[ $verbose == "true" ]]; then
        echo perl Source/functions/SMITO/search_pcr_primers_pe_pass2.pl $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m unaligned 2\> $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/search_pcr_primers_pe_pass2.log
    fi
    perl Source/functions/SMITO/search_pcr_primers_pe_pass2.pl $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m unaligned 2> $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/search_pcr_primers_pe_pass2.log

    for v in SNV{1..12}; do 
        for s in fw_1 rv_1; do 
            snv=${v}_${s}
            if [[ $verbose == "true" ]]; then
                echo mkdir -p $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/$snv/notfound_2
                echo mv $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/$snv/notfound_1.fq.gz $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/$snv/notfound_2/unaligned_1.fq.gz
                echo mv $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/$snv/notfound_2.fq.gz $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/$snv/notfound_2/unaligned_2.fq.gz
            fi
            mkdir -p $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/$snv/notfound_2
            mv $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/$snv/notfound_1.fq.gz $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/$snv/notfound_2/unaligned_1.fq.gz
            mv $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/$snv/notfound_2.fq.gz $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/$snv/notfound_2/unaligned_2.fq.gz
        done
    done
    if [[ $verbose == "true" ]]; then
        echo mkdir -p $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/notfound_1/notfound_2
        echo mv $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/notfound_1/notfound_1.fq.gz $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/notfound_1/notfound_2/unaligned_1.fq.gz
        echo mv $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/notfound_1/notfound_2.fq.gz $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/notfound_1/notfound_2/unaligned_2.fq.gz
    fi
    mkdir -p $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/notfound_1/notfound_2
    mv $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/notfound_1/notfound_1.fq.gz $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/notfound_1/notfound_2/unaligned_1.fq.gz
    mv $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/notfound_1/notfound_2.fq.gz $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/notfound_1/notfound_2/unaligned_2.fq.gz
}

## Search and count PCR forward and reverse primers for SNV1-12 successively, 
## works for single-end and paired-end sequencing, calls the worker routines 
## `primer_search_se` and `primer_search_pe`, respectively. 
## baseDir: the top-level user folder (default: ./Data)
## exptDir: the user folder that follows the PennSCAP-T structure (for the specification and more details see https://github.com/kimpenn/ngs-pipeline) 
## endType: the type of sequencing end, choose from "SE" (single-end) and "PE" (paired-end)
## sampleID: the library name
## bcidx: the mito-barcode indices, separated by space (default: '1 2 3 4 5 6 7 8 9 10')
## ncores: the number of mito-barcodes to search in parallel (default: 5)
function primer_stats {
    local baseDir="Data"
    local exptDir="E.smito"
    local endType="SE"
    local bcidx='1 2 3 4 5 6 7 8 9 10'
    local verbose=false
    local verboseFlag=""
    local ncores=5
    while [[ $# -gt 0 ]]; do
        case $1 in
            --baseDir)
                baseDir=$2
                shift; shift
                ;;
            --exptDir)
                exptDir=$2
                shift; shift
                ;;
            --endType)
                endType=$2
                shift; shift
                ;;
            --sampleID)
                sampleID=$2
                shift; shift
                ;;
            --bcidx)
                bcidx=$2
                shift; shift
                ;;
            --ncores)
                ncores=$2
                shift; shift
                ;;
            --verbose)
                verbose=true
                shift
                ;;
            *)
                echo "primer_stats [--baseDir Data] [--endType SE] [--exptDir E.smito] [--ncores 5] [--bcidx '1 2 3 4 5 6 7 8 9 10'] --sampleID L19R23P3 [--verbose]"
                return 1
                ;;
        esac
    done

    sampleName=Sample_$sampleID
    if [[ $verbose==true ]]; then
        verboseFlag="--verbose"
    fi

    if [[ $endType == "SE" ]]; then
        for i in $bcidx; do
            m=BC_M$i
            if [[ $verbose == "true" ]]; then
                echo perl Source/functions/SMITO/parse_cutdemux_log_se.pl $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/*.cutdemux.log \> $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/stats.tsv
            fi
            perl Source/functions/SMITO/parse_cutdemux_log_se.pl $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/*.cutdemux.log > $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/stats.tsv
            if [[ $verbose == "true" ]]; then
                echo gzip -cd $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/unaligned_1.fq.gz \| wc -l \| awk \'{ print \$1 / 4 }\' \> $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/stats_nreadpairs.txt
            fi
            gzip -cd $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/unaligned_1.fq.gz | wc -l | awk '{ print $1 / 4 }' > $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/stats_nreadpairs.txt
        done
        if [[ $verbose == "true" ]]; then
            echo gzip -cd $baseDir/$exptDir/analyzed/$sampleName/cutdemux/notfound/unaligned_1.fq.gz \| wc -l \| awk \'{ print \$1 / 4 }\' \> $baseDir/$exptDir/analyzed/$sampleName/cutdemux/notfound/stats_nreadpairs.txt
        fi
        gzip -cd $baseDir/$exptDir/analyzed/$sampleName/cutdemux/notfound/unaligned_1.fq.gz | wc -l | awk '{ print $1 / 4 }' > $baseDir/$exptDir/analyzed/$sampleName/cutdemux/notfound/stats_nreadpairs.txt

        echo $bcidx | tr ' ' "\n" | parallel -j $ncores "primer_search_se --baseDir $baseDir --exptDir $exptDir --sampleName $sampleName --mitoID BC_M{1} $verboseFlag" 
    else
        for i in $bcidx; do 
            m=BC_M$i
            if [[ $verbose == "true" ]]; then
                echo perl Source/functions/SMITO/parse_cutdemux_log_pe.pl $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/*.cutdemux.log \> $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/stats.tsv
            fi
            perl Source/functions/SMITO/parse_cutdemux_log_pe.pl $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/*.cutdemux.log > $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/stats.tsv

            if [[ $verbose == "true" ]]; then
                echo gzip -cd $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/unaligned_1.fq.gz \| wc -l \| awk \'{ print \$1 / 4 }\' \> $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/stats_nreadpairs.txt
            fi
            gzip -cd $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/unaligned_1.fq.gz | wc -l | awk '{ print $1 / 4 }' > $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/stats_nreadpairs.txt
        done

        if [[ $verbose == "true" ]]; then
            echo gzip -cd $baseDir/$exptDir/analyzed/$sampleName/cutdemux/notfound/unaligned_1.fq.gz \| wc -l \| awk \'{ print \$1 / 4 }\' \> $baseDir/$exptDir/analyzed/$sampleName/cutdemux/notfound/stats_nreadpairs.txt
        fi
        gzip -cd $baseDir/$exptDir/analyzed/$sampleName/cutdemux/notfound/unaligned_1.fq.gz | wc -l | awk '{ print $1 / 4 }' > $baseDir/$exptDir/analyzed/$sampleName/cutdemux/notfound/stats_nreadpairs.txt

        echo $bcidx | tr ' ' "\n" | parallel -j $ncores "primer_search_pe --baseDir $baseDir --exptDir $exptDir --sampleName $sampleName --mitoID BC_M{1} $verboseFlag" 

        for i in $bcidx; do
            m=BC_M$i
            for vi in SNV{1..12}; do
                for si in fw_1 rv_1; do
                    snvi=${vi}_${si}
                    for vj in SNV{1..12}; do
                        for sj in fw_2 rv_2; do
                            snvj=${vj}_${sj}
                            if [[ $verbose == "true" ]]; then
                                echo gzip -cd $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/$snvi/$snvj/unaligned_1.fq.gz \| wc -l \| awk \'{ print \$1 / 4 }\' \> $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/$snvi/$snvj/unaligned_1.nreads.txt
                            fi
                            gzip -cd $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/$snvi/$snvj/unaligned_1.fq.gz | wc -l | awk '{ print $1 / 4 }' > $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/$snvi/$snvj/unaligned_1.nreads.txt
                        done
                    done
                    if [[ $verbose == "true" ]]; then
                        echo gzip -cd $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/$snvi/notfound_2/unaligned_1.fq.gz \| wc -l \| awk \'{ print \$1 / 4 }\' \> $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/$snvi/notfound_2/unaligned_1.nreads.txt
                    fi
                    gzip -cd $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/$snvi/notfound_2/unaligned_1.fq.gz | wc -l | awk '{ print $1 / 4 }' > $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/$snvi/notfound_2/unaligned_1.nreads.txt
                done
            done
            for vj in SNV{1..12}; do
                for sj in fw_2 rv_2; do
                    snvj=${vj}_${sj}
                    if [[ $verbose == "true" ]]; then
                        echo gzip -cd $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/notfound_1/$snvj/unaligned_1.fq.gz \| wc -l \| awk \'{ print \$1 / 4 }\' \> $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/notfound_1/$snvj/unaligned_1.nreads.txt
                    fi
                    gzip -cd $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/notfound_1/$snvj/unaligned_1.fq.gz | wc -l | awk '{ print $1 / 4 }' > $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/notfound_1/$snvj/unaligned_1.nreads.txt
                done
            done
            if [[ $verbose == "true" ]]; then
                echo gzip -cd $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/notfound_1/notfound_2/unaligned_1.fq.gz \| wc -l \| awk \'{ print \$1 / 4 }\' \> $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/notfound_1/notfound_2/unaligned_1.nreads.txt
            fi
            gzip -cd $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/notfound_1/notfound_2/unaligned_1.fq.gz | wc -l | awk '{ print $1 / 4 }' > $baseDir/$exptDir/analyzed/$sampleName/cutdemux/$m/notfound_1/notfound_2/unaligned_1.nreads.txt
        done
    fi
}

## Use samtools and Picard to generate sequence alignment stats 
## baseDir: the top-level user folder (default: ./Data)
## exptDir: the user folder that follows the PennSCAP-T structure (for the specification and more details see https://github.com/kimpenn/ngs-pipeline) 
## endType: the type of sequencing end, choose from "SE" (single-end) and "PE" (paired-end)
## sampleID: the library name
## bcidx: the mito-barcode indices, separated by space (default: '1 2 3 4 5 6 7 8 9 10')
## ncores: the multithread parameter for Java GC (default: 5)
## picard: the location of the Picard JAR file (default: $HOME/Applications/Picard/2.17.0/picard.jar)
function sam_stats {
    local baseDir="Data"
    local exptDir="E.smito"
    local endType="SE"
    local ncores=5
    local bcidx='1 2 3 4 5 6 7 8 9 10'
    local verbose=false
    local picard=$HOME/Applications/Picard/2.17.0/picard.jar
    while [[ $# -gt 0 ]]; do
        case $1 in
            --baseDir)
                baseDir=$2
                shift; shift
                ;;
            --exptDir)
                exptDir=$2
                shift; shift
                ;;
            --endType)
                endType=$2
                shift; shift
                ;;
            --sampleID)
                sampleID=$2
                shift; shift
                ;;
            --bcidx)
                bcidx=$2
                shift; shift
                ;;
            --ncores)
                ncores=$2
                shift; shift
                ;;
            --picard)
                picard=$2
                shift; shift
                ;;
            --verbose)
                verbose=true
                shift
                ;;
            *)
                echo "sam_stats [--baseDir Data] [--endType SE] [--exptDir E.smito] [--bcidx '1 2 3 4 5 6 7 8 9 10'] --sampleID L19R23P3 [--ncores 5] [--picard $HOME/Applications/Picard/2.17.0/picard.jar] [--verbose]"
                return 1
                ;;
        esac
    done

    sampleName=Sample_$sampleID

    for i in $bcidx; do
        m=M$i
        if [[ $verbose == "true" ]]; then
            echo java -XX:ParallelGCThreads=$ncores -jar $picard MarkDuplicates INPUT=$baseDir/$exptDir/analyzed/$sampleName/${sampleName}_$m/star/${sampleName}_$m.star.posSorted.bam OUTPUT=$baseDir/$exptDir/analyzed/$sampleName/${sampleName}_$m/star/${sampleName}_$m.star.DupMarked.bam METRICS_FILE=$baseDir/$exptDir/analyzed/$sampleName/${sampleName}_$m/star/${sampleName}_$m.star.DupMarked.stats.txt REMOVE_DUPLICATES=false USE_JDK_DEFLATER=true USE_JDK_INFLATER=true 2\> $baseDir/$exptDir/analyzed/$sampleName/${sampleName}_$m/star/${sampleName}_$m.star.DupMarked.bam.log 
        fi
        java -XX:ParallelGCThreads=$ncores -jar $picard MarkDuplicates INPUT=$baseDir/$exptDir/analyzed/$sampleName/${sampleName}_$m/star/${sampleName}_$m.star.posSorted.bam OUTPUT=$baseDir/$exptDir/analyzed/$sampleName/${sampleName}_$m/star/${sampleName}_$m.star.DupMarked.bam METRICS_FILE=$baseDir/$exptDir/analyzed/$sampleName/${sampleName}_$m/star/${sampleName}_$m.star.DupMarked.stats.txt REMOVE_DUPLICATES=false USE_JDK_DEFLATER=true USE_JDK_INFLATER=true 2> $baseDir/$exptDir/analyzed/$sampleName/${sampleName}_$m/star/${sampleName}_$m.star.DupMarked.bam.log 
        if [[ $verbose == "true" ]]; then
            echo samtools flagstat $baseDir/$exptDir/analyzed/$sampleName/${sampleName}_$m/star/${sampleName}_$m.star.DupMarked.bam \> $baseDir/$exptDir/analyzed/$sampleName/${sampleName}_$m/star/${sampleName}_$m.star.DupMarked.flagstat.txt
        fi
        samtools flagstat $baseDir/$exptDir/analyzed/$sampleName/${sampleName}_$m/star/${sampleName}_$m.star.DupMarked.bam > $baseDir/$exptDir/analyzed/$sampleName/${sampleName}_$m/star/${sampleName}_$m.star.DupMarked.flagstat.txt
        if [[ $verbose == "true" ]]; then
            echo samtools view -h $baseDir/$exptDir/analyzed/$sampleName/${sampleName}_$m/star/${sampleName}_$m.star.DupMarked.bam \| perl -ne \'if \(/^@/\) { print } else { print if /NH:i:1\\t/ }\' \| samtools view -Sb - \> $baseDir/$exptDir/analyzed/$sampleName/${sampleName}_$m/star/${sampleName}_$m.star.unique.posSorted.bam
        fi
        samtools view -h $baseDir/$exptDir/analyzed/$sampleName/${sampleName}_$m/star/${sampleName}_$m.star.DupMarked.bam | perl -ne 'if (/^@/) { print } else { print if /NH:i:1\t/ }' | samtools view -Sb - > $baseDir/$exptDir/analyzed/$sampleName/${sampleName}_$m/star/${sampleName}_$m.star.unique.posSorted.bam
        if [[ $verbose == "true" ]]; then
            echo samtools index $baseDir/$exptDir/analyzed/$sampleName/${sampleName}_$m/star/${sampleName}_$m.star.unique.posSorted.bam
        fi
        samtools index $baseDir/$exptDir/analyzed/$sampleName/${sampleName}_$m/star/${sampleName}_$m.star.unique.posSorted.bam
        if [[ $verbose == "true" ]]; then
            echo samtools depth -d 0 -l 0 -q 0 -Q 0 -g 0x400 -J $baseDir/$exptDir/analyzed/$sampleName/${sampleName}_$m/star/${sampleName}_$m.star.unique.posSorted.bam \| gzip -c \> $baseDir/$exptDir/analyzed/$sampleName/${sampleName}_$m/star/${sampleName}_$m.star.unique.depth.tsv.gz
        fi
        samtools depth -d 0 -l 0 -q 0 -Q 0 -g 0x400 -J $baseDir/$exptDir/analyzed/$sampleName/${sampleName}_$m/star/${sampleName}_$m.star.unique.posSorted.bam | gzip -c > $baseDir/$exptDir/analyzed/$sampleName/${sampleName}_$m/star/${sampleName}_$m.star.unique.depth.tsv.gz

        ## Add an additional filter for reads that (1) start from the 5' end of the PCR forward primer, (2) have length >130bp, (3) on the forward strand 
        if [[ $verbose == "true" ]]; then
            echo perl Source/functions/SMITO/filter_sam.pl --infile $baseDir/$exptDir/analyzed/$sampleName/${sampleName}_$m/star/${sampleName}_$m.star.unique.posSorted.bam --outfile $baseDir/$exptDir/analyzed/$sampleName/${sampleName}_$m/star/${sampleName}_$m.star.strict.posSorted.bam --snvfile Data/snv_loci_v2.csv --minlen 135 --offup 1 --offdn 1 2\> $baseDir/$exptDir/analyzed/$sampleName/${sampleName}_$m/star/${sampleName}_$m.star.strict.posSorted.bam.log
        fi
        perl Source/functions/SMITO/filter_sam.pl --infile $baseDir/$exptDir/analyzed/$sampleName/${sampleName}_$m/star/${sampleName}_$m.star.unique.posSorted.bam --outfile $baseDir/$exptDir/analyzed/$sampleName/${sampleName}_$m/star/${sampleName}_$m.star.strict.posSorted.bam --snvfile Data/snv_loci_v2.csv --minlen 135 --offup 1 --offdn 1 2> $baseDir/$exptDir/analyzed/$sampleName/${sampleName}_$m/star/${sampleName}_$m.star.strict.posSorted.bam.log
        if [[ $verbose == "true" ]]; then
            echo samtools index $baseDir/$exptDir/analyzed/$sampleName/${sampleName}_$m/star/${sampleName}_$m.star.strict.posSorted.bam
        fi
        samtools index $baseDir/$exptDir/analyzed/$sampleName/${sampleName}_$m/star/${sampleName}_$m.star.strict.posSorted.bam
        if [[ $verbose == "true" ]]; then
            echo samtools depth -d 0 -l 0 -q 0 -Q 0 -g 0x400 -J $baseDir/$exptDir/analyzed/$sampleName/${sampleName}_$m/star/${sampleName}_$m.star.strict.posSorted.bam \| gzip -c \> $baseDir/$exptDir/analyzed/$sampleName/${sampleName}_$m/star/${sampleName}_$m.star.strict.depth.tsv.gz
        fi
        samtools depth -d 0 -l 0 -q 0 -Q 0 -g 0x400 -J $baseDir/$exptDir/analyzed/$sampleName/${sampleName}_$m/star/${sampleName}_$m.star.strict.posSorted.bam | gzip -c > $baseDir/$exptDir/analyzed/$sampleName/${sampleName}_$m/star/${sampleName}_$m.star.strict.depth.tsv.gz
    done
}

## Get coverage and mismatch call for each base position using `samtools mpileup`
## Here we don't distinguish SNV regions.
## baseDir: the top-level user folder (default: ./Data)
## exptDir: the user folder that follows the PennSCAP-T structure (for the specification and more details see https://github.com/kimpenn/ngs-pipeline) 
## subsam_n: the number of reads to subsample from the input BAM (useful if the input is too large; default: 500,000)
## subsam_l: the label for $subsam_n (default: "500k")
## sampleID: the library name
## genome: the location of the mito genome FASTA (default: /lab/repo/resources/src/mm10.mito/chrM.fa)
## bcidx: the mito-barcode indices, separated by space (default: '1 2 3 4 5 6 7 8 9 10')
function sam_mpileup_subsam_wholegenome {
    local baseDir="Data"
    local exptDir="E.smito"
    local verbose=false
    local subsam_n=500000
    local subsam_l='500k'
    local genome='/lab/repo/resources/src/mm10.mito/chrM.fa'
    local bcidx='1 2 3 4 5 6 7 8 9 10'
    while [[ $# -gt 0 ]]; do
        case $1 in
            --baseDir)
                baseDir=$2
                shift; shift
                ;;
            --exptDir)
                exptDir=$2
                shift; shift
                ;;
            --sampleID)
                sampleID=$2
                shift; shift
                ;;
            --bcidx)
                bcidx=$2
                shift; shift
                ;;
            --subsam_n)
                subsam_n=$2
                shift; shift
                ;;
            --subsam_l)
                subsam_l=$2
                shift; shift
                ;;
            --genome)
                genome=$2
                shift; shift
                ;;
            --verbose)
                verbose=true
                shift
                ;;
            *)
                echo "sam_mpileup_subsam_wholegenome [--baseDir Data] [--exptDir E.smito] [--subsam_n 500000] [--subsam_l 500k] [--genome /lab/repo/resources/src/mm10.mito/chrM.fa] [--bcidx '1 2 3 4 5 6 7 8 9 10'] --sampleID L19R23P3 [--verbose]"
                return 1
                ;;
        esac
    done

    sampleName=Sample_$sampleID
    for i in $bcidx; do
        m=M$i
        starDir=$baseDir/$exptDir/analyzed/$sampleName/${sampleName}_${m}/star
        snvDir=$starDir
        if [[ $verbose == "true" ]]; then
            echo samtools mpileup --excl-flags UNMAP,SECONDARY,QCFAIL --count-orphans --min-BQ 0 --min-MQ 0 --max-depth $subsam_n --reverse-del --fasta-ref $genome $snvDir/${sampleName}_${m}.star.unique.posSorted.bam \| gzip -c \> $snvDir/${sampleName}_${m}.star.unique.mpileup.sub$subsam_l.tsv.gz
        fi
        samtools mpileup --excl-flags UNMAP,SECONDARY,QCFAIL --count-orphans --min-BQ 0 --min-MQ 0 --max-depth $subsam_n --reverse-del --fasta-ref $genome $snvDir/${sampleName}_${m}.star.unique.posSorted.bam | gzip -c > $snvDir/${sampleName}_${m}.star.unique.mpileup.sub$subsam_l.tsv.gz

        if [[ $verbose == "true" ]]; then
            echo gzip -cd $snvDir/${sampleName}_${m}.star.unique.mpileup.sub$subsam_l.tsv.gz \| perl Source/functions/SMITO/rm_mpileup_indel.pl \| gzip -c \> $snvDir/${sampleName}_${m}.star.unique.mpileup.sub$subsam_l.noindel.tsv.gz
        fi
        gzip -cd $snvDir/${sampleName}_${m}.star.unique.mpileup.sub$subsam_l.tsv.gz | perl Source/functions/SMITO/rm_mpileup_indel.pl | gzip -c > $snvDir/${sampleName}_${m}.star.unique.mpileup.sub$subsam_l.noindel.tsv.gz

        if [[ $verbose == "true" ]]; then
            echo "gzip -cd $snvDir/${sampleName}_${m}.star.unique.mpileup.sub$subsam_l.tsv.gz | perl Source/functions/SMITO/get_mpileup_ins.pl | gzip -c > $snvDir/${sampleName}_${m}.star.unique.mpileup.sub$subsam_l.ins.tsv.gz"
        fi
        gzip -cd $snvDir/${sampleName}_${m}.star.unique.mpileup.sub$subsam_l.tsv.gz | perl Source/functions/SMITO/get_mpileup_ins.pl | gzip -c > $snvDir/${sampleName}_${m}.star.unique.mpileup.sub$subsam_l.ins.tsv.gz
    done
}

## Get coverage and mismatch call for each base position using `samtools mpileup` 
## for each SNV region; bases overflowing each SNV region will be ignored.
## baseDir: the top-level user folder (default: ./Data)
## exptDir: the user folder that follows the PennSCAP-T structure (for the specification and more details see https://github.com/kimpenn/ngs-pipeline) 
## subsam_n: the number of reads to subsample from the input BAM (useful if the input is too large; default: 500,000)
## subsam_l: the label for $subsam_n (default: "500k")
## sampleID: the library name
## genome: the location of mouse mito genome FASTA (default: /lab/repo/resources/src/mm10.mito/chrM.fa)
## bcidx: the mito-barcode indices, separated by space (default: '1 2 3 4 5 6 7 8 9 10')
function sam_mpileup_subsam {
    local baseDir="Data"
    local exptDir="E.smito"
    local verbose=false
    local subsam_n=500000
    local subsam_l='500k'
    local genome='/lab/repo/resources/src/mm10.mito/chrM.fa'
    local bcidx='1 2 3 4 5 6 7 8 9 10'
    local starbase="unique"
    while [[ $# -gt 0 ]]; do
        case $1 in
            --baseDir)
                baseDir=$2
                shift; shift
                ;;
            --exptDir)
                exptDir=$2
                shift; shift
                ;;
            --sampleID)
                sampleID=$2
                shift; shift
                ;;
            --bcidx)
                bcidx=$2
                shift; shift
                ;;
            --subsam_n)
                subsam_n=$2
                shift; shift
                ;;
            --subsam_l)
                subsam_l=$2
                shift; shift
                ;;
            --genome)
                genome=$2
                shift; shift
                ;;
            --starbase)
                starbase=$2
                shift; shift
                ;;
            --verbose)
                verbose=true
                shift
                ;;
            *)
                echo "sam_mpileup_subsam [--baseDir Data] [--exptDir E.smito] [--subsam_n 500000] [--subsam_l 500k] [--genome /lab/repo/resources/src/mm10.mito/chrM.fa] [--bcidx '1 2 3 4 5 6 7 8 9 10'] --sampleID L19R23P3 [--starbase unique] [--verbose]"
                return 1
                ;;
        esac
    done

    sampleName=Sample_$sampleID
    for i in $bcidx; do
        m=M$i
        starDir=$baseDir/$exptDir/analyzed/$sampleName/${sampleName}_${m}/star
        for snv in SNV{1..12}; do 
            range=$(awk -vFS=',' -vsnv=$snv '$1==snv { print $6 }' Data/snv_loci_v2.csv)
            snvDir=$starDir/$snv
            mkdir -p $snvDir
            if [[ $verbose == "true" ]]; then
                echo samtools view -h $starDir/${sampleName}_${m}.star.${starbase}.posSorted.bam chrM:$range \| samtools view -Sb - \> $snvDir/${sampleName}_${m}.star.${starbase}.posSorted.bam
            fi
            samtools view -h $starDir/${sampleName}_${m}.star.${starbase}.posSorted.bam chrM:$range | samtools view -Sb - > $snvDir/${sampleName}_${m}.star.${starbase}.posSorted.bam
            
            if [[ $verbose == "true" ]]; then
                echo samtools index $snvDir/${sampleName}_${m}.star.${starbase}.posSorted.bam
            fi
            samtools index $snvDir/${sampleName}_${m}.star.${starbase}.posSorted.bam

            if [[ $verbose == "true" ]]; then
                echo samtools mpileup --excl-flags UNMAP,SECONDARY,QCFAIL --count-orphans --min-BQ 0 --min-MQ 0 --max-depth $subsam_n --reverse-del --fasta-ref $genome $snvDir/${sampleName}_${m}.star.${starbase}.posSorted.bam --region chrM:$range \| gzip -c \> $snvDir/${sampleName}_${m}.star.${starbase}.mpileup.sub$subsam_l.tsv.gz #--output-extra FLAG,POS
            fi
            samtools mpileup --excl-flags UNMAP,SECONDARY,QCFAIL --count-orphans --min-BQ 0 --min-MQ 0 --max-depth $subsam_n --reverse-del --fasta-ref $genome $snvDir/${sampleName}_${m}.star.${starbase}.posSorted.bam --region chrM:$range | gzip -c > $snvDir/${sampleName}_${m}.star.${starbase}.mpileup.sub$subsam_l.tsv.gz #--output-extra FLAG,POS

            if [[ $verbose == "true" ]]; then
                echo gzip -cd $snvDir/${sampleName}_${m}.star.${starbase}.mpileup.sub$subsam_l.tsv.gz \| perl Source/functions/SMITO/rm_mpileup_indel.pl \| gzip -c \> $snvDir/${sampleName}_${m}.star.${starbase}.mpileup.sub$subsam_l.noindel.tsv.gz
            fi
            gzip -cd $snvDir/${sampleName}_${m}.star.${starbase}.mpileup.sub$subsam_l.tsv.gz | perl Source/functions/SMITO/rm_mpileup_indel.pl | gzip -c > $snvDir/${sampleName}_${m}.star.${starbase}.mpileup.sub$subsam_l.noindel.tsv.gz

            if [[ $verbose == "true" ]]; then
                echo "gzip -cd $snvDir/${sampleName}_${m}.star.${starbase}.mpileup.sub$subsam_l.tsv.gz | perl Source/functions/SMITO/get_mpileup_ins.pl | gzip -c > $snvDir/${sampleName}_${m}.star.${starbase}.mpileup.sub$subsam_l.ins.tsv.gz"
            fi
            gzip -cd $snvDir/${sampleName}_${m}.star.${starbase}.mpileup.sub$subsam_l.tsv.gz | perl Source/functions/SMITO/get_mpileup_ins.pl | gzip -c > $snvDir/${sampleName}_${m}.star.${starbase}.mpileup.sub$subsam_l.ins.tsv.gz
        done
    done 
}

## The overall pipeline 
## repoDir: the repo folder (default: /lab/repo)
## baseDir: the top-level user folder (default: ./Data)
## exptDir: the user folder that follows the PennSCAP-T structure (for the specification and more details see https://github.com/kimpenn/ngs-pipeline) 
## picard: the location of the Picard JAR file (default: $HOME/Applications/Picard/2.17.0/picard.jar)
## subsam_n: the number of reads to subsample from the input BAM (useful if the input is too large; default: 500,000)
## subsam_l: the label for $subsam_n (default: "500k")
## ncores: the number of cores running in parallel (default: 5)
## sampleID: the library name
## exptID: the library's run ID in /lab/repo
## genome: the location of the mito genome FASTA (default: /lab/repo/resources/src/mm10.mito/chrM.fa)
## bcidx: the mito-barcode indices, separated by space (default: '1 2 3 4 5 6 7 8 9 10')
function run_pipe {
    local repoDir="/lab/repo"
    local baseDir="Data"
    local exptDir="E.smito"
    local endType="SE"
    local ncores=5
    local picard=$HOME/Applications/Picard/2.17.0/picard.jar
    local subsam_n=500000
    local subsam_l='500k'
    local genome='/lab/repo/resources/src/mm10.mito/chrM.fa'
    local bcidx='1 2 3 4 5 6 7 8 9 10'
    local verbose=false
    local verboseFlag=''
    while [[ $# -gt 0 ]]; do
        case $1 in
            --repoDir)
                repoDir=$2
                shift; shift
                ;;
            --baseDir)
                baseDir=$2
                shift; shift
                ;;
            --endType)
                endType=$2
                shift; shift
                ;;
            --exptDir)
                exptDir=$2
                shift; shift
                ;;
            --exptID)
                exptID=$2
                shift; shift
                ;;
            --sampleID)
                sampleID=$2
                shift; shift
                ;;
            --bcidx)
                bcidx=$2
                shift; shift
                ;;
            --ncores)
                ncores=$2
                shift; shift
                ;;
            --picard)
                picard=$2
                shift; shift
                ;;
            --subsam_n)
                subsam_n=$2
                shift; shift
                ;;
            --subsam_l)
                subsam_l=$2
                shift; shift
                ;;
            --genome)
                genome=$2
                shift; shift
                ;;
            --verbose)
                verbose=true
                verboseFlag='--verbose'
                shift
                ;;
            *)
                echo "run_pipe [--repoDir /lab/repo] [--baseDir Data] [--exptDir E.smito] [--endType SE] [--bcidx '1 2 3 4 5 6 7 8 9 10'] --exptID 730 --sampleID L19R23P3 [--ncores 5] [--picard $HOME/Applications/Picard/2.17.0/picard.jar] [--genome=/lab/repo/resources/src/mm10.mito/chrM.fa] [--subsam_n 500000] [--subsam_l 500k] [--verbose]"
                return 1
                ;;
        esac
    done

    if [[ $verbose == "true" ]]; then
        echo init_repo --repoDir $repoDir --baseDir $baseDir --exptDir $exptDir --exptID $exptID --sampleID $sampleID --bcidx "$bcidx" $verboseFlag
    fi
    init_repo --repoDir $repoDir --baseDir $baseDir --exptDir $exptDir --exptID $exptID --sampleID $sampleID --bcidx "$bcidx" $verboseFlag

    if [[ $verbose == "true" ]]; then
        echo primer_stats --baseDir $baseDir --endType $endType --exptDir $exptDir --sampleID $sampleID --bcidx "$bcidx" --ncores $ncores $verboseFlag
    fi
    primer_stats --baseDir $baseDir --endType $endType --exptDir $exptDir --sampleID $sampleID --bcidx "$bcidx" --ncores $ncores $verboseFlag

    if [[ $verbose == "true" ]]; then
        echo sam_stats --baseDir $baseDir --endType $endType --exptDir $exptDir --sampleID $sampleID --bcidx "$bcidx" --ncores $ncores --picard $picard $verboseFlag
    fi
    sam_stats --baseDir $baseDir --endType $endType --exptDir $exptDir --sampleID $sampleID --bcidx "$bcidx" --ncores $ncores --picard $picard $verboseFlag

    if [[ $verbose == "true" ]]; then
        echo sam_mpileup_subsam --baseDir $baseDir --exptDir $exptDir --subsam_n $subsam_n --subsam_l $subsam_l --genome $genome --bcidx "$bcidx" --starbase strict --sampleID $sampleID $verboseFlag
    fi
    sam_mpileup_subsam --baseDir $baseDir --exptDir $exptDir --subsam_n $subsam_n --subsam_l $subsam_l --genome $genome --bcidx "$bcidx" --starbase strict --sampleID $sampleID $verboseFlag
}

## Get coverage and mismatch call for each base position using `samtools mpileup` 
## This is particular for BWA output.
## for each SNV region; bases overflowing each SNV region will be ignored.
## $1: the library name
## $2: the libary name and mito-barcode concatenated
function sam_mpileup_bwa_permito {
    local libraryID=$1
    local libraryMitoID=$2
    baseDir=Data
    exptDir=E.smito
    sampleName=Sample_$libraryID
    sampleMitoName=Sample_$libraryMitoID
    bwaDir=$baseDir/$exptDir/analyzed/$sampleName/$sampleMitoName/bwa
    subsam_n=500000
    subsam_l='500k'
    genome=Data/mm10.mito/chrM.fa
    for snv in SNV{1..12}; do 
        range=$(awk -vFS=',' -vsnv=$snv '$1==snv { print $6 }' Data/snv_loci_v2.csv)
        snvDir=$bwaDir/$snv; 
        mkdir -p $snvDir
        samtools view -h $bwaDir/$sampleMitoName.bwa.q60.posSorted.bam chrM:$range | samtools view -Sb - > $snvDir/$sampleMitoName.bwa.q60.posSorted.bam
        samtools index $snvDir/$sampleMitoName.bwa.q60.posSorted.bam
        samtools mpileup --excl-flags UNMAP,SECONDARY,QCFAIL --count-orphans --min-BQ 0 --min-MQ 0 --max-depth $subsam_n --reverse-del --fasta-ref $genome $snvDir/$sampleMitoName.bwa.q60.posSorted.bam --region chrM:$range | gzip -c > $snvDir/$sampleMitoName.bwa.q60.mpileup.sub$subsam_l.tsv.gz #--output-extra FLAG,POS
        gzip -cd $snvDir/$sampleMitoName.bwa.q60.mpileup.sub$subsam_l.tsv.gz | perl Source/functions/SMITO/rm_mpileup_indel.pl | gzip -c > $snvDir/$sampleMitoName.bwa.q60.mpileup.sub$subsam_l.noindel.tsv.gz
        gzip -cd $snvDir/$sampleMitoName.bwa.q60.mpileup.sub$subsam_l.tsv.gz | perl Source/functions/SMITO/get_mpileup_ins.pl | gzip -c > $snvDir/$sampleMitoName.bwa.q60.mpileup.sub$subsam_l.ins.tsv.gz
    done
}

export -f init_repo
export -f primer_search_se
export -f primer_search_pe
export -f primer_stats
export -f sam_stats
export -f sam_mpileup_subsam
export -f sam_mpileup_subsam_wholegenome 
export -f run_pipe
export -f sam_mpileup_bwa_permito
