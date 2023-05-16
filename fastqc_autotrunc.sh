#!/bin/bash
###########################################################################################
###########################################################################################
####                                                                                   ####
####                   AUTOMATIC TRUNCATION AND TRIMMING OF NGS DATA                   ####
####                                                                                   ####
###########################################################################################
###########################################################################################

################################ Pablo Catarecha, 2023 ####################################

#### This script performs FASTQC analysis on bulk or pre-processed fastq data and picks the
#### optimal truncation values to be used in downstream DADA2 denoising.
#### User-provided sequences must be filtered prior to running this script in order to get
#### rid of adapter sequences.
#### This script works on the folder that stores the reads that will be imported into the
#### QIIME2 pipeline, and can be placed as part of it, right before denoising.
#### Please bear in mind that there are some folders that must exist before running this.
#### Additionally, fastqc must be installed in your system for this script to work.

####--------------------------Setting working directories------------------------------####

#### These MUST EXIST before running this script.

#### Default working paths.
WD=/your_working_folder ### This is the working folder. Everything hangs from this.
RAW_DATA=$WD/where_reads_are ### This is where fastq reads are stored to be imported into QIIME2.

####----------------------------------App structure------------------------------------####

echo "####----------------Process start----------------####"
echo " "
echo "####------------Performing FASTQC quality control------------####"
echo " "

mkdir $WD/fastqc_output

FASTQC=$WD/fastqc_output

fastqc --extract -o $FASTQC $RAW_DATA/*

echo " "
echo "####------------Calculating truncation points on FASTQC data------------####"
echo " "

#### This script parses FASTQC data and automates the calculation for the average truncation point including all reads in the experiment.
#### The truncation point is calculated for each fastq file as the base where median quality falls below average quality.
#### This point roughly matches the inflection point at which the observed quality score starts to be too low when analyzing fastq raw reads.
#### It is advisable to check several truncation values, in order to maximize the feature count in downstream analysis.

# Initialize count variables.
R1SUM=0
R2SUM=0
COUNT=0
LENGTH_SUM=0
# Parse raw fastq reads.
for q in $(ls $RAW_DATA)
do
# Extract unique IDs for each read.
# This line assumes your read filenames have the structure "experiment_sample_ID_sequencerCode_R1|2_sequencerLane.fastq.gz." If this is not your name format, please adapt it to your actual data.
qID=$(cut -d "_" -f3 <<< "$q")
# Check if picked read is on the forward side, according to the name format above.
if [[ "$q" =~ '_R1_' ]]
then
# Parse FASTQC output from the same read ID and find the base where median quality drops below average quality.
((R1SUM+=$(sed -e '1,/>>Per/d' -e '/>>END_MODULE/,$d' $FASTQC/*$qID*_R1_*/fastqc_data.txt | cut -f1,2,3 | sed -e '1,10d' | awk '($3-$2)<0' | head -1 | cut -d " " -f1 | cut -d "-" -f1) ))
# Count the number of reads. Only forward reads are counted in paired reads.
((COUNT++))
# Check the length of the forward reads. This will allow to set a minimum overlap between forward and reverse reads.
((LENGTH_SUM+=$(cat $FASTQC/*$qID*_R1_*/fastqc_data.txt | grep -A1 '#Length' | tail -1 | cut -f1) ))
else
# If picked read is on the reverse side, just find truncation point.
((R2SUM+=$(sed -e '1,/>>Per/d' -e '/>>END_MODULE/,$d' $FASTQC/*$qID*_R2_*/fastqc_data.txt | cut -f1,2,3 | sed -e '1,10d' | awk '($3-$2)<0' | head -1 | cut -d " " -f1 | cut -d "-" -f1) ))
fi
done
# Average truncation point for forward reads.
R1TRUN=$(($R1SUM / $COUNT))
# Average truncation point for reverse reads.
R2TRUN=$(($R2SUM / $COUNT))
# Average read length calculated on forward reads. Ideally, forward and reverse reads should be of the same length.
LENGTH=$(($LENGTH_SUM / $COUNT))
# Check if truncated reads from the same sample ID have a 12 nt overlap at least. If not, extend truncation point of forward reads to that limit.
if [[ "$R1TRUN + $RTRUN" < "$LENGTH + 12" ]]
then
gap=$(( $LENGTH + 12 - $R1TRUN - $R2TRUN))
((R1TRUN+=$gap))
fi
echo "Truncation points are $R1TRUN for forward reads and $R2TRUN for reverse reads."

echo " "
echo "####------------Truncation complete------------####"
echo " "

#### Truncation values calculated above can be integrated into the denoise step with DADA2, using the QIIME2 pipeline. The line below is just an example.

#qiime dada2 denoise-paired --i-demultiplexed-seqs your_imported_artifact_from_raw_data_folder.qza --p-trunc-len-f $R1TRUN --p-trunc-len-r $R2TRUN --o-table your_denoised_feature_table.qza --o-representative-sequences your_denoised_sequences.qza --o-denoising-stats your_denoising_stats.qza

echo " "
echo "####------------Process ending sucessfully!------------####"
