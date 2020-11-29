# Script to size filter, then cross ref vs RefSeq TSS/Exons
# 2 steps: #1 filt size, #2 filt for refseq
# usage: sh filtsizeandrefseq.sh [file] [size]
#size="$2"
#infile="$1"
#infileName=${infile%.*}
#awk -v size="$size" 'BEGIN{OFS=FS="\t"}{if (($3-$2)>size) {print $1,$2,$3,($3-$2)}}' $infile | bedtools intersect -v -a - -b NCBI.RefSeq.TSS.Exons.merged.txt > $outName

#In the past, we've captured HMRs 25 bp or longer. Now, we've switched to 50, assuming enhacers have multiple TF motifs#
#Run as a bash script#

module restore anaconda #load the module 'anaconda' which has the Anaconda module saved
#conda init bash #initiliaze Conda so our shell will be "conda-aware"
source activate bedtools #bedtools originally installed using conda and saved as a conda environment 'bedtools'

##Set variables and paths##
NCBI='/data/hodges_lab/aganve/data/DNA_Methylation/7_tracks/filedrop_temp/NCBI.RefSeq.TSS.Exons.merged.txt' #BED file of RefSeq promoters and exons

#First, filter for HMRs that are 50 bp or larger and create a new (temporary) file with a fourth column that lists HMR size; hmr.tmp and v2.hmr.tmp are the same file with the exception that v2.hmr.tmp has a fourth column listing HMR size
#Next, cross reference to RefSeq promoters and exons to identify HMRs in -a (Alpha) that don't overlap with -b (Beta) and aren't promoters or exons using NCBI RefSeq
#Don't use -v to identify overlap between Alpha and TSS/Exons
#Second filter step is also done for -b (Beta)
awk -v size="50" 'BEGIN{OFS=FS="\t"}{if (($3-$2)>size) {print $1,$2,$3,($3-$2)}}' Alpha_WGBS_rep1_S3.hmr.tmp > Alpha_WGBS_rep1_S3_v2.hmr.tmp  | bedtools intersect -v -a Alpha_WGBS_rep1_S3_v2.hmr.tmp -b Beta_WGBS_rep1_S2_v2.hmr.tmp ${NCBI} > Alpha_WGBS_rep1_S3.minsize50.filtforrefseqTSSexons.txt
awk -v size="50" 'BEGIN{OFS=FS="\t"}{if (($3-$2)>size) {print $1,$2,$3,($3-$2)}}' Beta_WGBS_rep1_S2.hmr.tmp >  Beta_WGBS_rep1_S2_v2.hmr.tmp | bedtools intersect -v -a Beta_WGBS_rep1_S2_v2.hmr.tmp -b Alpha_WGBS_rep1_S3_v2.hmr.tmp ${NCBI} > Beta_WGBS_rep1_S2.minsize50.filtforrefseqTSSexons.txt

#If running in the command-line...
#bedtools intersect -v -a Alpha_WGBS_rep1_S3_v2.hmr.tmp -b Beta_WGBS_rep1_S2_v2.hmr.tmp NCBI.RefSeq.TSS.Exons.merged.txt > Alpha_WGBS_rep1_S3.minsize50.filtforrefseqTSSexons.txt
#bedtools intersect -v -a Beta_WGBS_rep1_S2_v2.hmr.tmp -b Alpha_WGBS_rep1_S3_v2.hmr.tmp NCBI.RefSeq.TSS.Exons.merged.txt > Beta_WGBS_rep1_S2.minsize50.filtforrefseqTSSexons.txt
