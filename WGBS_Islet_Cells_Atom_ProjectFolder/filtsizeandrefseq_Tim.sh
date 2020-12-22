# Script to size filter, then cross ref vs RefSeq TSS/Exons
# 2 steps: #1 filt size, #2 filt for refseq
# usage: sh filtsizeandrefseq.sh [file] [size]
size="$2"
infile="$1"
infileName=${infile%.*}
outName=$infileName.minsize$size.filtforrefseqTSSexons.txt
awk -v size="$size" 'BEGIN{OFS=FS="\t"}{if (($3-$2)>size) {print $1,$2,$3,($3-$2)}}' $infile | bedtools intersect -v -a - -b NCBI.RefSeq.TSS.Exons.merged.txt > $outName

size="$2"
infile="$1"
NCBI='/data/hodges_lab/aganve/data/DNA_Methylation/5_hmrs/NCBI.RefSeq.TSS.Exons.merged.txt'
outName=$infile.minsize$size.filtforrefseqTSSexons.txt
awk -v size="$size" 'BEGIN{OFS=FS="\t"}{if (($3-$2)>size) {print $1,$2,$3,($3-$2)}}' $infile.hmr.bed > $infile.hmr_v2.bed | bedtools intersect -v -a $infile.hmr_v2.bed -b ${NCBI} > $outName
