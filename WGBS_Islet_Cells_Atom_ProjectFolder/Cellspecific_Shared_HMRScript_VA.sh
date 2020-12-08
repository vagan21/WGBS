##bedtools intersect allows one to screen for overlaps between 2 sets of genomic features and/or identify features unique to each cell type##
module restore anaconda
source activate bedtools
cd /data/hodges_lab/aganve/data/DNA_Methylation/7_tracks

##Cell-specific HMRs for Alpha and Beta cells##
#'-v' only reports those entries in -a (Alpha) that have no overlap in -b (Beta); likewise do the same for -b (Beta) -- report HMRs that are beta-cell specific
bedtools intersect -v -a Alpha_WGBS_rep1_S3.minsize50.filtforrefseqTSSexons.txt -b Beta_WGBS_rep1_S2.minsize50.filtforrefseqTSSexons.txt > CT\_A_cellspecific.txt
bedtools intersect -v -a Beta_WGBS_rep1_S2.minsize50.filtforrefseqTSSexons.txt -b Alpha_WGBS_rep1_S3.minsize50.filtforrefseqTSSexons.txt > CT\_B_cellspecific.txt

##Shared HMRs across Alpha and Beta cells##
#'-u'(unique) reports the mere presence of any overlapping features between -a (Alpha) and -b (Beta), with the 'a' feature being reported once if any overlaps found in -b
#in other words, just report the fact at least one overlap was found in -b
bedtools intersect -u -a Alpha_WGBS_rep1_S3.minsize50.filtforrefseqTSSexons.txt -b Beta_WGBS_rep1_S2.minsize50.filtforrefseqTSSexons.txt > CT\_all_shared.txt

------------------------------------
##Filtering for size (>= 50 bp) but not TSS/Exons#
##Cell-specific HMRs for Alpha and Beta cells##
bedtools intersect -v -a Alpha_WGBS_rep1_S3_v2.hmr.tmp -b Beta_WGBS_rep1_S2_v2.hmr.tmp > CT\_A_cellspecific_nofilter.txt
bedtools intersect -v -a Beta_WGBS_rep1_S2_v2.hmr.tmp -b Alpha_WGBS_rep1_S3_v2.hmr.tmp > CT\_B_cellspecific_nofilter.txt

##Shared HMRs across Alpha and Beta cells with no filtering for TSS/Exons##
bedtools intersect -u -a Alpha_WGBS_rep1_S3_v2.hmr.tmp -b Beta_WGBS_rep1_S2_v2.hmr.tmp > CT\_all_shared_nofilter.txt

##Identify Shared HMRs that were originally not filtered for TSS/Exons and now cross reference to NCBI RefSeq filter to identify HMRs unique to shared HMRs#
bedtools intersect -v -a CT\_all_shared_nofilter.txt -b NCBI.RefSeq.TSS.Exons.merged.txt > CT\_all_shared.filtforrefseqTSSexons.txt

##Identify HMRs that overlap between shared HMRs and NCBI RefSeq TSS/Exons##
bedtools intersect -u -a CT\_all_shared_nofilter.txt -b NCBI.RefSeq.TSS.Exons.merged.txt > CT\_all_shared_nonunique_overlaps_with_refseq.txt
