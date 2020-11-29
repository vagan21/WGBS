##bedtools intersect allows one to screen for overlaps between 2 sets of genomic features and/or identify features unique to each cell type##
module restore anaconda
source activate bedtools
cd /data/hodges_lab/aganve/data/DNA_Methylation/7_tracks

##Cell-specific HMRs for Alpha and Beta cells##
#'-v' only reports those entries in -a (Alpha) that have no overlap in -b (Beta); likewise do the same for -b (Beta) -- report HMRs that are B-cell specific
bedtools intersect -v -a Alpha_WGBS_rep1_S3.minsize50.filtforrefseqTSSexons.txt -b Beta_WGBS_rep1_S2.minsize50.filtforrefseqTSSexons.txt > CT\_A_cellspecific.txt
bedtools intersect -v -a Beta_WGBS_rep1_S2.minsize50.filtforrefseqTSSexons.txt -b Alpha_WGBS_rep1_S3.minsize50.filtforrefseqTSSexons.txt > CT\_B_cellspecific.txt

##Shared HMRs across Alpha and Beta cells##
#'-u'reports overlaps between -a (Alpha) and -b (Beta)
bedtools intersect -u -a Alpha_WGBS_rep1_S3.minsize50.filtforrefseqTSSexons.txt -b Beta_WGBS_rep1_S2.minsize50.filtforrefseqTSSexons.txt > CT\_all_shared.txt

##No filtering for TSS/Exon/Size##
bedtools intersect -v -a Alpha_WGBS_rep1_S3.hmr -b Beta_WGBS_rep1_S2.hmr > CT\_A_cellspecific_nofilter.txt
bedtools intersect -v -a Beta_WGBS_rep1_S2.hmr -b Alpha_WGBS_rep1_S3.hmr > CT\_B_cellspecific_nofilter.txt

##Shared HMRs across Alpha and Beta cells with no filtering for TSS/Exon/Size##
bedtools intersect -u -a Alpha_WGBS_rep1_S3.hmr -b Beta_WGBS_rep1_S2.hmr > CT\_all_shared_nofilter.txt

##Identify Shared HMRs that aren't TSSs or Exons##
bedtools intersect -v -a CT\_all_shared_nofilter.txt -b /data/hodges_lab/aganve/data/DNA_Methylation/7_tracks/NCBI.RefSeq.TSS.Exons.merged.txt > CT\_all_shared.filtforrefseqTSSexons.txt
