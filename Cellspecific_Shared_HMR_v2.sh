##bedtools intersect allows one to screen for overlaps between 2 sets of genomic features and/or identify features unique to each cell type##
module restore anaconda
source activate bedtools
cd /data/hodges_lab/aganve/data/DNA_Methylation/7_tracks

##Considering the top BP for genes putatively regulated by beta-cell specific HMRs is the B lymphocyte signaling pathway, let's identify shared HMRs as well as those that are not shared between these cell types##
##Identify overlaps between beta-cell specific HMRs and B cell HMRs##
##30041 B cell HMRs; 14522 beta cell HMRs##
##306 overlaps##
#'-u'reports overlaps between -a (Alpha) and -b (Beta) --> writes original 'a' entry once if any overlaps found in 'b'
bedtools intersect -u -a Bcell.minsize50.filterforrefseqTSSexons.txt -b Beta_WGBS_rep1_S2.minsize50.filtforrefseqTSSexons.txt > CT\_all_shared_Bcell_Beta.txt

##Identify non-overlaps between beta-cell specific HMRs and B cell HMRs (i.e., beta-cell specific HMRs)##
##After removing HMRs that are B-cell specific, giving B cells were a top BP hit, will we see HMRs that are specific to beta cells?##
#29735 B-cell specific HMRs when compared to Beta cells; 14216 Beta-cell specific HMRs when compared to B cells##
#'-v' only reports those entries in -a (B cell) that have no overlap in -b (Beta)
bedtools intersect -v -a Bcell.minsize50.filterforrefseqTSSexons.txt -b Beta_WGBS_rep1_S2.minsize50.filtforrefseqTSSexons.txt > CT\_Bcell_specific.txt
bedtools intersect -v -a Beta_WGBS_rep1_S2.minsize50.filtforrefseqTSSexons.txt -b Bcell.minsize50.filterforrefseqTSSexons.txt > CT\_Beta_cell_specific.txt

##Use bedtools intersect to identify shared intervals between B cells and beta cells without -u flag##
##313 overlaps##
bedtools intersect -a Bcell.minsize50.filterforrefseqTSSexons.txt -b Beta_WGBS_rep1_S2.minsize50.filtforrefseqTSSexons.txt > CT\_all_shared_Bcell_Beta_v2.txt

##What happens if we subtract B lymphocyte genomic features from beta-cell genomic features? Will we get genomic features specific to beta cells?##
##14490 beta cell HMRs left after subtracting B cell HMRs##
#bedtools substract removes each overlapping interval in 'b' (i.e., B cell) from 'a' (i.e., Beta), leaving us with beta-cell specific HMRs (potentially?) if that wasn't identified in the previous step...
substractBed -a Beta_WGBS_rep1_S2.minsize50.filtforrefseqTSSexons.txt -b Bcell.minsize50.filterforrefseqTSSexons.txt >CT\_Beta_cell_specific_bedtoolsubstract.txt
