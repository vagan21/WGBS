#!/bin/bash
# annotate_HMRs.sh


# USAGE: annotate_HMRs.sh [CELLTYPE of interest] [CELLTYPE BED] [OTHER CELLTYPES BED] [INVERSE REGIONS BED FOR CLUSTERING] [OUTPUT DIRECTORY]
# [CELLTYPE]: Name of celltype (i.e. "Bcell")
# [CELLTYPE BED]: BED FILE of HMRs in cell type of interest
# [OTHER CELLTYPES BED]: BED FILE of HMRs in other celltypes (i.e. for B cells: Liver, Heart, and H1ESCs)
# [BLACKLIST REGIONS FOR CLUSTER]: BED FILE of regions that are boundaries for clustering: (i.e. NCBI.RefSeq.TSS.Exons.txt)
# [OUTPUT DIRECTORY]: The directory for outputing files

# NOTE WARNING!: IF ANY FILE IS NULL, REPLACE WITH blank.txt, an empty text file

# sh auto_annotate_hmr_internalClusters.sh Bcell Bcell.HMR.min50.txt Tcell.HMR.min50.txt Macrophage.HMR.min50.txt Rvent.Liver.H1ESC.min50.txt NCBI.RefSeq.TSS.Exons.merged.txt /Users/scottt7/Desktop/HMR_BEDfiles

# OUTPUT (11): Clustered/Unclustered files further delineated by cell-type specificity - up to four different levels of celltype membership binning
# (1) CT_unclustered
# (2) CT_unclustered_cellspecific
# (3) CT_unclustered_shared

# (4) CT_clustered_endToEnd
# (5) CT_clustered_individualHMRs
# (6) CT_clustered_individualHMRs_cellspecific
# (7) CT_clustered_individualHMRs_shared


######################       ISSUE REGIONS       ############################
# chr10	22623401	22623588 - Says Lineage, but one of those with HMRs spannign TSS/Exon
# Remake "Other" file:
#			cat Rvent.HMR.min50.txt Liver.HMR.min50.txt H1ESC.HMR.min50.txt | bedtools sort -i > Rvent.Liver.H1ESC.min50.txt

################################################################################
# Expects 5 cell type BED files as input; Declare as variables
CT=$1
CTFile=$2
OthFile=$3
invFile=$4
outDir=$5
TSSExonRef="/scratch/scottt7/paper_wases/HMRs/NCBI.RefSeq.TSS.Exons.merged.only24chrom.txt"


# Echo the variable assignments for good measure
echo "Cell-type identified: ${CT}"
echo "Cell-type file identified: ${CTFile}"
echo "Other cell-type file identified: ${OthFile}"
echo "Clustering Regions Whitelist identified: ${invFile}"
echo "Output Directory identified: ${outDir}"
echo $"\n"

################################################################################
# Find Unclustered HMRs

##### ((1)) CT_unclustered
##### 
# Filt for TSS/Exons, merge, and filter for unmerged ones 
bedtools intersect -v -a $CTFile -b $TSSExonRef | bedtools merge -c 2 -o count -d 6000 -i - | awk 'BEGIN{FS=OFS="\t"}{if ($4<2) print}' -  | awk 'BEGIN{OFS=FS="\t"}{print $1,$2,$3}' - > $outDir/$CT\_unclustered.txt
echo "Found unclustered HMRs."


##### ((2)) CT_unclustered_cellspecific
##### 
# Inverse intersect (1) with comparison file to find (a) ct-specific unclustered HMRs
bedtools intersect -v -a $outDir/$CT\_unclustered.txt -b $OthFile | awk 'BEGIN{OFS=FS="\t"}{print $1,$2,$3}' - > $outDir/$CT\_unclustered_cellspecific.txt
echo "Found cell-specific unclustered HMRs."


##### ((3) CT_unclustered_shared
##### 
# Intersect (1) with comparison file to find (b) ct-shared unclustered HMRs
bedtools intersect -u -a $outDir/$CT\_unclustered.txt -b $OthFile | awk 'BEGIN{OFS=FS="\t"}{print $1,$2,$3}' - > $outDir/$CT\_unclustered_shared.txt
echo "Found shared unclustered HMRs."


################################################################################
# Find Clustered HMRs

##### ((4)) CT_internalClusters
##### 
# internal Clusters that don't cross at TSS/Exon, that also compose 3+ HMRs
# Strategy: 
#		- (pre-a) Create a whitelist from the TSS/Exon blacklist
# 		- (a) Find all inverse regions (Whitelist regions) that contain 3+ HMRs (potential clusters)
# 		- (b) Make BED that has two BED coordinates per line - Left: Whitelist Region; Right: HMR region
#		- (c) Merge using this altered BED-BED file and a custom awk script instead of BEDtools Merge

# (pre-a) Turn blacklist into whitelist
# To download the hg19 chr size file: 
# mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from hg19.chromInfo"  > hg19.genome
# awk 'BEGIN{OFS=FS="\t"}{print $1,$2,$2}' /scratch/scottt7/paper_wases/HMRs/hg19.genome | bedtools sort -i - > /scratch/scottt7/paper_wases/HMRs/hg19.genome_sorted.txt
bedtools sort -i $TSSExonRef | bedtools complement -i - -g /scratch/scottt7/paper_wases/HMRs/hg19.chrom.sizes.twoCol.sorted > $CT\_whitelist.txt
echo "Made whitelist."


# (a) Find whitelist regions with 3+ HMRs
bedtools intersect -c -F 1.0 -a $CT\_whitelist.txt -b $CTFile | awk 'BEGIN{OFS=FS="\t"}{if ($4>2) print}' - > $outDir/$CT\_whitelistRegions_contain3ormoreHMRs.txt
echo "Found whitelist regions that contain 3+ HMRs."

# (b)Find what HMRs are in these to create a double BED file in a parsable list for custom merging: [1-3]: BED coordinates of Whitelist region [4]: #HMR in Whitelist region, [5-7]: BED coordinates of HMRs in Whitelist region
# Diagnostic: bedtools intersect -loj -F 1.0 -a Bcell_whitelistRegions_contain3ormoreHMRs.txt -b Bcell.HMR.min50.txt | awk 'BEGIN{OFS=FS="\t"}{print}' - | bedtools sort -i - | head
# awk 'BEGIN{OFS=FS="\t"}{print $5,$6,$7,$1,$2,$3,$4}' temp.txt | head

bedtools intersect -loj -F 1.0 -a $outDir/$CT\_whitelistRegions_contain3ormoreHMRs.txt -b $CTFile | awk 'BEGIN{OFS=FS="\t"}{print $5,$6,$7,$1,$2,$3,$4}' - | bedtools sort -i - > $outDir/$CT\_whitelistRegions_contain3ormoreHMRs_individualHMRs_doubleBED_rev.txt
echo "Created double BED for whitelist regions that contain 3+ HMRs."

# (c) Custom merge within these regions
awk -v whitelistStart=1 -v hmrsCount=1 -v clusterChr=1 -v clusterStart=1 -v clusterEnd=2 'BEGIN{OFS=FS="\t";dist=6000;whitelistStart=1;hmrsCount=1;clusterChr=1;clusterStart=1;clusterEnd=2;} {
	# Check if we are in the same WhitelistBoundary
	if ($5!=whitelistStart) {
		print clusterChr,clusterStart,clusterEnd,hmrsCount;
		whitelistStart=$5;
		hmrsCount=1;
		clusterChr=$1;
		clusterStart=$2;
		clusterEnd=$3;
	} else {
		# case if (1) we are staying within a Whitelist boundary, but (2) the previous HMR was >6000bp away
		if (($2-clusterEnd)>6000) {
			print clusterChr,clusterStart,clusterEnd,hmrsCount;
			clusterStart=$2;
			clusterEnd=$3;
			hmrsCount=1;
		} else {
			clusterEnd=$3;
			hmrsCount+=1;
		}
	}
}' $outDir/$CT\_whitelistRegions_contain3ormoreHMRs_individualHMRs_doubleBED_rev.txt | awk 'BEGIN{OFS=FS="\t"}{if ($4>2) print $0}' - > $outDir/$CT\_internalClusters.txt
echo "Found internalClusters containing 3+ HMRs."

#########################
#########################
#########################
# What does this code do? 
awk 'BEGIN{OFS=FS="\t"}{if ($1 !~ /\_/)print $1,$2,$3}' $outDir/$CT\_internalClusters.txt > $outDir/$CT\_internalClusters_BED.txt
echo $"\n"


##### ((5)) internalClusters_individualHMRs 
#####
bedtools intersect -u -a $CTFile -b $outDir/$CT\_internalClusters_BED.txt > $outDir/$CT\_internalClusters_individualHMRs.txt
echo "Found the HMRs composing internalClusters containing 3+ HMRs."

##### ((6)) internalClusters_individualHMRs_cellspecific
##### 
# Inverse intersect (1) with our comparison file to find (a) ct-specific and (b) at least partially shared regions
bedtools intersect -v -a $outDir/$CT\_internalClusters_individualHMRs.txt -b $OthFile > $outDir/$CT\_internalClusters_individualHMRs_cellspecific.txt
echo "Found cell-specific clustered HMRs."


##### ((7) CT_internalClusters_individualHMRs_sharedBeyondLineage
##### 
#	Intersect CT with OtherFile to find shared
bedtools intersect -u -a $outDir/$CT\_internalClusters_individualHMRs.txt -b $OthFile > $outDir/$CT\_internalClusters_individualHMRs_shared.txt
echo "Found shared clustered HMRs."







# test Other: cat H1ESC.HMR.min50.txt Liver.HMR.min50.txt Rvent.HMR.min50.txt | bedtools sort -i - | bedtools intersect -v -a - -b NCBI.RefSeq.TSS.Exons.merged.txt | bedtools merge -i - > H1ESC.Liver.Rvent.HMR.min50.filtForRefSeqTSSExons.txt
# Before:    61384 H1ESC.Liver.Rvent.HMR.min50.filtForRefSeqTSSExons.txt
# After not merging:   63434 H1ESC.Liver.Rvent.HMR.min50.filtForRefSeqTSSExons.txt
# | awk '$1 ~ /chr9/' - | cat
# > H1ESC.Liver.Rvent.HMR.min50.filtForRefSeqTSSExons.txt
# test: 
# sh auto_annotate_hmr_internalClusters.sh Bcell Bcell.HMR.min50.txt Tcell.HMR.min50.txt Macrophage.HMR.min50.txt H1ESC.Liver.Rvent.HMR.min50.filtForRefSeqTSSExons.txt blank.txt /Users/scottt7/Desktop/HMR_BEDfiles

# wc -l Bcell_unclustered_sharedBeyondLineage.txt Bcell_unclustered_sharedLineage.txt Bcell_unclustered_sharedSubLineage.txt Bcell_unclustered_sharedAtAll.txt 
# wc -l 

# wc -l Bcell_internalClusters.txt Bcell_internalClusters_individualHMRs.txt Bcell_internalClusters_individualHMRs_sharedSubLineage.txt Bcell_internalClusters_individualHMRs_sharedLineage.txt Bcell_internalClusters_individualHMRs_sharedBeyondLineage.txt Bcell_internalClusters_individualHMRs_sharedAtAll.txt Bcell_internalClusters_individualHMRs_cellspecific.txt




# Bcell_internalClusters_individualHMRs_cellspecific.txt Bcell_internalClusters_individualHMRs_sharedSubLineage.txt Bcell_internalClusters_individualHMRs_sharedLineage.txt Bcell_internalClusters_individualHMRs_sharedBeyondLineage.txt Bcell_internalClusters_individualHMRs_sharedAtAll.txt 





