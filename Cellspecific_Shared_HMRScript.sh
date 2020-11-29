

##### Cell specific 
bedtools intersect -v -a CT_ofInterest -b Other_CT > CT\_all_cellspecific.txt


##### ((1)) ALL_shared
bedtools intersect -u -a CT_ofInterest -b $Other_CT > CT\_all_shared.txt