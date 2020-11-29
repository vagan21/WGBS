# make_anti.sh
# Usage: sh make_anti.sh [other cell type HMR files]
# Can be any length
Anchor_CT=$1
shift
Other_Files="$@"

cat $Other_Files | bedtools sort -i - | bedtools merge -i - > anti_$Anchor_CT.txt
