import sys
import csv
import pandas as pd
from matplotlib import pyplot as plt
from math import log10
from glob import glob

#negative control is bedtools shuffle
#positive control is TSS 
enrich_dir = sys.argv[1]
#shuffle_control_dir = sys.argv[2]
tissue = sys.argv[2]

pdx_tissue_dict = {
    "Bcell": "Cells_EBV-transformed_lymphocytes",
    "Liver": "Liver",
    "fHeart": "Heart_Left_Ventricle",
    "fSpinal": "Brain_Spinal_cord_cervical_c-1",
    "Adrenal": "Adrenal_Gland",
    "by_hmr_type": "abhls"}

pdx_tissue = pdx_tissue_dict[tissue]

enrich_files = glob("%s/%s/*" % (enrich_dir, tissue))
#shuffle_control_files = glob("%s/%s/*" % (shuffle_control_dir, tissue))
#positive_control_file = "/data/hodges_lab/doranser/results/enrichment_hmr_predixcan/hmr_controls/refseqGenes_TSS_1000_2000_Whole_Blood_enrich.txt"

enrich_df = pd.DataFrame(columns=["hmr_type", "fold_change", "p_value"])

for file in enrich_files:
    filename = file.split("/")[-1]
    hmr_type = filename.split("_%s" % pdx_tissue)[0]
    print(file)
    print(filename)
    print(hmr_type)
    with open(file) as f:
        records = csv.reader(f, delimiter="\t")
        next(records)
        next(records)
        for record in records:
            enrich_df = enrich_df.append({
                "hmr_type": hmr_type,
                "fold_change": float(record[3]),
                "p_value": float(record[6])
                }, ignore_index=True)

label_dict = {
    "all": "All",
    "all_cellspecific": "All CS",
    "all_shared": "All shared",
    "unclustered_cellspecific": "Unclustered CS",
    "unclustered": "Unclustered",
    "unclustered_shared": "Unclustered shared",
    "internalClusters": "Clusters",
#    "clustered_cellspecific": "Clusters indiv. HMRs CS",
    "internalClusters_thatContainCellSpecific": "Clusters containing CS",
    "internalClusters_containsCS": "Clusters containing CS",
    "internalClusters_individualHMRs": "Clusters indiv. HMRs",
    "internalClusters_individualHMRs_cellspecific": "Clusters indiv. HMRs CS",
    "internalClusters_individualHMRs_shared": "Clusters indiv. HMRs shared"
    }
print(enrich_df)

enrich_df["neg_log10_p"] = enrich_df.apply(lambda row: log10(row.p_value) * -1, axis=1)
enrich_df["plot_label"] = enrich_df.apply(lambda row: label_dict[row.hmr_type], axis=1)
rounded_df = enrich_df.round(decimals=4)

print(enrich_df)

sample_sizes = {
    "Bcell": 147,
    "Liver": 208,
    "fHeart": 386,
    "fSpinal": 126,
    "Adrenal": 233,
    }

groups = rounded_df.groupby("hmr_type")
for name, group in groups:
    plt.plot(group["fold_change"], group["neg_log10_p"], marker = "o", linestyle="", label=name, markersize=4.5)
plt.gca().set_xlim(left=0.6)
plt.gca().set_ylim(bottom=0)
plt.legend(sorted(rounded_df["plot_label"]), loc="lower left")
if tissue == "by_hmr_type":
    plt.title("HMR enrichment for PrediXcan DB SNPs across tissues")
else:
    plt.title("%s HMR enrichment for PrediXcan DB SNPs (GTEx n=%s)" % (tissue, sample_sizes[tissue]))
plt.xlabel("Enrichment (fold change)")
plt.ylabel("-log(p)")
plt.tight_layout()
plt.savefig("/data/hodges_lab/doranser/results/enrichment_hmr_predixcan/enrichment_plots/%s_HMR_pdx_enrichment_plot.pdf" % tissue)

with open("/data/hodges_lab/doranser/results/enrichment_hmr_predixcan/enrichment_tables/%s_HMR_pdx_enrichment_tbl.txt" % tissue, "w") as f:
    f.write(rounded_df.__repr__())
