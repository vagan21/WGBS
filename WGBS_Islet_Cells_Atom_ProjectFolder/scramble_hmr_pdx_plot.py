import csv
import sys
import re
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import colors
from glob import glob

scrambled_dir = sys.argv[1]
scrambled_files = glob("%s/*/*" % scrambled_dir)

matched_dir = sys.argv[2]

pdx_tissues = set()
hmr_subsets = set()
raw_hmr_tissues = set()
hmr_tissues = set()

pdx_tissue_dict = {
    "Bcell": "Cells EBV-transformed lymphocytes",
    "Liver": "Liver",
    "fHeart": "Heart Left Ventricle",
    "fSpinal": "Brain Spinal cord cervical c-1",
    "Tcell": "Cells EBV-transformed lymphocytes",
    "Adrenal": "Adrenal Gland"
    }

hmr_tissue_dict = {
    "Adrenal internalClusters individualHMRs cellspecific": "Adrenal clusters indiv. HMRs cellspecific",
    "Adrenal internalClusters containsCS": "Adrenal clusters thatContainCellSpecific",
    "Liver internalClusters containsCS": "Liver clusters thatContainCellSpecific",
    "Liver internalClusters individualHMRs cellspecific": "Liver clusters indiv. HMRs cellspecific",
    "Liver clustered cellspecific": "Liver clusters indiv. HMRs cellspecific",
    "Bcell clustered cellspecific": "Bcell clusters indiv. HMRs cellspecific",
    "Bcell internalClusters individualHMRs cellspecific": "Bcell clusters indiv. HMRs cellspecific",
    "Bcell internalClusters thatContainCellSpecific": "Bcell clusters thatContainCellSpecific",
    "fHeart internalClusters containsCS": "fHeart clusters thatContainCellSpecific",
    "fHeart internalClusters individualHMRs cellspecific": "fHeart clusters indiv. HMRs cellspecific",
    "fSpinal internalClusters containsCS": "fSpinal clusters thatContainCellSpecific",
    "fSpinal internalClusters individualHMRs cellspecific": "fSpinal clusters indiv. HMRs cellspecific"
    }

for f in scrambled_files:
    name = f.split("/")[-1].split(".")[0]
    name_split = re.split('_(?=[A-Z])', name)
    raw_hmr_tissues.add(f.split("/")[-1].split("_")[0])
    hmr_tissues.add(f.split("/")[-1].split("_")[0])
    hmr_filelabel = name_split[0].replace("_", " ")
    hmr_subset = hmr_tissue_dict.setdefault(hmr_filelabel, hmr_filelabel)
    hmr_subsets.add(hmr_subset)
    pdx_tissues.add("_".join(name_split[1:]).split("_pdx")[0].replace("_", " "))

pdx_tissues = sorted(list(pdx_tissues))
hmr_subsets = sorted(list(hmr_subsets))
hmr_tissues = sorted(list(hmr_tissues))

enrich_df = pd.DataFrame(columns=pdx_tissues, index = hmr_subsets)
sig_df = pd.DataFrame(columns=pdx_tissues, index = hmr_subsets)

cs_hmr_subsets = []
for subset in hmr_subsets:
    if "cellspecific" in subset.lower():
        cs_hmr_subsets.append(subset)
cs_hmr_subsets.sort()
print(cs_hmr_subsets)

cs_enrich_df = pd.DataFrame(columns=pdx_tissues, index=cs_hmr_subsets)
cs_sig_df = pd.DataFrame(columns=pdx_tissues, index=cs_hmr_subsets)

for f in scrambled_files:
    name = f.split("/")[-1].split(".")[0]
    name_split = re.split('_(?=[A-Z])', name)
    hmr_tissue = "_".join(f.split(("/")[-1].split("_")[0]))
    hmr_filelabel = name_split[0].replace("_", " ")
    hmr_subset = hmr_tissue_dict.setdefault(hmr_filelabel, hmr_filelabel)
    pdx_tissue = "_".join(name_split[1:]).split("_pdx")[0].replace("_", " ")
    with open(f) as file:
        records = csv.reader(file, delimiter="\t")
        next(records)
        next(records)
        for record in records:
            enrich_df.at[hmr_subset, pdx_tissue] = float(record[3])
            sig_df.at[hmr_subset, pdx_tissue] = float(record[6])
            if hmr_subset in cs_hmr_subsets:
                cs_enrich_df.at[hmr_subset, pdx_tissue] = float(record[3])
                cs_sig_df.at[hmr_subset, pdx_tissue] = float(record[6])

for hmr_tissue in raw_hmr_tissues:
    pdx_tissue = pdx_tissue_dict[hmr_tissue]
    for f in glob("%s/%s/*" % (matched_dir, hmr_tissue)):
        name = f.split("/")[-1].split(".")[0]
        name_split = re.split('_(?=[A-Z])', name)
        hmr_filelabel = name_split[0].replace("_", " ")
        hmr_subset = hmr_tissue_dict.setdefault(hmr_filelabel, hmr_filelabel)
        print(hmr_subset)
        with open(f) as file:
            records = csv.reader(file, delimiter="\t")
            next(records)
            next(records)
            for record in records:
                enrich_df.at[hmr_subset, pdx_tissue] = float(record[3])
                sig_df.at[hmr_subset, pdx_tissue] = float(record[6])
                if hmr_subset in cs_hmr_subsets:
#                    print(hmr_subset, pdx_tissue)
                    cs_enrich_df.at[hmr_subset, pdx_tissue] = float(record[3])
       	       	    cs_sig_df.at[hmr_subset, pdx_tissue] = float(record[6])

cols_order = ["Adrenal Gland", "Cells EBV-transformed lymphocytes", "Liver", "Heart Left Ventricle", "Brain Spinal cord cervical c-1"]

enrich_df.sort_index(axis=0, ascending=False, inplace=True)
cs_enrich_df.sort_index(axis=0, ascending=False, inplace=True)

rounded_df = enrich_df.round(4)
rounded_df = rounded_df.fillna(0)
rounded_df = rounded_df[cols_order]

cs_rounded_df = cs_enrich_df.round(4)
cs_rounded_df = cs_rounded_df.fillna(0)
cs_rounded_df = cs_rounded_df[cols_order]

#norm = colors.DivergingNorm(vcenter=1)

plt.figure(figsize=(8, 11))
plt.pcolor(rounded_df, cmap=plt.get_cmap("seismic"), vmin =0.17, vmax=2.1)
#plt.colorbar()
plt.yticks(np.arange(0.5, len(rounded_df.index), 1), rounded_df.index)
plt.xticks(np.arange(0.5, len(rounded_df.columns), 1), rounded_df.columns, rotation="vertical")
plt.title("Enrichment of PrediXcan SNPs in HMRs")
plt.colorbar()
plt.tight_layout()
plt.savefig("/data/hodges_lab/doranser/results/enrichment_hmr_predixcan/enrichment_plots/scrambled_hmr_pdx_plot.pdf")

plt.figure(figsize=(7, 5))
plt.pcolor(cs_rounded_df, cmap=plt.get_cmap("YlOrRd"), vmin=0.9, vmax=2.1)
plt.yticks(np.arange(0.5, len(cs_rounded_df.index), 1), cs_rounded_df.index, fontsize=8)
plt.xticks(np.arange(0.5, len(cs_rounded_df.columns), 1), cs_rounded_df.columns, rotation=40, ha="right", fontsize=8)
plt.title("Enrichment of HMRs for PrediXcan SNPs")
plt.colorbar()
plt.tight_layout()
plt.savefig("/data/hodges_lab/doranser/results/enrichment_hmr_predixcan/enrichment_plots/scrambled_hmr_pdx_cs_only_ylorred_plot.pdf")
for col in rounded_df:
    for i, row_value in rounded_df[col].iteritems():
        if sig_df[col][i] < 0.05:
            rounded_df[col][i] = "%s*" % rounded_df[col][i]
rounded_df.to_csv("/data/hodges_lab/doranser/results/enrichment_hmr_predixcan/enrichment_tables/scrambled_hmr_pdx_table.txt", sep="\t")

for col in cs_rounded_df:
    for i, row_value in cs_rounded_df[col].iteritems():
        if cs_sig_df[col][i] < 0.05:
            cs_rounded_df[col][i] = "%s*" % cs_rounded_df[col][i]
cs_rounded_df.to_csv("/data/hodges_lab/doranser/results/enrichment_hmr_predixcan/enrichment_tables/scrambled_hmr_pdx_cs_only_table.txt", sep="\t")

#with open("/data/hodges_lab/doranser/results/enrichment_hmr_predixcan/enrichment_tables/scrambled_hmr_pdx_table.txt", "w") as f:
#    f.write(rounded_df.__repr__())
