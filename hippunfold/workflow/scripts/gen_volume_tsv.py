import nibabel as nib
import numpy as np
import pandas as pd

lookup_df = pd.read_table(snakemake.input.lookup_tsv, index_col="index")

# get indices and names from lookup table
indices = lookup_df.index.to_list()
names = lookup_df.abbreviation.to_list()
hemis = ["L", "R"]

# collect output rows
rows = []

for in_img, hemi in zip(snakemake.input.segs, hemis):
    img_nib = nib.load(in_img)
    img = img_nib.get_fdata()
    zooms = img_nib.header.get_zooms()

    # voxel size in mm^3
    voxel_mm3 = np.prod(zooms)

    new_entry = {
        "subject": "sub-{subject}".format(subject=snakemake.wildcards["subject"]),
        "hemi": hemi,
    }

    for index, name in zip(indices, names):
        # add volume as value, name as key
        new_entry[name] = np.sum(img == index) * voxel_mm3

    rows.append(new_entry)

# create dataframe from collected rows
df = pd.DataFrame(rows, columns=["subject", "hemi"] + names)

df.to_csv(snakemake.output.tsv, sep="\t", index=False)
