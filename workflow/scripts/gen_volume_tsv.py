import nibabel as nib
import pandas as pd
import numpy as np

lookup_df = pd.read_table(snakemake.input.lookup_tsv,index_col='index')

#get indices and names from lookup table
indices = lookup_df.index.to_list()
names = lookup_df.abbreviation.to_list()
hemis = ['L','R']

#create the output dataframe

df = pd.DataFrame(columns = ['subject','hemi'] + names)


for in_img,hemi in zip(snakemake.input.segs,hemis):

    img_nib = nib.load(in_img)
    img = img_nib.get_fdata()
    zooms = img_nib.header.get_zooms()

    #voxel size in mm^3
    voxel_mm3 = np.prod(zooms)

    new_entry = {'subject': 'sub-{subject}'.format(subject = snakemake.wildcards['subject']),'hemi': hemi}
    for index,name in zip(indices,names):
        # add volume as value, name as key
        new_entry[name] = np.sum(img==index) * voxel_mm3

    #now create a dataframe from it
    df = df.append(new_entry, ignore_index=True)
   
df.to_csv(snakemake.output.tsv,sep='\t', index=False)
