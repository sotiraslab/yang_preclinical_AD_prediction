from glob import glob
import os
import shutil

import pandas as pd


def find_png(base_dir):
    
    png_files = glob(os.path.join(base_dir, "**/*.png"), recursive = True)
    png_df = pd.DataFrame(png_files, columns = ["path"])
    png_df["id"] = png_df.loc[:,"path"].str.replace(base_dir + "/", "").str.split("/", expand = True).iloc[:,0]
    return png_df

def move_png(png_df, odir, group_col = "id", file_col = "path"):

    if not os.path.exists(odir):
        os.makedirs(odir, exist_ok = True)

    # Create group-specific directories and copy files
    for _, row in png_df.iterrows():

        group_id = row[group_col]
        file_path = row[file_col]

        # Create the group's folder if it doesn't exist
        group_folder = os.path.join(odir, group_id)
        os.makedirs(group_folder, exist_ok=True)

        # Copy the file to the group's folder
        if os.path.exists(file_path):
            shutil.copy(file_path, group_folder)
        else:
            print(f"File not found: {file_path}")

base_dir = "/ceph/chpc/shared/aristeidis_sotiras_group/b.y.yang_scratch/pet_analysis/nipype"
png_df = find_png(base_dir)
move_png(png_df, "/ceph/chpc/shared/aristeidis_sotiras_group/b.y.yang_scratch/pet_analysis/qc")
