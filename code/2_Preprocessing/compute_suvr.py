# ==============================================================================

# Name: compute_suvr.py
# Author: Murat Bilgel, Braden Yang
# Description: this workflow is for computing SUVR for amyloid PET scans

# Notes
# - this script relies on code from Murat Bilgel's PET-MRI processing pipeline;
#   to gain access to his python package, please contact Murat
# - we are using the recommended time windows to compute an average PET image
#   - PIB: 30-60 min post injection
#   - FBP: 50-70 min post injection
# - the input CSV should contain the following:
#   - path to PET image
#   - path to PET JSON
#   - path to T1 image
#   - path to MUSE ROI mask
#   - path to brain mask

# TODO:
# - check that CSV has all necessary columns

# ==============================================================================

# =================
INTERACTIVE = False
# =================

# =================
from time import time
start_time = time()
# =================

import argparse
import json
import os
import sys

# import pandas as pd
from nipype import config, logging

sys.path.append("/home/b.y.yang/sotiraslab/BradenADLongitudinalPrediction/code")
from project_utils.pet.helper import *

# ===========================
# ===== PARSE ARGUMENTS =====
# ===========================

# parse arguments
parser = argparse.ArgumentParser(
    description="run Murat Bilgel's (modified) PET-MRI processing pipeline to compute regional MUSE SUVRs",
    formatter_class=argparse.RawDescriptionHelpFormatter
)

data_input = parser.add_mutually_exclusive_group()
data_input.add_argument("-d", "--data_table", help="path to spreadsheet containing paths to necessary files; if specified, then multithreaded processing will be used to process multiple subjects")
data_input.add_argument("-i", "--id", type = str, help="subject ID; if specified, then single threaded processing will be performed on a single subject")

# NOTE: these arguments will be ignored if --data_table is specified
parser.add_argument("--petpath", type = str, default = None, help="path to PET NIFTI file")
parser.add_argument("--petjsonfile", type = str, default = None, help="path to PET JSON file")
parser.add_argument("--musemaskpath", type = str, default = None, help="path to brain mask NIFTI file")
parser.add_argument("--musemriskullpath", type = str, default = None, help="path to T1 NIFTI file (with skull)")
parser.add_argument("--muselabelpath", type = str, default = None, help="path to MUSE ROI label NIFTI file")
parser.add_argument("--manual_coreg_mat", type = str, default = None, help="if PET-MRI coregistration has been manually corrected, input path to the manual affine matrix here")

parser.add_argument("-o", "--output_dir", help="output directory")
parser.add_argument("--n_procs", type=int, default=8, help="number of processes for multithreaded processing; ignored if single subject is specified")
parser.add_argument("--t_start", type=float, help="start time for averaging PET image")
parser.add_argument("--t_end", type=float, help="end time for averaging PET image")
parser.add_argument("--realign_t_start", type=float, default=None, help="(optional) start time for computing an average image as \
                    reference for realignment")
parser.add_argument("--realign_t_end", type=float, default=None, help="(optional) end time for computing an average image as \
                    reference for realignment")
parser.add_argument("--reg_smooth_fwhm", type=float, default=4, help="(optional) FWHM of smoothing operation if you want to \
                    smooth before coregistration or realignment; note that this does *not* smooth the final SUVR image, \
                    it only creates a smoothed intermediate image to help with registration steps")
parser.add_argument("--iterative_smoothing", action="store_true", help="apply iterative smoothing to SUVR image")
parser.add_argument("--pet_is_3d", action="store_true", help="specify if PET image is 3D and has already been averaged")
parser.add_argument("--pvc", action="store_true", help="apply PVC correction")
parser.add_argument("--modify_json_func", type=str, default=None, help="name of python function to modify PET JSON files")
parser.add_argument("--overwrite", action="store_true", help="run pipeline and overwrite existing outputs anyways even if SUVR CSV already exists")
parser.add_argument("--create_pet_json", action="store_true", help="create a dummy JSON header file if one doesn't exist for the PET you are trying to process; assumed that the input --petjsonfile (or corresponding column in --data_table) contains the path to the dummy JSON to create")
parser.add_argument("--create_pet_json_n_frames", type=int, default=4, help="number of frames for dummy JSON header file")
parser.add_argument("--create_pet_json_frame_duration_min", type=float, default=5, help="duration of each frame for dummy JSON header file")

if INTERACTIVE:
    args = parser.parse_args(args = [])
    # define preset arguments
    args.data_table = "/home/b.y.yang/sotiraslab/BradenADLongitudinalPrediction/data/subj/proc/petanalysis.csv"
    args.output_dir = "/scratch/b.y.yang/pet_analysis_test"
    args.n_procs = 1
    args.t_start = 0
    args.t_end = 99999
    args.pet_is_3d = True
else:
    args = parser.parse_args()

# handle output directory
if not args.id is None:
    output_dir = os.path.join(args.output_dir, args.id)  # create subdirectory for single subject
else:
    output_dir = args.output_dir
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# handle modify_json_func argument
if not args.modify_json_func is None:
    if args.modify_json_func in globals():
        modify_json_func = globals()[args.modify_json_func]  # Get function reference
    else:
        print(f"Error: Function '{args.modify_json_func}' not found.")
else:
    modify_json_func = None

# print args
print("")
print("+++++++++++++++++++++++++++++++")
print("ARGUMENTS")
print("---------")
for arg, value in vars(args).items():
    print(f"{arg}: {value}")
print("+++++++++++++++++++++++++++++++")
print("")

# check if outputs already generated
suvr_odir = os.path.join(output_dir, "wf/suvr_csv/suvr_csv")
if not args.overwrite and os.path.exists(suvr_odir) and any(file.endswith('.csv') for file in os.listdir(suvr_odir)):
    print("SUVR CSV already generated for this subject; exiting")
    sys.exit(0)

# =====================
# ===== CONFIGURE =====
# =====================

# config.enable_debug_mode()
config.update_config({'execution': {'stop_on_first_crash': True}})
logging.update_logging(config)

# list of subjects to manually exclude
# TODO: move elsewhere (e.g in config or params file)
subj_exclude = [
    "OASIS_OAS30053_64.37_PIB",  # missing frame timing info (could guess)
    "OASIS_OAS30143_68.58_FBP",  # missing frame timing info (could guess)
    "OASIS_OAS30321_73.99_FBP",  # missing frame timing info (could guess)
    "OASIS_OAS31381_71.27_PIB",  # only 3 frames out of 26,
    "OASIS_OAS30253_65.82_PIB",  # only 1 frame out of 26
    "A4_B54288722_71.89_FBP",    # only 3 frames out of 4
    "A4_B96962266_85.28_FBP",    # only 1 frame out of 4
    "A4_B65206013_75.96_FBP",    # only 1 frame out of 4
    "A4_B30242622_73.32_FBP",    # only 1 frame out of 4
    "A4_B18964625_73.21_FBP",    # weird error (ValueError: substring not found, SPM step)
    "A4_B99716486_71.14_FBP",    # 3 out of 4 frames, but modify JSON func is adding 4 frames to timing file (easy fix)
]

# =====================
# ===== LOAD DATA =====
# =====================

if not args.data_table is None:

    # load data table
    data_table, missing_files = load_data_table(args.data_table)
    data_table["manual_coreg_mat"] = None  # add column for manual coregistration matrix

    # save data_table and missing_files to spreadsheet
    data_table.to_csv(os.path.join(output_dir,'data_table.csv'))
    missing_files.to_csv(os.path.join(output_dir,'missing_data.csv'))

    # manually exclude subjects
    data_table = data_table[~data_table["id"].isin(subj_exclude)]

else:

    # manually exclude subjects
    if args.id in subj_exclude:
        raise ValueError(f"Subject {args.id} is in the list of subjects to exclude from PET analysis")

    # create data table with single subject
    args_dict = {
        "id": args.id,
        "petpath": args.petpath,
        "petjsonfile": args.petjsonfile,
        "musemaskpath": args.musemaskpath,
        "musemriskullpath": args.musemriskullpath,
        "muselabelpath": args.muselabelpath,
        "manual_coreg_mat": args.manual_coreg_mat,
    }
    data_table = pd.Series(args_dict)
    with open(os.path.join(output_dir, 'args.json'), 'w') as f:
        json.dump(args_dict, f)

# ========================
# ===== RUN PIPELINE =====
# ========================

# create JSON files
if args.create_pet_json:
    def map_create_pet_json(petjsonfile):
        create_pet_json(
            args.create_pet_json_n_frames,
            args.create_pet_json_frame_duration_min,
            os.path.dirname(petjsonfile),
            os.path.basename(petjsonfile).replace(".json", "")
        )
    if isinstance(data_table, pd.DataFrame):
        data_table["petjsonfile"].map(map_create_pet_json)
    else:
        map_create_pet_json(data_table["petjsonfile"])

wf = build_suvr_pipeline(
    data_table,
    output_dir,
    t_start = args.t_start,
    t_end = args.t_end,
    realign_t_start = args.realign_t_start,
    realign_t_end = args.realign_t_end,
    reg_smooth_fwhm = args.reg_smooth_fwhm,
    pet_is_4D = not args.pet_is_3d,
    modify_json_func = modify_json_func,
    iterative_smoothing = args.iterative_smoothing,
)

# run pipeline
wf.write_graph('pipeline.dot', graph2use='colored', simple_form=True)
if not args.data_table is None:
    result = wf.run('MultiProc', plugin_args={'n_procs': args.n_procs})
else:
    result = wf.run()


print("")
print("+++++++++++++++++++++++++++++++")
print("Elapsed time: %f minutes" % ((time() - start_time) / 60))
print("+++++++++++++++++++++++++++++++")
print("")
