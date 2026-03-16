# helper functions for PET SUVR pipeline

import json
import os
import shutil
import sys
from glob import glob

import numpy as np
import pandas as pd
from nipype.interfaces import fsl
from nipype.interfaces.base import BaseInterface, File, TraitedSpec, traits
from nipype.interfaces.fsl import BinaryMaths
from nipype.interfaces.io import DataSink
from nipype.interfaces.utility import Function, IdentityInterface
from nipype.pipeline.engine import JoinNode, Node, Workflow
from nipype.utils.filemanip import split_filename

sys.path.insert(0,"/home/b.y.yang/other_repos/pet-analysis")
from pet_analysis.ratio.suvr import create_pipeline
from pet_analysis.util.nipype_misc_funcs import ConcatenateSpreadsheets


def check_file_exists(df, columns_to_check):
    exists_df = df[columns_to_check].map(os.path.exists).add_suffix(".exists")
    exists_mask = exists_df.all(axis=1)
    return exists_mask, exists_df

def modify_json_add_FrameTimesStart(d):

    """
    add a "FrameTimesStart" field, which is the cummulative sum of
    FrameDuration list; always starts at 0

    use for the following datasets:
    - A4
        - note: for certain scans, there is no "FrameDuration" field; in
        these cases, assume 4 frames of 5-min (300s) each
    - MCSA
    """

    if not "FrameTimesStart" in d:

        # if FrameTimesStart doesn't exist, assume 4 frames of 5-min each
        if not "FrameDuration" in d:
            d["FrameDuration"] = [300, 300, 300, 300]
        
        # make FrameTimesStart
        d["FrameTimesStart"] = np.cumsum(
            np.concatenate([(0,), np.array(d["FrameDuration"])[:-1]])
        ).tolist()
    
    return d

def modify_json_change_FrameTimes_to_Time(d):

    """
    change "FrameTimes" to "Time"

    use for the following datasets:
    - OASIS
    """

    if "FrameTimes" in d:
        d["Time"] = d["FrameTimes"]
        del d["FrameTimes"]
    return d

def create_pet_json(
    n_frames: int,
    frame_duration_min: float,
    odir: str,
    prefix: str = "pet"
):
    
    frame_duration = np.repeat(frame_duration_min, n_frames)
    frame_times_start = np.cumsum(np.concatenate([(0,), frame_duration[:-1]]))
    
    d = {
        "FrameDuration": frame_duration.tolist(),
        "FrameTimesStart": frame_times_start.tolist(),
    }
    
    if not os.path.exists(odir): os.makedirs(odir, exist_ok = True)
    json_path = os.path.join(odir, f"{prefix}.json")
    with open(json_path, "w") as fp: json.dump(d, fp)
    return json_path

class ModifyJsonInputSpec(TraitedSpec):
    in_json = File(exists=True, desc="input JSON file", position=0, mandatory=True)
    func = traits.Callable(desc="function to modify JSON file", postion=1, mandatory=True)

class ModifyJsonOutputSpec(TraitedSpec):
    out_json = File(exists=True, desc="output JSON file")

class ModifyJson(BaseInterface):
    input_spec = ModifyJsonInputSpec
    output_spec = ModifyJsonOutputSpec

    def _run_interface(self, runtime):

        # if d is None or (isinstance(d, float) and np.isnan(d)):

        #     # if no JSON given, we'll assume that the input function
        #     # creates a new JSON
        #     d = self.inputs.func()

        # else:

        # load json
        with open(self.inputs.in_json, "r") as fpin:
            d = json.load(fpin)

        # if JSON is a list, take the first element
        if isinstance(d, (list, tuple)):
            d = d[0]
        
        # modify json using input function
        d = self.inputs.func(d)

        # save modified json
        with open(self._list_outputs()["out_json"], "w") as fpout:
            json.dump(d, fpout)
        return runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        fname = self.inputs.in_json
        _, base, ext = split_filename(fname)
        outputs["out_json"] = os.path.abspath(base + "_modified" + ext)
        return outputs

def create_muse_qc_wf(odir, name="muse_qc_wf"):

    def plot_roi_wrapper(roi, bg, subj_id, suffix, roi_idx = None):

        import os

        from nibabel import Nifti1Image, load
        from nilearn.plotting import plot_roi
        from numpy import isin, where

        if not roi_idx is None:
            roi = load(roi)
            roi_data = roi.get_fdata()
            mask = isin(roi_data, roi_idx)   # MUSE ROI indices for left and right cerebellum exterior
            masked_data = where(mask, roi_data, 0)
            roi = Nifti1Image(masked_data, affine = roi.affine, header = roi.header)
        
        output_file = os.path.abspath(subj_id + f"_{suffix}.png")
        plot_roi(
            roi_img = roi,
            bg_img = bg,
            display_mode = "mosaic",
            cut_coords = 6,
            threshold = 0.0001,  # this is to set background to black
            output_file = output_file
        )
        return output_file

    inputspec = Node(IdentityInterface(fields=['in_mri','in_label','mri_to_mni_mat','subj_id']), name="inputspec")
    template = fsl.Info.standard_image('MNI152_T1_1mm.nii.gz')
    transform_mri_node = Node(fsl.ApplyXFM(apply_xfm=True, reference = template), name='transform_mri')
    transform_label_node = Node(fsl.ApplyXFM(apply_xfm=True, reference = template, interp = "nearestneighbour"), name='transform_label')

    muse_roi_idx = [4,11,23,30,31,32,36,37,38,39,47,48,49,50,51,52,55,56,57,58,59,60,71,72,73,75,76,100,101,102,103,104,105,106,107,108,109,112,113,114,115,116,117,118,119,120,121,122,123,124,125,128,129,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207]
    plot_muse_node = Node(Function(input_names=["roi","bg","subj_id","suffix","roi_idx"], output_names=["out_png"], function=plot_roi_wrapper), name='plot_muse')
    plot_muse_node.inputs.suffix = "MUSE_QC"
    plot_muse_node.inputs.roi_idx = muse_roi_idx

    cerebellar_gm_roi_idx = [38, 39]
    plot_cerebellar_gm_node = Node(Function(input_names=["roi","bg","subj_id","suffix","roi_idx"], output_names=["out_png"], function=plot_roi_wrapper), name='plot_cerebellar_gm')
    plot_cerebellar_gm_node.inputs.suffix = "MUSE_QC_cerebellar_gm"
    plot_cerebellar_gm_node.inputs.roi_idx = cerebellar_gm_roi_idx
    
    datasink = Node(DataSink(), name = "datasink")
    datasink.inputs.base_directory = odir

    muse_qc_wf = Workflow(name = name, base_dir = odir)
    muse_qc_wf.connect([
        (inputspec, transform_mri_node, [('in_mri','in_file')]),
        (inputspec, transform_mri_node, [('mri_to_mni_mat','in_matrix_file')]),
        (inputspec, transform_label_node, [('in_label','in_file')]),
        (inputspec, transform_label_node, [('mri_to_mni_mat','in_matrix_file')]),
        
        (transform_mri_node, plot_muse_node, [('out_file','bg')]),
        (transform_label_node, plot_muse_node, [('out_file','roi')]),
        (inputspec, plot_muse_node, [('subj_id','subj_id')]),
        (plot_muse_node, datasink, [('out_png','qc_muse')]),

        (transform_mri_node, plot_cerebellar_gm_node, [('out_file','bg')]),
        (transform_label_node, plot_cerebellar_gm_node, [('out_file','roi')]),
        (inputspec, plot_cerebellar_gm_node, [('subj_id','subj_id')]),
        (plot_cerebellar_gm_node, datasink, [('out_png','qc_cerebellar_gm')]),
    ])

    return muse_qc_wf

def create_overlay_snapshot_wf(odir, name="overlay_wf"):

    def plot_roi_wrapper(roi, bg, subj_id):

        import os

        from nilearn.plotting import plot_roi

        output_file = os.path.abspath(subj_id + "_MUSE_QC.png")
        plot_roi(
            roi_img = roi,
            bg_img = bg,
            display_mode = "mosaic",
            cut_coords = 5,
            threshold = 0.0001,  # this is to set background to black
            output_file = output_file
        )
        return output_file
    
    inputspec = Node(IdentityInterface(fields=['background','mask','subj_id']), name="inputspec")
    plot_overlay_node = Node(Function(input_names=["roi","bg","subj_id"], output_names=["out_png"], function=plot_roi_wrapper), name='plot_overlay')
    datasink = Node(DataSink(), name = "datasink")
    datasink.inputs.base_directory = odir

    overlay_wf = Workflow(name = name, base_dir = odir)
    overlay_wf.connect([
        (inputspec, plot_overlay_node, [('background','bg'), ("mask", "roi"), ("subj_id", "subj_id")]),
        (plot_overlay_node, datasink, [('out_png',name)])
    ])

    return overlay_wf

def build_suvr_pipeline(
    data,
    output_dir,
    t_start,
    t_end,
    realign_t_start=None,
    realign_t_end=None,
    reg_smooth_fwhm = 4,
    pet_is_4D = True,
    label_schema = '/home/b.y.yang/other_repos/pac-pet/code/pipeline/MUSE_label_schema_for_PiB.json',
    ref_reg = 'cerebellar GM',
    modify_json_func = None,
    iterative_smoothing = False,
):
    
    """
    Wrapper for creating entire SUVR nipype pipeline
    """

    # pipeline params
    realign_params = {
        'quality': 1,
        'separation': 2,
        'fwhm': 7,
        'register_to_mean': False,
        'interp': 2,
        'write_interp': 4,
        'out_prefix': 'r',
        'dyn_mean_wts': 'frameduration',
        'unpadsize': 0,
        'padsize': 0,
    }
    # handle realign params
    if not realign_t_start is None:
        realign_params["target_t_start"] = realign_t_start
    if not realign_t_end is None:
        realign_params["target_t_end"] = realign_t_end

    registration_params = {
        'cost': 'mutualinfo',
        'dof': 6,
        'smooth_fwhm': reg_smooth_fwhm,
    }

    # handle input data
    # - if data is a pd.DataFrame, set up workflow as a multithreaded iterable
    # - if data is a pd.Series or dictionary, set up workflow as a single process
    if isinstance(data, pd.DataFrame):

        iterable_subj = True

        # get params from input data
        id_list = data['id'].values.tolist()
        musemask_list = data['musemaskpath'].values.tolist()
        musemriskull_list = data['musemriskullpath'].values.tolist()
        muselabel_list = data['muselabelpath'].values.tolist()
        pib_list = data['petpath'].values.tolist()
        pibtiming_list = data['petjsonfile'].values.tolist()
        manualcoregmat_list = data['manual_coreg_mat'].values.tolist()

        # placeholder Node to enable iteration over scans
        infosource = Node(IdentityInterface(fields=['idvi']),
                        iterables=('idvi', id_list),
                        name="infosource")
        
    elif isinstance(data, (pd.Series, dict)):

        iterable_subj = False

        # get params from input data
        id_list = [data["id"]]
        musemask_list = [data['musemaskpath']]
        musemriskull_list = [data['musemriskullpath']]
        muselabel_list = [data['muselabelpath']]
        pib_list = [data['petpath']]
        pibtiming_list = [data['petjsonfile']]
        manualcoregmat_list = [data['manual_coreg_mat']]

        # placeholder Node to enable iteration over scans
        infosource = Node(IdentityInterface(fields=['idvi']),
                        name="infosource")
        infosource.inputs.idvi = id_list[0]
    
    else:
        
        raise ValueError("input data must be a pd.DataFrame, pd.Series, or dictionary")
    
    # convert param list to dict
    pib_dict = dict(zip(id_list, pib_list))
    pibtiming_dict = dict(zip(id_list, pibtiming_list))
    musemask_dict = dict(zip(id_list, musemask_list))
    musemriskull_dict = dict(zip(id_list, musemriskull_list))
    muselabel_dict = dict(zip(id_list, muselabel_list))
    manualcoregmat_dict = dict(zip(id_list, manualcoregmat_list))

    def get_value(key, dict):
        return dict[key]
    getpib = Node(Function(input_names=['key','dict'], output_names=['value'],
                        function=get_value), name='getpib')
    getpib.inputs.dict = pib_dict

    # get full path to the txt file listing the duration of each PIB time frame
    #  - number of rows must be the same as the number of PIB time frames,
    # with each row listing the time in minutes
    getpibtiming = getpib.clone(name='getpibtiming')
    getpibtiming.inputs.dict = pibtiming_dict

    # get full path to brainmask corresponding to id from spreadsheet,
    # in same space as MUSE labels
    getmusemask = getpib.clone(name='getmusemask')
    getmusemask.inputs.dict = musemask_dict

    # get full path to MRI with skull corresponding to idvi from spreadsheet,
    # in same space as MUSE labels
    getmusemriskull = getpib.clone(name='getmusemriskull')
    getmusemriskull.inputs.dict = musemriskull_dict

    # get full path to MUSE label image corresponding to idvi from spreadsheet,
    # in same space as MRI
    getmuselabel = getpib.clone(name='getmuselabel')
    getmuselabel.inputs.dict = muselabel_dict

    getmanualcoregmat = getpib.clone(name='getmanualcoregmat')
    getmanualcoregmat.inputs.dict = manualcoregmat_dict

    # TODO: include option to perform PVC correction
    pvc_params = None

    # params for iterative smoothing of SUVR image
    if iterative_smoothing:
        smoothing_params = {
            "target_fwhm": (10,10,10),
            "stepsize": 0.5,
            "tolerance": 0.5
        }
    else:
        smoothing_params = None

    # node to use a binary brain mask to skull-strip MRI (doesn't actually create a brain mask)
    skullstrip = Node(BinaryMaths(operation='mul'), name='skullstrip')

    # get main PET SUVR workflow
    pet_pipeline = create_pipeline(label_schema=label_schema,
                                ref=ref_reg,
                                pet_is_4D=pet_is_4D,
                                use_skullstripped_mri_for_coreg=False,
                                realign_params=realign_params,
                                registration_params=registration_params,
                                pvc_params=pvc_params,
                                t_start_coreg=t_start, t_end_coreg=t_end,
                                t_start_suvr=t_start, t_end_suvr=t_end,
                                smoothing_params=smoothing_params,
                                n_procs=1, name='pet_workflow')

    if iterable_subj:
        suvr_csv = JoinNode(ConcatenateSpreadsheets(outputname='ROI_SUVR'),
                            joinsource='infosource', joinfield=['sheetlist'],
                            unique=True, name='suvr_csv')

    # BYY: create MUSE qc snapshots
    muse_qc_wf = create_muse_qc_wf(odir=output_dir)

    # create and connect workflow
    wf = Workflow(name="wf", base_dir=output_dir)

    if not modify_json_func is None:
        modify_json = Node(ModifyJson(func=modify_json_func), name="modify_json")
        wf.connect([
            (infosource, getpibtiming, [('idvi','key')]),
            (getpibtiming, modify_json, [("value","in_json")]),
            (modify_json, pet_pipeline, [("out_json","inputspec.pettiming")])
        ])
    else:
        wf.connect([
            (infosource, getpibtiming, [('idvi','key')]),
            (getpibtiming, pet_pipeline, [('value','inputspec.pettiming')]),
        ])

    wf.connect([
        (infosource, getpib, [('idvi','key')]),
        (infosource, getmusemask, [('idvi','key')]),
        (infosource, getmusemriskull, [('idvi','key')]),
        (infosource, getmuselabel, [('idvi','key')]),
        
        #(infosource, getdeform, [('idvi','key')]),

        (getmusemriskull, skullstrip, [('value','in_file')]),
        (getmusemask, skullstrip, [('value','operand_file')]),

        (infosource, pet_pipeline, [('idvi','inputspec.id')]),
        (getpib, pet_pipeline, [('value','inputspec.pet')]),
        (skullstrip, pet_pipeline, [('out_file','inputspec.mri')]),
        (getmusemriskull, pet_pipeline, [('value','inputspec.mriskull')]),
        (getmuselabel, pet_pipeline, [('value','inputspec.label')]),
        #(getdeform, pet_pipeline, [('value','inputspec.mri_to_template_transforms')]),

        (pet_pipeline, muse_qc_wf, [("preprocess_mriskull_wf.outputspec.img", "inputspec.in_mri")]),  # MUSE QC
        (pet_pipeline, muse_qc_wf, [("preprocess_label_wf.outputspec.img", "inputspec.in_label")]),
        (pet_pipeline, muse_qc_wf, [("affine_MNI_wf.mri_to_mni.out_matrix_file", "inputspec.mri_to_mni_mat")]),
        (infosource, muse_qc_wf, [("idvi", "inputspec.subj_id")]),
    ])

    if not data['manual_coreg_mat'] is None:
        wf.connect([
            (infosource, getmanualcoregmat, [('idvi','key')]),
            (getmanualcoregmat, pet_pipeline, [('value','coreg_wf.inputspec.in_matrix_file')]),
        ])

    if iterable_subj:
        wf.connect([(pet_pipeline, suvr_csv, [('outputspec.suvr_csv','sheetlist')]),])
    else:
        suvr_csv_datasink = Node(DataSink(), name="suvr_csv")
        wf.connect([(pet_pipeline, suvr_csv_datasink, [('outputspec.suvr_csv','suvr_csv')]),])

    return wf

def find_png(dir, idx = 7):
    png_files = glob(os.path.join(dir, "**/*.png"), recursive=True)

    temp_df = pd.DataFrame([item.split("/") for item in png_files])
    png_df = temp_df.iloc[:,idx].str.replace("_idvi_", "")
    return pd.concat([png_df, pd.Series(png_files)], axis = 1, ignore_index = True)

def move_png(png_df, odir):

    if not os.path.exists(odir):
        os.makedirs(odir, exist_ok = True)

    # Create group-specific directories and copy files
    for _, row in png_df.iterrows():

        group_id = row.iloc[0]
        file_path = row.iloc[1]

        # Create the group's folder if it doesn't exist
        group_folder = os.path.join(odir, group_id)
        os.makedirs(group_folder, exist_ok=True)

        # Copy the file to the group's folder
        if os.path.exists(file_path):
            shutil.copy(file_path, group_folder)
        else:
            print(f"File not found: {file_path}")
