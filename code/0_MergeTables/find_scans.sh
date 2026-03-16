#!/bin/bash

# ===========================
# ===== PARSE ARGUMENTS =====
# ===========================

usage() {
    echo "Usage: $0 [-o <str>] [-h]"
    echo "Options:"
    echo "  -o <str>: path to output directory; default is PWD"
    echo "  -h: display this help message"
    echo ""
    exit 1
}

# define defaults
odir=$(pwd)

# read arguments
while getopts o:h arg
do
	case $arg in
    o)  odir=${OPTARG};;
    h)  usage;;
	?)	echo ""
		echo "Unknown arguments passed; exiting."
		echo ""
		usage;
	esac
done

if [ ! -d "$odir" ]; then
    mkdir -p "$odir"
fi

# ================
# ===== ADNI =====
# ================

adni_dir="/ceph/chpc/shared/aristeidis_sotiras_group/aris_data/ADNI/bids/rawdata"

# amyloid PET
if [[ ! -f ${odir}/adni_amyloid.txt ]]; then
    echo "+++ finding ADNI amyloid PET scans +++"
    find \
        ${adni_dir} \
        -mindepth 4 -maxdepth 4 \
        -type f \
        \( -name "sub-*_ses-*_trc-AV45*.nii.gz" -o -name "sub-*_ses-*_trc-FBB*.nii.gz" \) \
        > ${odir}/adni_amyloid.txt
fi

# T1
if [[ ! -f ${odir}/adni_t1.txt ]]; then
    echo "+++ finding ADNI T1 MRI scans +++"
    find \
        ${adni_dir} \
        -mindepth 4 -maxdepth 4 \
        -type f \
        -name "sub-*_ses-*_acq-*_proc-gradbias_T1w.nii.gz" \
        > ${odir}/adni_t1.txt
fi

# =================
# ===== OASIS =====
# =================

oasis_dir="/ceph/chpc/rcif_datasets/oasis/OASIS3"

# amyloid PET
if [[ ! -f ${odir}/oasis_amyloid.txt ]]; then
    echo "+++ finding OASIS amyloid PET scans +++"
    # find \
    #     ${oasis_dir} \
    #     -mindepth 4 -maxdepth 4 \
    #     -type f \
    #     \( -name "sub-*_ses-*_acq-AV45*.nii.gz" -o -name "sub-*_ses-*_acq-PIB*.nii.gz" \) \
    #     > ${odir}/oasis_amyloid.txt
    find \
        ${oasis_dir} \
        -mindepth 5 -maxdepth 5 \
        -type f \
        \( -name "sub-*_ses*_acq-AV45*.nii.gz" -o -name "sub-*_ses*_acq-PIB*.nii.gz" \) \
        > ${odir}/oasis_amyloid.txt
fi

# T1
if [[ ! -f ${odir}/oasis_t1.txt ]]; then
    echo "+++ finding OASIS T1 MRI scans +++"
    # find \
    #     ${oasis_dir} \
    #     -mindepth 4 -maxdepth 4 \
    #     -type f \
    #     -name "sub-*_ses-*_*T1w*.nii.gz" \
    #     > ${odir}/oasis_t1.txt
    find \
        ${oasis_dir} \
        -mindepth 5 -maxdepth 5 \
        -type f \
        -name "sub-*_ses*_*T1w*.nii.gz" \
        > ${odir}/oasis_t1.txt
fi

# tau PET
oasis_tau_dir1="/ceph/chpc/rcif_datasets/oasis/OASIS3AV1451"
oasis_tau_dir2="/ceph/chpc/rcif_datasets/oasis/OASIS3AV1451L"
if [[ ! -f ${odir}/oasis_tau.txt ]]; then
    echo "+++ finding OASIS tau PET scans +++"
    find \
        ${oasis_tau_dir1} ${oasis_tau_dir2} \
        -mindepth 5 -maxdepth 5 \
        -type f \
        -name "sub-*_ses*_*AV1451*.nii.gz" \
        > ${odir}/oasis_tau.txt
fi

# ================
# ===== MCSA =====
# ================

mcsa_dir="/ceph/chpc/shared/aristeidis_sotiras_group/aris_data/MCSA/images/bids"

# amyloid PET
if [[ ! -f ${odir}/mcsa_amyloid.txt ]]; then
    echo "+++ finding MCSA amyloid PET scans +++"
    find \
        ${mcsa_dir} \
        -mindepth 4 -maxdepth 4 \
        -type f \
        \( -name "sub-*_ses-*_trc-PIB*_pet.nii.gz" \) \
        > ${odir}/mcsa_amyloid.txt
fi

# T1
if [[ ! -f ${odir}/mcsa_t1.txt ]]; then
    echo "+++ finding MCSA T1 MRI scans +++"
    find \
        ${mcsa_dir} \
        -mindepth 4 -maxdepth 4 \
        -type f \
        -name "sub-*_ses-*_T1w.nii.gz" \
        > ${odir}/mcsa_t1.txt
fi

# ===============
# ===== PAC =====
# ===============

pac_pet_dir="/ceph/chpc/shared/aristeidis_sotiras_group/aris_data/PAC/images/nifti/pet/unprocessed"
pac_mri_dir="/ceph/chpc/shared/aristeidis_sotiras_group/aris_data/PAC/images/nifti/mri/rawdata"

# amyloid PET
if [[ ! -f ${odir}/pac_amyloid.txt ]]; then
    echo "+++ finding PAC amyloid PET scans +++"
    find \
        ${pac_pet_dir} \
        -mindepth 2 -maxdepth 3 \
        -type f \
        \( -name "*.nii*" \) \
        > ${odir}/pac_amyloid.txt
fi

# T1
if [[ ! -f ${odir}/pac_t1.txt ]]; then
    echo "+++ finding PAC T1 MRI scans +++"
    find \
        ${pac_mri_dir} \
        -mindepth 3 -maxdepth 3 \
        -type f \
        \( -name "*_T1_LPS.nii.gz" \) \
        > ${odir}/pac_t1.txt
fi

# MUSE segmentation files
if [[ ! -f ${odir}/pac_muse.txt ]]; then
    echo "+++ finding PAC MUSE segmentations +++"
    find \
        ${pac_mri_dir} \
        -mindepth 3 -maxdepth 3 \
        -type f \
        \( -name "*_fastbc_muse.nii.gz" \) \
        > ${odir}/pac_muse.txt
fi

# DLICV masks
if [[ ! -f ${odir}/pac_dlicv.txt ]]; then
    echo "+++ finding PAC DLICV masks +++"
    find \
        ${pac_mri_dir} \
        -mindepth 3 -maxdepth 3 \
        -type f \
        \( -name "*_dlicvmask.nii.gz" \) \
        > ${odir}/pac_dlicv.txt
fi
