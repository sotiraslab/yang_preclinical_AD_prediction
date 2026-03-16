#!/bin/bash

# ===========================
# ===== PARSE ARGUMENTS =====
# ===========================

usage() {
    echo "Usage: $0 [-i <str>] [-s <str>] [-h]"
    echo "Options:"
    echo "  -i <str>: path to input T1"
    echo "  -s <str>: scan ID"
	echo "  -o <str>: output directory"
	echo "  -p <str>: path to MUSE repository (optional)"
	echo "  -c <int>: number of cores for multithreaded process"
    echo "  -h: Display this help message"
	echo ""
    exit 1
}

# define defaults
proj_dir="/ceph/chpc/shared/aristeidis_sotiras_group/b.y.yang_scratch/UPenn-sMRI-pipeline/sMRI"
num_core="4"

# read arguments
while getopts i:s:o:p:c:h arg
do
	case $arg in
    i)  t1_path=${OPTARG};;
	s)	scan_id=${OPTARG};;
	o)  odir=${OPTARG};;
	p)  proj_dir=${OPTARG};;
	c)  num_core=${OPTARG};;
    h)  usage;;
	?)	echo ""
		echo "Unknown arguments passed; exiting."
		echo ""
		usage;
	esac
done

# checks
if [[ ! -f $t1_path ]]; then
	echo "Path to input T1 does not exist; exiting."
	exit 1
fi

# =================
# ===== SETUP =====
# =================

module load apptainer

temp_dir=$(mktemp -d /scratch/${USER}/tmp.XXXXXXX)
export TMPDIR=${temp_dir}

PID=$$
t1_tmp=${TMPDIR}/${scan_id}_T1.nii.gz

# # get params from param txt file
# param_path="/home/b.y.yang/sotiraslab/BradenADLongitudinalPrediction/data/subj/proc/muse.txt"
# param_arr=($(cat ${param_path} | awk -v line=${idx} '{if (NR == line) print $0}'))

# define variables
# PID=$$
# scan_id=${param_arr[0]}
# idvi=${param_arr[1]}
# t1_path=${param_arr[2]}
# t1_tmp=${TMPDIR}/${scan_id}_T1.nii.gz

# Paths
# The ${PROJ} path should be provided as absolute path.
PROJ=${proj_dir}
DEST=${odir}
SCR=${PROJ}/sMRI_ProcessingPipeline
confFile=${SCR}/configs/config_3D.sh
log_dir="${DEST}/logs"

# print

echo "+++ Running MUSE pipeline on ${scan_id} +++"
echo "+++ scan path: ${t1_path} +++"
echo ""
echo "+++ PROJ: ${PROJ} +++"
echo "+++ DEST: ${DEST} +++"
echo "+++ SCR: ${SCR} +++"
echo "+++ config: ${confFile} +++"
echo "+++ logs: ${log_dir} +++"

### Make dir
mkdir -p ${DEST}
mkdir -p ${log_dir}

if [ ! -d ${DEST}/Protocols/RAVENS/${scan_id} ]; then  # check for RAVENS output; if doesn't already exist, run MUSE pipeline

	### copy T1 to tmp dir
	cp ${t1_path} ${t1_tmp}

	### Run ProcessingPipeline_subject_3D
	echo -e "\nRunning commands on          : `hostname`"
	echo -e "Start time                   : `date +%F-%H:%M:%S` \n"

	### This ensure the $TMPDIR is always set
	if [ ! -d "$TMPDIR" ]
	then
		export TMPDIR=${PWD}/TMPDIR
	fi
	echo "TMPDIR=$TMPDIR"

	echo -e "\n"
	set -x

	apptainer \
		run \
		-B ${TMPDIR},${PROJ},${DEST} \
		${PROJ}/Container/cbica-muse-pipeline_1.0.0.sif \
		${PROJ}/sMRI_ProcessingPipeline/Scripts/ProcessingPipeline_subject_3D.sh \
		-ID ${scan_id} \
		-T1 ${t1_tmp} \
		-dest ${DEST} \
		-MT ${num_core} \
		-config ${confFile} \
		|& tee ${log_dir}/ProcessingPipeline_subject_3D.sh-${PID}.log

	set +x
	echo -e "\n"

else

	echo "RAVENS outputs exist; MUSE pipeline not run, exiting."

fi

module unload apptainer
rm -rf $temp_dir
