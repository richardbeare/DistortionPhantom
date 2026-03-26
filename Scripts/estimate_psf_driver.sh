#!/bin/bash

export TD=$(mktemp -d)
trap 'rm -rf $TD; exit 0' 0 1 2 3 14 15
# Iterative registration and psf estimation
# First stage registration is linear, to the ideal image
# subsequent stages are nonlinear and to the blurred ideal
# image, with blurring based on the previous PSF estimation

#IDEALORIG=/group/deve2/ACTIVE/data/richard.beare/REVAMP/revamp_malawi/DistortionPhantom/Scripts/ImportantImages/final_reorient2.nii.gz
IDEALORIG=/group/deve2/ACTIVE/data/richard.beare/REVAMP/revamp_malawi/DistortionPhantom/Scripts/ImportantImages/final_reorient_copper.nii.gz
AXIAL=${1}
PSFD=${2}

RESULTFILE=${2}/warped_$(basename ${AXIAL})
[ -e ${RESULTFILE} ] && { echo "already run" ; exit; }

mkdir -p ${PSFD}

mkdir -p ${TD}/ANTS/

# Initialize
# flirt -ref ${IDEAL} -in ${AXIAL} -o ${TD}/ax_aligned.nii.gz -omat ${TD}/x.mat
# OBSERVED=${TD}/ax_aligned.nii.gz

# Externsion - blur the final with the default blurring parameters to make registration work b etter
spack unload -a
source ${HOME}/miniconda3/bin/activate psf_estimation

II=$(basename $IDEALORIG)
IDEAL=${TD}/${II}

python ./apply_default_psf.py -p ${IDEALORIG} -a ${AXIAL} -o ${IDEAL} 
conda deactivate

TARG=${TD}/IT_0
mkdir -p ${TARG}/ANTS/
spack load ants@2.5.1


antsRegistration --verbose 1 --dimensionality 3 --float 0 --collapse-output-transforms 1 \
    --output [ ${TARG}/ANTS/ants,${TARG}/ANTS/antsWarped.nii.gz,${TARG}/ANTS/antsInverseWarped.nii.gz ] --interpolation Linear \
    --use-histogram-matching 0 --winsorize-image-intensities [ 0.005,0.995 ] --initial-moving-transform [ ${IDEAL},${AXIAL},0 ] \
    --transform Rigid[ 0.1 ] --metric MI[ ${IDEAL},${AXIAL},1,32,Regular,0.5 ] --convergence [ 1000x500x250x0,1e-6,10 ] \
    --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0 \
    --transform Affine[ 0.1 ] --metric MI[ ${IDEAL},${AXIAL},1,32,Regular,0.5 ] --convergence [ 1000x500x250x0,1e-6,10 ] \
    --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0 \
    --transform SyN[ 0.1,1,0 ] --metric CC[ ${IDEAL},${AXIAL},1,1] \
    --convergence [ 100x70x50x20,1e-6,10 ] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0
    #--metric MI[ ${IDEAL},${AXIAL},1,32] \

# we warp the real scan to the blurred ideal one,
# but use the unblurred ideal in deconvolution
spack unload -a

OBSERVED=${TARG}/ANTS/antsWarped.nii.gz
W=${TARG}/ANTS/ants1Warp.nii.gz
WINV=${TARG}/ANTS/ants1InverseWarp.nii.gz
#TARG=${TD}/IT_${iteration}

source ${HOME}/miniconda3/bin/activate psf_estimation
#conda activate psf_estimation
mkdir -p ${TARG}/ANTS/
python ./estimate_psf_2.py -p ${IDEALORIG} -r ${OBSERVED} -a ${AXIAL} -o ${TARG}

psfout=psf_$(basename ${AXIAL})
mv ${TARG}/psf.nii.gz $2/${psfout}
mv ${OBSERVED} ${RESULTFILE}

BB=$(basename ${AXIAL})
BB=${BB/.nii.gz/}
mv ${W} ${2}/${BB}_Warp.nii.gz
mv ${WINV} ${2}/${BB}_InverseWarp.nii.gz

# iterative re-estimation seems to go nuts.
# Try a brute force optimization

# for iteration in {1..3} ; do
 
#     #/usr/local/ANTS/master/bin/antsRegistrationSyNQuick.sh -d 3 -f ${TARG}/blurred.nii.gz -m $AXIAL -o ${TD}/ANTS/ants
#     # This is a modified commandline from above, with the initialisation changed to geometric
#     /usr/local/ANTS/master/bin//antsRegistration --verbose 1 --dimensionality 3 --float 0 --collapse-output-transforms 1 \
#     --output [ ${TARG}/ANTS/ants,${TARG}/ANTS/antsWarped.nii.gz,${TARG}/ANTS/antsInverseWarped.nii.gz ] --interpolation Linear \
#     --use-histogram-matching 0 --winsorize-image-intensities [ 0.005,0.995 ] --initial-moving-transform [ ${TARG}/blurred.nii.gz,${AXIAL},0 ] \
#     --transform Rigid[ 0.1 ] --metric MI[ ${TARG}/blurred.nii.gz,${AXIAL},1,32,Regular,0.25 ] --convergence [ 1000x500x250x0,1e-6,10 ] \
#     --shrink-factors 12x8x4x2 --smoothing-sigmas 4x3x2x1vox --transform Affine[ 0.1 ] \
#     --metric MI[ ${TARG}/blurred.nii.gz,${AXIAL},1,32,Regular,0.25 ] --convergence [ 1000x500x250x0,1e-6,10 ] \
#     --shrink-factors 12x8x4x2 --smoothing-sigmas 4x3x2x1vox --transform SyN[ 0.1,3,0 ]\
#     --metric MI[ ${TARG}/blurred.nii.gz,${AXIAL},1,32] \
#     --convergence [ 100x100x70x50x0,1e-6,10 ] --shrink-factors 10x6x4x2x1 --smoothing-sigmas 5x3x2x1x0vox

#     OBSERVED=${TARG}/ANTS/antsWarped.nii.gz
# done
