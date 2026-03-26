#!/usr/bin/env python

import numpy as np
import argparse
import sys
from scipy.optimize import minimize
import scipy

def is_readable_file(parser, arg):
    try:
        f = open(arg, 'r')
        f.close()
    except Exception:
        raise argparse.ArgumentTypeError("{0} does not exist or is not readable".format(arg))

    return(arg)


parser = argparse.ArgumentParser(description="PSF default blur")

parser.add_argument("-p", "--phantom",
                    type=lambda x: is_readable_file(parser, x),
                    required=True,
                    help="Synthetic phantom nifti file")

parser.add_argument("-a", "--axial",
                    type=lambda x: is_readable_file(parser, x),
                    required=True,
                    help="axial scan (nifti file)")

#parser.add_argument("-c", "--coronal",
#                    type=lambda x: is_readable_file(parser, x),
#                    required=True,
#                    help="coronal scan (nifti file)")

# parser.add_argument("-s", "--sagittal",
#                     type=lambda x: is_readable_file(parser, x),
#                     required=True,
#                     help="coronal scan (nifti file)")


parser.add_argument("-o", "--output",
                    required=True,
                    help="Output prefix (for optimal estimated scans)")


args = parser.parse_args()

import SimpleITK as sitk

def getStandardPSF(Im):
    spacing = np.array(Im.GetSpacing())
    # different parameters for the through plane direction
    thruplane = np.argmax(spacing)
    cv = (1.2 * spacing)**2/(8*np.log(2))
    cv[thruplane] = spacing[thruplane]**2/(8*np.log(2))
    sig = np.sqrt(cv)
    return sig

def getModPSF(Im, Fast=False):
    # pass 1
    axO = np.array([0.55, 0.55, 4.31])
    axF = np.array([0.68, 0.68, 3.81])
    corO=np.array([0.5, 3.4, 0.5])
    corF=np.array([0.69, 3.11, 0.69])
    sagO=np.array([4.24, 0.5, 0.5])
    sagF=np.array([4.93, 0.7, 0.7])
    # pass 2
    axO = np.array([0.54, 0.54, 4.25])
    axF = np.array([0.7, 0.7, 3.43])
    corO=np.array([0.5, 2.86, 0.5])
    corF=np.array([0.69, 2.62, 0.69])
    sagO=np.array([3.99, 0.5, 0.5])
    sagF=np.array([5.19, 0.77, 0.77])

    if Fast:
        ax=axF
        cor=corF
        sag=sagF
    else:
        ax=axO
        cor=corO
        sag=sagO

    spacing = np.array(Im.GetSpacing())
    # different parameters for the through plane direction
    thruplane = np.argmax(spacing)
    if thruplane == 2:
        cv = ax
    elif thruplane == 1:
        cv = cor
    else:
        cv = sag
    sig = np.sqrt(cv)
    return sig

import os.path
def main():
    # Load everything
    ideal = sitk.ReadImage(args.phantom, sitk.sitkFloat32)
    axial = sitk.ReadImage(args.axial)
    # scale the inputs
    ideal = sitk.Cast(sitk.Normalize(ideal), sitk.sitkFloat32)

    # is it a fast sequence
    bn = os.path.basename(args.axial)
    if 'Fast' in bn:
        FST = True
        print("Fast sequence")
    else:
        FST = False
        print("Original sequence")
    # get the default blurring
    defaultPSF = getStandardPSF(axial)
    #defaultPSF = getModPSF(axial, FST)
    print(defaultPSF)
    idealsmoothed = sitk.SmoothingRecursiveGaussian(ideal, defaultPSF)
    sitk.WriteImage(idealsmoothed, args.output)

    return
    #coronal = sitk.ReadImage(args.coronal)
    #sagittal = sitk.ReadImage(args.sagittal)

    # Initalise the PSF - note that we need to estimate one per direction
    #axpsf = mkKernel(axial, ideal)
    # setup the nparray
    #npkern = setupNPKern(axpsf)
    #print(axpsf)
    #x = minimize(testKernel, x0 = npkern.ravel(), args=(axpsf, ideal, axialresampled, npkern.shape))
    #breakpoint()
main()
