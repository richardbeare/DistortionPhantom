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


parser = argparse.ArgumentParser(description="PSF estimation")

parser.add_argument("-p", "--phantom",
                    type=lambda x: is_readable_file(parser, x),
                    required=True,
                    help="Synthetic phantom nifti file")

parser.add_argument("-a", "--axial",
                    type=lambda x: is_readable_file(parser, x),
                    required=True,
                    help="axial scan (nifti file)")

parser.add_argument("-r", "--axialresampled",
                    type=lambda x: is_readable_file(parser, x),
                    required=True,
                    help="axial scan aligned with phantom and resampled (nifti file)")

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
    return cv

# Lets assume the kernel is symmetric - therefore we have fewer parameters to estimate
# Kernel needs to be an image, which we pass to the ConvolutionImageFilter
# We can create it by filtering an image with a single bright spot
# Set up a default kernel
def mkKernel(Im, target):
    psf = getStandardPSF(Im)

    # create the kernel in the ideal image space
    SampleIm = sitk.Cast(target, sitk.sitkFloat32)
    SampleIm = 0*SampleIm
    sz = SampleIm.GetSize()
    middle = np.array(sz)/2
    print(psf)
    SampleIm.SetPixel(int(middle[0]), int(middle[1]), int(middle[2]), 100)
    smoothed = sitk.SmoothingRecursiveGaussian(SampleIm, [np.sqrt(psf[0]), np.sqrt(psf[1]), np.sqrt(psf[2])])
    xx = int(middle[0])
    yy = int(middle[1])
    zz = int(middle[2])

    # how many voxels in each direction
    fs = 10

    smoothed = smoothed[(xx-fs):(xx+fs+1), (yy-fs):(yy+fs+1), (zz-fs):(zz+fs+1)]
    return smoothed
    sitk.WriteImage(smoothed, "smoothed.nii.gz")
    print(middle)

def setupNPKern(imkern):
    sz = imkern.GetSize()
    middle= np.int_(np.floor(np.array(sz, dtype=np.int32)/2))
    print(middle)
    return sitk.GetArrayFromImage(imkern[0:(middle[0]+1), 0:(middle[1]+1), 0:(middle[2]+1)])

def npkernToImage(npkern, image):
    sz = npkern.shape
    npkern1 = npkern[::-1,:,:]
    npkern1 = npkern1[1:,:,:]
    b = np.concatenate((npkern, npkern1), axis=0)
    sz = b.shape
    npkern2 = b[:,::-1,:]
    npkern2 = npkern2[:,1:,:]
    b = np.concatenate((b, npkern2), axis=1)
    sz = b.shape
    npkern3 = b[:,:,::-1]
    npkern3 = npkern3[:,:,1:]
    b = np.concatenate((b, npkern3), axis=2)
    bi = sitk.GetImageFromArray(b)
    bi.CopyInformation(image)
    bi = sitk.Cast(bi, sitk.sitkFloat32)
    return bi


# will want to call this from scipy.optimize
def testKernel(optkern, imkernel, ideal, comparison, kernshape):
    # steps
    # 1. convert optkern to imkernel - imkernel is a template.
    #    optkern is a np.array containing 1/8 of the total
    # 2. Smooth the ideal image
    # 3. resample it to the comparison
    # 4. sum of squared differences (need to mess around with normalization etc)
    # 5. return
    # 
    # perhaps reregister?? Hope not 
    print("start cost estimate")
    oo = optkern
    oo = np.reshape(oo, kernshape)
    localkernel = npkernToImage(oo, imkernel)
    synthlowres = sitk.Convolution(ideal, localkernel)
    sitk.WriteImage(synthlowres, "sr.nii.gz")
    breakpoint()
    #synthlowres = sitk.Resample(synthlowres, referenceImage = comparison)
    difference = synthlowres - comparison
    difference = difference * difference
    score = sitk.GetArrayViewFromImage(difference)
    score = score.sum()
    print("Complete cost estimate " + str(score))
    return score

def estPSF_fourier(ideal, observed):
    idealF = scipy.fft.rfftn(sitk.GetArrayViewFromImage(ideal))
    observedF = scipy.fft.rfftn(sitk.GetArrayViewFromImage(observed))
    FF = observedF/idealF
    psf = scipy.fft.irfftn(FF)
    res = sitk.GetImageFromArray(np.real(psf))
    res.CopyInformation(ideal)
    res = sitk.FFTShift(res)
    # Apply the filter to the ideal
    convIdeal = scipy.fft.irfftn(idealF * FF)
    convIdeal = sitk.GetImageFromArray(np.real(convIdeal))
    convIdeal.CopyInformation(ideal)
    return { "PSF": res, "blurred" : convIdeal}

def main():
    # Load everything
    ideal = sitk.ReadImage(args.phantom, sitk.sitkFloat32)
    axial = sitk.ReadImage(args.axial)
    axialresampled = sitk.ReadImage(args.axialresampled)
    # scale the inputs
    ideal = sitk.Cast(sitk.Normalize(ideal), sitk.sitkFloat32)

    # test with the noise free data
    #idealsmoothed = sitk.SmoothingRecursiveGaussian(ideal, (3,3,3))
    #sitk.WriteImage(idealsmoothed, "sm.nii.gz")
    #PSF = estPSF_fourier(ideal, idealsmoothed)
    #sitk.WriteImage(PSF, args.output)
    #return
    axialresampled = sitk.Cast(sitk.Normalize(axialresampled), sitk.sitkFloat32)
    PSF = estPSF_fourier(ideal, axialresampled)
    sitk.WriteImage(PSF["PSF"], args.output + "/psf.nii.gz")
    sitk.WriteImage(PSF["blurred"], args.output + "/blurred.nii.gz")
    return
    #coronal = sitk.ReadImage(args.coronal)
    #sagittal = sitk.ReadImage(args.sagittal)

    # Initalise the PSF - note that we need to estimate one per direction
    axpsf = mkKernel(axial, ideal)
    # setup the nparray
    npkern = setupNPKern(axpsf)
    #print(axpsf)
    x = minimize(testKernel, x0 = npkern.ravel(), args=(axpsf, ideal, axialresampled, npkern.shape))
    breakpoint()
main()