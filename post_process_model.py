#!/usr/bin/env python3

import numpy as np
import argparse
import sys

def is_readable_file(parser, arg):
    try:
        f = open(arg, 'r')
        f.close()
    except Exception:
        raise argparse.ArgumentTypeError("{0} does not exist or is not readable".format(arg))

    return(arg)



parser = argparse.ArgumentParser(description="Phantom model hacking")

parser.add_argument("-p", "--phantom",
                    type=lambda x: is_readable_file(parser, x),
                    required=True,
                    help="Synthetic phantom nifti file")

parser.add_argument("-m", "--mask",
                    type=lambda x: is_readable_file(parser, x),
                    required=True,
                    help="mask blocking the phantom outlets so a fill will work")

parser.add_argument("-s", "--coppersulphate", 
                    required=False,
                    type=float, 
                    default=100,
                    help="Brightness for solution")

parser.add_argument("-b", "--body", 
                    required=False,
                    type=float, 
                    default=10,
                    help="Brightness for plastic body")

parser.add_argument("-o", "--output",
                    required=True,
                    help="Output nifti file")


args = parser.parse_args()


import SimpleITK as sitk

def main():
    phantom = sitk.ReadImage(args.phantom)
    plugs = sitk.ReadImage(args.mask)

    pluggedA = sitk.Maximum(phantom, plugs*255) > 0  
    # pad boundary just in case
    sitk.WriteImage(pluggedA, "pluggedorig.nii.gz")
    
    #plugged = sitk.BinaryMorphologicalClosing(pluggedA, (3,3,3))
    #sitk.WriteImage(plugged - pluggedA, "diff.nii.gz")
    #plugged = sitk.ConstantPad(plugged, (1,1,1), (1,1,1))
    filled = sitk.BinaryFillhole(pluggedA, foregroundValue = 1)

    cuso4 = sitk.Cast(filled - pluggedA, sitk.sitkFloat32) * args.coppersulphate
    plastic = sitk.Cast(phantom > 0, sitk.sitkFloat32) * args.body

    result = cuso4 + plastic
    sitk.WriteImage(result, args.output)

main()