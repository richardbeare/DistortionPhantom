# conda activate PySuperRes2

import SimpleITK as sitk
import numpy as np
P = sitk.ReadImage("/group/deve2/ACTIVE/data/richard.beare/REVAMP/GHOST/data/T2_phantom.nii.gz", sitk.sitkUInt16);
ignore = sitk.ReadImage("/group/deve2/ACTIVE/data/richard.beare/REVAMP/GHOST/data/T2_checkermask.nii.gz", sitk.sitkUInt16);

# make edges in the bright side of a boundary

perode = sitk.GrayscaleErode(P, (1,1,1))
pdilate = sitk.GrayscaleDilate(P, (1,1,1))

edges = ((pdilate - perode) > 100) * (ignore == 0)
sitk.WriteImage(edges, "/group/deve2/ACTIVE/data/richard.beare/REVAMP/GHOST/data/T2_edges.nii.gz")

# to determine phantom volume
k = sitk.GrayscaleMorphologicalOpening(P > 5, [3,3,3], sitk.sitkBox)

k1 = sitk.GrayscaleDilate(k, [3,3,3], sitk.sitkBox)
k2 = sitk.BinaryFillhole(k1)
k3 = sitk.GrayscaleErode(k2, [3,3,3], sitk.sitkBox)

ls = sitk.LabelStatisticsImageFilter()
ls.Execute(k3,k3)
ls.GetCount(1)*np.array(k3.GetSpacing()).prod()
#np.float64(1987495.4363859568) = 1.987 liters
