# conda activate PySuperRes2

import SimpleITK as sitk

P = sitk.ReadImage("ImportantImages/final_reorient2.nii.gz", sitk.sitkUInt16);

# version without the plastic
PP = sitk.Mask(P, P==4000, 0, 0)
sitk.WriteImage(PP, "ImportantImages/final_reorient_copper.nii.gz")

# make edges in the bright side of a boundary

perode = sitk.GrayscaleErode(P, (1,1,1))
edges = ((P - perode) > 0) * (P > 2000)
sitk.WriteImage(edges, "ImportantImages/final_reorient2_edges.nii.gz")

# Edges aren't good for marking reliable warp information - just use the spheres
# Actually - this is wrong - lining up the warp field and the original shows that
# the warp information does come from the edges - the low res makes it appear to
# be the middle of the spheres if you don't have a reference.
# Does mean that we need to be careful about which edges we use though

# probably could have done this by finding holes in the copper mask....
manualmarkers = sitk.ReadImage("ImportantImages/manual_spheres.nii.gz", sitk.sitkUInt8)

plastic = P == 700
copper = P == 4000

copperfilled = sitk.GrayscaleMorphologicalClosing(copper, (3,3,3))
copperfilled = sitk.BinaryFillhole(copperfilled, False, 1)
copperedge = sitk.GrayscaleDilate(copperfilled, (1,1,1)) - sitk.GrayscaleErode(copperfilled, (1,1,1))
popen = sitk.GrayscaleMorphologicalOpening(plastic, (5,5,5))
# Now keep the blobs with a marker
popen2 = sitk.ReconstructionByDilation((manualmarkers*popen), popen)
edges2 = (sitk.GrayscaleDilate(popen2, (2,2,2)) - popen2) * copper
# place some of the edge inside the sphere
ee = sitk.GrayscaleDilate(edges2, (2,2,2)) * popen2
edges2 = sitk.Maximum(ee, edges2)
edges2 = sitk.Maximum(copperedge, edges2)

sitk.WriteImage(edges2, "ImportantImages/final_reorient2_edges.nii.gz")

