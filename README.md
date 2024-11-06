# Creating images of the distortion phantom model

For estimation of point spread function, distortion characterization etc

## Building

This uses ITK and tclap - the fetching and build is described in setup.sh

Clone this repo, change into the repo, run setup.sh

```
git clone https://github.com/richardbeare/DistortionPhantom

cd DistortionPhantom

bash ./setup.sh

```

This will produce an executable called `BuildLoadMesh/loadMesh`

This will convert an STL file to a 3D image. In order to do this
with the Calibre distortion phantom the 3D STEP model is required (supplied
by Calibre). I converted it to STL using FreeCad.

Conversion is then performed using a command like:

```
BuildLoadMesh/loadMesh -v 0.5 -b 5 -o output.nii.gz -i path/to/phantom.stl
```

The `v` and `b` flags are voxel buffer sizes in mm, with buffer referring to the 
spacing added to the mesh bounding box.

There will be a couple more steps to perform to get contrasts set up, which I'll
document once I've seen some real scans.

## Problems
The mesh rendering tool isn't perfect. Sometimes it leaves blank lines through the
solid or fills isolate lines in the empty space - I've no idea why.

These are minor, but annoying, because they may stop the fill process from
working. I've included corrections in my plug mask - this is manual.

Finding the problems is a pain - do a tophat and review.