#!/bin/bash


if [ ! -d BuildITK ] ; then 
	( 
git clone https://github.com/InsightSoftwareConsortium/ITK.git && \
   mkdir -p BuildITK && \
   cd BuildITK && \
   cmake -DModule_IOMeshSTL=ON -DModule_MeshToPolyData=ON ../ITK && \
   make -j4
)

fi

if [ ! -e Code/tclap ] ; then
(
cd Code && git clone https://github.com/mirror/tclap
)
fi

(
mkdir -p BuildLoadMesh && cd BuildLoadMesh && cmake -DITK_DIR=../BuildITK/ ../Code && make
)
