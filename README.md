<p align="center">
<img width="250" height="154" src="Common/doc/logoSU2small.png">
</p>


# SU2-Efficient-Output-Control 
SU2 version 7.3.1 is modified to efficiently perform large unsteady simulations/optimizations compatible with OpenCFS H5 output format.

SU2-HDF5-Output uses the C++ interface HighFive (https://github.com/BlueBrain/HighFive) to write the mesh and the solutions of all timesteps into a single HDF5 file in parallel. The output file is fully compatible with OpenCFS (https://opencfs.gitlab.io/userdocu/Installation/ParaView/) and can be opened with the OpenCFS Paraview Plugin. Furthermore, for aeroacoustic applications, the source term computed in SU2 can be directly written to the OpenCFS native input file circumventing any unnecessary format conversions.


Installation
------------
The code has been compiled and tested with gcc version 11 and openmpi version 3.1.2.
You can use meson for configuration (e.g. `./meson.py -Denable-autodiff=true -Denable-directdiff=true build --prefix=YourPath`). Then the compilation is done using ninja (e.g. `./ninja -v -C build/ install`). For more information read https://su2code.github.io/docs_v7/Build-SU2-Linux-MacOS/#configuration-and-compilation


Usage
------------
To use the OpenCFS H5 output format, set `PARAVIEW_CFS_H5` as the output type in the config file e.g. `OUTPUT_FILES=(RESTART,PARAVIEW_CFS_H5)`.
To reduce the number of output files for unsteady optimizations, run the `shape_optimization.py` script with the option `--minimal True`. This option should be used only with `PARAVIEW_CFS_H5` output type. It removes all unnecessary files preceding the latest two design steps except the H5 solutions.

Precompiled Paraview binaries (V5.9) with the CFS reader plugin is placed inside 'Paraview_CFS'. To build the plugin for other versions of Paraview see: https://opencfs.gitlab.io/userdocu/Installation/ParaView/


Acknowledgment
--------------
We gratefully acknowledge the support from 'KONWIHR'.

KONWIHR [https://www.konwihr.de/] (Competence Network for Scientific High Performance Computing in Bavaria).
