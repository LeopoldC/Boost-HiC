# Boost-HiC

Software requiered
=================
Boost-HiC current implementation is in python 2.7 and need the current list of package :
-h5py
-numpy
-copy
-sklearn
-scipy
-skimage

Input :
=================
Boost-HiC use HiC-Pro output format (described in : http://nservant.github.io/HiC-Pro/MANUAL.html#browsing-the-results ) for raw contact map.

The contact map is stored in tab separated file as :
bin_i / bin_j / counts_ij
Only no zero values are stored. Contact map are symmetric

The bin are described in a separated bed file which give the genomic coordinate of each bin.

In a first step, the contact map are convert in hdf5 by the pipeline.

Output :
=================
-The raw matrix without filtered bin by the procedure.
-The BoostHiC matrix, filtered two.
-pos out : the position of every filtered bin in myarray.hdf5

How to use it
=================

python main.py "Boost-HiC" myarray.hdf5 myarray_filtered.hdf5 myarray_boosted.hdf5 pos_outname.txt

