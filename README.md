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

How to use it
=================
The script is made to be most easy to use ad possible, just eddit the file main.py and run it as :
python main.py 

Important parameter are :
-bedfilename : bed file of genomic coordinate of each bin
-Operation : 'Boost' or 'Sample'
-repositoryout : where output file are stored

-resolution : default is 10000 basepair, change it if you need
-alpha : the alpha used to Boost-HiC

Some usefull tool are available in HiCtools.py if you need to made your own script.

Output :
=================
Every ouput contact map are stored in hdf5 with the key 'data' in the hdf5 dict. 

-inputmat.hdf5 : The contact map that is load in .matrix at the begin, just convert in hdf5
-inputmat_filtered.hdf5 : The original contact map with some filtered bin
-filteredbin.txt : list of raw/col that are filtered or not in the contact map as a boolean list.
-boostedmat.hdf5 : The input mat improved by BoostHiC procedure.


