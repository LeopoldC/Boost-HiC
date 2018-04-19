# Boost-HiC

How to use it
=================
python main.py "Boost-HiC" myarray.hdf5 myarray_filtered.hdf5 myarray_boosted.hdf5 pos_outname.txt

Input :
=================
Actually, the code work for matrix saved in hdf5 format who contain a numpy array. Your contact do not have received any change (no normalisation).

Output :
=================
-The raw matrix without filtered bin by the procedure.
-The BoostHiC matrix, filtered two.
-pos out : the position of every filtered bin in myarray.hdf5
