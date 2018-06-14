import h5py
import numpy as np
import sys

#my own toolkit
import HiCutils
import utils
import convert

### YOU ARE SUPPOSED TO ONLY MODIFY VALUE HERE ###
#input file
bedfilename='/users/invites/carron/Documents/Boost-HiC/test_dataset/rep1_10000_abs.bed'
matrixfilename='/users/invites/carron/Documents/Boost-HiC/test_dataset/chr16ES_10000.matrix'
Operation='Sample'
repositoryout='/run/media/carron/0ac0fffa-3350-431d-b1d1-865f8a21db21/data/Hi-C/Mouse/boyan/test/'

#default parameter
resolution=10000 #default : 10kb
achr="chr16"
alpha=0.2
species="mouse"
chrlist=utils.dictchr[species]
###


def BoostHiC(amat):
	normmat=HiCutils.SCN(np.copy(amat))
	FFmat=np.power(HiCutils.fastFloyd(1/np.power(normmat.copy(),alpha)),-1/alpha) #to dist, FF, to contact in one line
	boostedmat=HiCutils.adjustPdS(normmat,FFmat)
	return boostedmat

def Sample(amat,repositoryout):
	percentofsample=[0.1,1.,10.]
	for j in percentofsample:
		print("Value of sample",j)
		chrmat_s=np.copy(amat)
		chrmat=HiCutils.downsample_basic(chrmat_s,j)
		fh5 = h5py.File(repositoryout+"inputmat_sampleat_"+str(j)+".hdf5", "w")
		fh5['data'] = chrmat
		fh5.close()



### CODE EXECUTION ###

# load the data
print("LOADING MATRIX")
D=convert.loadabsdatafile(bedfilename)
beginfend=D[achr][0]
endfend=D[achr][1]
print("Data fend :",beginfend,endfend)
basemat=convert.loadmatrixselected(matrixfilename,beginfend,endfend)

#matrix filtering
print("FILTERING")
pos_out=HiCutils.get_outliers(basemat)
basematfilter=basemat[np.ix_(~pos_out, ~pos_out)]
basematfilter=np.copy(basematfilter)
#basematfilter=basematfilter[0:1000,0:1000]
print(len(basemat),len(basematfilter))
fh5 = h5py.File(repositoryout+"inputmat.hdf5", "w")
fh5['data'] = basemat
fh5.close()
fh5 = h5py.File(repositoryout+"inputmat_filtered.hdf5", "w")
fh5['data']=basematfilter
fh5.close()
utils.savematrixasfilelist3(pos_out,repositoryout+"filteredbin.txt")

if Operation=="Boost":
	print("Boost Hic")
	boosted=BoostHiC(basematfilter)
	#save
	fh5 = h5py.File(repositoryout+"boostedmat.hdf5", "w")
	fh5['data']=boosted
	fh5.close()
elif Operation=="Sample":
	print("SAMPLING")
	Sample(basematfilter,repositoryout)




