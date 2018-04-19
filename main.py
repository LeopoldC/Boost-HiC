import h5py
import numpy as np
import sys

import HiCutils
import utils

if sys.argv[1]=="Boost-HiC":
	basefile=sys.argv[2]
	basematfiltered=sys.argv[3]
	boostedmat=sys.argv[4]
	posoutname=sys.argv[5]
	chrlist=utils.dictchr[species]
	#load 
	fh5 = h5py.File(basefile, "r")
	chrmat=np.array(fh5['data'])
	fh5.close()
	#filt
	pos_out=HiCutils.get_outliers(chrmat)
	chrmat=chrmat[np.ix_(~pos_out, ~pos_out)]
	#boost
	newchrmat,basicmat,alpha=HiCutils.findalpha(chrmat,0.1,0.01)
	#adjust
	boosted=HiCutils.adjustPdS(np.copy(HiCutils.SCN(basicmat)),newchrmat)
	#save
	print("alpha trouve:",alpha)
	fh5 = h5py.File(basematfiltered, "w")
	fh5['data']=basicmat
	fh5.close()
	fh5 = h5py.File(boostedmat, "w")
	fh5['data']=boosted
	fh5.close()
	utils.savematrixasfilelist3(pos_out,posoutname)
elif sys.argv[1]=="Sample":
	baserep=sys.argv[2]
	species=sys.argv[3]
	repositoryout=sys.argv[4]
	resolution=float(sys.argv[5])
	chrlist=utils.dictchr[species]
	for i in chrlist:
		print("i: ",i)
		filenamein=baserep+i+".hdf5"
		fh5 = h5py.File(filenamein, "r")
		chrmat=np.array(fh5['data'])
		fh5.close()
		if resolution>10000:
			print('on va biner de 10000 a ',resolution)
			chrmat=HiCutils.binamatrixin2d(chrmat,10000,resolution)
		pos_out=HiCutils.get_outliers(chrmat)
		chrmat=chrmat[np.ix_(~pos_out, ~pos_out)]
		chrmatC=np.copy(chrmat)
		percentofsample=[0.1,1.,10.,100.]
		for j in percentofsample:
			print("valeur d'ech",j)
			chrmat_s=np.copy(chrmatC)
			chrmat=HiCutils.downsample_basic(chrmat_s,j)
			fh5 = h5py.File(repositoryout+i+"_basic_"+str(j)+".hdf5", "w")
			fh5['data'] = chrmat
			fh5.close()
		Vout=HiCutils.convertposoutseginmatlab(pos_out)
		utils.savematrixasfilelist3(Vout,repositoryout+i+"seg.txt")
		utils.savematrixasfilelist3(pos_out,repositoryout+i+"posout.txt")

