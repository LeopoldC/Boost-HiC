#Leopold Carron
#december 2016
#CC-By-SA
#now in py 3


import sys
import os.path as op
import utils
import numpy as np
import scipy.io as scpio
from scipy import sparse
import h5py


def loadabsdatafile(filein):
	"""
	in a _abs.bedfile
	out : a dict with chr=[begin,end] in the bed, as finaly in the matrix
	"""
	fin=open(filein,"r")
	d={}
	d["Total"]=0
	l=fin.readline()
	while l:
		ls=l.split()
		if d.__contains__(ls[0]):
			d[ls[0]][1]=float(ls[3])
		else: #init
			d[ls[0]]=[float(ls[3]),float(ls[3])]
		d["Total"]+=1
		l=fin.readline()
	fin.close()
	return d

def loadmatrix(filein,sizemat):
	"""
	in : a matrix file, is size
	out : the matrix generated
	"""
	print(sizemat)
	mat=np.zeros((sizemat,sizemat))
	fin=open(filein,"r")
	l=fin.readline()
	while l:
		ls=l.split()
		i=int(float(ls[0])-1)
		j=int(float(ls[1])-1)
		v=float(ls[2])
		mat[i,j]=v
		mat[j,i]=v
		l=fin.readline()
	fin.close()
	print("Number of contact in map :",np.sum(np.triu(mat)))
	return mat

def loadmatrixselected(filein,B,E):
	"""
	-in : a matrix file, is size
	-out : the matrix generated
	-i-B>=0 j-E>=0
	"""
	#E-=1
	B-=1 #file start at one
	sizemat=int(E-B)
	print("Matrix size :",sizemat)
	mat=np.zeros((sizemat,sizemat))
	fin=open(filein,"r")
	l=fin.readline()
	while l:
		ls=l.split()
		i=int(float(ls[0])-B)
		j=int(float(ls[1])-B)
		if i>=0 and j>=0 and i<sizemat and j<sizemat:
			#print(ls[0],ls[1],i,j,B,E)
			v=float(ls[2])
			mat[i,j]=v
			mat[j,i]=v
		l=fin.readline()
	fin.close()
	print("Number of contact in the map :",np.sum(np.triu(mat)))
	return mat

