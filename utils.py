#2017
#CC-By-SA
#Carron Leopold
#generic toolbox contain every usefull function that i need

import sys
import numpy as np
import numpy.random as npr
import os.path as op

###Common variable for every scripts
listchrhuman=["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM"]
listchrmouse=["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr16","chr15","chr17","chr18","chr19","chrX"]

#all species that can be generated
dictchr={"human":listchrhuman,
"mouse":listchrmouse}


###tools
#no common tools neded actually in this setup!

####LOAD
def loadchrsizedict(path,resolution):
	"""
	in : path to chromosomesize file from UCSC assembled genome
	out : a dict with size of each chr, containt both: TotalSize HicUsedTotalSize HicChrBegin
	"""
	f=open(path,"r")
	d={}
	l=f.readline()
	d["TotalSize"]=0
	d["HicUsedTotalSize"]=0
	while l:
		#print(l.strip('\n'))
		ls=l.split()
		d[ls[0]]=int(ls[1])
		d["TotalSize"]+=int(ls[1])
		F1=ls[0].__contains__("random")==False
		F2=ls[0].__contains__("Un")==False
		F3=str(ls[0])!="chrM"
		#F3=True
		F4=ls[0].__contains__("hap")==False
		#print(ls[0],ls[1],F1,F2,F3)
		if F1 & F2 & F3 & F4: #!!!! CARE !!!! All that part from the dict is generated for a specific resolution
			d[ls[0]]=np.float(ls[1])
			d["HicChrBegin"+ls[0]]=d["HicUsedTotalSize"]
			d["HicUsedTotalSize"]+=np.ceil(np.float(ls[1])/resolution)
		l=f.readline()
	f.close()
	return d

def LoadChrCentrom(filenamein):
	"""
	in : a centromeric file from ucsc
	out : a dict d[achr]:pos_of_centrom
	"""
	f=open(filenamein,'r')
	dout=dict()
	l=f.readline()
	while l:
		ls=l.split("\t")
		dout[ls[0]]=float(ls[1])
		l=f.readline()
	f.close()
	return dout

def loadfilelist(path):
	"""
	in : a file as a list of element
	out : return the list
	"""
	f=open(path,"r")
	thelist=list()
	l=f.readline()
	while l:
		thelist.append(l.strip("\n"))
		l=f.readline()
	f.close()
	return thelist

def loadfilefloatlist(path):
	"""
	in : a file as a list of element
	out : return the list
	"""
	f=open(path,"r")
	thelist=list()
	l=f.readline()
	while l:
		thelist.append(float(l.strip("\n")))
		l=f.readline()
	f.close()
	return thelist

def loadfileintlist(path):
	"""
	in : a file as a list of element
	out : return the list
	"""
	f=open(path,"r")
	thelist=list()
	l=f.readline()
	while l:
		thelist.append(int(float(l.strip("\n"))))
		l=f.readline()
	f.close()
	return thelist

def loadfiledict(path):
	"""
	in : a file to convert in dict
	out : a dict as d[ls[0]]=float(ls[1])
	"""
	f=open(path,"r")
	d=dict()
	l=f.readline()
	while l:
		ls=l.split()
		d[ls[0]]=float(ls[1])	
		l=f.readline()
	f.close()
	return d

def loadstrfiledict(path):
	"""
	in : a file to convert in dict
	out : a dict as d[ls[0]]=str(ls[1])
	"""
	f=open(path,"r")
	d=dict()
	l=f.readline()
	while l:
		ls=l.split()
		d[ls[0]]=ls[1]
		l=f.readline()
	f.close()
	return d


####SAVE

def savelistofset(setlist,nameout):
	"""
	in : a list of set
	out : set tabulate, each line a set
	"""
	out=open(nameout,"w")
	for i in setlist:
		astr=str()
		for j in i:
			astr+=j+"\t"
		astr=astr.strip("\t")+"\n" #faster to write
		out.write(astr)	
	out.close()

def savematrixasfilelist(matrix,fileout,loop):
	"""
	in : a numpy matrix, name how to vectorise it, savediagornot
	out : a vector file two have the distribution of the matrix contain
	"""
	fout=open(fileout,'w')
	size=matrix.shape
	i=0
	j=0
	if loop:
		while i<size[0]:
			while j<size[1]: #safest call
				if matrix[i,j]>0:
					fout.write(str(matrix[i,j])+"\n")			
				j+=1
			j=0
			i+=1
	else:
		while i<size[0]:
			while j<size[1]: #safest call
				if i!=j and matrix[i,j]>0:
					fout.write(str(np.float(matrix[i,j]))+"\n")			
				j+=1
			j=0
			i+=1
	fout.close()

def savematrixasfilelist2(matrix,fileout):
	"""
	in : a numpy matrix, name how to vectorise it, savediagornot
	out : a vector file two have the distribution of the matrix contain
	null value is destroy
	"""
	fout=open(fileout,'w')
	matrix=matrix[matrix.nonzero()]
	for i in matrix:
		fout.write(str(i)+"\n")			
	fout.close()

def savematrixasfilelist3(matrix,fileout):
	"""
	same as savematrixasfilelist2 without destroy
	"""
	fout=open(fileout,'w')
	for i in matrix:
		fout.write(str(i)+"\n")			
	fout.close()

def savenumberoflinkfrommataslist(matrix,fileout):
	"""
	in : a numpy matrix,name how to save it
	out : a vector for label number of node saving
	"""
	out=open(fileout,'w')
	i=0
	L=matrix.shape[0]
	while i<L:
		s=np.sum(matrix[i,:])
		out.write(str(s)+"\n")
		i+=1			
	out.close()

def savealist(alist,filename):
	"""
	in a list , name to save it
	save a list in a text file 1 element by line
	"""
	out=open(filename,'w')
	for i in alist:
		out.write(str(i)+"\n") #if i is numeric
	out.close()

def givematrixsize(sizedict,chrlist,achr,chrindex):
	"""
	out : end and begin of matrix pos
	"""	
	Li=len(chrlist)
	#print(achr,chrindex,Li)
	if chrindex==(Li-1): 
		print("fin de chromosome")
		return sizedict["HicUsedTotalSize"],sizedict["HicChrBegin"+achr]
	else:
		return sizedict["HicChrBegin"+chrlist[chrindex+1]],sizedict["HicChrBegin"+achr]


