#python 3
#2017
#CC-By-SA
#Carron Leopold & Vincent Matthys
#specific utils for hic need

import sys
import numpy as np

from copy import deepcopy
from sklearn.neighbors import KernelDensity
from scipy.stats import *
from skimage.measure import *

import utils


#############################################################################
def windowssize2(n):
	"""
	size of windows in contact probability generator
	little bit tricky for complexity optimisation
	VdT : index of every bin in Vs
	Vs: size of every probability blow
	"""
	Vs=np.ones(n)
	#generate Vs
	i=1
	while i<n:
		Vs[i]=Vs[i-1]*1.01
		i+=1
	Vs=np.floor(Vs) #not the reel one at this step
	k=0
	i=0
	while i<n and k<n:
		Vs[i]=int(Vs[i])
		k+=int(Vs[i])
		i+=1
	#adjust
	Vs=Vs[0:i]
	Lj=len(Vs)
	Z=k-n
	Vs[Lj-1]=Vs[Lj-1]-Z
	Vs[Lj-2]=Vs[Lj-1]+Vs[Lj-2]
	Lj-=1
	i-=1
	Vs=Vs[0:i]
	print(Lj,n,i,Z,k,np.sum(Vs))
	#
	VdT=np.ones(n)
	i=0
	k=0
	while i<Lj:
		j=0
		while j<Vs[i]:
			VdT[k]=i
			j+=1
			k+=1
		i+=1	
	return VdT,Vs
	
def contactprobability(amat,Vs):
	"""
	return contact probability of an array
	"""
	L=len(Vs)
        matsize=np.shape(amat)
        probability=np.zeros(L)
        i=0
	K=0
        while i<L:
		j=0
		avec=np.array(list())
		while j<Vs[i]:
			avec=np.append(avec,np.diag(amat,k=K))
			j+=1
			K+=1
		probability[i]=np.mean(avec)
		i+=1	
        return probability

def adjustPdS(normmat,boostmat):
	"""
	Adjusting boostmap contact probability with normmat one
	"""
	matsize=np.shape(normmat)
	print(matsize[0])
	VdT,Vs=windowssize2(matsize[0])
	PC_base=contactprobability(normmat,Vs)
	PC_FF=contactprobability(boostmat,Vs)
	returnmat=np.copy(boostmat)
	#adjustement : a l'air ok!
	i=0
        j=0
        while i<matsize[0]: 
		adjust=PC_base[int(VdT[i])]/PC_FF[int(VdT[i])]
                while j<matsize[1]:
                        returnmat[j-i,j]=boostmat[j-i,j]*(adjust)
                        returnmat[j,j-i]=boostmat[j,j-i]*(adjust)
                        j+=1
                i+=1
                j=i
	return returnmat

######


def fastFloyd(contact):    
	n = contact.shape[0]    
	shortest = contact    
	for k in range(n):        
		i2k = np.tile(shortest[k,:], (n, 1))        
		k2j = np.tile(shortest[:, k], (n, 1)).T        
		shortest = np.minimum(shortest, i2k + k2j)    
	return shortest

def boost(normmat,alpha):
	"""
	Boost a mat
	"""
        matsize=np.shape(normmat)
        FFmat=np.power(fastFloyd(1/np.power(normmat.copy(),alpha)),-1/alpha)
	test=np.absolute(normmat-FFmat)
	nbo10=len(np.where(test>0.000000000000001)[0])
	returnmat=np.copy(FFmat)
        return returnmat,nbo10


### set of sub for the boost-hic algo

def makedmat():
	"""
	sub of find alpha
	"""
	d=dict()
	d['E']=list()
	d['A']=list()
	d['B']=list()
	d['C']=list()
	return d

def initdmat(d,nr,boostmat,actualalpha,step):
	"""
	sub of find alpha
	"""
	d['E']=[nr,boostmat,actualalpha-3*step]
	d['A']=[nr,boostmat,actualalpha-2*step]
	d['B']=[nr,boostmat,actualalpha-step]
	d['C']=[nr,boostmat,actualalpha]
	return d

def movestepfromd(d,nr,boostmat,actualalpha):
	"""
	sub of find alpha
	"""
	d['E']=d['A']
	d['A']=d['B']
	d['B']=d['C']
	d['C']=[nr,boostmat,actualalpha]
	return d

def continuefromd(d):
	"""
	sub of find alpha
	"""
	Dea=d['A'][0]-d['E'][0]
	Dab=d['B'][0]-d['A'][0]
	Dbc=d['C'][0]-d['B'][0]
	print("D1 : ",Dea,Dab,Dbc,Dbc-Dab,Dab-Dea)
	if (Dbc-Dab)>(Dab-Dea):
		return True
	return False

def printthedict(d):
	"""
	sub of find alpha
	"""
	print(d['E'][0],d['A'][0],d['B'][0],d['C'][0])

###

def findalpha(amat,froma,basestep):
        #init param
        alpha=froma
	step=basestep
	StepDict=makedmat()
	normmat=np.copy(amat)
	normmat=SCN(normmat)
	print("matrix size : ",normmat.shape)
	#first boost init
	newmat,nr=boost(normmat,alpha)
	StepDict=initdmat(StepDict,nr,newmat,alpha,step)
	alpha+=step
	nrnew=nr
	while nrnew==nr:
		newmat,nrnew=boost(normmat,alpha)
		StepDict=movestepfromd(StepDict,nrnew,newmat,alpha)
		alpha+=step
	return StepDict['B'][1],amat,StepDict['B'][2]


                       
#############################################################################

def binamatrixin2d(anumpyarray,resolutionfrom,resolutionto):
	"""
	in : A numpy array , number of bin in raw and in col
	out : the matrix binned
	"""
	convertionfactor=np.ceil(resolutionto/resolutionfrom)
	s=anumpyarray.shape
	print("dimension de la matrice:",s)
	#has to be identical as result in other function like chrsizedict)
	newsizei=np.ceil(s[0]*resolutionfrom/resolutionto)
	newsizej=np.ceil(s[1]*resolutionfrom/resolutionto)
	print(newsizei,newsizej)
	newarray=np.zeros((int(newsizei),int(newsizej)))
	print("taille de la matrice appres rescale :",newarray.shape)
	i=0
	j=0
	while i<newsizei:
		while j<newsizej:
			ifrom=int(i*convertionfactor)
			ito=int((i+1)*convertionfactor)
			jfrom=int(j*convertionfactor)
			jto=int((j+1)*convertionfactor)
			if i==newsizei-1:
				asum=np.sum(anumpyarray[ifrom:,jfrom:jto])
			elif j==newsizej-1:
				asum=np.sum(anumpyarray[ifrom:ito,jfrom:])
			elif i==newsizei-1 and j==newsizej-1:
				asum=np.sum(anumpyarray[ifrom:,jfrom:])
			else:
				asum=np.sum(anumpyarray[ifrom:ito,jfrom:jto])
			newarray[i,j]=asum
			#newarray[j,i]=asum
			j+=1
		i+=1
		j=0
	return newarray

#############################################################################
#Filter

def kde_sklearn(reads, bandwidth = 1000, **kwargs):
	"""
	Kernel Density Estimation with Scikit-learn :
	Estimate the density from the reads distribution
	Entry :
		- reads, as the number of reads per bin
		- bandwith : width of the gaussien parameter (can be modified)
	Output :
		- returns the density as a matrix for each point from 0 to max(reads)
	"""
	x_grid = np.linspace(0, np.max(reads), int(np.max(reads)) + 1)
	kde_skl = KernelDensity(kernel = 'gaussian', bandwidth=bandwidth, **kwargs)
	kde_skl.fit(reads)
	# score_samples() returns the log-likelihood of the samples
	log_pdf = kde_skl.score_samples(x_grid[:, np.newaxis])
	res = np.exp(log_pdf) / np.exp(log_pdf).sum()
	return res

def kde_outliers(reads, threshold = [0.05, 0.95], w = 1000):
	"""
	outliers as given by the kernel density estimator
	Entry :
		- reads : the number of reads per bin
		- threshold : list of [threshold min, threshold max].
	Below threshold min and above threshold max,
	bins are considered as outliers
	If threshold = float between 0 and 1, consider the threshold x and 1-x
	- w : bandwith of the gaussian window (as needed in kde_sklearn)
	Output :
		- positions of outliers : array of bool where True = outlier, False = to be kept
		- limits of outliers as a tuple (limit_below, limite_above)
		- density as estimated by kde_sklearn function
	"""
	# Before everything, 0 bins have to be cleaned to evaluate density properly
	# We will then work on reads[reads > 0]
	cleaned = reads[reads > 0]

	# We need first the density estimation
	density = kde_sklearn(cleaned[:, None], bandwidth = w)

	# Then we calculate the cumulative sum
	cum_sum = np.cumsum(density)

	# We can know find the limits
	if type(threshold) == float and threshold < 1 and threshold > 0 :
		threshold = [threshold, 1 - threshold]
	limit_below = np.abs(cum_sum - threshold[0]).argmin()
	limit_above = np.abs(cum_sum - threshold[1]).argmin()

	# Now use reads again to add 0 bins to the outliers
	return ((reads < limit_below) + (reads > limit_above)),(limit_below, limit_above), density


def get_outliers(mat, threshold = [0.05, 0.995]):
	"""
	Given a mat return outliers bin
	"""
	# Position of outliers : limitis and density
	#print(mat.sum(axis = 0))
	pos_out, lims, density = kde_outliers(mat.sum(axis = 0), threshold = threshold, w = 2000)
	return pos_out


def clean_out(data, pos_out):
	"""
	Given the positions of outliers, clean the data. This mean, where
	put to 0 all reads where pos_out == True
	Entry :
		- data : raw data of contact map
		- pos_out : positions of outliers bins, numpy.ndarray of one dimension
	containing booleans indicating for each position True if this bin is
	a outlier, and False if this bin is not.
	"""
	data[pos_out] = 0
	data[:, pos_out] = 0
	return data



#############################################################################
def divide_diag(data, mean = True):
    """
    Divide by the sum of diag if mean = False
    Divide by the mean of diag if mean = True
    """
    # Dimension of data
    dim = data.shape[0]

    for k in range(-dim + 1, dim):
        diagonal = np.diag(data, k)
        if mean == False :
            div = float(np.maximum(1, diagonal.sum()))
        else :
            # For a given diagonal, take the mean of the non-null element
            ## To estimate if the number of contact is above or below the
            ## mean of contact for this genomic distance
            div = np.nan_to_num(float(diagonal[diagonal != 0].mean()))
            if div == 0:
                div = 1
        data += np.diagflat(-diagonal + diagonal / div, k)
    return data

###############################################################################

def downsample_basic(contact, k):
	"""
	-make binomial sampling on each contact
	-k has to be between 1 and 100
	"""
	B = np.zeros(contact.shape)
	S = contact.sum()
	i, j = 0, 0
	L = contact.shape[0]
	while i < L:
		while j < L:
			if contact[i, j] != 0:
				B[i, j] = npr.binomial(contact[i, j], k*1.0/100.)
			j += 1
		i += 1
		j = i
	B = B + B.T
	S = B.sum()
	return B

def downsample_dicho(mat, k):
	matA=downsample_basic(mat, k)
	matB=mat-matA
	print("test unit:",matA.shape,matB.shape,np.sum(np.sum(matA)),np.sum(np.sum(matB)))
	return matA,matB


def SCN(D, max_iter = 10, mean = True):    
	# Iteration over max_iter    
	for i in range(max_iter):        
		D /= np.maximum(1, D.sum(axis = 0))       
		D /= np.maximum(1, D.sum(axis = 1)[:, None])    
		# To make matrix symetric again   
	return (D + D.T)/2 

def SCNandOE(D, max_iter = 10, after = 1, mean = True):
    """
    SCN method to normalize contact map.
    Entry :
        - Matrix  of contact map
        - max_iteration of algorithm
        - after : 1 if diagonal is divided after iterations, 0 to divide for each iteration
        - mean : if set to true, each diagonal is divided by its mean. If set to false, each diagonal is divided by its sum
    """
    # Cpy of D
    data = D
    # Iteration over max_iter
    for i in range(max_iter):
        data /= np.maximum(1, data.sum(axis = 0))
        data /= np.maximum(1, data.sum(axis = 1)[:, None])
        if after == 0:
            divide_diag(data, mean = mean)
    if after == 1:
        divide_diag(data, mean = mean)
    # To make matrix symetric again
    return (data + data.T)/2


def observed_expected(Hicmat):
	"""
	Run observed_expected on a Hicmat, return it
	"""
	OE=deepcopy(Hicmat)
	i=0
	j=0
	L=len(Hicmat)
	while j<L:
		thediag=np.diag(Hicmat,k=j)
		mtg=np.mean(thediag)
		while i<(L-j):
			v=Hicmat[i,i+j]/mtg
			OE[i,i+j]=v
			OE[i+j,i]=v
			i+=1
		i=0
		j+=1
	return OE
			

