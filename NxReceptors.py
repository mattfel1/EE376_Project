import numpy as np
from operator import mul    # or mul=lambda x,y:x*y
from fractions import Fraction
# import plotly.plotly as py
# import plotly.graph_objs as go
import matplotlib.pyplot as plt
from matplotlib import cm as CM
import scipy.interpolate
from textwrap import wrap
import matplotlib.pylab as pylab
params = {'legend.fontsize': 'x-large',
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)
## Paper: https://arxiv.org/pdf/1604.03508.pdf

withInteraction = 0# 0 = first plot in paper, 1 = second plot in paper
show45 = True

# Sweeping parameters
pstep = 1.0/(2**8)

n_rec = 2 # Nmuber of receptors
n = n_rec + 1 # Number of states
aH = 10 # Rate constant for U -> B in high concentration
aL = 1 # Rate constant for U -> B in low concentration
B = 20 # Rate constant for B -> U, independent of concentration
a = np.zeros(n_rec) # State transition probability in state k, not given X
p = np.zeros(n_rec) # Probability of X=H when in state k (apriori distr)
A = np.zeros(n) # Eigenvector components
A[0] = B**(n_rec)

def nCk(n,k): 
  return int( reduce(mul, (Fraction(n-i, i+1) for i in range(k)), 1) )

def phi(p):
	if (p == 0):
		return 0
	else:
		return -p*np.log(p)

def pi(k,A,Z):
	return A[k] / Z

def PHI(k):
	return (n_rec-k) * psi(p[k])

def PHI2(k):
	return psi(p[k])

def psi(p):
	return phi(p*aH+(1-p)*aL) - (p*phi(aH) + (1-p)*phi(aL))

def IXY(A,Z):
	I = 0
	for k in range(0,n_rec):
		if (withInteraction == 0):
			I = I + pi(k,A,Z) * PHI(k)
		else:
			I = I + pi(k,A,Z) * PHI2(k)			
	if p.sum() > 1:
		return 0
	else:
		return I
	# return I


dim = int(1/pstep)
p_ranges = np.zeros((n_rec,dim))
for k in range(0,n_rec):
	p_ranges[k,:] = np.arange(0,1,pstep)
P = np.meshgrid(*p_ranges)
numpoints = P[0].flatten().size
I = np.zeros(numpoints)
for dp in np.arange(0,numpoints):
	# Set up p, a, and A vector for this datapoint
	for k in range(0,n_rec):
		p[k] = P[k].flatten()[dp]
		a[k] = aL + (aH-aL)*p[k]
	for i in range(1,n):
		prod = 1
		for ii in range(0,i):
			prod = prod * a[ii]
		if (withInteraction == 0):
			A[i] = nCk(n_rec,i)*B**(n_rec-i)*prod
		else:
			A[i] = B**(n_rec-i)*prod

	Z = A.sum()

	# Copmute MI for this point
	I[dp] = IXY(A,Z)

maxp = [I.argmax()/dim, I.argmax()%dim]

if (n_rec == 2):	
	if (withInteraction == 0):
		summary = "Mutual Information W/O Interaction                                      Max MI = %f @ (%f, %f) for %d receptors ( ~%f per receptor )" % (I.max(), pstep*maxp[0], pstep*maxp[1], n_rec, I.max()/n_rec)
	else:
		summary = "Mutual Information With Interaction                                      Max MI = %f @ (%f, %f) for %d receptors ( ~%f per receptor )" % (I.max(), pstep*maxp[0], pstep*maxp[1], n_rec, I.max()/n_rec)

	print summary

	x,y = np.meshgrid(p_ranges[0], p_ranges[1])
	z=I.reshape(dim,dim)

	if show45:
		for i in range(0,dim):
			for j in range(0,dim):
				if (i == j):
					z[i,j] = 0
				dist_from_max = np.sqrt((i - maxp[0])**2 + (j - maxp[1])**2)
				radius = 7
				if (dist_from_max > radius - 1 and dist_from_max < radius + 1):
					z[i,j] = I.max()/4


	heatmap = plt.imshow(z, cmap='hot', interpolation='nearest', origin='lower', extent=[0, 1, 0, 1])
	ax = plt.gca
	cbar = plt.colorbar(heatmap)
	plt.title("\n".join(wrap(summary, 60)), fontsize=25)
	plt.xlabel('p0', fontsize=20)
	plt.ylabel('p1', fontsize=20)
	plt.show()

else:
	print "Max MI is %f for %d receptors ( ~%f per receptor )" % (I.max(), n_rec, I.max()/n_rec)

# # Greedy plotly way
# data = [
#     go.Heatmap(
#         z=I.reshape(dim,dim),
#         x=p_ranges[0],
#         y=p_ranges[1],
#         zsmooth="best",
#         zmin=0.5,
#         zmax=I.max()
#     )
# ]
# py.iplot(data, filename='basic-heatmap')


