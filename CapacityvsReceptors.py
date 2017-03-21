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

numreceptors = np.arange(2,8,1)
C = np.zeros(numreceptors.size)
withInteraction = 1# 0 = first plot in paper, 1 = second plot in paper
show45 = True
for n_rec_id in range(0, numreceptors.size):

	# Sweeping parameters
	pstep = 1.0/(2**3)

	n_rec = numreceptors[n_rec_id] # Nmuber of receptors
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

	print "Max MI is %f for %d receptors ( ~%f per receptor )" % (I.max(), n_rec, I.max()/n_rec)
	C[n_rec_id] = I.max()

e = plt.figure(1)
plt.plot(numreceptors, C)
plt.title("Max Capacity vs Number of Receptors (Interaction = %d)\n" % withInteraction)
plt.xlabel("Number of receptors")
plt.ylabel("Capacity (nats)")

print C
# Without inter
# e = plt.figure(1)
# C = [3.573502, 5.348241 , 7.080807 , 8.740483 , 10.303278, 12.063548]
# plt.plot(numreceptors, C)
# plt.title("Max Capacity vs Number of Receptors (Interaction = %d)\n" % withInteraction)
# plt.xlabel("Number of receptors")
# plt.ylabel("Capacity (nats)")

# With inter
# e = plt.figure(1)
# C =  [2.095334, 2.354862, 2.454205, 2.484172, 2.514171, 2.524171]
# plt.plot(numreceptors, C)
# plt.title("Max Capacity vs Number of Receptors (Interaction = %d)\n" % withInteraction)
# plt.xlabel("Number of receptors")
# plt.ylabel("Capacity (nats)")

plt.show()
