import numpy as np
from operator import mul    # or mul=lambda x,y:x*y
from fractions import Fraction
# import plotly.plotly as py
# import plotly.graph_objs as go
import matplotlib.pyplot as plt
from matplotlib import cm as CM
from matplotlib.pyplot import plot, draw, show
import scipy.interpolate
from textwrap import wrap
import matplotlib.pylab as pylab
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (15, 5),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)

def getW(T01, T10):
	wstep = 0.001
	for w in np.arange(wstep,10,wstep):
		last_pt = (w-wstep)**(T01 + T10) - 2*(w-wstep)**(T01 + T10 - 1) + (w-wstep)**(T01 + T10 - 2)
		this_pt = w**(T01 + T10) - 2*w**(T01 + T10 - 1) + w**(T01 + T10 - 2)
		if (last_pt <= 1 and this_pt >= 1):
			return w
	return 0

def getC(w):
	return np.log2(w)

## Paper: https://arxiv.org/pdf/1105.1969.pdf
# Setup channel params
print "Generating Transfer Function..."
D = .282 # Water = 0.282
R = 1
F = 1.3
xx = np.arange(0.1, 10, .1)
pi = np.pi
C = np.zeros(xx.size)
for x_id in range(0, xx.size):
	x = xx[x_id]
	step = 0.005
	t = np.arange(.0001,100,step)
	gxt = np.zeros(t.size)
	for i in range(0, t.size):
		gxt[i] = 1/(4*pi*D*t[i]) * np.exp(-x**2/(4*D*t[i])) 

	print "Transfer Function Generated!"

	print "Generating Transition Times..."
	# Setup sweep params
	directions = [01, 10]
	power = F
	start_time = 0
	xfer_time = np.zeros(2) + t[t.size-1]
	for d in range(0,2):
		transition = directions[d]
		S = 20

		# Configure step function
		T_change = start_time
		r = np.zeros(t.size)
		for ii in range(0, t.size):
			if (transition == 01):
				if (t[ii] < T_change):
					r[ii] = 0
				else:
					r[ii] = F
			else:
				if (t[ii] < T_change):
					r[ii] = F
				else:
					r[ii] = 0

		# Get channel response
		c = np.convolve(gxt, r)[0:t.size]

		if (transition == 01):	
			for j in range(1, c.size):
				if (c[j-1] <= 2*S and c[j] >= 2*S):
					xfer_time[d] = t[j] - T_change
					start_time = t[j] - T_change
					break
		else:
			for j in range(1, c.size):
				if (c[j-1] >= S and c[j] <= S):
					xfer_time[d] = t[j] - T_change
					break
		if (transition == 01):
			series = "T_01"
		else:
			series = "T_10"
		print "Generated series %s!" % (series)

	print "Computing Capacity..."
	# Find capacity for T01 T10
	W = getW(xfer_time[0], xfer_time[1])
	C[x_id] = getC(W)
	print "Capacity Computed!"



e = plt.figure(1)
plt.plot(xx, C)
plt.title("Capacity vs distance (Power = %f (units))\n" % F)
plt.xlabel("Distance (units)")
plt.ylabel("Capacity (nats)")


# Wait
show()