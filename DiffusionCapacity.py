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

def getW(T01, T10):
	wstep = 0.01
	for w in np.arange(wstep,1.5,wstep):
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
pi = np.pi
x = 1
step = 0.01
t = np.arange(.01,100,step)
gxt = np.zeros(t.size)
for i in range(0, t.size):
	gxt[i] = 1/(4*pi*D*t[i]) * np.exp(-x**2/(4*D*t[i])) 

print "Transfer Function Generated!"

print "Generating Transition Times..."
# Setup sweep params
directions = [01, 10]
powers = np.arange(.4,8,.1)
start_time = np.zeros(powers.size)
xfer_time = np.zeros((2, powers.size)) + t[t.size-1]
for d in range(0,2):
	transition = directions[d]
	S = 20

	for i in range(0,powers.size):
		T_change = start_time[i]
		r = np.zeros(t.size)
		F = powers[i]
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

		c = np.convolve(gxt, r)[0:t.size]
		if (transition == 01):	
			for j in range(1, c.size):
				if (c[j-1] <= 2*S and c[j] >= 2*S):
					xfer_time[d,i] = t[j] - T_change
					start_time[i] = t[j] - T_change
					break
		else:
			for j in range(1, c.size):
				if (c[j-1] >= S and c[j] <= S):
					xfer_time[d,i] = t[j] - T_change
					break
	if (transition == 01):
		series = "T_01"
	else:
		series = "T_10"
	print "Generated series %s!" % (series)

print "Computing Capacities..."
# Find capacity vs T01 T10
C = np.zeros(powers.size)
for i in range(0, powers.size):
	W = getW(xfer_time[0,i], xfer_time[1,i])
	C[i] = getC(W)
print "Capacites Computed!"

m = C.argmax()
print "Maximum Capacity = %f (F = %f, T_01 = %f, T_10 = %f" % (C.max(), powers[m], xfer_time[0,m], xfer_time[1,m])

f = plt.figure(1)
plt.plot(t[0:(int)(10/step)],gxt[0:(int)(10/step)])
plt.title("Transfer Function of Diffusion Channel")
plt.xlabel("time (s)")
plt.ylabel("Diffusion Concentration")
# f.show()

g = plt.figure(2)
plt.plot(powers, xfer_time[0,:])
plt.title("Transition time for 01 vs power")
plt.xlabel("F (units)")
plt.ylabel("T_01 (s)")
# g.show()

h = plt.figure(3)
plt.plot(powers, xfer_time[1,:])
plt.title("Transition time for 10 vs power")
plt.xlabel("F (units)")
plt.ylabel("T_01 (s)")
# h.show()

e = plt.figure(4)
plt.plot(powers, C)
plt.title("Capacity vs power")
plt.xlabel("F (units)")
plt.ylabel("Capacity (nats)")

# Wait
show()