import numpy as np
from operator import mul    # or mul=lambda x,y:x*y
from fractions import Fraction
# import plotly.plotly as py
# import plotly.graph_objs as go
import matplotlib.pyplot as plt
from matplotlib import cm as CM
import scipy.interpolate
from textwrap import wrap

## Paper: https://arxiv.org/pdf/1105.1969.pdf
# Setup channel parmas
D = .282 # Water = 0.282
R = 1
pi = np.pi
x = 1
step = 0.01
t = np.arange(.01,100,step)
gxt = np.zeros(t.size)
for i in range(0, t.size):
	gxt[i] = 1/(4*pi*D*t[i]) * np.exp(-x**2/(4*D*t[i])) 

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



# plt.plot(t[0:(int)(10/step)],gxt[0:(int)(10/step)])
# plt.title("Transfer Function of Diffusion Channel")
# plt.show()

plt.plot(powers, xfer_time[0,:])
plt.title("Transition time for 01 vs power")
plt.show()

plt.plot(powers, xfer_time[1,:])
plt.title("Transition time for 10 vs power")
plt.show()
