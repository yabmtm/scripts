#!/bin/python

# script used to plot integrin/cilengitide free energy data
# for acs 2016. - Matt Hurley

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# 1.75 data
r175, = plt.plot(0.4698, 3.752, 'ro', label='(NMeR)GDfV : 1.75') # R
g175, = plt.plot(1.722, -0.662, 'bo', label='R(NMeG)DfV : 1.75') # G
d175, = plt.plot(3.225, 4.357, 'go', label='RG(NMeD)fV : 1.75') # D
f175, = plt.plot(3.771, 7.817, 'yo', label='RGD(NMeF)V : 1.75') # F
v175, = plt.plot(-0.871, -0.698, 'ko', label='RGDf(NMeV) : 1.75') # V

# 2.00 data
r200, = plt.plot(0.4698, 3.865, 'r^', label='(NMeR)GDfV : 2.00') # R
g200, = plt.plot(1.722, 1.474, 'b^', label='R(NMeG)DfV : 2.00') # G
d200, = plt.plot(3.225, 3.548, 'g^', label='RG(NMeD)fV : 2.00') # D
f200, = plt.plot(3.771, 0.598, 'y^', label='RGD(NMef)V : 2.00') # F
v200, = plt.plot(-0.871, -0.750, 'k^', label='RGDf(NMeV) : 2.00') # V

# old data generated by Amanda
plt.plot(0.4698, 1.001, 'rs') # R
plt.plot(1.722, -0.951, 'bs') # G
# plt.plot(3.225, -10.589, 'gs') # D (again exempt)
plt.plot(3.771, 0.562, 'ys') # F
plt.plot(-0.871, 0.700, 'ks') # V

# other D data exempt from data to focus on values that hit closer to experiment
# plt.plot(3.225, 12.309, 'gD') # 1.50
# plt.plot(3.225, 12.230, 'g*') # 1.66
# plt.plot(3.225, 16.732, 'g+') # 1.71

# y = x line for reference
plt.plot([-1, 0, 1, 2, 3, 4, 5, 6, 7], [-1, 0, 1, 2, 3, 4, 5, 6, 7] )

# axis labels
plt.title('Comparison of simulation and experimental data.')
plt.xlabel('Experimental Data (kcal/mole)')
plt.ylabel('Simulation Data (kcal/mole)')

# legend
# plt.legend(numpoints=1)

# save plot as a png file
plt.savefig('1.png')


