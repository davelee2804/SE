#!/usr/bin/env python

import math, sys, numpy, matplotlib.pyplot as plt

L = 1.0e+6
H = 5.0e+2
U = 1.0
g = 0.031
tau = 0.08
rho = 1.022e+3
nu = 3.0e+2
gamma = 5.0e-8
dt = 0.0092592592592592587

file = open( sys.argv[1], 'r' )

size = int(sys.argv[2])

time = numpy.zeros( size )
KE = numpy.zeros( size )
PE = numpy.zeros( size )
ws = numpy.zeros( size )
fr = numpy.zeros( size )

for i in range( 0, size, 1 ):
	line = file.readline()
	l2 = line.split()
	time[i] = float(l2[0])*dt
	KE[i] = float(l2[2])
	PE[i] = float(l2[3])
	ws[i] = float(l2[4])
	fr[i] = float(l2[5])

time *= (L/U)
time /= (60.0*60.0*24.0*365.25)
KE *= U*U*H*rho*L*2.0*L
PE *= U*U*H*rho*L*2.0*L
fr *= U*U*U*H*rho*2.0*L
ws *= U*U*U*H*rho*2.0*L

'''
fig = plt.figure()
ax = fig.add_subplot(111)
plt.plot( time, KE )
plt.plot( time, PE )
leg = ax.legend( ('KE', 'PE'), 'upper left' )
ax.set_xlabel( 'time (years)' )
ax.set_ylabel( 'energy (Joules)' )
plt.show()
'''
fig = plt.figure()
ax = fig.add_subplot(111)
plt.plot( time, ws )
plt.plot( time, fr )
leg = ax.legend( ('wind-stress', 'bottom-friction') )
ax.set_xlabel( 'time (days)' )
ax.set_ylabel( 'power (Joules/second)' )
plt.show()
