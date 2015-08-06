#!/usr/bin/env python

import math, sys, numpy, matplotlib.pyplot as plt

L = 1.0e+6
H = 1.0e+3
U = 1.0
g = 0.02
tau = 0.1
rho = 1.0e+3
nu = 3.3e+2
gamma = 0.0
dt = 0.0037037037037037038

file = open( sys.argv[1], 'r' )

size = int(sys.argv[2])

time = numpy.zeros( size )
KE1 = numpy.zeros( size )
KE2 = numpy.zeros( size )
PE = numpy.zeros( size )
ps1 = numpy.zeros( size )
ps2 = numpy.zeros( size )
ws = numpy.zeros( size )

for i in range( 0, size, 1 ):
	line = file.readline()
	l2 = line.split()
	time[i] = float(l2[0])*dt
	KE1[i] = float(l2[2])
	KE2[i] = float(l2[3])
	PE[i] = float(l2[4])
	ps1[i] = float(l2[5])
	ps2[i] = float(l2[6])
	ws[i] = float(l2[7])

time *= (L/U)
time /= (60.0*60.0*24.0)
KE1 *= U*U*H*rho*L*L
KE2 *= U*U*H*rho*L*L
PE *= U*U*H*rho*L*L
ps1 *= rho*H*U*U*U*L
ps2 *= rho*H*U*U*U*L
ws *= rho*H*U*U*U*L

fig = plt.figure()
ax = fig.add_subplot(111)
plt.plot( time, KE1 )
plt.plot( time, KE2 )
plt.plot( time, PE )
leg = ax.legend( ('KE-1', 'KE-2', 'PE'), 'upper left' )
ax.set_xlabel( 'time (days)' )
ax.set_ylabel( 'energy/unit area (J/m^2)' )
plt.show()
'''
'''
fig = plt.figure()
ax = fig.add_subplot(111)
plt.plot( time, ps1 )
plt.plot( time, ps2 )
plt.plot( time, ws )
leg = ax.legend( ('pressure-1', 'pressure-2', 'wind-stress') )
ax.set_xlabel( 'time (days)' )
ax.set_ylabel( 'power/unit area (J/s.m^2)' )
plt.show()
