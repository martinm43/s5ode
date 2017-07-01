# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 20:49:39 2017

@author: Original as pNYU Physics as shown at http://www.physics.nyu.edu/pine/pymanual/html/chap9/chap9_scipy.html
"""

import numpy as np
from math import pi
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def f(y, t, params):
    theta, omega = y      # unpack current values of y
    I, D, M = params  # unpack parameters
    derivs = [omega,      # list of dy/dt=f functions
             -D*omega**2/I + M*np.sin(theta)/I]
    return derivs

# Parameters
I = 2.0      # Moment of inertia
D = 1.5      # Viscous/fluid damping factor
M = 0.65     # Magnetic force factor.

# Initial values
theta0 = pi/2     # initial angular displacement
omega0 = 0.0     # initial angular velocity

# Bundle parameters for ODE solver
params = [I, D, M]

# Bundle initial conditions for ODE solver
y0 = [theta0, omega0]

# Make time array for solution
tStop = 200.
tInc = 0.05
t = np.arange(0., tStop, tInc)

# Call the ODE solver
psoln = odeint(f, y0, t, args=(params,))

# Plot results
fig = plt.figure(1, figsize=(8,8))

# Plot theta as a function of time
ax1 = fig.add_subplot(311)
ax1.plot(t, psoln[:,0])
ax1.set_xlabel('time')
ax1.set_ylabel('theta')

# Plot omega as a function of time
ax2 = fig.add_subplot(312)
ax2.plot(t, psoln[:,1])
ax2.set_xlabel('time')
ax2.set_ylabel('omega')

# Plot omega vs theta
ax3 = fig.add_subplot(313)
twopi = 2.0*np.pi
ax3.plot(psoln[:,0]%twopi, psoln[:,1], '.', ms=1)
ax3.set_xlabel('theta')
ax3.set_ylabel('omega')
ax3.set_xlim(0., twopi)

plt.tight_layout()
plt.show()