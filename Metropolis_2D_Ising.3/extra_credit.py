from pylab import *
from joblib import Parallel, delayed
import sys
from scipy.optimize import curve_fit

def monte_carlo():
	for time_step in range(N):
		i, j = int(random.random()*Nx), int(random.random()*Nx)
		Metropolis(i, j)

def Metropolis(i, j):
	spin_old = spins[i][j]
	spin_new = - 1 * spin_old
	up = spins[(i + 1) % Nx][j]
	down = spins[(i - 1) % Nx][j]
	left = spins[i][(j - 1) % Nx]
	right = spins[i][(j + 1) % Nx]
	neighbor = up + down + left + right
	del_E = - h * (spinNEW - spinOLD) - J * neighs * (spinNEW - spinOLD)
	if(random.random()< exp(- deltaE / T)):
		spins[i][j] = spin_new



# Parameters #
##
##
# Temperature
T = float(sys.argv[1])
# Exchange coupling
J = float(sys.argv[2])
# External field
h = float(sys.argv[3])
# External field divided by J
#hinvJ = h/J
# Size of lattice, Nx x Nx
Nx = 20
# Total number of spins
N = Nx*Nx
# Number of sweeps
Nsweep = float(sys.argv[4])
# Averaging frequency
Nsamp = 10


