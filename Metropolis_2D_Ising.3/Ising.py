from pylab import *
from numpy import random
import sys
from scipy.optimize import curve_fit

def Monte_Carlo():

 # Preform N random spin flips
 for step in range(N):
  # Draw random spin
  ii=int(random.random()*Nx)
  jj=int(random.random()*Nx)
  # attempt to flip
  Metropolis(ii,jj)

def Metropolis(ii,jj):

  # Store old spin
  spinOLD=spins[ii][jj]
  # New spin flips sign
  spinNEW=-1*spinOLD

  # Neighbors on a square lattice
  # Enforcing periodic boundary conditions
  up = spins[(ii+1)%Nx][jj]
  down = spins[(ii-1)%Nx][jj]
  left = spins[ii][(jj-1)%Nx]
  right = spins[ii][(jj+1)%Nx]

  neighs= up+down+left+right

  # Compute change in energy
  deltaE = - h * (spinNEW - spinOLD) - J * neighs * (spinNEW - spinOLD)########INSERT CODE#########
  
  # Metropolis Acceptance Criteria
  # Draw random number btw 0:1
  # If random number is  < change in Boltzmann dist
  # accept, else continue
  if(random.random()< exp(- deltaE / T)): ######INSERT CODE######
   # If move accepted, change spin
   spins[ii][jj]=spinNEW 


# Parameters #
##
##
# Temperature
T = float(sys.argv[1])
# Exchange coupling
J = 1.
# External field
h_collection = arange(-0.005, 0.006, 0.001)
# External field divided by J
#hinvJ = h/J
# Size of lattice, Nx x Nx
Nx = 15
# Total number of spins
N = Nx*Nx
# Number of sweeps
Nsweep = 3000
# Averaging frequency
Nsamp = 10
# Initialize flag
StartFlag = 0
# If 1 make plot of time versus M/N
plotflag = 0
# If 1 make image of system
imflag = 0


# For Ising model in 2d
# If StartFlag = 1 start all down
if StartFlag:
 spins=zeros([Nx,Nx])-1
else:
 # otherwise start in random arangement
 temp=random.randint(2,size=(N))*2-1
 spins=reshape(temp,(Nx,Nx))

# Container for magnetization and fluctuations
netMag=array([]) #instantaneous total magnetization
avgMag=array([]) #running average of instantaneous magnetization per spin, for plotting only
corr={}
m_series = []
m_sd_series = []
with open("magT%.2f.dat"%(T), "w") as file:
  for h in h_collection:
    print "MC simulation"
    print "T = %.2f h = %.2f J = %.2f" %(T,h,J)
    print "Nx = %d by Nx = %d" %(Nx,Nx)
    print "Nsweeps = %d and Nsamp = %d" %(Nsweep,Nsamp)
    print
    print "Step \t M \t\t <M>"
    print "%d \t %.4f \t %.4f" %(0, float(sum(spins))/N, sum(spins)/N)

    file.write("MC simulation\n")
    file.write("T = %.2f h = %.2f J = %.2f\n" %(T,h,J))
    file.write("Nx = %d by Nx = %d\n" %(Nx,Nx))
    file.write("Nsweeps = %d and Nsamp = %d\n" %(Nsweep,Nsamp))
    file.write("\n")
    file.write("Step \t M \t\t <M>\n")
    file.write("%d \t %.4f \t %.4f\n" %(0, sum(spins)/N, sum(spins)/N))

  # Preform Nsweeps through the lattice
    for step in range(Nsweep): 
      # Perform N random spin flips
      Monte_Carlo()
      # Every Nsamp steps, compute magnetization
      # and its running average
      if(step%Nsamp==0 and step>0):
        tempmag=sum(spins)
        netMag=append(netMag,tempmag)
        avgMag=append(avgMag,mean(netMag)/N)
        print "%d \t %.4f \t %.4f" %(step, float(tempmag)/N, mean(netMag)/N)
        file.write("%d \t %.4f \t %.4f\n" %(step, float(tempmag)/N, mean(netMag)/N))
        # Compute correlation
        # for i1 in range(Nx/2):
        #   for j1 in range(Nx/2):
        #     for i2 in range(i1,Nx/2):
        #       for j2 in range(j1,Nx/2):
        #         sisj=float(spins[i1][j1])*float(spins[i2][j2])
        #         si=float(spins[i1][j1])
        #         sj=float(spins[i2][j2])
        #         dist = sqrt(((i2-i1)%(Nx))**2+((j2-j1)%(Nx))**2)
        #         if dist not in corr:
        #           corr[dist]=[]
        #         corr[dist].append([sisj,si,sj])
      
    # rij=[]
    # cij=[]
    # for dist in corr:
    #   sisjls=[elem[0] for elem in corr[dist]]
    #   sils=[elem[1] for elem in corr[dist]]
    #   sjls=[elem[2] for elem in corr[dist]]
    #   rij.append(dist)
    #   cij.append(mean(sisjls)-mean(sils)*mean(sjls))

    
    varMag = mean([(mag - mean(netMag)/N)**2 for mag in avgMag]) # mean(netMag)/N is our best estimation to ensemble average <M>
    errorMag = sqrt(varMag)/(Nsweep/Nsamp)
    m_series.append(mean(netMag)/N)
    m_sd_series.append(varMag)
    print "<M> = %.4f \t with sampling error = %.8f \t and Var(M) = %.8f " %(mean(netMag)/N, errorMag, varMag)
    file.write("<M> = %.4f \t with sampling error = %.8f \t and Var(M) = %.8f \n" %(mean(netMag)/N, errorMag, varMag))


# with open("corrtable.dat","a") as file:
#   file.write("%.2f \n"%T)
#   for item in rij:
#     file.write("%f,"%item)
#   file.write("\n")
#   for item in cij:
#     file.write("%f,"%item)
#   file.write("\n")

# def func(x, a, b, c):
#   return a * exp(-b * x) + c

# popt, pcov = curve_fit(func, rij, cij)
# x = arange(0, 15, 0.1)
# y = func(x, *popt)
# plot(x, y)
# scatter(rij, cij)
# title("T=%s" % T)
# xlabel(r'$\vec{r_{ij}}$')
# ylabel(r'$c_{ij}$')
# savefig("T=%s.png" % T)
# show()

if plotflag:
 plot(range(Nsweep/Nsamp-1),avgMag,'-b')
 ylim([-1,1])
 xlabel(r'MC time',size=20)
 ylabel(r'$<m>$',size=20)
 savefig('MagnetizationT%.1f.png'%(T))
 cla()

if imflag:
 imshow(spins,cmap='binary',interpolation='nearest')
 xlabel(r'$x$',size=20)
 ylabel(r'$y$',size=20)
 savefig('Ising-ConfigT%.1f.png'%(T))
 cla()


