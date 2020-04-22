class diatomic:
    def __init__(self,bounds,dx,atomLocations,wellWidth,depth2 = 0):

        from numpy import arange, zeros_like
        #Here we define our problem, basically BC's and IC's
        self.bounds = bounds
        self.dx = dx
        self.atomLocations = atomLocations
        self.nAtoms = len(atomLocations)
        self.wellWidth = wellWidth
        self.x = arange(bounds[0], bounds[1]+dx, dx)
        self.V = zeros_like(self.x)

        self.V[(self.x > atomLocations[0] - self.wellWidth/2) & (self.x < atomLocations[0] + self.wellWidth/2)] = -20
        self.V[(self.x > atomLocations[1] - self.wellWidth/2) & (self.x < atomLocations[1] + self.wellWidth/2)] = -20 + depth2
        

    def plotV(self):
        #Plotting the potential energy
        from matplotlib.pyplot import plot,show
        plot(self.x,self.V)
        show()

    def sOrbital(self,atomLoc,x):
        #The evaluates the wave function
        from numpy import sqrt, exp, abs
        return 1/sqrt(0.57) * exp(-abs(x - atomLoc)/0.57)

    def plotBasis(self,atomLoc):
        #Plots the function we calculate in sOrbital
        from matplotlib.pyplot import plot, show
        plot(self.x,self.sOrbital(atomLoc,self.x))
        show()

    def DDbasis(self,atomLoc,x):
        from numpy import sqrt, exp, abs
        #returns the derivative of the sOrbital Equation
        #Useful for calculating kinetic energy PDF
        return 1/0.57**2 * 1/sqrt(0.57) * exp(-abs(x - atomLoc)/0.57)

    def HamiltonianMatrix(self):
        #We will build our eval/evec problem
        from numpy import zeros
        from scipy.integrate import trapz
        self.H = zeros([self.nAtoms,self.nAtoms])
        for iatom,atomOne in enumerate(self.atomLocations):#Iterate over all atoms
            for jatom,atomTwo in enumerate(self.atomLocations):
                self.H[iatom,jatom] = trapz(self.sOrbital(atomOne,self.x) * self.V * self.sOrbital(atomTwo,self.x), self.x) + self.KE(atomOne, atomTwo)#Make matrix out of integral and then sum of potential and kinetic energy terms

    def KE(self,i,j):
        from scipy.integrate import trapz
        return -trapz(self.sOrbital(i,self.x)*self.DDbasis(j,self.x), self.x)#calculate kinetic energy

    def solveProblem(self):
        from numpy.linalg import eig
        #use nunpy to solve our eignevalue problem
        vals,self.vecs = eig(self.H)
        self.key = sorted(range(len(vals)), key = lambda x:vals[x])#sort the indeces according to the eigenvals
        self.vals = sorted(vals)
        return (self.vecs[1,0]**2/self.vecs[0,0]**2, self.vecs[1,1]**2/self.vecs[0,1]**2)

    def plotSolution(self,n):
        from numpy import zeros_like,array
        from matplotlib.pyplot import plot,ylim,show
        solVec = self.vecs[:,self.key[n]]
        gsWaveFunc = zeros_like(self.x)
        for iatom,atom in enumerate(self.atomLocations):
            gsWaveFunc += solVec[iatom] * self.sOrbital(atom,self.x)#wavefunction is sum of individual solutions
        plot(gsWaveFunc)#plot wavefunction and square of wavefunction
        plot(gsWaveFunc**2, '--')
        ylim(-1,1)
        show()

dx = 0.001
bounds = [0,10]
atomLocations = [4,5]
wellWidth = 1
from numpy import arange
bonding = []
antibonding = []
for i in arange(0,20,0.1):
    myMolecule = diatomic(bounds,dx,atomLocations,wellWidth,depth2 = i)
    myMolecule.HamiltonianMatrix()
    new_bonding, new_antibonding = myMolecule.solveProblem()
    bonding.append(new_bonding)
    antibonding.append(new_antibonding)


from matplotlib.pyplot import plot, show
plot(bonding)
plot(antibonding)
show()