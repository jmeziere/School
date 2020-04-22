class diatomic:
    def __init__(self,dx,nAtoms,wellWidth):

        from numpy import arange, zeros_like
        #Here we define our problem, basically BC's and IC's
        self.dx = dx
        self.atomLocations = arange(1,2 * nAtoms,2)
        self.nAtoms = nAtoms
        self.wellWidth = wellWidth
        self.x = arange(0, 2*nAtoms+dx, dx)
        self.V = zeros_like(self.x)
        
        for atom in self.atomLocations:
            self.V[(self.x > atom - self.wellWidth/2) & (self.x < atom + self.wellWidth/2)] = -20

    def plotV(self):
        #Plotting the potential energy
        from matplotlib.pyplot import plot,show
        plot(self.x,self.V)
        show()

    def sOrbital(self,atomLoc,x):
        #The evaluates the wave function
        from numpy import sqrt, exp, abs
        return 1/sqrt(0.57) * exp(-abs(x - atomLoc)/0.57)

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
        return self.vals

    def plotSolution(self,n):
        from numpy import zeros_like,array
        from matplotlib.pyplot import plot,ylim,show
        solVec = self.vecs[:,self.key[n]]
        gsWaveFunc = zeros_like(self.x)
        for iatom,atom in enumerate(self.atomLocations):
            gsWaveFunc += solVec[iatom] * self.sOrbital(atom,self.x)#wavefunction is sum of individual solutions
        plot(gsWaveFunc)#plot wavefunction and square of wavefunction
        plot(gsWaveFunc**2, '--')
        show()

dx = 0.001
wellWidth = 1
from matplotlib.pyplot import plot, show
for natoms in [5,10,20,50]:
    myMolecule = diatomic(dx,natoms,wellWidth)
    myMolecule.HamiltonianMatrix()
    evals = myMolecule.solveProblem()
    myMolecule.plotSolution(0)
    myMolecule.plotSolution(1)
    temp_eval = evals.copy()
    evals.reverse()
    evals.extend(temp_eval)
    plot(evals,'.')
    show()
