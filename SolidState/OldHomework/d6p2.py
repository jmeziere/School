class diatomic:
    def __init__(self,bounds,dx,atomLocations,wellWidth,nBasis):

        from numpy import arange, zeros_like
        #Here we define our problem, basically BC's and IC's
        self.bounds = bounds
        self.dx = dx
        self.nBasis = nBasis
        self.wellWidth = wellWidth
        self.x = arange(bounds[0], bounds[1]+dx, dx)
        self.a = (bounds[1] - bounds[0])
        self.V = zeros_like(self.x)

        for iatom,atom in enumerate(atomLocations):
            self.V[(self.x > atom - self.wellWidth/2) & (self.x < atom + self.wellWidth/2)] = -20

    def basis(self,n,x):
        #The evaluates the wave function
        from numpy import sqrt, sin, pi
        return sqrt(2/self.a) * sin(n * pi * x / self.a)

    def DDbasis(self,n,x):
        from numpy import sqrt, sin, pi
        #returns the derivative of the sOrbital Equation
        #Useful for calculating kinetic energy PDF
        return -sqrt(2/self.a) * n**2 * pi**2 * sin(n * pi * x / self.a) / self.a**2

    def HamiltonianMatrix(self):
        #We will build our eval/evec problem
        from numpy import zeros
        from scipy.integrate import trapz
        self.H = zeros([self.nBasis,self.nBasis])
        for i in range(1,self.nBasis+1):#Iterate over all atoms
            for j in range(1,self.nBasis+1):
                self.H[i-1,j-1] = trapz(self.basis(i,self.x) * self.V * self.basis(j,self.x), self.x) + self.KE(i, j)#Make matrix out of integral and then sum of potential and kinetic energy terms

    def KE(self,i,j):
        from scipy.integrate import trapz
        return -trapz(self.basis(i,self.x)*self.DDbasis(j,self.x), self.x)#calculate kinetic energy

    def solveProblem(self):
        from numpy.linalg import eig
        from numpy import array
        #use nunpy to solve our eignevalue problem
        vals,self.vecs = eig(self.H)
        self.key = sorted(range(len(vals)), key = lambda x:vals[x])#sort the indeces according to the eigenvals
        self.vals = array(sorted(vals))
        return(self.vals)

    def plotSolution(self,n):
        from numpy import zeros_like
        from matplotlib.pyplot import plot,ylim,show
        solVec = self.vecs[:,self.key[n]]
        gsWaveFunc = zeros_like(self.x)

        for i in range(1,self.nBasis+1):
            new_Func = solVec[i-1] * self.basis(i,self.x)
            gsWaveFunc += new_Func #wavefunction is sum of individual solutions
        plot(self.x,gsWaveFunc)#plot wavefunction and square of wavefunction
        plot(self.x,gsWaveFunc**2, '--')
        show()

dx = 0.001
bounds = [0,10]
atomLocations = [3,6]
wellWidth = 1
#nBasis = 100

evals = []
for nBasis in range(10,101,10):
    myMolecule = diatomic(bounds,dx,atomLocations,wellWidth,nBasis)
    myMolecule.HamiltonianMatrix()
    evals.append(myMolecule.solveProblem())
    #myMolecule.plotSolution(0)
    #myMolecule.plotSolution(1)
    print(nBasis)

error = []
for j in evals:
    from numpy import abs, sum
    error.append(sum(abs(j - (evals[-1])[:len(j)])))

from matplotlib.pyplot import show, plot
plot(range(10,101,10),error,'.')
show()
