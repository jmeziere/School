class diatomic:
    def __init__(self,bounds,dx,atomLocations,wellWidth,nBasis,V):

        from numpy import arange, zeros_like
        #Here we define our problem, basically BC's and IC's
        self.bounds = bounds
        self.dx = dx
        self.atomLocations = atomLocations
        self.wellWidth = wellWidth
        self.nBasis = nBasis
        self.x = arange(bounds[0], bounds[1]+dx, dx)
        self.a = bounds[1] - bounds[0]
        self.V = zeros_like(self.x)

        for iatom,atom in enumerate(self.atomLocations):
            self.V[(self.x > atom - self.wellWidth/2) & (self.x < atom + self.wellWidth/2)] = V[iatom]

    def basis(self,n,k,x):
        #The evaluates the wave function
        from numpy import sqrt, exp, pi
        return 1/sqrt(self.a) * exp(complex(0,1) * (2 * pi * n * x/self.a + k * x))

    def DDbasis(self,n,k,x):
        from numpy import sqrt, exp, pi
        #returns the derivative of the basis function
        #Useful for calculating kinetic energy PDF
        return -1/sqrt(self.a) * (2*pi*n/self.a+k)**2 *exp(complex(0,1) * (2 * pi * n * x/self.a + k * x))

    def HamiltonianMatrix(self,k):
        #We will build our eval/evec problem
        from numpy import zeros

        self.k = k
        self.H = zeros([self.nBasis+1,self.nBasis+1],dtype = complex)
        for iindex,ibasis in enumerate(range(-self.nBasis//2, self.nBasis//2 + 1)):#Iterate over all atoms
            for jindex,jbasis in enumerate(range(-self.nBasis//2, self.nBasis//2 + 1)):
                self.H[iindex,jindex] = self.dx * sum(self.basis(ibasis,k,self.x).conjugate() * self.V * self.basis(jbasis,k,self.x)) + self.KE(ibasis,jbasis)
                #Make matrix out of integral and then sum of potential and kinetic energy terms
        print(self.H)

    def KE(self,i,j):
        return -self.dx * sum(self.basis(i,self.k,self.x).conjugate() * self.DDbasis(j,self.k,self.x))
        #calculate kinetic energy

    def solveProblem(self):
        from numpy.linalg import eig
        #use nunpy to solve our eignevalue problem
        vals,self.vecs = eig(self.H)
        self.key = sorted(range(len(vals)), key = lambda x:vals[x])#sort the indeces according to the eigenvals
        self.vals = sorted(vals)

    def plotSolution(self,n):
        from numpy import zeros_like,linspace
        from matplotlib.pyplot import plot,figure,show
        
        plotx,dx = linspace(0,2 * self.a, 1000,retstep = True)
        solVec = self.vecs[:,self.key[n]]
        gsWaveFunc = zeros_like(plotx,dtype = complex)
        for iindex,ibasis in enumerate(range(-self.nBasis//2, self.nBasis//2 + 1)):
            gsWaveFunc += solVec[iindex] * self.basis(ibasis,self.k,plotx)
            #wavefunction is sum of individual solutions
        
        normalization = dx * sum(gsWaveFunc.conjugate() * gsWaveFunc)
        figure(1)
        plot(plotx,gsWaveFunc.conjugate() * gsWaveFunc / normalization,'r--')#plot wavefunction and square of wavefunction
        show()

from numpy import pi
nBasis = 5
dx = 0.001
bounds = [0,3]
atomLocations = [1.5]
wellWidth = 1
V = [-20]
myMolecule = diatomic(bounds,dx,atomLocations,wellWidth,nBasis,V)
myMolecule.HamiltonianMatrix(pi/3/2)
myMolecule.solveProblem()
for i in range(50):
    myMolecule.plotSolution(i)
