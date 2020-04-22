class diatomic:
    def __init__(self,bounds,dx,atomLocations,wellWidth,nBasis,V):

        from numpy import arange, full_like
        #Here we define our problem, basically BC's and IC's
        self.bounds = bounds
        self.dx = dx
        self.atomLocations = atomLocations
        self.wellWidth = wellWidth
        self.nBasis = nBasis
        self.x = arange(bounds[0], bounds[1]+dx, dx)
        self.a = bounds[1] - bounds[0]
        self.V = full_like(self.x,-10000)

        for iatom,atom in enumerate(self.atomLocations):
            self.V[(self.x < atom - 0.01) | (self.x > atom + 0.01)] = -1/(self.x[(self.x < atom - 0.01) | (self.x > atom + 0.01)] - atom)**2

        from matplotlib.pyplot import plot, show
        plot(self.V)
        show()

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

    def KE(self,i,j):
        return -self.dx * sum(self.basis(i,self.k,self.x).conjugate() * self.DDbasis(j,self.k,self.x))
        #calculate kinetic energy

    def solveProblem(self):
        from numpy.linalg import eig
        #use nunpy to solve our eignevalue problem
        vals,self.vecs = eig(self.H)
        self.key = sorted(range(len(vals)), key = lambda x:vals[x])#sort the indeces according to the eigenvals
        self.vals = sorted(vals)
        return self.vals

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

    def plotBand(self):
        from numpy import arange, pi, full_like
        all_evals = []
        all_kvals = []
        for i in arange(-pi/self.a,pi/self.a,0.1):
            self.HamiltonianMatrix(i)
            evals = self.solveProblem()
            kvals = full_like(evals,i)
            all_evals.append(evals)
            all_kvals.append(kvals)
        
        from matplotlib.pyplot import plot, show
        plot(all_kvals, all_evals, 'bo',markersize = 0.5)
        show()


from numpy import pi
nBasis = 10
dx = 0.001
bounds = [0,3]
atomLocations = [1.5]
wellWidth = 1
V = [-20]
myMolecule = diatomic(bounds,dx,atomLocations,wellWidth,nBasis,V)
myMolecule.HamiltonianMatrix(pi/3/2)
myMolecule.solveProblem()
myMolecule.plotBand()
for i in range(2):
    myMolecule.plotSolution(i)
