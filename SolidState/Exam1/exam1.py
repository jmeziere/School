print("This file contains problem 2\n")
print("For problem 2c I used the following values:\np=20\nl=1\na=10\nnBasis=100\n")

class schrodinger:
    def __init__(self,p,l,a,dx,nBasis):
        from numpy import arange, zeros_like
        #Setup for the problem            
        self.a = a
        self.dx = dx
        self.nBasis = nBasis
        self.x = arange(0,a+dx,dx)
        self.V = zeros_like(self.x)
        for i in range(len(self.x)):
            if self.x[i] > a/2-l and self.x[i] < a/2+l:
                 self.V[i] = p

    def basis(self,basIndex):
        from numpy import sqrt, sin, pi
        #Setup our basis functions
        return sqrt(2/self.a) * sin(basIndex * pi * self.x / self.a)

    def DDbasis(self,basIndex):
        from numpy import sqrt, sin, pi
        #Take the second derivative for our basis functions
        return -sqrt(2/self.a) * (basIndex * pi / self.a)**2 * sin(basIndex * pi * self.x / self.a)

    def makeMatrix(self):
        from numpy import zeros
        #Make our Hamiltonian matrix
        self.H = zeros((self.nBasis,self.nBasis))

        for i in range(self.nBasis):
            for j in range(self.nBasis):
                self.H[i,j] = self.dx * sum(self.basis(i+1) * self.V * self.basis(j+1)) -self.dx * sum(self.basis(i+1) * self.DDbasis(j+1))

    def solveMatrix(self):
        from numpy.linalg import eig
        #Solve the eigenvalue problem
        evals, self.evecs = eig(self.H)
        self.key = sorted(range(len(evals)),key = lambda k:evals[k])
        self.evals = sorted(evals)
        
    def plotSolution(self,n):
        #Plot our solutions
        from matplotlib.pyplot import plot, show

        solVec = self.evecs[:,self.key[n]]
        solution = 0
        for i in range(len(solVec)):
            solution += solVec[i] * self.basis(i+1)
        plot(self.x,solution)
        show()


a = 10
dx = 0.001
nBasis = 100
l = 1
p = 20

mySchrodinger = schrodinger(p,l,a,dx,nBasis)
mySchrodinger.makeMatrix()
mySchrodinger.solveMatrix()
mySchrodinger.plot(0)
mySchrodinger.plot(1)


