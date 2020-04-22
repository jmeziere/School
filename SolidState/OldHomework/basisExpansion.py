class basisExpansion:
    def __init__(self,atomicSpacing,dx,atomicLocation,wellWidth):
        from numpy import arange, zeros_like

        self.dx = dx
        self.a = atomicSpacing
        self.atomicLacation = atomicLocation
        self.x = arange(0,3 * (atomicSpacing+1),dx)
        self.wellWidth = wellWidth

        self.V = zeros_like(self.x)
        for iatom,atom in enumerate(arange(atomicLocation,3*atomicLocation,atomicSpacing)):
            self.V[(self.x > atom - self.wellWidth/2) & (self.x < atom + self.wellWidth/2)] = -20
    
    def basis(self,n,l,atomLoc,x):
        from numpy import sqrt, exp
        lookup = {'1s':1/sqrt(0.15) * exp(-abs(x - atomLoc)/0.15),
                  '2s':1/(2*sqrt(2)*0.15**(3/2))*(2-abs(x-atomLoc)/0.15)*exp(-abs(x-atomLoc)/0.15)}

        return lookup[str(n) + l]

    def plotBasis(self,n,l,atomLoc):
        from matplotlib.pyplot import plot, show
        from numpy import linspace

        x = linspace(-self.a,2*self.a,1000)
        plot(x,self.basis(n,l,atomLoc,x))
        show()

spacing = 5
location = 3
dx = .001
wellWidth = 1
myBasisExp = basisExpansion(spacing,dx,location,wellWidth)
myBasisExp.plotBasis(2,'s',location)
