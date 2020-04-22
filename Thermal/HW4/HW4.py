#3.24 Use a computer to study the entropy, temperature, and heat capacity of an Einstein solid, as follows, Let the solid contain 50 oscillators, and from 0 to 100 units of energy. Make a table, analogous to Table 3.2, in which each row represents a different value for the energy. Use separate columns for the energy, multiplicity, entropy, temperature, and heat capacity. To calculate the temperature, evaluate \Delta U / \Delta S for two nearby rows in the table. The heat capacity can be omputed in a similar way.

from numpy import arange, log, zeros_like, insert
from scipy.special import factorial
from matplotlib.pyplot import plot, show

q = arange(0,101,1)
N = 50
omega = factorial(q + N - 1)/(factorial(q) * factorial(N - 1)) 
S = log(omega)
temp = (2/ (S[2:] - S[:-2]))
temp = insert(temp,0,0)
capacity = 2 / (50 * (temp[2:] - temp[:-2]))

plot(temp[1:-1], capacity)
show()

