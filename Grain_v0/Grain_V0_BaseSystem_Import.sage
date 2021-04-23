import scipy.special
from itertools import combinations 
import timeit

start = timeit.timeit()

#define the paramaters

#how many clocks are required
no_clocks = 86000

#the size of the register
state_size = 80

#how many output bits are needed
no_output_bits=86000

#degree of the system of equations to be solved
degree = 3

#the number of monomials of degree 1->max degree
#the variables 'no_monomials' will be the number of linearisation variables
#'scipy.special.comb' provides the binomial coefficient
no_monomials = 0
for i in range(1,degree+1):
    no_monomials = no_monomials + binomial(state_size,i)

#define polynomial ring in which the the INITIAL STATE, OUTPUT and LINEARISATION variables will be defined
R=BooleanPolynomialRing(state_size+no_output_bits+no_monomials,'x')
#Inject these variables for use
R.inject_variables()

with open("Grain_V0_BaseSystem.txt","r") as myfile:
    EQ_1 = eval(myfile.read().replace('\n', ''))

print(EQ_1[100])