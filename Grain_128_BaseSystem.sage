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

#define initial state
s = [0]*state_size
for i in range(state_size):
    s[i] = eval("x" + str(i))

# print(s)
#define output
Y = [0]*no_output_bits
for i in range(no_output_bits):
    Y[i] = eval("x" + str(state_size+i))

# print(Y)
#We can consider 's' as a register that will be shifted for a specified number of times. Thus, there will exist more entries than simply those defined in 's'.
#We therefore initialise a vector 'add' that will be appended to s to hold all the update bits (one vector for equations with respect to keystream bits S and one with respect the the linear approximation$
add=[0]*no_clocks
S = s+add
Left = s+add
Right = s+add

#Generate equations for update using output bits
# S[15] = S[0]+S[7]
for i in range(len(s),no_clocks):
     S[i] = (S[i-80]+S[i-67]+S[i-57]+S[i-42]+S[i-29]+S[i-18])
     
#Generate equations for update using linear approximation
for i in range(len(s),no_clocks):
    #s1s2z=s1s2s9
    #Here we are going to build two systems; one containing z, the other not
    # Gz[i] = (S[i-5]+S[i-14]*S[i-9]+S[i-12]*S[i-9]+S[i-9]*S[i-5])*Y[i-15]
    #G[i]= (S[i-5]+S[i-14]*S[i-9]+S[i-12]*S[i-9]+S[i-9]*S[i-5])*(S[i-12]+S[i-14]*S[i-5]+S[i-9]*S[i-5]+S[i-14]*S[i-12]*S[i-9]+S[i-14]*S[i-9]*S[i-5])
    #Left[i] = S[i-(state_size-3)]*S[i-(state_size-25)]*S[i-(state_size-46)]*S[i-(state_size-64)]+S[i-(state_size-3)]*S[i-(state_size-25)]*S[i-(state_size-46)]+S[i-(state_size-3)]*S[i-(state_size-46)]*S[i-(state_size-64)]+S[i-(state_size-1)]*S[i-(state_size-3)]*S[i-(state_size-46)]+S[i-(state_size-25)]*S[i-(state_size-46)]*S[i-(state_size-64)]+S[i-(state_size-25)]*S[i-(state_size-46)]*S[i-(state_size-64)]+S[i-(state_size-1)]*S[i-(state_size-25)]*S[i-(state_size-46)]+S[i-(state_size-1)]*S[i-(state_size-46)]*S[i-(state_size-64)]
    Right[i] = import scipy.special
from itertools import combinations 
import timeit

start = timeit.timeit()

#define the paramaters

#how many clocks are required
no_clocks = 275585000

#the size of the register
state_size = 128

#how many output bits are needed
no_output_bits=275585000

#degree of the system of equations to be solved
degree = 5

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

#define initial state
s = [0]*state_size
for i in range(state_size):
    s[i] = eval("x" + str(i))

# print(s)
#define output
Y = [0]*no_output_bits
for i in range(no_output_bits):
    Y[i] = eval("x" + str(state_size+i))

# print(Y)
#We can consider 's' as a register that will be shifted for a specified number of times. Thus, there will exist more entries than simply those defined in 's'.
#We therefore initialise a vector 'add' that will be appended to s to hold all the update bits (one vector for equations with respect to keystream bits S and one with respect the the linear approximation$
add=[0]*no_clocks
S = s+add
Left = s+add
Right = s+add

#Generate equations for update using output bits
# S[15] = S[0]+S[7]
for i in range(len(s),no_clocks):
     S[i] = ((S[i-32]+S[i-47]+S[i-58]+S[i-90]+S[i-121]+S[i-128])
     
#Generate equations for update using linear approximation
for i in range(len(s),no_clocks):
    #s1s2z=s1s2s9
    #Here we are going to build two systems; one containing z, the other not
    # Gz[i] = (S[i-5]+S[i-14]*S[i-9]+S[i-12]*S[i-9]+S[i-9]*S[i-5])*Y[i-15]
    #G[i]= (S[i-5]+S[i-14]*S[i-9]+S[i-12]*S[i-9]+S[i-9]*S[i-5])*(S[i-12]+S[i-14]*S[i-5]+S[i-9]*S[i-5]+S[i-14]*S[i-12]*S[i-9]+S[i-14]*S[i-9]*S[i-5])
    #Left[i] = S[i-(state_size-3)]*S[i-(state_size-25)]*S[i-(state_size-46)]*S[i-(state_size-64)]+S[i-(state_size-3)]*S[i-(state_size-25)]*S[i-(state_size-46)]+S[i-(state_size-3)]*S[i-(state_size-46)]*S[i-(state_size-64)]+S[i-(state_size-1)]*S[i-(state_size-3)]*S[i-(state_size-46)]+S[i-(state_size-25)]*S[i-(state_size-46)]*S[i-(state_size-64)]+S[i-(state_size-25)]*S[i-(state_size-46)]*S[i-(state_size-64)]+S[i-(state_size-1)]*S[i-(state_size-25)]*S[i-(state_size-46)]+S[i-(state_size-1)]*S[i-(state_size-46)]*S[i-(state_size-64)]
    Right[i] = (S[i-(state_size-8)]+1)*(S[i-(state_size-42)]+1)*(S[i-(state_size-95)]+1)*Y[i-state_size]+S[i-(state_size-2)]*S[i-(state_size-8)]*S[i-(state_size-42 )]+ S[i-(state_size-2)]*S[i-(state_size-8)]*S[i-(state_size-95 )]+ S[i-(state_size-2)]*S[i-(state_size-8 )]+ S[i-(state_size-2)]*S[i-(state_size-42)]*S[i-(state_size-95 )]+ S[i-(state_size-2)]*S[i-(state_size-42 )]+ S[i-(state_size-2)]*S[i-(state_size-95 )]+ S[i-(state_size-2 )]+ S[i-(state_size-15)]*S[i-(state_size-8)]*S[i-(state_size-42 )]+ S[i-(state_size-15)]*S[i-(state_size-8)]*S[i-(state_size-95 )]+ S[i-(state_size-15)]*S[i-(state_size-8 )]+ S[i-(state_size-15)]*S[i-(state_size-42)]*S[i-(state_size-95 )]+ S[i-(state_size-15)]*S[i-(state_size-42 )]+ S[i-(state_size-15)]*S[i-(state_size-95 )]+ S[i-(state_size-15 )]+ S[i-(state_size-36)]*S[i-(state_size-8)]*S[i-(state_size-42 )]+ S[i-(state_size-36)]*S[i-(state_size-8)]*S[i-(state_size-95 )]+ S[i-(state_size-36)]*S[i-(state_size-8 )]+ S[i-(state_size-36)]*S[i-(state_size-42)]*S[i-(state_size-95 )]+ S[i-(state_size-36)]*S[i-(state_size-42 )]+ S[i-(state_size-36)]*S[i-(state_size-95 )]+ S[i-(state_size-36 )]+ S[i-(state_size-45)]*S[i-(state_size-8)]*S[i-(state_size-42 )]+ S[i-(state_size-45)]*S[i-(state_size-8)]*S[i-(state_size-95 )]+ S[i-(state_size-45)]*S[i-(state_size-8 )]+ S[i-(state_size-45)]*S[i-(state_size-42)]*S[i-(state_size-95 )]+ S[i-(state_size-45)]*S[i-(state_size-42 )]+ S[i-(state_size-45)]*S[i-(state_size-95 )]+ S[i-(state_size-45 )]+ S[i-(state_size-64)]*S[i-(state_size-8)]*S[i-(state_size-42 )]+ S[i-(state_size-64)]*S[i-(state_size-8)]*S[i-(state_size-95 )]+ S[i-(state_size-64)]*S[i-(state_size-8 )]+ S[i-(state_size-64)]*S[i-(state_size-42)]*S[i-(state_size-95 )]+ S[i-(state_size-64)]*S[i-(state_size-42 )]+ S[i-(state_size-64)]*S[i-(state_size-95 )]+ S[i-(state_size-64 )]+ S[i-(state_size-73)]*S[i-(state_size-8)]*S[i-(state_size-42 )]+ S[i-(state_size-73)]*S[i-(state_size-8)]*S[i-(state_size-95 )]+ S[i-(state_size-73)]*S[i-(state_size-8 )]+ S[i-(state_size-73)]*S[i-(state_size-42)]*S[i-(state_size-95 )]+ S[i-(state_size-73)]*S[i-(state_size-42 )]+ S[i-(state_size-73)]*S[i-(state_size-95 )]+ S[i-(state_size-73 )]+ S[i-(state_size-89)]*S[i-(state_size-8)]*S[i-(state_size-42 )]+ S[i-(state_size-89)]*S[i-(state_size-8)]*S[i-(state_size-95 )]+ S[i-(state_size-89)]*S[i-(state_size-8 )]+ S[i-(state_size-89)]*S[i-(state_size-42)]*S[i-(state_size-95 )]+ S[i-(state_size-89)]*S[i-(state_size-42 )]+ S[i-(state_size-89)]*S[i-(state_size-95 )]+ S[i-(state_size-89 )]+ S[i-(state_size-93)]*S[i-(state_size-8)]*S[i-(state_size-42 )]+ S[i-(state_size-93)]*S[i-(state_size-8)]*S[i-(state_size-95 )]+ S[i-(state_size-93)]*S[i-(state_size-8 )]+ S[i-(state_size-93)]*S[i-(state_size-42)]*S[i-(state_size-95 )]+ S[i-(state_size-93)]*S[i-(state_size-42 )]+ S[i-(state_size-93)]*S[i-(state_size-95 )]+ S[i-(state_size-93 )]+ S[i-(state_size-8)]*S[i-(state_size-13)]*S[i-(state_size-20 )]+ S[i-(state_size-8)]*S[i-(state_size-60)]*S[i-(state_size-79 )]+ S[i-(state_size-13)]*S[i-(state_size-20)]*S[i-(state_size-42 )]+ S[i-(state_size-13)]*S[i-(state_size-20)]*S[i-(state_size-95 )]+ S[i-(state_size-13)]*S[i-(state_size-20 )]+ S[i-(state_size-42)]*S[i-(state_size-60)]*S[i-(state_size-79 )]+ S[i-(state_size-60)]*S[i-(state_size-79)]*S[i-(state_size-95 )]+ S[i-(state_size-60)]*S[i-(state_size-79)]
    #for now we are just going to focus on the system without z. This is the system whose degree we want to reduce.

#create a matrix that will hold the final equations (the sum of those given with respect to the keystream bits + those given by the linear approximation)
EQ = list(var('eq%d' %i) for i in (0..no_clocks))

# produce final system
for i in range(len(s),no_clocks):
    EQ[i] = (Right[i])

# get rid of redundant equations
EQ_1=EQ[len(s):no_clocks]

# original_output = sys.stdout
# fp = open("Grain_V0_BaseSystem.txt","w")
# sys.stdout = fp

# print(EQ_1) 
Alpha = []
check = 0
for i in Alpha:
    check = check+EQ_1[i]
print(check)
    #for now we are just going to focus on the system without z. This is the system whose degree we want to reduce.