import scipy.special
from itertools import combinations 
import timeit
import random

start = timeit.timeit()
#define the paramaters

#how many clocks are required
no_clocks = 5000

#the size of the register
state_size = 15

#how many output bits are needed
no_output_bits=5000

#degree of the system of equations to be solved
degree = 4

#the number of monomials of degree 1->max degree
#the variables 'no_monomials' will be the number of linearisation variables
#'scipy.special.comb' provides the binomial coefficient


#Calculate the number of monomials there will be in the system
#This will be used for the linearisation process
no_monomials = 0
for i in range(1,degree+1):
    no_monomials = no_monomials + binomial(state_size,i)

#define polynomial ring in which the the INITIAL STATE and LINEARISATION variables will be defined
#Inject these variables for use
R=BooleanPolynomialRing(state_size+no_output_bits+no_monomials,'x')
R.inject_variables()

#define initial state variables
s = [0]*state_size
for i in range(state_size):
    s[i] = eval("x" + str(i))

#define output
Y = [0]*no_output_bits
for i in range(no_output_bits):
    Y[i] = eval("x" + str(state_size+i))
    
with open("Better_n_15_System.txt", "r") as f:
    for line in f:
        Final_System = eval(line.strip())
print('done')

# #Randomly load the NFSR and LFSR
# #Print the loaded states so the initial state can be viewed on execution
# s_out = [0]*state_size
# b_out = [0]*state_size

# for i in range(state_size):
#     s_out[i] = random.randint(0,1)
#     b_out[i] = random.randint(0,1)

# print(s_out)
# print(b_out)


# #Initialise state sequence arrays for the NFSR and the LFSR
# #Produce output (assumed to be known to the attack)
# #Update the state of each register
# #Output -> (s3s4)*b1+(s1+s2)
# #NFSR update -> B_out[i] = 1+b2b3
# #LFSR update ->  S_out[i] = s0+s3
# add=[0]*no_clocks
# S_out = s_out+add
# B_out = b_out+add

# O = [0]*(no_clocks-state_size)
# for i in range(state_size,no_clocks):
#     O[i-state_size] =mod(S_out[i-(state_size-1)]+S_out[i-(state_size-2)]+S_out[i-(state_size-3)]*S_out[i-(state_size-4)]*B_out[i-(state_size-1)],2)
#     B_out[i] = mod(1+B_out[i-(state_size-2)]*B_out[i-(state_size-3)],2)
#     S_out[i] = mod((S_out[i-5]+S_out[i-2]),2)

# for i in range(len(Final_System)):
#     for j in range(len(O)):
#         m = eval("x" + str(state_size+j))
#         n = eval('O['+str(j)+']')
#         Final_System[i]=Final_System[i].subs({m:n})


