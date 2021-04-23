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

#Randomly load the NFSR and LFSR
#Print the loaded states so the initial state can be viewed on execution
s_out = [0]*state_size
b_out = [0]*state_size

for i in range(state_size):
    s_out[i] = random.randint(0,1)
    b_out[i] = random.randint(0,1)

print(s_out)
print(b_out)


#Initialise state sequence arrays for the NFSR and the LFSR
#Produce output (assumed to be known to the attack)
#Update the state of each register
#Output -> (s13+s1s9+s4s9+s9s13)*b4+(s4+s1s13+s9s13+s1s4s9+s1s9s13)
#NFSR update -> B_out[i] = s0 + b0 + b13 + b2b14 + b3b5
#LFSR update ->  S_out[i] = s0+s7
add=[0]*no_clocks
S_out = s_out+add
B_out = b_out+add

O = [0]*(no_clocks-state_size)
for i in range(state_size,no_clocks):
    O[i-state_size] =mod((S_out[i-(state_size-13)]+S_out[i-(state_size-1)]*S_out[i-(state_size-9)]+S_out[i-(state_size-4)]*S_out[i-(state_size-9)]+S_out[i-(state_size-9)]*S_out[i-(state_size-13)])*B_out[i-(state_size-4)]+(S_out[i-(state_size-4)]+S_out[i-(state_size-1)]*S_out[i-(state_size-13)]+S_out[i-(state_size-9)]*S_out[i-(state_size-13)]+S_out[i-(state_size-1)]*S_out[i-(state_size-4)]*S_out[i-(state_size-9)]+S_out[i-(state_size-1)]*S_out[i-(state_size-9)]*S_out[i-(state_size-13)]),2)
    B_out[i] = mod(S_out[i-(state_size-0 )]+ B_out[i-(state_size-0 )]+ B_out[i-(state_size-13 )]+ B_out[i-(state_size-2)]*B_out[i-(state_size-14 )]+ B_out[i-(state_size-3)]*B_out[i-(state_size-5 )],2)
    S_out[i] = mod((S_out[i-15]+S_out[i-8]),2)


#Calculate the number of monomials there will be in the system
#This will be used for the linearisation process
no_monomials = 0
for i in range(1,degree+1):
    no_monomials = no_monomials + binomial(state_size,i)

#define polynomial ring in which the the INITIAL STATE and LINEARISATION variables will be defined
#Inject these variables for use
R=BooleanPolynomialRing(state_size+no_monomials,'x')
R.inject_variables()

#define initial state variables
s = [0]*state_size
for i in range(state_size):
    s[i] = eval("x" + str(i))

#Initialise arrays for LFSR state sequence and system of equations to be solved
add=[0]*no_clocks
S = s+add
Final_System = [0]*(no_clocks-state_size)

for j in range(len(Final_System)):
#     Final_System[i] = S[(i+1)]*S[(i+3 )]+S[(i+1 )]+S[(i+2)]*S[(i+3 )]+S[(i+2)]+(S[(i+3)]+1)*O[i]
    Final_System[j] = Final_System[j]+S[(state_size+j)-(state_size-4)]*S[(state_size+j)-(state_size-1)]*S[(state_size+j)-(state_size-13)]*S[(state_size+j)-(state_size-9 )]+ S[(state_size+j)-(state_size-4)]*S[(state_size+j)-(state_size-1)]*S[(state_size+j)-(state_size-9 )]+ S[(state_size+j)-(state_size-4)]*S[(state_size+j)-(state_size-13)]*S[(state_size+j)-(state_size-9 )]+ S[(state_size+j)-(state_size-4)]*S[(state_size+j)-(state_size-13 )]+ S[(state_size+j)-(state_size-4)]*S[(state_size+j)-(state_size-9 )]+ S[(state_size+j)-(state_size-4 )]+ S[(state_size+j)-(state_size-1)]*S[(state_size+j)-(state_size-13)]*S[(state_size+j)-(state_size-9 )]+ S[(state_size+j)-(state_size-1)]*S[(state_size+j)-(state_size-13 )]+ S[(state_size+j)-(state_size-1)]*S[(state_size+j)-(state_size-9 )]+ S[(state_size+j)-(state_size-1)]+(S[(state_size+j)-(state_size-13)]+S[(state_size+j)-(state_size-1)]*S[(state_size+j)-(state_size-9)]+S[(state_size+j)-(state_size-4)]*S[(state_size+j)-(state_size-9)]+S[(state_size+j)-(state_size-9)]*S[(state_size+j)-(state_size-13)]+1)*O[j]
    S[(j+state_size)] = (S[(j+state_size)-15]+S[(j+state_size)-8])

# print(Final_System)
    
U = [0]*no_monomials
for i in range(no_monomials):
    U[i] = eval("x" + str(state_size+i))

monomials =[1]*1
for i in range(1,len(U)):
    C = Combinations(s,i)
    no = C.cardinality()
    monomials_temp = [1]*no
    for i in range(no):
        for x in C[i]:
            monomials_temp[i]=monomials_temp[i]*x
    monomials = monomials+monomials_temp

#initialise the dictionary
d = {}
#fill the dictionary with the substituion variables
for i in range(1,len(U)+1):
    d[monomials[i]] = U[i-1]

Final_System_SubWithoutY = [0]*len(Final_System)
Final_System_SubWithY = [0]*len(Final_System)
Final_System_Sub = [0]*len(Final_System)

#For each equation in the system
for i in range(len(Final_System)):

    #NOTE:if you sub for the output bits first, you then only need to do the first line of code as there will no longer be output variables
    #make substitution for monomials without output bits
    Final_System_Sub[i] = sum(Final_System[i].monomial_coefficient(m) * v for m,v in d.items())

F = GF(2)
M = matrix(F,len(Final_System_Sub),len(U))
#Fill coefficient matrix with 1s in posisitons corresponding to initial
# state bits and odd keystream bits
for j in range(len(Final_System_Sub)):
    for i in range(len(U)):
        M[j,i] = Final_System_Sub[j].monomial_coefficient(U[i])
w= vector([0]*len(Final_System_Sub))
# w = [0]*len(Final_System)

# M = transpose(M)

for i in range(len(Final_System)):
#     if (EQ_1[i].monomial_coefficient(1)==1):
    if (Final_System[i].monomial_coefficient(1)==1):      
        w[i]=1
# print(w.parent())
print(M.solve_right(w))
A = M.right_kernel()
# # print(A.dimension())
print(A)
for i in range(A.dimension()):
    print(A[i])
stop = timeit.timeit()

print(stop-start)