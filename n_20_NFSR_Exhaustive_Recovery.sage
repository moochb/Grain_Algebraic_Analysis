import scipy.special
from itertools import combinations 
import time
import random
import itertools



#define the paramaters

#how many clocks are required
no_clocks = 15000

#the size of the register
state_size = 20

#how many output bits are needed
no_output_bits=15000

#degree of the system of equations to be solved
degree = 2

#the number of monomials of degree 1->max degree
#the variables 'no_monomials' will be the number of linearisation variables
#'scipy.special.comb' provides the binomial coefficient

s_out = [0]*state_size
b_out = [0]*state_size

for i in range(state_size):
    s_out[i] = random.randint(0,1)
    b_out[i] = random.randint(0,1)
    
print(s_out)
print(b_out)

add=[0]*no_clocks
S_out = s_out+add
B_out = b_out+add



O = [0]*(no_clocks-state_size)
for i in range(state_size,no_clocks):
    O[i-state_size] =mod(S_out[i-(state_size-4)]*S_out[i-(state_size-1)]*S_out[i-(state_size-9 )]+ S_out[i-(state_size-4)]*S_out[i-(state_size-9)]*B_out[i-(state_size-12 )]+ S_out[i-(state_size-4 )]+ S_out[i-(state_size-1)]*S_out[i-(state_size-13)]*S_out[i-(state_size-9 )]+ S_out[i-(state_size-1)]*S_out[i-(state_size-9)]*B_out[i-(state_size-12 )]+ S_out[i-(state_size-1 )]+ S_out[i-(state_size-13)]*S_out[i-(state_size-9)]*B_out[i-(state_size-12 )]+ S_out[i-(state_size-13)]*S_out[i-(state_size-9 )]+ S_out[i-(state_size-13)]*B_out[i-(state_size-12 )]+ S_out[i-(state_size-13 )]+ B_out[i-(state_size-12)],2)
    B_out[i] = mod(S_out[i-(state_size-0 )]+ B_out[i-(state_size-0 )]+ B_out[i-(state_size-13 )]+ B_out[i-(state_size-19 )]+ B_out[i-(state_size-15 )]+ B_out[i-(state_size-2)]*B_out[i-(state_size-15 )]+ B_out[i-(state_size-3)]*B_out[i-(state_size-5 )]+ B_out[i-(state_size-7)]*B_out[i-(state_size-8 )]+ B_out[i-(state_size-14)]*B_out[i-(state_size-19 )]+ B_out[i-(state_size-13)]*B_out[i-(state_size-6)]*B_out[i-(state_size-17)]*B_out[i-(state_size-18 )]+ B_out[i-(state_size-10)]*B_out[i-(state_size-11)]*B_out[i-(state_size-12 )],2)
    S_out[i] = mod((S_out[i-20]+S_out[i-9]+S_out[i-5]+S_out[i-3]),2)

#number of stages recoverd by partial recovery
m=8

start = time.perf_counter()
#set up set up candidates for rest of state
partial_recovered = list(itertools.product([0, 1], repeat=state_size-m))

for j in range(len(partial_recovered)):
    
    #initialise candidate state and output sequence
    test_state = [0]*state_size
    O_test= [0]*(no_clocks-state_size)
    partial_recovered[j]=list(partial_recovered[j])
    test_state = b_out[0:m]+partial_recovered[j]
    TEST_STATE = test_state+add
#     if test_state==b_out:
#         print('yes')
#         print(test_state)
    
    #update and produce output
    for i in range(state_size,no_clocks):
        O_test[i-state_size] =mod(S_out[i-(state_size-4)]*S_out[i-(state_size-1)]*S_out[i-(state_size-9 )]+ S_out[i-(state_size-4)]*S_out[i-(state_size-9)]*TEST_STATE[i-(state_size-12 )]+ S_out[i-(state_size-4 )]+ S_out[i-(state_size-1)]*S_out[i-(state_size-13)]*S_out[i-(state_size-9 )]+ S_out[i-(state_size-1)]*S_out[i-(state_size-9)]*TEST_STATE[i-(state_size-12 )]+ S_out[i-(state_size-1 )]+ S_out[i-(state_size-13)]*S_out[i-(state_size-9)]*TEST_STATE[i-(state_size-12 )]+ S_out[i-(state_size-13)]*S_out[i-(state_size-9 )]+ S_out[i-(state_size-13)]*TEST_STATE[i-(state_size-12 )]+ S_out[i-(state_size-13 )]+ TEST_STATE[i-(state_size-12)],2)    
        TEST_STATE[i] = mod(S_out[i-(state_size-0 )]+ TEST_STATE[i-(state_size-0 )]+ TEST_STATE[i-(state_size-13 )]+ TEST_STATE[i-(state_size-19 )]+ TEST_STATE[i-(state_size-15 )]+ TEST_STATE[i-(state_size-2)]*TEST_STATE[i-(state_size-15 )]+ TEST_STATE[i-(state_size-3)]*TEST_STATE[i-(state_size-5 )]+ TEST_STATE[i-(state_size-7)]*TEST_STATE[i-(state_size-8 )]+ TEST_STATE[i-(state_size-14)]*TEST_STATE[i-(state_size-19 )]+ TEST_STATE[i-(state_size-13)]*TEST_STATE[i-(state_size-6)]*TEST_STATE[i-(state_size-17)]*TEST_STATE[i-(state_size-18 )]+ TEST_STATE[i-(state_size-10)]*TEST_STATE[i-(state_size-11)]*TEST_STATE[i-(state_size-12 )],2)    
    if O==O_test:
        print("yes",test_state)
    
stop = time.perf_counter()
print("online",stop-start)