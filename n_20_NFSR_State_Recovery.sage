## This code simulates the recovery of the NFSR in a "Grain-like" structure, assuming both the LFSR and the output is
##known.
## A few things need to be specified:
## - the number of trials to be performed (no_trials)
## - the number of clocks to do each trial (no_clocks). Note that the number of output bits produced is no_clocks+state_size
## 

import scipy.special
from itertools import combinations 
import timeit
import random

#how many trials
no_trials  = 10000

#the size of the register
state_size = 20

#how many clocks (number of output bits = no_clocks + state_size)
no_clocks = 500

#define array for number of bits recovered each trial
count = [0]*no_trials


#define polynomial ring in which the the NFSR STATE BITS will be defined
R=BooleanPolynomialRing(state_size+no_clocks,'b')

#Inject these variables for use
R.inject_variables()

#repeat for the number of trials
for l in range(no_trials):
    
    
#     start = timeit.timeit()

    #define initial states for the NFSR and the LFSR and initalise them "randomly" using the randint function.
    s_out = [0]*state_size
    b_out = [0]*state_size
    for i in range(state_size):
        s_out[i] = random.randint(0,1)
        b_out[i] = random.randint(0,1)

    #define arrays to hold the state sequences for the LFSR and the NFSR
    add=[0]*no_clocks
    S_out = s_out+add
    B_out = b_out+add


    #define array to hold the output sequence
    O = [0]*(no_clocks-state_size)
    
    #NOTE: THIS HAS TO BE DEFINED FOR EACH UNIQUE STRUCTURE
    #THIS IS THE FEEDBACK FUNCTIONS TO THE LFSR/NFSR AND THE OUTPUT FUNCTION.
    #The next block of code will produce the output sequence that will be visable to the attacker.

    #produce output and update each register for the perscribed number of clocks
    for i in range(len(s_out),no_clocks):
        O[i-state_size] =mod((S_out[i-(state_size-13)]+S_out[i-(state_size-1)]*S_out[i-(state_size-9)]+S_out[i-(state_size-4)]*S_out[i-(state_size-9)]+S_out[i-(state_size-9)]*S_out[i-(state_size-13)])*B_out[i-(state_size-4)]+(S_out[i-(state_size-4)]+S_out[i-(state_size-1)]*S_out[i-(state_size-13)]+S_out[i-(state_size-9)]*S_out[i-(state_size-13)]+S_out[i-(state_size-1)]*S_out[i-(state_size-4)]*S_out[i-(state_size-9)]+S_out[i-(state_size-1)]*S_out[i-(state_size-9)]*S_out[i-(state_size-13)]),2)
        B_out[i] = mod(S_out[i-(state_size-0 )]+ B_out[i-(state_size-0 )]+ B_out[i-(state_size-13 )]+ B_out[i-(state_size-19 )]+ B_out[i-(state_size-15 )]+ B_out[i-(state_size-2)]*B_out[i-(state_size-15 )]+ B_out[i-(state_size-3)]*B_out[i-(state_size-5 )]+ B_out[i-(state_size-7)]*B_out[i-(state_size-8 )]+ B_out[i-(state_size-14)]*B_out[i-(state_size-19 )]+ B_out[i-(state_size-13)]*B_out[i-(state_size-6)]*B_out[i-(state_size-17)]*B_out[i-(state_size-18 )]+ B_out[i-(state_size-10)]*B_out[i-(state_size-11)]*B_out[i-(state_size-12 )],2)
        S_out[i] = mod((S_out[i-20]+S_out[i-9]+S_out[i-5]+S_out[i-3]),2)

    #we are assuming that the attacked knows the LFSR contents and the output, but not the NFSR contents
    #this next block of code builds equations that relate the NFSR state bits to the known LFSR state bits and output bits

    #define and initialise an array to hold the symbolic NFSR state bits
    #Note that we are not using the NFSR update here and are simply just assigning a new variable to each state bit
    B = [0]*no_clocks
    for i in range(no_clocks):
        B[i] = eval("b" + str(i))
    
    #define an array for the system of equations 
    System = [0]*no_clocks

    #NOTE: THIS HAS TO BE DEFINED FOR EACH UNIQUE STRUCTURE
    #THIS IS THE SYMBOLIC FEEDBACK FUNCTIONS TO THE LFSR/NFSR.
    for i in range(len(s_out),no_clocks):
        System[i] = (S_out[i-(state_size-13)]+S_out[i-(state_size-1)]*S_out[i-(state_size-9)]+S_out[i-(state_size-4)]*S_out[i-(state_size-9)]+S_out[i-(state_size-9)]*S_out[i-(state_size-13)])*B[i-(state_size-4)]+(S_out[i-(state_size-4)]+S_out[i-(state_size-1)]*S_out[i-(state_size-13)]+S_out[i-(state_size-9)]*S_out[i-(state_size-13)]+S_out[i-(state_size-1)]*S_out[i-(state_size-4)]*S_out[i-(state_size-9)]+S_out[i-(state_size-1)]*S_out[i-(state_size-9)]*S_out[i-(state_size-13)])+ O[i-state_size] 
#         print(System[i])
    
    #define array to count number of NFSR recovered for CURRENT trial
    count_this_trial = [0]*len(System) 
    
    #We are going to run through sliding windows of length 128 to see which window recovers the most state bits.
    #all equations with a constant monomial of 1 are recoverable
    #if the degree of such an eqution is 2, we recover 2 bits, otherwise, we recover 1.
    for j in range(len(System)-40):
        for i in range(20+j,40+j):
            if System[i].monomial_coefficient(1)==1:
                if(System[i].degree())==2:
                    count_this_trial[j] = count_this_trial[j]+2
                else:
                    count_this_trial[j] = count_this_trial[j]+1
    
    #store in count the maximum number of bits recovered across all of the tested windows.
    count[l]=count[l]+max(count_this_trial)
    
# print(int(sum(count)/len(count)))

print(count)

#build an array that holds all of the possible bits recoverd and the count for how many were recovered (a histogram)
count1 = [0]*state_size
for i in range(len(count)):
    for j in range(state_size):
        if(count[i]==j):
            count1[j] = count1[j]+1

for i in range(state_size):
    if(count1[i]):
        print(i,count1[i])