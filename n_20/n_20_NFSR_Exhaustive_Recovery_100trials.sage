#This code exhuastively searchs a set of candidates for the NFSR in a Grain-like structure.
#This code assumes that the LFSR has been recovered and that the NFSR state has been partially recovered
#You need to specify:
    #how many bits were partially recovered (the value 'm')
    #how many trials are going to be done
    #the amount of keystream to be used to check
    #the state size of the register
    #the update and output functions for the entire Grain-like structure
#Here, we use the registers and the functions used in the paper "Algebraic Attacks on Grain-like Keystream Generators" in the simulations section.


import scipy.special
from itertools import combinations 
import time
import random
import itertools

#initialise a total time counter
total = 0

#number of stages recoverd by partial recovery
m=8

#define number of trials to perform
no_trials = 100

#For the defined number of trials we are going to:
    #randomly define an initial state and produce output
    #produce all possible NFSR candidates (following partial recovery of $m$ bits)
    #produce output for each candidate (using correct LFSR state)
    #check to see if output of candidate matches correct output
#NOTE: we only test against 100 bits of output. Assuming the output function is 1-1
#after 20 bits, every output sequence should be unique
for l in range(no_trials):

    #how many clocks are required
    no_clocks = 100

    #the size of the register
    state_size = 20

    #randomly generate initial states for the LFSR and the NFSR
    #print initial states
    #produce output
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


    #starting timing for each trial
    start = time.perf_counter()

    #for each possible candidate (candidates will be of length state_size-m)
     for j in range(2^(state_size-m)):

        #initialise candidate state and output sequence
        test_state = [0]*state_size
        O_test= [0]*(no_clocks-state_size)

        #produce binary candidate
        a = '{:012b}'.format(j)
        partial_recovered = [int(i) for i in a]
        test_state = b_out[0:m]+partial_recovered

        #initialise state sequence for candidate
        TEST_STATE = test_state+add

        #Produce output and update state
        for i in range(state_size,no_clocks):
            O_test[i-state_size] =mod(S_out[i-(state_size-4)]*S_out[i-(state_size-1)]*S_out[i-(state_size-9 )]+ S_out[i-(state_size-4)]*S_out[i-(state_size-9)]*TEST_STATE[i-(state_size-12 )]+ S_out[i-(state_size-4 )]+ S_out[i-(state_size-1)]*S_out[i-(state_size-13)]*S_out[i-(state_size-9 )]+ S_out[i-(state_size-1)]*S_out[i-(state_size-9)]*TEST_STATE[i-(state_size-12 )]+ S_out[i-(state_size-1 )]+ S_out[i-(state_size-13)]*S_out[i-(state_size-9)]*TEST_STATE[i-(state_size-12 )]+ S_out[i-(state_size-13)]*S_out[i-(state_size-9 )]+ S_out[i-(state_size-13)]*TEST_STATE[i-(state_size-12 )]+ S_out[i-(state_size-13 )]+ TEST_STATE[i-(state_size-12)],2)    
            TEST_STATE[i] = mod(S_out[i-(state_size-0 )]+ TEST_STATE[i-(state_size-0 )]+ TEST_STATE[i-(state_size-13 )]+ TEST_STATE[i-(state_size-19 )]+ TEST_STATE[i-(state_size-15 )]+ TEST_STATE[i-(state_size-2)]*TEST_STATE[i-(state_size-15 )]+ TEST_STATE[i-(state_size-3)]*TEST_STATE[i-(state_size-5 )]+ TEST_STATE[i-(state_size-7)]*TEST_STATE[i-(state_size-8 )]+ TEST_STATE[i-(state_size-14)]*TEST_STATE[i-(state_size-19 )]+ TEST_STATE[i-(state_size-13)]*TEST_STATE[i-(state_size-6)]*TEST_STATE[i-(state_size-17)]*TEST_STATE[i-(state_size-18 )]+ TEST_STATE[i-(state_size-10)]*TEST_STATE[i-(state_size-11)]*TEST_STATE[i-(state_size-12 )],2)    
        
        #test whether output sequence produced by the candidate matches the correct one
        if O[0:no_clocks-state_size]==O_test[0:no_clocks-state_size]:
            #if it does calculate the time taken to find the solution
            stop = time.perf_counter()
            #add to running total and move on to next trial
            total = total+(stop-start)
            break

#once all trials are finished, print the average time taken to find correct NFSR initial state
print(total/l)

