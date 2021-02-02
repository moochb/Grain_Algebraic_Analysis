import scipy.special
from itertools import combinations 
import timeit

start = timeit.timeit()
#how many clocks are required
no_clocks = 1000

#the size of the register
state_size = 15

#how many output bits are needed
no_output_bits=1000

#degree of the system of equations to be solved
degree = 4
s = [0,0,1,0,1,1,0,1,0,0,1,1,1,0,0]
add=[0]*no_clocks
S = s+add
G = s+add



LFSR_output = [0]*no_clocks
for i in range(len(s),no_clocks):
    LFSR_output[i-state_size] = (S[i-5]+S[i-14]*S[i-9]+S[i-12]*S[i-9]+S[i-9]*S[i-5])*(S[i-12]+S[i-14]*S[i-5]+S[i-9]*S[i-5]+S[i-14]*S[i-12]*S[i-9]+S[i-14]*S[i-9]*S[i-5])
    S[i] = (S[i-15]+S[i-8])
    
print(berlekamp_massey(LFSR_output))