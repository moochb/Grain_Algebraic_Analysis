import scipy.special
from itertools import combinations 
import timeit

start = timeit.timeit()
#how many clocks are required
no_clocks = 102100

#the size of the register
state_size = 40

s = [0,0,1,0,1,1,0,1,0,0,1,1,1,0,0,0,0,1,0,1,1,0,1,0,0,1,1,1,0,0,0,0,1,0,1,1,0,1,0,0,1,1,1,0,0,0,0,1,0,1,1,0,1,0]
len(s)
# add=[0]*no_clocks
# S = s+add
# G = s+add



# Left = [0]*no_clocks
# for i in range(len(s),no_clocks):
#     Left[i-state_size] = mod(S[i-(state_size-3)]*S[i-(state_size-25)]*S[i-(state_size-46)]+S[i-(state_size-1)]*S[i-(state_size-3)]*S[i-(state_size-46)]+S[i-(state_size-1)]*S[i-(state_size-25)]*S[i-(state_size-46)]+S[i-(state_size-1)]*S[i-(state_size-46)]*S[i-(state_size-64)],2)
#     S[i] = mod((S[i-80]+S[i-67]+S[i-57]+S[i-42]+S[i-29]+S[i-18]),2)
    
# print(berlekamp_massey(LFSR_output))

#right = mod(S[i-(state_size-64)]+S[i-(state_size-3)]*S[i-(state_size-46)]+S[i-(state_size-25)]*S[i-(state_size-46)]+S[i-(state_size-46)]*S[i-(state_size-64)]+S[i-(state_size-1)]*S[i-(state_size-64)]+S[i-(state_size-3)]*S[i-(state_size-46)]+S[i-(state_size-25)]*S[i-(state_size-64)],2)