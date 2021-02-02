import scipy.special
from itertools import combinations 
import timeit

start = timeit.timeit()
#how many clocks are required
no_clocks = 275585000

#the size of the register
state_size = 128

s = [0,0,1,0,1,1,0,1,0,0,1,1,1,0,0,0,0,1,0,1,1,0,1,0,0,1,1,1,0,0,0,0,1,0,1,1,0,1,0,0,1,1,1,0,0,0,0,1,0,1,1,0,1,0,0,1,1,1,0,0,0,0,1,0,1,1,0,1,0,0,1,1,1,0,0,1,1,1,1,0,0,0,1,0,1,1,0,1,0,0,0,0,1,0,1,1,0,1,0,0,0,1,1,0,1,0,0,1,1,1,1,0,0,0,1,0,1,1,0,1,1,0,0,1,1,1,0,0]
len(s)
add=[0]*no_clocks
S = s+add
G = s+add



Left = [0]*no_clocks
for i in range(len(s),no_clocks):
    Left[i-state_size] = mod(S[i-(state_size-2)]*S[i-(state_size-8)]*S[i-(state_size-42)]*S[i-(state_size-95 )]+ S[i-(state_size-15)]*S[i-(state_size-8)]*S[i-(state_size-42)]*S[i-(state_size-95 )]+S[i-(state_size-36)]*S[i-(state_size-8)]*S[i-(state_size-42)]*S[i-(state_size-95 )]+ S[i-(state_size-45)]*S[i-(state_size-8)]*S[i-(state_size-42)]*S[i-(state_size-95 )]+ S[i-(state_size-64)]*S[i-(state_size-8)]*S[i-(state_size-42)]*S[i-(state_size-95 )]+ S[i-(state_size-73)]*S[i-(state_size-8)]*S[i-(state_size-42)]*S[i-(state_size-95 )]+ S[i-(state_size-89)]*S[i-(state_size-8)]*S[i-(state_size-42)]*S[i-(state_size-95 )]+ S[i-(state_size-93)]*S[i-(state_size-8)]*S[i-(state_size-42)]*S[i-(state_size-95 )]+ S[i-(state_size-8)]*S[i-(state_size-13)]*S[i-(state_size-20)]*S[i-(state_size-42)]*S[i-(state_size-95 )]+ S[i-(state_size-8)]*S[i-(state_size-13)]*S[i-(state_size-20)]*S[i-(state_size-42 )]+ S[i-(state_size-8)]*S[i-(state_size-13)]*S[i-(state_size-20)]*S[i-(state_size-95 )]+  S[i-(state_size-8)]*S[i-(state_size-42)]*S[i-(state_size-60)]*S[i-(state_size-79)]*S[i-(state_size-95 )]+ S[i-(state_size-8)]*S[i-(state_size-42)]*S[i-(state_size-60)]*S[i-(state_size-79 )]+ S[i-(state_size-8)]*S[i-(state_size-60)]*S[i-(state_size-79)]*S[i-(state_size-95 )]+ S[i-(state_size-13)]*S[i-(state_size-20)]*S[i-(state_size-42)]*S[i-(state_size-95 )]+ S[i-(state_size-42)]*S[i-(state_size-60)]*S[i-(state_size-79)]*S[i-(state_size-95)],2)
    S[i] = mod((S[i-32]+S[i-47]+S[i-58]+S[i-90]+S[i-121]+S[i-128]),2)
    
print(berlekamp_massey(Left))