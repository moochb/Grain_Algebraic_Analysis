A = matrix([0,0,0,0])
dec = 0
for i in range(A.ncols()):
    dec = dec+A[:,i]*2^(A.ncols()-1-i)
dec