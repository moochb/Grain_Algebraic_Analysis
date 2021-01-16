n = 17
#define polynomial ring in which the initial state variables and the keystream bits exist
R.<s0,s1,s2,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,u0,u1,u2,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,u13,u14,u15,u16,u17,u18,u19,u20> = BooleanPolynomialRing(45)
#create a matrix 's' to hold the initial state variables
s =[s0,s1,s2]
s_0_length=len(s)

#create a matrix to hold all the output variables
Y = [y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20]


#We can consider 's' as a register that will be shifted for a specified number of times. Thus, there will exist more entries than simply those defined in 's'.
#We therefore initialise a vector 'add' that will be appended to s to hold all the update bits (one vector for equations with respect to keystream bits S and one with respect the the linear approximation$
add=list(0for i in (1..1000))
S = s+add
G = s+add

#Generate equations for update using output bits
# S[15] = S[0]+S[7]
for i in range(len(s),n):
     S[i] = (S[i-3]+S[i-2])


#Generate equations for update using linear approximation
for i in range(len(s),n):
    #(s0s1s2+1)z0+(s0s1s2+1)s9=0
    G[i] = (S[i-1]+1)*Y[i-3]+S[i-2]*(S[i-1]+1)

    # print(i,G[i])
    #create a matrix that will hold the final equations (the sum of those given with respect to the keystream bits + those given by the linear approximation)
EQ = list(var('eq%d' %i) for i in (0..1000))


# Combine both forms of equations
# print(S)
# print(G)
for i in range(len(s),n):
    EQ[i] = (G[i])

EQ_1=EQ[len(s):n]

######################################################
#Now substitute in output bits
# S0 = [0,0,1,0,1,1,0,1,0,0,1,1,1,0,0]
# B0 = [1,1,0,0,1,0,1,1,1,0,0,0,1,0,1]

O=[1,0,1,1,0,1,0,1,0,0,1,1,0,0,1,0,1,1]
for i in range(len(EQ_1)):
    EQ_1[i]=EQ_1[i].subs(y0=O[0],y1=O[1],y2=O[2],y3=O[3],y4=O[4],y5=O[5],y6=O[6],y7=O[7],y8=O[8],y9=O[9],y10=O[10])

############################
#Now linearise
# EQ_1
#Define linear substitution variables
U = [u0,u1,u2,u3,u4,u5]

#create monomial list and count the number of combinations
C = Combinations(s)
no = C.cardinality()

#Create each of the monomials
#NOTE: The upper bound in the value of monomials is sum_{k = 0 to d}(n choose k) (number of monomials of degree up to and including d in the n variables)
#There will be len(U)+1 monominals because they include the monomial '1' here
monomials =[1]*(len(U)+1)
for i in range(len(U)+1):
    for x in C[i]:
        monomials[i]=monomials[i]*x
monomials
#initialise the dictionary
d = {}
#fill the dictionary with the substituion variables
for i in range(1,len(U)+1):
    d[monomials[i]] = U[i-1]
d
EQ_SubWithoutY = [0]*len(EQ_1)
EQ_SubWithY = [0]*len(EQ_1)
EQ_Sub = [0]*len(EQ_1)
#For each equation in the system
for i in range(len(EQ_1)):

    #NOTE:if you sub for the output bits first, you then only need to do the first line of code as there will no longer be output variables
    #make substitution for monomials without output bits
    EQ_Sub[i] = sum(EQ_1[i].monomial_coefficient(m) * v for m,v in d.items())

#     #make substitution for monomials with output bits
#     for j in range(len(W)):
#         EQ_SubWithY[i]= sum(EQ_1[i].monomial_coefficient(m*Y[j]) * v for m,v in d.iteritems())
#         EQ_Sub[i] = EQ_Sub[i]+EQ_SubWithY[i]*Y[j]
#We are interested in the coefficients of the variables
#Initialise coefficient matrix
F = GF(2)
M = matrix(F,len(EQ_Sub),len(U))
#Fill coefficient matrix with 1s in posisitons corresponding to initial
# state bits and odd keystream bits
for j in range(len(EQ_Sub)):
    for i in range(len(U)):
        M[j,i] = EQ_Sub[j].monomial_coefficient(U[i])
w= vector(F,[0]*len(EQ_Sub))

for i in range(len(EQ_Sub)):
    if (EQ_1[i].monomial_coefficient(1)==1):
        w[i]=1
# M
print(M.solve_right(w))
# print(w)
print(n)