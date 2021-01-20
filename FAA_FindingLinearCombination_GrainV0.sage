import scipy.special
from itertools import combinations 
import timeit

start = timeit.timeit()

#define the paramaters

#how many clocks are required
no_clocks = 86000

#the size of the register
state_size = 80

#how many output bits are needed
no_output_bits=86000

#degree of the system of equations to be solved
degree = 3

#the number of monomials of degree 1->max degree
#the variables 'no_monomials' will be the number of linearisation variables
#'scipy.special.comb' provides the binomial coefficient
no_monomials = 0
for i in range(1,degree+1):
    no_monomials = no_monomials + binomial(state_size,i)

#define polynomial ring in which the the INITIAL STATE, OUTPUT and LINEARISATION variables will be defined
R=BooleanPolynomialRing(state_size+no_output_bits+no_monomials,'x')
#Inject these variables for use
R.inject_variables()

#define initial state
s = [0]*state_size
for i in range(state_size):
    s[i] = eval("x" + str(i))

# print(s)
#define output
Y = [0]*no_output_bits
for i in range(no_output_bits):
    Y[i] = eval("x" + str(state_size+i))

# print(Y)
#We can consider 's' as a register that will be shifted for a specified number of times. Thus, there will exist more entries than simply those defined in 's'.
#We therefore initialise a vector 'add' that will be appended to s to hold all the update bits (one vector for equations with respect to keystream bits S and one with respect the the linear approximation$
add=[0]*no_clocks
S = s+add
Gz = s+add
G = s+add

#Generate equations for update using output bits
# S[15] = S[0]+S[7]
for i in range(len(s),no_clocks):
     S[i] = (S[i-80]+S[i-67]+S[i-57]+S[i-42]+S[i-29]+S[i-18])
     
#Generate equations for update using linear approximation
for i in range(len(s),no_clocks):
    #s1s2z=s1s2s9
    #Here we are going to build two systems; one containing z, the other not
    # Gz[i] = (S[i-5]+S[i-14]*S[i-9]+S[i-12]*S[i-9]+S[i-9]*S[i-5])*Y[i-15]
    #G[i]= (S[i-5]+S[i-14]*S[i-9]+S[i-12]*S[i-9]+S[i-9]*S[i-5])*(S[i-12]+S[i-14]*S[i-5]+S[i-9]*S[i-5]+S[i-14]*S[i-12]*S[i-9]+S[i-14]*S[i-9]*S[i-5])
    G[i] = (S[i-16]+S[i-77]*S[i-34]+S[i-55]*S[i-34]+S[i-34]*S[i-16])*(S[i-79]+S[i-55]+S[i-77]*S[i-16]+S[i-34]*S[i-16]+S[i-77]*S[i-55]*S[i-34]+S[i-77]*S[i-55]*S[i-16]) 
#for now we are just going to focus on the system without z. This is the system whose degree we want to reduce.

#create a matrix that will hold the final equations (the sum of those given with respect to the keystream bits + those given by the linear approximation)
EQ = list(var('eq%d' %i) for i in (0..no_clocks))

# produce final system
for i in range(len(s),no_clocks):
    EQ[i] = (G[i])

# get rid of redundant equations
EQ_1=EQ[len(s):no_clocks]

############################
#Now linearise

#Define linear substitution variables
U = [0]*no_monomials
for i in range(no_monomials):
    U[i] = eval("x" + str(state_size+no_output_bits+i))
# print(len(U))
# print(s)
#create monomial list and count the number of combinations
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

EQ_SubWithoutY = [0]*len(EQ_1)
EQ_SubWithY = [0]*len(EQ_1)
EQ_Sub = [0]*len(EQ_1)

#For each equation in the system
for i in range(len(EQ_1)):

    #NOTE:if you sub for the output bits first, you then only need to do the first line of code as there will no longer be output variables
    #make substitution for monomials without output bits
    EQ_Sub[i] = sum(EQ_1[i].monomial_coefficient(m) * v for m,v in d.items())

#make substitution for monomials with output bits
#     for j in range(len(W)):
#         EQ_SubWithY[i]= sum(EQ_1[i].monomial_coefficient(m*Y[j]) * v for m,v in d.iteritems())
#         EQ_Sub[i] = EQ_Sub[i]+EQ_SubWithY[i]*Y[j]
#We are interested in the coefficients of the variables
#Initialise coefficient matrix
F = GF(2)
M = matrix(F,len(EQ_Sub),len(U)-state_size)
#Fill coefficient matrix with 1s in posisitons corresponding to initial
# state bits and odd keystream bits


#we are looking to get rid of all higher order monomials (these monomials were linearised to the variables x26 and beyond, so we are looking to eliminate these variables.)
dontgetrid = 0
for i in range(1,3):
    dontgetrid = dontgetrid+binomial(state_size,i)
    #dontgetrid=120

for j in range(len(EQ_Sub)):
    #we are going to run through all of the variables beyond the ones we dont want to get rid of
    for i in range(len(U)-dontgetrid):
        M[j,i] = EQ_Sub[j].monomial_coefficient(U[i+dontgetrid])
w= vector(F,[0]*len(EQ_Sub))

#In order to find which rows are linearly dependent, we transpose and find a basis for the column space
M_1 = transpose(M)
# print(N_1)

# R = M_1.echelon_form()
# R = R%2
# R = R.rref()
# R = R%2
# # print(R[0:81])

# # EQ[80]+EQ[80+23]+EQ[80+38]+EQ[80+51]+EQ[80+62]+EQ[80+67]+EQ[80+80]
# for i in range(l):
#     print(i+1,sum(R[:,i]))
R = M_1.echelon_form()
R = R%2
R = R.rref()
R = R%2
# print(R[0:81])

#find the first column of the coefficient matrix that is linearly dependent
# col = 0
# for i in range(1,R.ncols()):
#     if sum(R[:,i])!=1:
#         col = i
#         print(col)
#         print(R[:,col])
#         break
col = binomial(state_size,degree)
#We are going to print the linear combination to a file for good keeping 
#save the original standard output as a reference       
original_stdout = sys.stdout

#create a file called stuff.txt with writing privledges 
fil = open("stuff.txt","w")

#set standard output to the file stuff.txt
sys.stdout = fil

#print the equation number that is depenedent
print("EQ_[{}]+".format(col),end='')

#for our linearly dependent row, print the combination of equations that produces it
for j in range(R.nrows()):
            if R[j,col]==1:
                print("EQ_[{}]+".format(j), end='')

sys.stdout = original_stdout
