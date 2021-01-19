import scipy.special
from itertools import combinations 
import timeit

start = timeit.timeit()

#define the paramaters

#how many clocks are required
no_clocks = 2000

#the size of the register
state_size = 15

#how many output bits are needed
no_output_bits=2000

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
G = s+add

#Generate equations for update using output bits
# S[15] = S[0]+S[7]
for i in range(len(s),no_clocks):
     S[i] = (S[i-15]+S[i-8])
     
#Generate equations for update using linear approximation
for i in range(len(s),no_clocks):
    #s1s2z=s1s2s9
    #Here we are going to build two systems; one containing z, the other not
    G[i] = (S[i-5]+S[i-14]*S[i-9]+S[i-12]*S[i-9]+S[i-9]*S[i-5])*Y[i-15]+(S[i-5]+S[i-14]*S[i-9]+S[i-12]*S[i-9]+S[i-9]*S[i-5])*(S[i-12]+S[i-14]*S[i-5]+S[i-9]*S[i-5]+S[i-14]*S[i-12]*S[i-9]+S[i-14]*S[i-9]*S[i-5])

#for now we are just going to focus on the system without z. This is the system whose degree we want to reduce.

#create a matrix that will hold the final equations (the sum of those given with respect to the keystream bits + those given by the linear approximation)
EQ = list(var('eq%d' %i) for i in (0..no_clocks))

# produce final system
for i in range(len(s),no_clocks):
    EQ[i] = (G[i])

# get rid of redundant equations
EQ_=EQ[len(s):no_clocks]

EQ_1 = EQ_

#From the code in "FAA_FindingLinearCombination_Generic_LFSR_NFSR_Combiner_AlgebraicAttackn=15_Homeogeneous"
#we see that the linear combination that will remove all higher order terms beyond 2 is
#     EQ_[0]+EQ_[5]+EQ_[7]+EQ_[17]+EQ_[22]+EQ_[24]+EQ_[28]+EQ_[29]+EQ_[30]+EQ_[31]+EQ_[33]+EQ_[34]+EQ_[36]+EQ_[37]+EQ_[38]+EQ_[39]+EQ_[41]+EQ_[42]+EQ_[43]+EQ_[52]+EQ_[54]+EQ_[56]+EQ_[57]+EQ_[58]+EQ_[63]+EQ_[65]+EQ_[66]+EQ_[68]+EQ_[69]+EQ_[71]+EQ_[72]+EQ_[74]+EQ_[78]+EQ_[79]+EQ_[81]+EQ_[82]+EQ_[85]+EQ_[86]+EQ_[87]+EQ_[88]+EQ_[89]+EQ_[90]+EQ_[93]+EQ_[94]+EQ_[97]+EQ_[98]+EQ_[99]+EQ_[100]+EQ_[106]+EQ_[107]+EQ_[110]+EQ_[112]+EQ_[115]+EQ_[116]+EQ_[117]+EQ_[118]+EQ_[123]+EQ_[125]+EQ_[127]+EQ_[131]+EQ_[132]+EQ_[138]+EQ_[139]+EQ_[141]+EQ_[142]+EQ_[143]+EQ_[144]+EQ_[145]+EQ_[147]+EQ_[150]+EQ_[151]+EQ_[154]+EQ_[160]+EQ_[162]+EQ_[168]+EQ_[169]+EQ_[171]+EQ_[172]+EQ_[173]+EQ_[175]+EQ_[177]+EQ_[182]+EQ_[183]+EQ_[187]+EQ_[190]+EQ_[192]+EQ_[197]+EQ_[199]+EQ_[201]+EQ_[203]+EQ_[204]+EQ_[206]+EQ_[207]+EQ_[208]+EQ_[209]+EQ_[212]+EQ_[214]+EQ_[215]+EQ_[219]+EQ_[220]+EQ_[221]+EQ_[222]+EQ_[224]+EQ_[225]+EQ_[226]+EQ_[228]+EQ_[231]+EQ_[233]+EQ_[234]+EQ_[237]+EQ_[238]+EQ_[241]+EQ_[242]+EQ_[244]+EQ_[246]+EQ_[248]+EQ_[249]+EQ_[250]+EQ_[254]+EQ_[255]+EQ_[259]+EQ_[260]+EQ_[261]+EQ_[263]+EQ_[264]+EQ_[265]+EQ_[266]+EQ_[268]+EQ_[271]+EQ_[273]+EQ_[274]+EQ_[277]+EQ_[278]+EQ_[279]+EQ_[286]+EQ_[288]+EQ_[289]+EQ_[293]+EQ_[295]+EQ_[296]+EQ_[297]+EQ_[298]+EQ_[300]+EQ_[301]+EQ_[302]+EQ_[306]+EQ_[308]+EQ_[313]+EQ_[314]+EQ_[319]+EQ_[320]+EQ_[322]+EQ_[323]+EQ_[324]+EQ_[325]+EQ_[326]+EQ_[328]+EQ_[329]+EQ_[333]+EQ_[335]+EQ_[337]+EQ_[338]+EQ_[339]+EQ_[340]+EQ_[341]+EQ_[342]+EQ_[343]+EQ_[344]+EQ_[345]+EQ_[349]+EQ_[351]+EQ_[355]+EQ_[356]+EQ_[357]+EQ_[358]+EQ_[359]+EQ_[360]+EQ_[362]+EQ_[363]+EQ_[364]+EQ_[366]+EQ_[367]+EQ_[369]+EQ_[371]+EQ_[373]+EQ_[374]+EQ_[375]+EQ_[378]+EQ_[379]+EQ_[381]+EQ_[382]+EQ_[384]+EQ_[386]+EQ_[387]+EQ_[388]+EQ_[389]+EQ_[390]+EQ_[393]+EQ_[394]+EQ_[396]+EQ_[398]+EQ_[400]+EQ_[401]+EQ_[405]+EQ_[406]+EQ_[407]+EQ_[412]+EQ_[414]+EQ_[415]+EQ_[417]+EQ_[418]+EQ_[420]+EQ_[423]+EQ_[425]+EQ_[429]+EQ_[431]+EQ_[432]+EQ_[434]+EQ_[437]+EQ_[439]+EQ_[447]+EQ_[450]+EQ_[455]
for i in range(len(EQ_)-456):
    EQ_1[i] = EQ_[0+i]+EQ_[5+i]+EQ_[7+i]+EQ_[17+i]+EQ_[22+i]+EQ_[24+i]+EQ_[28+i]+EQ_[29+i]+EQ_[30+i]+EQ_[31+i]+EQ_[33+i]+EQ_[34+i]+EQ_[36+i]+EQ_[37+i]+EQ_[38+i]+EQ_[39+i]+EQ_[41+i]+EQ_[42+i]+EQ_[43+i]+EQ_[52+i]+EQ_[54+i]+EQ_[56+i]+EQ_[57+i]+EQ_[58+i]+EQ_[63+i]+EQ_[65+i]+EQ_[66+i]+EQ_[68+i]+EQ_[69+i]+EQ_[71+i]+EQ_[72+i]+EQ_[74+i]+EQ_[78+i]+EQ_[79+i]+EQ_[81+i]+EQ_[82+i]+EQ_[85+i]+EQ_[86+i]+EQ_[87+i]+EQ_[88+i]+EQ_[89+i]+EQ_[90+i]+EQ_[93+i]+EQ_[94+i]+EQ_[97+i]+EQ_[98+i]+EQ_[99+i]+EQ_[100+i]+EQ_[106+i]+EQ_[107+i]+EQ_[110+i]+EQ_[112+i]+EQ_[115+i]+EQ_[116+i]+EQ_[117+i]+EQ_[118+i]+EQ_[123+i]+EQ_[125+i]+EQ_[127+i]+EQ_[131+i]+EQ_[132+i]+EQ_[138+i]+EQ_[139+i]+EQ_[141+i]+EQ_[142+i]+EQ_[143+i]+EQ_[144+i]+EQ_[145+i]+EQ_[147+i]+EQ_[150+i]+EQ_[151+i]+EQ_[154+i]+EQ_[160+i]+EQ_[162+i]+EQ_[168+i]+EQ_[169+i]+EQ_[171+i]+EQ_[172+i]+EQ_[173+i]+EQ_[175+i]+EQ_[177+i]+EQ_[182+i]+EQ_[183+i]+EQ_[187+i]+EQ_[190+i]+EQ_[192+i]+EQ_[197+i]+EQ_[199+i]+EQ_[201+i]+EQ_[203+i]+EQ_[204+i]+EQ_[206+i]+EQ_[207+i]+EQ_[208+i]+EQ_[209+i]+EQ_[212+i]+EQ_[214+i]+EQ_[215+i]+EQ_[219+i]+EQ_[220+i]+EQ_[221+i]+EQ_[222+i]+EQ_[224+i]+EQ_[225+i]+EQ_[226+i]+EQ_[228+i]+EQ_[231+i]+EQ_[233+i]+EQ_[234+i]+EQ_[237+i]+EQ_[238+i]+EQ_[241+i]+EQ_[242+i]+EQ_[244+i]+EQ_[246+i]+EQ_[248+i]+EQ_[249+i]+EQ_[250+i]+EQ_[254+i]+EQ_[255+i]+EQ_[259+i]+EQ_[260+i]+EQ_[261+i]+EQ_[263+i]+EQ_[264+i]+EQ_[265+i]+EQ_[266+i]+EQ_[268+i]+EQ_[271+i]+EQ_[273+i]+EQ_[274+i]+EQ_[277+i]+EQ_[278+i]+EQ_[279+i]+EQ_[286+i]+EQ_[288+i]+EQ_[289+i]+EQ_[293+i]+EQ_[295+i]+EQ_[296+i]+EQ_[297+i]+EQ_[298+i]+EQ_[300+i]+EQ_[301+i]+EQ_[302+i]+EQ_[306+i]+EQ_[308+i]+EQ_[313+i]+EQ_[314+i]+EQ_[319+i]+EQ_[320+i]+EQ_[322+i]+EQ_[323+i]+EQ_[324+i]+EQ_[325+i]+EQ_[326+i]+EQ_[328+i]+EQ_[329+i]+EQ_[333+i]+EQ_[335+i]+EQ_[337+i]+EQ_[338+i]+EQ_[339+i]+EQ_[340+i]+EQ_[341+i]+EQ_[342+i]+EQ_[343+i]+EQ_[344+i]+EQ_[345+i]+EQ_[349+i]+EQ_[351+i]+EQ_[355+i]+EQ_[356+i]+EQ_[357+i]+EQ_[358+i]+EQ_[359+i]+EQ_[360+i]+EQ_[362+i]+EQ_[363+i]+EQ_[364+i]+EQ_[366+i]+EQ_[367+i]+EQ_[369+i]+EQ_[371+i]+EQ_[373+i]+EQ_[374+i]+EQ_[375+i]+EQ_[378+i]+EQ_[379+i]+EQ_[381+i]+EQ_[382+i]+EQ_[384+i]+EQ_[386+i]+EQ_[387+i]+EQ_[388+i]+EQ_[389+i]+EQ_[390+i]+EQ_[393+i]+EQ_[394+i]+EQ_[396+i]+EQ_[398+i]+EQ_[400+i]+EQ_[401+i]+EQ_[405+i]+EQ_[406+i]+EQ_[407+i]+EQ_[412+i]+EQ_[414+i]+EQ_[415+i]+EQ_[417+i]+EQ_[418+i]+EQ_[420+i]+EQ_[423+i]+EQ_[425+i]+EQ_[429+i]+EQ_[431+i]+EQ_[432+i]+EQ_[434+i]+EQ_[437+i]+EQ_[439+i]+EQ_[447+i]+EQ_[450+i]+EQ_[455+i]
# EQ_[0]+EQ_[5]+EQ_[7]+EQ_[17]+EQ_[22]+EQ_[24]+EQ_[28]+EQ_[29]+EQ_[30]+EQ_[31]+EQ_[33]+EQ_[34]+EQ_[36]+EQ_[37]+EQ_[38]+EQ_[39]+EQ_[41]+EQ_[42]+EQ_[43]+EQ_[52]+EQ_[54]+EQ_[56]+EQ_[57]+EQ_[58]+EQ_[63]+EQ_[65]+EQ_[66]+EQ_[68]+EQ_[69]+EQ_[71]+EQ_[72]+EQ_[74]+EQ_[78]+EQ_[79]+EQ_[81]+EQ_[82]+EQ_[85]+EQ_[86]+EQ_[87]+EQ_[88]+EQ_[89]+EQ_[90]+EQ_[93]+EQ_[94]+EQ_[97]+EQ_[98]+EQ_[99]+EQ_[100]+EQ_[106]+EQ_[107]+EQ_[110]+EQ_[112]+EQ_[115]+EQ_[116]+EQ_[117]+EQ_[118]+EQ_[123]+EQ_[125]+EQ_[127]+EQ_[131]+EQ_[132]+EQ_[138]+EQ_[139]+EQ_[141]+EQ_[142]+EQ_[143]+EQ_[144]+EQ_[145]+EQ_[147]+EQ_[150]+EQ_[151]+EQ_[154]+EQ_[160]+EQ_[162]+EQ_[168]+EQ_[169]+EQ_[171]+EQ_[172]+EQ_[173]+EQ_[175]+EQ_[177]+EQ_[182]+EQ_[183]+EQ_[187]+EQ_[190]+EQ_[192]+EQ_[197]+EQ_[199]+EQ_[201]+EQ_[203]+EQ_[204]+EQ_[206]+EQ_[207]+EQ_[208]+EQ_[209]+EQ_[212]+EQ_[214]+EQ_[215]+EQ_[219]+EQ_[220]+EQ_[221]+EQ_[222]+EQ_[224]+EQ_[225]+EQ_[226]+EQ_[228]+EQ_[231]+EQ_[233]+EQ_[234]+EQ_[237]+EQ_[238]+EQ_[241]+EQ_[242]+EQ_[244]+EQ_[246]+EQ_[248]+EQ_[249]+EQ_[250]+EQ_[254]+EQ_[255]+EQ_[259]+EQ_[260]+EQ_[261]+EQ_[263]+EQ_[264]+EQ_[265]+EQ_[266]+EQ_[268]+EQ_[271]+EQ_[273]+EQ_[274]+EQ_[277]+EQ_[278]+EQ_[279]+EQ_[286]+EQ_[288]+EQ_[289]+EQ_[293]+EQ_[295]+EQ_[296]+EQ_[297]+EQ_[298]+EQ_[300]+EQ_[301]+EQ_[302]+EQ_[306]+EQ_[308]+EQ_[313]+EQ_[314]+EQ_[319]+EQ_[320]+EQ_[322]+EQ_[323]+EQ_[324]+EQ_[325]+EQ_[326]+EQ_[328]+EQ_[329]+EQ_[333]+EQ_[335]+EQ_[337]+EQ_[338]+EQ_[339]+EQ_[340]+EQ_[341]+EQ_[342]+EQ_[343]+EQ_[344]+EQ_[345]+EQ_[349]+EQ_[351]+EQ_[355]+EQ_[356]+EQ_[357]+EQ_[358]+EQ_[359]+EQ_[360]+EQ_[362]+EQ_[363]+EQ_[364]+EQ_[366]+EQ_[367]+EQ_[369]+EQ_[371]+EQ_[373]+EQ_[374]+EQ_[375]+EQ_[378]+EQ_[379]+EQ_[381]+EQ_[382]+EQ_[384]+EQ_[386]+EQ_[387]+EQ_[388]+EQ_[389]+EQ_[390]+EQ_[393]+EQ_[394]+EQ_[396]+EQ_[398]+EQ_[400]+EQ_[401]+EQ_[405]+EQ_[406]+EQ_[407]+EQ_[412]+EQ_[414]+EQ_[415]+EQ_[417]+EQ_[418]+EQ_[420]+EQ_[423]+EQ_[425]+EQ_[429]+EQ_[431]+EQ_[432]+EQ_[434]+EQ_[437]+EQ_[439]+EQ_[447]+EQ_[450]+EQ_[455]

O = [0,0,1,1,0,1,0,1,1,1,0,0,0,0,0,0,0,1,0,0,0,1,1,0,1,0,0,1,1,0,0,1,1,1,0,0,1,1,1,0,1,1,0,1,1,0,0,1,0,0,0,1,0,0,0,0,1,0,1,1,1,1,0,0,0,0,0,1,1,1,1,1,1,0,0,1,1,1,0,0,1,1,1,0,1,0,0,0,0,0,1,1,1,1,1,1,0,1,0,0,0,1,1,0,1,1,1,1,0,0,0,0,1,0,1,0,1,1,0,0,0,1,0,0,0,1,1,1,1,1,1,0,1,0,1,0,0,0,0,1,0,1,0,1,1,0,1,0,1,1,1,0,1,0,1,1,1,1,0,1,1,0,1,0,1,0,0,0,0,0,0,0,1,1,1,0,0,1,0,1,0,1,0,1,1,1,1,0,0,1,1,0,1,1,0,1,1,1,1,1,0,1,0,0,1,0,1,0,0,0,1,1,0,0,1,0,0,1,1,0,0,1,1,1,0,0,1,1,1,1,0,1,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,1,1,0,0,0,1,1,0,1,1,0,1,0,0,1,1,1,1,0,0,1,1,1,0,1,0,1,1,1,1,1,0,0,1,1,0,1,0,1,0,1,0,1,1,1,1,1,1,1,0,1,1,1,0,0,1,1,0,1,0,1,1,0,0,0,1,1,0,0,0,0,0,0,0,0,1,1,0,0,1,0,0,0,1,1,0,0,0,0,1,0,0,1,0,0,0,1,1,1,1,0,1,1,0,0,0,1,0,0,1,1,1,0,0,1,1,0,0,1,1,1,0,0,0,0,1,1,1,1,1,1,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,1,0,0,0,1,0,0,1,1,0,0,0,1,1,1,0,1,0,1,0,0,1,0,1,1,0,1,1,0,1,1,0,0,0,1,1,0,0,0,1,0,0,0,0,1,1,0,0,1,1,1,0,0,1,0,1,1,0,0,1,0,0,0,0,0,0,1,1,1,0,0,0,0,1,0,1,1,1,0,1,0,1,0,1,0,0,0,0,0,0,0,0,0,1,1,1,1,0,1,0,0,1,1,0,1,0,1,1,0,0,0,1,1,1,0,1,0,1,1,0,1,0,1,1,0,1,0,1,0,1,0,0,1,1,1,0,0,1,0,0,1,0,1,0,0,1,1,1,1,0,0,1,1,1,1,1,1,1,0,1,0,0,0,1,0,0,1,1,0,1,1,1,1,1,0,1,1,0,0,0,1,0,1,1,0,1,0,0,0,0,0,1,1,0,0,1,0,0,1,1,0,1,0,1,1,0,0,0,1,1,0,0,0,1,0,1,0,1,1,0,1,0,0,1,1,1,1,0,0,0,0,0,0,1,0,1,1,0,0,0,1,0,0,1,0,1,0,1,1,0,0,0,0,1,1,1,0,0,1,1,0,0,0,1,1,0,1,0,1,1,0,0,1,1,0,0,0,0,1,1,0,1,1,0,0,0,1,0,1,1,1,1,0,0,0,1,0,0,0,0,1,0,1,1,0,1,1,0,1,1,1,1,0,1,0,1,1,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,1,0,0,0,0,0,0,1,1,0,0,0,1,1,1,0,1,0,1,1,1,1,0,0,1,0,1,0,1,0,0,1,0,1,0,1,0,1,1,1,0,1,1,1,0,0,0,1,0,0,1,1,0,1,0,0,0,1,0,1,0,1,1,0,1,1,0,1,1,0,1,1,1,0,0,0,1,0,1,1,1,1,1,0,1,0,1,0,1,1,0,1,1,0,1,0,0,0,0,0,0,1,0,0,1,0,1,1,0,0,0,0,1,1,0,1,1,1,0,1,0,1,1,0,1,0,0,1,1,1,1,0,0,1,1,0,0,0,1,1,1,1,1,0,1,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,1,1,0,1,0,0,1,0,1,0,1,0,0,1,0,0,1,1,0,0,1,1,1,0,0,1,0,0,0,1,0,1,0,1,0,0,1,1,0,0,1,0,1,0,0,0,0,1,1,0,0,1,1,1,1,0,1,1,1,1,0,0,1,1,1,1,0,1,0,1,1,0,1,1,1,0,1,0,1,0,0,0,1,1,0,1,0,1,0,1,0,1,0,0,1,1,1,1,1,0,0,1,0,1,0,1,0,1,1,0,0,0,0,0,1,1,1,0,0,1,1,1,0,0,1,1,1,0,0,1,0,0,0,0,1,0,0,1,0,0,0,0,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,0,1,0,0,1,1,1,1,0,0,0,0,1,1,1,0,0,1,1,1,1,1,0,0,1,0,0,0,0,0,0,0,1,0,0,1,0,0,1,1,0,1,1,1,1,0,0,1,1,0,0,1,1,0,1,0,0,1,1,1,0,0,1,0,0,0,1,1,1,0,0,1,1,0,0,1,0,1,1,1,0,0,0,1,0,0,1,0,0,1,1,1,1,1,0,0,0,0,1,0,1,1,0,0,1,0,0,1,1,1,1,0,0,1,0,0,0,1,0,0,1,1,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,1,0,1,1,1,1,1,1,1,0,1,0,0,0,0,1,1,0,1,1,0,0,1,1,0,0,1,1,1,0,0,1,1,0,0,0,0,0,0,1,0,1,0,0,1,1,1,0,0,1,0,0,0,0,1,1,1,1,0,0,1,0,1,1,0,1,1,1,1,0,0,1,1,1,1,1,1,1,1,1,1,0,1,1,1,0,1,0,1,1,1,0,0,1,0,1,1,0,0,1,0,1,1,1,0,1,0,1,1,1,0,0,1,0,1,0,1,1,0,1,0,1,0,1,0,1,0,0,1,0,0,0,1,0,1,1,0,0,1,1,0,0,1,1,0,0,1,0,1,0,1,1,1,0,1,1,1,1,0,1,1,0,1,0,1,0,1,1,1,0,1,0,1,1,0,1,0,1,1,1,1,1,1,0,0,0,0,0,0,1,0,0,0,0,1,0,1,1,1,1,0,1,1,1,0,1,0,1,1,1,0,1,0,0,1,0,1,0,0,0,0,0,1,1,0,0,1,1,1,1,1,0,1,0,0,1,0,0,1,1,0,0,1,0,1,0,1,1,0,0,0,1,0,0,1,0,0,1,0,0,0,1,0,1,1,0,1,1,0,0,1,0,1,1,1,1,1,1,0,0,0,1,0,1,0,1,0,1,1,0,0,1,1,0,1,1,0,0,1,1,0,0,0,0,0,0,0,1,0,1,0,0,1,0,0,1,1,1,0,0,0,0,1,0,1,0,0,1,1,0,1,1,1,1,0,0,1,0,0,0,1,1,1,1,1,0,0,1,1,1,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1,0,0,1,0,1,0,0,1,0,0,1,0,0,1,1,0,1,0,1,0,1,1,0,1,1,1,1,0,0,1,1,0,0,0,1,0,1,1,0,1,1,0,0,0,1,1,1,1,1,0,0,0,1,0,0,0,0,1,1,0,1,0,1,1,1,1,1,0,0,1,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,1,0,0,0,1,0,0,0,0,0,0,0,1,1,0,1,1,1,1,1,0,0,1,1,1,0,0,0,1,0,1,0,1,1,0,0,0,1,1,1,1,1,1,1,0,0,0,1,0,1,0,1,0,1,1,1,0,0,0,0,0,1,0,1,0,0,0,1,1,0,0,0,0,0,0,0,1,0,1,1,0,0,0,1,0,1,1,1,1,1,1,0,0,0,0,1,0,0,1,0,0,1,1,0,1,0,1,1,0,0,1,1,0,1,0,1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,1,0,1,0,0,0,0,1,1,1,1,1,0,0,1,0,1,1,0,0,1,1,0,1,0,0,1,0,1,0,0,1,1,1,1,1,0,0,1,0,1,0,0,0,1,0,0,1,0,0,1,0,1,1,0,0,1,0,1,1,1,0,1,1,1,0,0,1,0,0,0,1,0,0,0,1,0,1,0,1,1,1,0,1,0,0,1,1,1,0,0,0,1,0,1,1,0,1,0,1,0,0,1,1,1,1,1,0,0,1,0,0,1,0,1,0,0,0,1,0,1,0,0,1,1,1,0,0,0,1,1,1,1,1,1,1,1,0,1,0,1,0,1,0,1,0,0,1,1,1,1,0,1,0,0,0,0,1,1,1,0,0,0,1,0,0,1,1,1,1,0,0,1,1,1,0,1,0,1,0,1,0,0,0,1,0,1,0,0,0,0,0,0,1,0,0,0,1,1,1,0,0,1,1,0,1,1,0,0,1]

#Here we are saving ourselves some time by creating the strings for the variables and the array posisiton to subsitute on the fly
for i in range(len(EQ_1)):
    for j in range(len(O)):
        m = eval("x" + str(state_size+j))
        n = eval('O['+str(j)+']')
        EQ_1[i]=EQ_1[i].subs({m:n})

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

# for i in range(len(EQ_Sub)):
# #     if (EQ_1[i].monomial_coefficient(1)==1):
#     if (EQ_1[i].monomial_coefficient(1)==1):      
#         w[i]=1
        
# print(M.solve_right(w))
A = M.right_kernel()
print(A[1])

end = timeit.timeit()
print(end - start)