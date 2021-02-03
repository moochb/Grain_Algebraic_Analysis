import scipy.special
from itertools import combinations 
import timeit

start = timeit.timeit()

#define the paramaters

#how many clocks are required
no_clocks = 5000

#the size of the register
state_size = 15

#how many output bits are needed
no_output_bits=5000

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
Left = [0]*no_clocks
Right = [0]*no_clocks

#Generate equations for update using output bits
# S[15] = S[0]+S[7]
for i in range(len(s),no_clocks):
     S[i] = (S[i-15]+S[i-8])

Alpha = [575 ,  571 ,  570 ,  567 ,  566 ,  563 ,  559 ,  557 ,  556 ,  555 ,  554 ,  553 ,  552 ,  550 ,  549 ,  548 ,  547 ,  543 ,  538 ,  537 ,  534 ,  532 ,  530 ,  526 ,  525 ,  523 ,  521 ,  520 ,  519 ,  517 ,  516 ,  513 ,  511 ,  508 ,  507 ,  506 ,  505 ,  503 ,  501 ,  500 ,  499 ,  495 ,  494 ,  487 ,  485 ,  482 ,  481 ,  480 ,  479 ,  476 ,  475 ,  474 ,  473 ,  469 ,  468 ,  467 ,  466 ,  463 ,  462 ,  461 ,  458 ,  457 ,  455 ,  454 ,  452 ,  449 ,  447 ,  445 ,  443 ,  439 ,  437 ,  435 ,  434 ,  431 ,  427 ,  426 ,  423 ,  418 ,  417 ,  406 ,  405 ,  404 ,  398 ,  397 ,  396 ,  395 ,  393 ,  392 ,  387 ,  385 ,  382 ,  380 ,  379 ,  377 ,  376 ,  375 ,  373 ,  372 ,  371 ,  370 ,  368 ,  367 ,  366 ,  364 ,  361 ,  359 ,  358 ,  357 ,  356 ,  355 ,  354 ,  353 ,  352 ,  351 ,  350 ,  348 ,  347 ,  346 ,  344 ,  343 ,  342 ,  341 ,  339 ,  338 ,  335 ,  334 ,  333 ,  332 ,  331 ,  326 ,  325 ,  324 ,  322 ,  321 ,  320 ,  317 ,  316 ,  315 ,  312 ,  309 ,  307 ,  306 ,  303 ,  301 ,  298 ,  287 ,  286 ,  284 ,  280 ,  277 ,  276 ,  275 ,  274 ,  270 ,  269 ,  268 ,  266 ,  251 ,  248 ,  247 ,  245 ,  242 ,  241 ,  238 ,  235 ,  233 ,  232 ,  228 ,  224 ,  222 ,  221 ,  220 ,  219 ,  218 ,  214 ,  208 ,  207 ,  205 ,  204 ,  202 ,  200 ,  197 ,  195 ,  194 ,  191 ,  189 ,  188 ,  187 ,  186 ,  184 ,  182 ,  179 ,  176 ,  174 ,  173 ,  172 ,  170 ,  169 ,  163 ,  162 ,  161 ,  159 ,  157 ,  155 ,  152 ,  151 ,  148 ,  147 ,  144 ,  141 ,  140 ,  139 ,  138 ,  136 ,  135 ,  133 ,  131 ,  130 ,  129 ,  127 ,  125 ,  123 ,  122 ,  121 ,  120 ,  119 ,  115 ,  112 ,  107 ,  105 ,  102 ,  101 ,  100 ,  97 ,  96 ,  95 ,  92 ,  88 ,  85 ,  82 ,  81 ,  80 ,  70 ,  69 ,  68 ,  67 ,  66 ,  62 ,  60 ,  57 ,  56 ,  51 ,  47 ,  44 ,  41 ,  40 ,  39 ,  36 ,  35 ,  34 ,  33 ,  29 ,  28 ,  26 ,  24 ,  23 ,  22 ,  17 ,  16 ,  11 ,  7 ,  5 , 0]
for j in range(no_clocks-state_size):
    for i in Alpha:
    #s1s2z=s1s2s9
    #Here we are going to build two systems; one containing z, the other not
           #Left[j]  = Left[j]+(S[(i+state_size)-(state_size-3)]*S[(i+state_size)-(state_size-1)]*S[(i+state_size)-(state_size-6 )]+ S[(i+state_size)-(state_size-1)]*S[(i+state_size)-(state_size-10)]*S[(i+state_size)-(state_size-6)])
        Left [j] = Left[j]+S[(i+state_size)-(state_size-3)]*S[(i+state_size)-(state_size-1)]*S[(i+state_size)-(state_size-6 )]+ S[(i+state_size)-(state_size-1)]*S[(i+state_size)-(state_size-10)]*S[(i+state_size)-(state_size-6)]+(S[(i+state_size)-(state_size-1)]*S[(i+state_size)-(state_size-10 )]+ S[(i+state_size)-(state_size-10)]*S[(i+state_size)-(state_size-6 )]+S[(i+state_size)-(state_size-3)])*Y[(i+state_size)-state_size]

