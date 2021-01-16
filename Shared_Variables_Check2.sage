import pickle
from sage.misc.persist import SagePickler
from sage.misc.persist import SageUnpickler
# import sage.misc.persist
R = BooleanPolynomialRing(5,'x')
R.inject_variables()
fp = open("hey.txt","r")
# hey = pickle.load(fp)
SageUnpickler.loads(fp)