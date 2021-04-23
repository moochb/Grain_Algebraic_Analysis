import pickle
from sage.misc.persist import SagePickler
R = BooleanPolynomialRing(5,'x')
R.inject_variables()
fp = open("hey.txt","w")
hey = SagePickler.dumps(x0)
pickle.dump(hey,fp)