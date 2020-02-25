import numpy as np
from icet.fitting import Optimizer
from icet.tools import Constraints


np.random.seed(42)
n, m = 20, 10

A = np.random.random((n, m))
y = np.random.random(n)

# constraint sum eci[inds] = 0
inds1 = [1, 3, 4, 5]
inds2 = [2, 6, 7, 8]
M = np.zeros((2, m))
M[0, inds1] = 1
M[1, inds2] = 1

c = Constraints(m)
c.add_constraint(M)

Ac = c.transform(A)
opt = Optimizer((Ac, y), fit_method='ridge')
opt.train()

parameters = c.inverse_transform(opt.parameters)
sum_1 = parameters[inds1].sum()
sum_2 = parameters[inds2].sum()
print('constraints 1, ', sum_1)
print('constraints 2, ', sum_2)

assert abs(sum_1) < 1e-12
assert abs(sum_2) < 1e-12
