from scipy.special import betainc, beta
import numpy as np

n_x = 10

n_a = 50
da = 17.4355

n_b = 50
db = 33.98305

result = np.zeros(shape=(n_x * n_a * n_b, 4), dtype=np.float32)
row = 0
for ix in range(1, n_x + 1):
    x = ix / n_x
    for ia in range(1, n_a + 1):
        a = ia * da
        for ib in range(1, n_b + 1):
            b = ib * db
            res = beta(a, b) * betainc(a, b, x)

            result[row, 0] = x
            result[row, 1] = a
            result[row, 2] = b
            result[row, 3] = res
            row += 1

np.savetxt("incbeta_table.csv", result, delimiter=", ", fmt='%1.16e')
