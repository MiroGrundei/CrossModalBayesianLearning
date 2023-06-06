from create_regressors_per_block import create_reg_cd
import numpy as np

# model = 'LHtoLH'
model = 'CItoAR'
n = 4
# base = 10 ** -n
# taus = np.logspace(1, 0, num=100, base=base)
exclude = []
SJs = [i for i in range(1, 30) if i not in exclude]
for s in SJs:
    sub = 'sub-{:02d}'.format(s)
    # taus = [0.1155, 0.0630, 0.0433, 0.0224, 0.0136, 0.0069, 0.0046, 0.0034, 0.0023, 0.0017, 0.0014, 0.0009, 0.0007, 0.0001]
    # taus = [0.001, 0.005, 0.01, 0.05]
    # taus = np.logspace(-3, -1, 100)
    taus = [0.0000, 0.001, 0.0015, 0.002, 0.003, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2]
    for tau in taus:
        tau = round(tau, n)
        create_reg_cd(tau, model, sub)