import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.special import exp1  # Exponential integral W(u)
with open("w1005-in.dat", "r") as f:
    lines = f.readlines()
T = float(lines[0].split()[1])  # hk
S = float(lines[1].split()[1])  # storativity
# Fixed parameters
Q = -50.85
r = 3.35
t1 = 3.88 / 24

# Time array
time_days = np.logspace(np.log10(0.000931), np.log10(0.18), 2731)
def theis_superposition(Q, T, S, r, t_array, t1):
    s_total = []
    epsilon = 1e-10  # small number to avoid division by zero or tiny times

    for t in t_array:
        # Early pumping period
        u1 = (r**2 * S) / (4 * T * max(t, epsilon))
        s1 = exp1(u1)

        if t <= t1:  # during pumping
            s = (Q / (4 * np.pi * T)) * s1
        else:  # after pumping stops (rebound)
            u2 = (r**2 * S) / (4 * T * max(t - t1, epsilon))
            s2 = exp1(u2)
            s = (Q / (4 * np.pi * T)) * (s1 - s2)

        s_total.append(s)

    return np.array(s_total)

drawdown = theis_superposition(Q, T, S, r, time_days, t1)
df_model = pd.DataFrame({"Drawdown": drawdown})
df_model.to_csv("w1005-out.dat", index=False, header=False)