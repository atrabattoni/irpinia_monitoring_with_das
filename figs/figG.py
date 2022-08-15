import numpy as np
import matplotlib.pyplot as plt

# Parameters
Pkd = 3.8e-9  # J.m^-3.s^-1
wlen = 6  # s
rho = 2720  # kg.m^-3
mu = 3.3e10  # Pa
R = 1e3 * np.arange(5, 59)  # m
ta = 1e6 * np.array([0.03, 0.1, 0.4])  # Pa

# Calculations
csm = np.sqrt(mu / rho)
Ekd = Pkd * wlen
Erd = 8 * np.pi * csm * R**2 * Ekd
M0d = mu / ta[:, np.newaxis] * Erd
MWd = (np.log10(M0d) - 9.1) / 1.5

# Plot
mm = 1 / 25.4
plt.figure(dpi=300, figsize=(95 * mm, 65 * mm), constrained_layout=True)
plt.style.use("figure.mplstyle")
styles = ["k:", "k-", "k:"]
for k in range(len(ta)):
    plt.plot(R / 1e3, MWd[k], styles[k])
    plt.text(R[-1] / 1e3, MWd[k][-1] - 0.29,
             f"$\\tau_a$={ta[k]/1e6:.2f} MPa", ha="right")
plt.xlim(0, 60)
plt.ylim(0, 3)
plt.xlabel("Distance [km]")
plt.ylabel("Magnitude")
plt.savefig("figG.png")
