import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# noise statistics
q10 = 6.980e-10
q50 = 1.252e-09
q90 = 4.842e-09

# estimated power density
dataset = pd.read_csv("../data/power_densities.csv")
signal_kinetic = dataset["signal_kinetic"]
signal_elastic = dataset["signal_elastic"]
noise_kinetic = dataset["noise_kinetic"]
noise_elastic = dataset["noise_elastic"]

# signal to noise ratio
snr = (signal_elastic - noise_elastic) / noise_elastic
mask = snr > 2

# fit line on detected events
p = np.polyfit(np.log(signal_kinetic[mask]), np.log(
    signal_elastic[mask]) - np.log(signal_kinetic[mask]), 0)
x = np.array([1e-9, 1e-5])
y = x * np.exp(p)

# plot
mm = 1 / 25.4
plt.style.use("figure.mplstyle")
plt.figure(dpi=300, figsize=(95 * mm, 80 * mm), constrained_layout=True)
plt.loglog(x, y, "--", color="black", label="Linear fit")
plt.axhline(q10, lw=0.75, c="black", ls=":", alpha=0.5, label="DAS noise 10th and 90th quantile")
plt.axhline(q50, lw=0.75, c="black", ls="-", alpha=0.5, label="DAS noise median")
plt.axhline(q90, lw=0.75, c="black", ls=":", alpha=0.5)
plt.loglog(signal_kinetic[~mask], signal_elastic[~mask], "X", c="C3", label="DAS SNR < 2")
plt.loglog(signal_kinetic[mask], signal_elastic[mask], "o", c="C2", label="DAS SNR > 2")
plt.legend()
plt.ylabel("Elastic power density (DAS) [$J m^{-3} s^{-1}$]")
plt.xlabel("Kinetic power density (COL3) [$J m^{-3} s^{-1}$]")
plt.savefig("fig9.png")
plt.close()
