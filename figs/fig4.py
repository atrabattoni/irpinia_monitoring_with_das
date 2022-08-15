import numpy as np
import matplotlib.pyplot as plt
resolution = np.finfo("float").resolution


def plot_pattern(ax, theta, r, label=None, **kwargs):
    positive = (r > -resolution)
    negative = (r < +resolution)
    rp = +r.copy()
    rp[~positive] = np.nan
    rn = -r.copy()
    rn[~negative] = np.nan
    ax.plot(theta, rp, ls="-", label=label, **kwargs)
    ax.plot(theta, rn, ls="--", **kwargs)


theta = np.linspace(0, 2*np.pi, 361)
velocity_P = np.cos(theta)
velocity_S = np.sin(theta)
strain_rate_P = np.cos(theta) * velocity_P 
strain_rate_S = np.cos(theta) * velocity_S 

plt.style.use("figure.mplstyle")
mm = 1/25.4
fig, axes = plt.subplots(ncols=2, figsize=(190*mm, 60*mm), dpi=300,
    subplot_kw={'projection': 'polar'}, constrained_layout=True)
ax = axes[0]
plot_pattern(ax, theta, velocity_P, c="C0", label=r"$v_x^P$")
plot_pattern(ax, theta, strain_rate_P, c="C1", label=r"$\dot{\varepsilon}_{xx}^{P}$")

ax = axes[1]
plot_pattern(ax, theta, velocity_S, c="C0", label=r"$v_x^S$")
plot_pattern(ax, theta, strain_rate_S, c="C1", label=r"$\dot{\varepsilon}_{xx}^S$")

type = "P"
for ax in axes:
    ax.plot([0, np.pi], [1, 1], c="black")
    ax.set_rmin(0)
    ax.set_rmax(1.01)
    ax.legend(loc="lower center", title=f"{type}-waves")
    ax.axis("off")
    type="S"

plt.savefig("fig4.png")
