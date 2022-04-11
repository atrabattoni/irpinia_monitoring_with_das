import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from utils import bandpass, fk_analysis, v2sr

freqmin = 1.0
freqmax = 20.0
h1 = 9.6
cs1 = 60
cs2 = 120
cp1 = 500
cp2 = 1500
shot_offset = 199
ng = 12

cases = ["DAS", "SIM"]
xarrs = {}
fks = {}

for case in cases:
    fname = f"../data/shot.{case}.nc"
    xarr = xr.open_dataarray(fname)
    ds = (xarr['offset'][1] - xarr['offset'][0]).values
    sampling_rate = 1/(xarr['time'][1] - xarr['time'][0]).values
    if case == "SIM":
        xarr["offset"] = xarr["offset"] + shot_offset
        xarr = v2sr(xarr, ng, ds)
        xarr = 4 * xarr
    f, k, fk = fk_analysis(xarr.values/1e9, 1/sampling_rate, ds,
                           window="hamming")
    xarr = xr.apply_ufunc(bandpass, xarr, kwargs=dict(
        freqmin=freqmin, freqmax=freqmax, df=sampling_rate, axis=0))
    xarrs[case] = xarr
    fks[case] = (f, k, fk)

plt.style.use("figure.mplstyle")
mm = 1/25.4
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(190 * mm, 130 * mm),
                         dpi=300, constrained_layout=True)

axes[0, 0].sharex(axes[0, 1])
axes[0, 0].sharey(axes[0, 1])
axes[1, 0].sharex(axes[1, 1])
axes[1, 0].sharey(axes[1, 1])

for index, case in enumerate(cases):
    xarr = xarrs[case]
    ax = axes[0, index]
    vlim = 1e4
    pcm = ax.pcolormesh(
        xarr.offset.values, xarr.time.values, xarr.values,
        rasterized=True, norm=colors.SymLogNorm(linthresh=100, vmin=-vlim, vmax=vlim))
    fig.colorbar(
        pcm, ax=ax, label=r"Strain rate $[\rm{n}m/\rm{m}]$", pad=0.02, aspect=30)
    s = xarr.offset - shot_offset
    t_p = np.abs(s) / cp1
    t_roll = np.abs(s) / cs1
    t_leak = t_roll / 2
    t_head = np.abs(s) / cs2 + 2 * h1 * np.sqrt(cs2**2 - cs1**2) / (cs1 * cs2)
    t_head[t_head > t_roll] = np.nan
    ax.plot(xarr.offset, t_p, c="C3", ls="--")
    ax.plot(xarr.offset, t_roll, c="C3", ls="--")
    ax.plot(xarr.offset, t_leak, c="C3", ls="--")
    ax.plot(xarr.offset, t_head, c="C3", ls="--")
    ax.set_xlim(50, 350)
    ax.set_xlabel(r"Offset $[\rm{m}]$")
    if index == 0:
        ax.set_ylim(3.0, -0.5)
        ax.set_ylabel(r"Time $[\rm{s}]$")

for index, case in enumerate(cases):
    f, k, fk = fks[case]
    ax = axes[1, index]
    pcm = ax.pcolormesh(k, f, fk, rasterized=True, vmin=np.quantile(
        fk, 0.01), vmax=np.max(fk))
    fig.colorbar(pcm, ax=ax, pad=0.02, aspect=30, label=r"PSD $[\rm{dB}]$")
    ax.set_xlabel(r"Wavenumber $[\rm{m}^{-1}]$")
    if index == 0:
        ax.set_ylim(0, 20)
        ax.set_ylabel(r"Frequency $[\rm{Hz}]$")
        ax.set_xlim(k[0], k[-1])

fig.savefig("fig4.jpg")
plt.close(fig)
