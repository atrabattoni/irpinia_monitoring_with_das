import string
import matplotlib.pyplot as plt
import numpy as np
import obspy
import scipy.signal as sp
import xarray as xr
from matplotlib.ticker import MultipleLocator
from matplotlib.offsetbox import AnchoredText
from utils import fk_analysis

# parameters
event = "20211003.042205"
before = 5.0
after = 45.0
limit_values = [50, 351, 900]
limit_names = ["A", "B", "C"]
limits = dict(zip(limit_names, limit_values))

# load COL3 data
tr = obspy.read("../data/" + event + "*.sac")[0]
sensibility = tr.stats.sac["user2"] / tr.stats.sac["user3"]
x_col = sensibility * tr.data
x_col = np.diff(x_col) * tr.stats.sampling_rate
t_col = (1e9 * tr.times("timestamp")).astype("datetime64[ns]")
t_col = t_col[:-1]

# get picked phases
ptoa = np.datetime64(tr.stats.starttime + tr.stats.sac["a"])
stoa = np.datetime64(tr.stats.starttime + tr.stats.sac["t0"])
p2s = tr.stats.sac["t0"] - tr.stats.sac["a"]

starttime = ptoa - np.timedelta64(round(1000 * before), "ms")
endtime = ptoa + np.timedelta64(round(1000 * after), "ms")

# load DAS data
xarr = xr.open_dataarray("../data/" + event + ".DAS.nc")
time_step = (xarr["time"][1] - xarr["time"][0]).values / np.timedelta64(1, "s")
sampling_rate = 1 / time_step
xarr = xarr.sel(
    time=slice(starttime - np.timedelta64(1, 's'),
               endtime + np.timedelta64(1, 's')),
    offset=slice(limits["B"], limits["C"]))
s_das, t_das, x_das = xarr.offset.values, xarr.time.values, xarr.values

# origin times
t0_col = (t_col[0] - ptoa) / np.timedelta64(1, 's')
t0_das = (t_das[0] - ptoa) / np.timedelta64(1, 's')

# spectrograms
f_col, t_col, S_col = sp.spectrogram(
    x_col, fs=tr.stats.sampling_rate, nperseg=256, noverlap=256-16)
f_das, t_das, S_das = sp.spectrogram(
    x_das/1e9, fs=sampling_rate, axis=0, nperseg=384, noverlap=384-16)
S_das = np.mean(S_das, axis=1)
S_col = 10 * np.log10(S_col)
# S_col = 10 * np.log10(S_col*f_col[:, None])
S_das = 10 * np.log10(S_das)

t_col += t0_col
t_das += t0_das

# FKs
f, k, fk = fk_analysis(x_das/1e9, 1/sampling_rate, 2.4,
                       window="hamming")  # ("tukey", 0.25))

# plot
plt.style.use("figure.mplstyle")
mm = 1 / 25.4
fig = plt.figure(constrained_layout=True, figsize=(190 * mm, 85 * mm), dpi=300)
subfigs = fig.subfigures(1, 2)

axes = subfigs[0].subplots(nrows=2, sharex=True, sharey=True)
ax = axes[0]
pcm = ax.pcolormesh(t_col, f_col, S_col, vmin=-180, vmax=-70, rasterized=True)
fig.colorbar(pcm, ax=ax, pad=0.02,
             label=r"PSD $[\rm{dB}\ \rm{re}\ \rm{m}^2\rm{s}^{-4}\rm{Hz}^{-1}]$")
ax.set_ylabel(r"Frequency $[\rm{Hz}]$")
ax.axvline(0, c='C3', ls=":")
ax.axvline(p2s, c='C3', ls=":")
ax.plot(label="COL3")

ax = axes[1]
pcm = ax.pcolormesh(t_das, f_das, S_das, vmin=-170, vmax=-110, rasterized=True)
fig.colorbar(pcm, ax=ax, pad=0.02,
             label=r"PSD $[\rm{dB}\ \rm{re}\ \rm{strain}^2\rm{s}^{-2}\rm{Hz}^{-1}]$")
ax.set_ylabel(r"Frequency $[\rm{Hz}]$")
ax.axvline(0, c='C3', ls=":")
ax.axvline(p2s, c='C3', ls=":")
ax.plot(label="DAS")

ax.set_xlabel(r"Time $[\rm{s}]$")
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.set_xlim(-before, after)
ax.set_ylim(0, 50)

ax = subfigs[1].subplots()
pcm = ax.pcolormesh(k, f, fk, vmin=-160, vmax=-100, rasterized=True)
fig.colorbar(pcm, ax=ax, pad=0.02, aspect=40,
             label=r"PSD $[\rm{dB}\ \rm{re}\ \rm{strain}^2\rm{s}^{-2}\rm{Hz}^{-1}\rm{m}^{-1}]$")
clist = [60, 130]
for c in clist:
    ax.plot([0, k[-1]], [0, k[-1]*c], c="C3", ls=":")
    ax.plot([0, k[0]], [0, abs(k[-1])*c], c="C3", ls=":")
ax.set_ylabel(r"Frequency $[\rm{Hz}]$")
ax.set_xlabel(r"Wavenumber $[\rm{m}^{-1}]$")
ax.set_ylim(0, 20)

letters = iter(string.ascii_letters)
for ax in [*axes, ax]:
    at = AnchoredText(
        f"{next(letters)})", loc="upper left", frameon=False, pad=0., borderpad=0.25,
        prop=dict(fontsize=14, weight="bold", color="white"))
    ax.add_artist(at)

boxes = iter(["COL3", "DAS"])
for ax in axes:
    at = AnchoredText(next(boxes), loc="upper right", pad=0.2, borderpad=0.6)
    at.patch.set_boxstyle("round4")
    ax.add_artist(at)

fig.savefig("fig3.png")
plt.close(fig)
