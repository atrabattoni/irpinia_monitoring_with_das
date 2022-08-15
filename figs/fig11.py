import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as sp
import xarray as xr
from obspy import read

from utils import bandpass, get_optimal

# parameters
event = "20210909.032630"
freqmin = 1
freqmax = 20
offset_lim = (60, 100)
cmin = 60
cmax = 300
cstep = 10

# load COL3 data
st = read(f"../data/{event}.COL3.C03.IX.CN.sac")
tr = st[0]
sensibility = tr.stats.sac["user2"] / tr.stats.sac["user3"]
tr.data = tr.data * sensibility
tr.detrend()
tr.taper(max_percentage=0.05)
tr.filter("bandpass", freqmin=freqmin, freqmax=freqmax)

# get picked phases
starttime = np.datetime64(tr.stats.starttime + tr.stats.sac.a - 5)
endtime = np.datetime64(tr.stats.starttime + tr.stats.sac.t0 + 10)
ptime = np.datetime64(tr.stats.starttime + tr.stats.sac.a)
stime = np.datetime64(tr.stats.starttime + tr.stats.sac.t0)

# load DAS data
xarr = xr.open_dataarray(f"../data/{event}.DAS.nc")
xarr = xarr.sel(
    time=slice(starttime, endtime),
    offset=slice(*offset_lim),
)
dt = (xarr["time"][1] - xarr["time"][0]).values / np.timedelta64(1, "s")
fs = 1 / dt
xarr.data = bandpass(xarr, freqmin, freqmax, fs, axis=0)

# beamform
copt, stack = get_optimal(xarr, cmin, cmax, cstep)
seismo_enveloppe = np.abs(sp.hilbert(tr.data))
coherent_enveloppe = np.abs(sp.hilbert(stack))
incoherent_enveloppe = np.abs(sp.hilbert(xarr, axis=0)).mean(axis=1)
tdas = xarr["time"].values
tseis = (1e6*tr.times("timestamp")).astype("datetime64[us]")

# plot
plt.style.use("figure.mplstyle")
mm = 1 / 25.4
fig, axes = plt.subplots(
    figsize=(95 * mm, 95*mm),
    ncols=4, sharey=True, constrained_layout=True, dpi=300)

ax = axes[0]
ax.set_title("COL3[Z]\n")
ax.fill_betweenx(tseis, -seismo_enveloppe * 1e6,
                 seismo_enveloppe * 1e6, color="dimgray")
ax.set_xlabel("Velocity\n[um/s]")
ax.set_ylabel("Time")
ax.set_xlim(-1, 1)
ax.tick_params(axis='y', labelrotation=90)

ax = axes[1]
ax.set_title("DAS\nincoherent")
ax.fill_betweenx(tdas, -incoherent_enveloppe,
                 incoherent_enveloppe, color="dimgray")
ax.set_xlabel("Strain-rate\n[nm/m]")
ax.set_xlim(-200, 200)
ax = axes[2]
ax.set_title("DAS\ncoherent")
ax.fill_betweenx(tdas, -coherent_enveloppe,
                 coherent_enveloppe, color="dimgray")
ax.set_xlabel("Strain-rate\n[nm/m]")
ax.set_xlim(-200, 200)
ax = axes[3]
xarr.plot(ax=ax, rasterized=True, vmin=-200, vmax=200, yincrease=False)
ax.set_title("DAS\n")
ax.set_ylim(endtime, starttime)
ax.set_xlim(60, 90)
ax.set_xticks([60, 100])
ax.set_xlabel("Offset\n[m]")
ax.set_ylabel("")

for ax in axes:
    ax.axhline(ptime, lw=1, c="C3")
    ax.axhline(stime, lw=1, c="C3")

fig.savefig("fig11.png")
plt.close(fig)
