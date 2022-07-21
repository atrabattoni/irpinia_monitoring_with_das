from matplotlib.colors import SymLogNorm
import matplotlib.pyplot as plt
import numpy as np
import obspy
from obspy.geodetics.base import locations2degrees, degrees2kilometers
import scipy.signal as sp
import xarray as xr

from utils import bandpass

# event I
event = "20211003.042205"
ylim = 100

# # event II
# event = "20210911.024800"
# ylim = 20

# parameters
clist = [60, 130]
freqmin = 1.0
freqmax = 20.0
before = 2.0
after = 20.0
limit_values = [50, 351, 900]
limit_names = ["A", "B", "C"]
limits = dict(zip(limit_names, limit_values))

# load COL3 data
tr = obspy.read("../data/" + event + "*.sac")[0]
sac = tr.stats.sac
sensibility = sac["user2"] / sac["user3"]
x_col = sensibility * tr.data
t_col = (1e9 * tr.times("timestamp")).astype("datetime64[ns]")
x_col = bandpass(x_col, freqmin, freqmax, tr.stats.sampling_rate, axis=0)

# get picked phases
ptoa = np.datetime64(tr.stats.starttime + sac["a"])
stoa = np.datetime64(tr.stats.starttime + sac["t0"])
p2s = sac["t0"] - sac["a"]
starttime = ptoa - np.timedelta64(round(1000 * before), "ms")
endtime = ptoa + np.timedelta64(round(1000 * after), "ms")

# load DAS data
xarr = xr.open_dataarray("../data/" + event + ".DAS.nc")
xarr = xarr.sel(time=slice(starttime, endtime))
time_step = (xarr["time"][1] - xarr["time"][0]).values / np.timedelta64(1, "s")
sampling_rate = 1 / time_step
xarr = xr.apply_ufunc(bandpass, xarr, kwargs=dict(
    freqmin=freqmin, freqmax=freqmax, df=sampling_rate, axis=0))
s_das, t_das, x_das = xarr.offset.values, xarr.time.values, xarr.values
x_das = sp.decimate(x_das, 4, axis=0)
t_das = t_das[::4]

# convert in local time
t_col = (t_col - ptoa) / np.timedelta64(1, 's')
t_das = (t_das - ptoa) / np.timedelta64(1, 's')

# estimate time of arrival on the cable
fiber = xr.open_dataset("../data/fiber_dgnss.nc")
fiber = fiber.interp_like(xarr)
degrees = locations2degrees(
    sac["evla"], sac["evlo"], fiber["latitude"], fiber["longitude"])
distance = np.sqrt(degrees2kilometers(degrees)**2 + (sac["evdp"] + 0.4)**2)
evd = np.sqrt(sac["dist"]**2 + (sac["evdp"] + sac["stel"]/1000)**2)
cp = evd / (sac["a"] - sac["o"])
cs = evd / (sac["t0"] - sac["o"])
pfront = -(evd - distance) / cp
sfront = p2s - (evd - distance) / cs

# plot
plt.style.use("figure.mplstyle")
mm = 1/25.4
fig, axes = plt.subplots(
    ncols=2, sharey=True, figsize=(95 * mm, 150 * mm), dpi=300,
    constrained_layout=True, gridspec_kw=dict(
        width_ratios=(1, 6)))

ax = axes[0]
ax.set_title("COL3 [Z]", pad=12)
ax.plot(1e6*x_col, t_col, c="black", lw=1)
ax.axhline(0, color="C3", ls=":")
ax.axhline(p2s, color="C3", ls=":")
ax.set_xlim(-ylim, ylim)
ax.set_xlabel(r"V $[\rm{µ}m/\rm{s}]$")
ax.set_ylabel(r"Time $[\rm{s}]$")
ax.set_ylim(after, -before)

ax = axes[1]
ax.set_title("DAS", pad=12)
vlim = np.max(np.abs(x_das))
pcm = ax.pcolormesh(
    s_das, t_das, x_das,
    rasterized=True, norm=SymLogNorm(vlim/50, vmin=-vlim, vmax=vlim))
fig.colorbar(pcm, ax=ax, label=r"Strain rate $[\rm{nm}/\rm{m}/\rm{s}]$",  location='bottom',
             pad=0.02, aspect=30)
ax.plot(fiber["offset"], pfront, color="C3", ls=":")
ax.plot(fiber["offset"], sfront, color="C3", ls=":")
for name, value in limits.items():
    ax.axvline(value, c="black", lw=0.5)
    ax.annotate(name, (value, 1.01), xycoords=(
        "data", "axes fraction"), horizontalalignment="center")
x0 = limits['A']
y0 = -before
for c in clist:
    dx = 150
    dy = dx / c
    ax.plot([x0, x0 + dx], [y0, y0 + dy], color="C3")
    ax.annotate(f"{c} m/s", (x0 + dx, y0 + dy), color="C3", weight="bold")
ax.set_xlabel(r"Offset $[\rm{m}]$")
ax.yaxis.set_tick_params(left=False, labelleft=False)

fig.savefig("fig2.jpg")
plt.close(fig)
