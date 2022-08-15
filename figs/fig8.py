import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import obspy
import xarray as xr

from utils import bandpass, v2sr

# parameters
event = "20211003.042205"
freqmin = 1.0
freqmax = 20.0
ng = 12
time_offset = np.timedelta64(500, "ms")
before = 1.0
after = 3.5
limit_values = [50, 351, 900]
limit_names = ["A", "B", "C"]
limits = dict(zip(limit_names, limit_values))

cases = ["model", "bedrock"]
velocities = {}
strain_rates = {}

# process simluated data
for case in cases:
    fname = f"../data/event.{case}.SIM.nc"
    xarr = xr.open_dataarray(fname)
    ds = (xarr['offset'][1] - xarr['offset'][0]).values
    sampling_rate = 1/(xarr['time'][1] - xarr['time'][0]).values
    xarr["offset"] = xarr["offset"] + 150
    velocities[case] = xarr
    xarr = v2sr(xarr, ng, ds)
    xarr = xr.apply_ufunc(bandpass, xarr, kwargs=dict(
        freqmin=freqmin, freqmax=freqmax, df=sampling_rate, axis=0))
    strain_rates[case] = xarr

# get bedrock amplitude value
strain_rate_scale = np.max(
    np.abs(strain_rates["bedrock"].isel(offset=slice(ng, -ng)).values))
velocity_scale = np.max(
    np.abs(velocities["bedrock"].isel(offset=slice(ng, -ng)).values))

# get picked phases
tr = obspy.read("../data/" + event + "*.sac")[0]
ptoa = np.datetime64(tr.stats.starttime + tr.stats.sac["a"]) + time_offset
stoa = np.datetime64(tr.stats.starttime + tr.stats.sac["t0"]) + time_offset
toa = ptoa
starttime = toa - np.timedelta64(round(1000 * before), "ms")
endtime = toa + np.timedelta64(round(1000 * after), "ms")

# load DAS data
xarr = xr.open_dataarray("../data/" + event + ".DAS.nc")
xarr = xarr.sel(
    time=slice(starttime, endtime),
    offset=slice(limits["A"], limits["B"]))
time_step = (xarr["time"][1] - xarr["time"][0]).values / np.timedelta64(1, "s")
sampling_rate = 1 / time_step
ds = (xarr["offset"][1] - xarr["offset"][0]).values
xarr = xr.apply_ufunc(bandpass, xarr, kwargs=dict(
    freqmin=freqmin, freqmax=freqmax, df=sampling_rate, axis=0))
xarr["time"] = (xarr["time"] - toa) / np.timedelta64(1, 's')
strain_rates["event"] = xarr

# plot
plt.style.use("figure.mplstyle")
mm = 1/25.4
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(190 * mm, 65 * mm),
                         dpi=300, constrained_layout=True, sharey=True)

case = "model"
xarr = velocities[case] / velocity_scale

ax = axes[0]
vlim = np.max(np.abs(xarr.values))
linthresh = vlim / 20
pcm = ax.pcolormesh(
    xarr.offset.values, xarr.time.values, xarr.values,
    rasterized=True, norm=colors.SymLogNorm(linthresh=linthresh, vmin=-vlim, vmax=vlim))
fig.colorbar(
    pcm, ax=ax, label=r"Normalized velocity", pad=0.02, aspect=30, orientation="horizontal")
ax.set_xlabel(r"Offset $[\rm{m}]$")
ax.set_ylabel(r"Time $[\rm{s}]$")

xarr = strain_rates[case] / strain_rate_scale

ax = axes[1]
vlim = np.max(np.abs(xarr.values))
linthresh = vlim / 20
pcm = ax.pcolormesh(
    xarr.offset.values, xarr.time.values, xarr.values,
    rasterized=True, norm=colors.SymLogNorm(linthresh=linthresh, vmin=-vlim, vmax=vlim))
fig.colorbar(
    pcm, ax=ax, label=r"Normalized strain rate", pad=0.02, aspect=30, orientation="horizontal")
ax.set_xlabel(r"Offset $[\rm{m}]$")

case = "event"
xarr = strain_rates[case]

ax = axes[2]
vlim = np.max(np.abs(xarr.values)) / 4
linthresh = vlim / 20
pcm = ax.pcolormesh(
    xarr.offset.values, xarr.time.values + 0.5, xarr.values,
    rasterized=True, norm=colors.SymLogNorm(linthresh=linthresh, vmin=-vlim, vmax=vlim))
fig.colorbar(
    pcm, ax=ax, label=r"Strain rate $[\rm{nm}.\rm{m}^{-1}.\rm{s}^{-1}]$", pad=0.02, aspect=30, orientation="horizontal")
ax.set_xlabel(r"Offset $[\rm{m}]$")

cmap = cm.get_cmap('viridis')
rgba = cmap(0.5)
annotations = ["a)", "b)", "c)"]
boxes = ["Model", "Model", "DAS"]
for k, ax in enumerate(axes):
    ax.text(0.02, 0.98, annotations[k], transform=ax.transAxes,
            fontweight='bold', fontsize=14, va='top', ha="left", color="white")
    ax.text(0.96, 0.96, boxes[k], transform=ax.transAxes,
            va='top', ha="right",
            bbox=dict(boxstyle="round4", fc="white", ec="black"))
    ax.set_facecolor(rgba)
    ax.set_ylim(3.0, -0.5)

fig.savefig("fig8.png")
plt.close(fig)
