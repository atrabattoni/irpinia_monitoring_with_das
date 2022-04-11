import numpy as np
import scipy.signal as sp
from matplotlib.dates import date2num
from scipy.signal import iirfilter, sosfilt, zpk2sos


def bandpass(data, freqmin, freqmax, df, corners=4, zerophase=False, axis=-1):
    fe = 0.5 * df
    low = freqmin / fe
    high = freqmax / fe
    z, p, k = iirfilter(corners, [low, high], btype='band',
                        ftype='butter', output='zpk')
    sos = zpk2sos(z, p, k)
    if zerophase:
        firstpass = sosfilt(sos, data, axis=axis)
        return sosfilt(sos, firstpass[::-1], axis=axis)[::-1]
    else:
        return sosfilt(sos, data, axis=axis)


def fk_analysis(x, dt=1.0, ds=1.0, window="hann", dbscale=True, return_onesided=True):
    """Compute the FK analysis of x with expected dims (time, offset)"""
    # Tapering
    wt = sp.get_window(window, x.shape[0])
    ws = sp.get_window(window, x.shape[1])
    w = wt[:, np.newaxis] * ws[np.newaxis, :]
    out = w * x
    # Scaling
    scale = np.sqrt(dt * ds / np.sum(w**2))
    out *= scale
    # FFT
    f = np.fft.fftshift(np.fft.fftfreq(out.shape[0], dt))
    k = np.fft.fftshift(np.fft.fftfreq(out.shape[1], ds))
    out = np.fft.fftshift(np.fft.fft2(out))
    out = out[:, ::-1]
    if return_onesided:
        mask = (f >= 0)
        out = out[mask, :]
        f = f[mask]
    if dbscale:
        out = 20 * np.log10(np.abs(out))
    return f, k, out


def v2sr(xarr, ng, ds):
    """Convert velocity into strain rate"""
    data = (np.roll(xarr.data, -ng//2, axis=1) -
            np.roll(xarr.data, ng//2, axis=1)) / (ng * ds)
    return xarr.copy(data=data)


def beamform(xarr, c):
    dt = (xarr["time"][1] - xarr["time"][0]).values / np.timedelta64(1, "s")
    fs = 1 / dt
    lags = (xarr["offset"] - xarr["offset"].mean("offset")).values / c * fs
    lags = np.rint(lags).astype("int")
    stack = np.zeros(xarr.shape[0])
    for k in range(xarr.shape[1]):
        signal = xarr[:, k].values
        signal = np.roll(signal, -lags[k])
        stack += signal
    stack = stack / xarr.shape[1]
    return stack


def get_optimal(xarr, cmin, cmax, cstep):
    cgrid = np.arange(cmin, cmax, cstep)
    score = np.zeros(cgrid.size)
    for k, c in enumerate(cgrid):
        stack = beamform(xarr, c)
        score[k] = get_score(stack)
    copt = cgrid[np.argmax(score)]
    stack = beamform(xarr, copt)
    return copt, stack


def get_score(stack):
    return np.sum(stack**2)
