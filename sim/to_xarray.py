from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


def get_metadata():
    metadata_path = Path("OUTPUT_FILES") / \
        "for_information_STATIONS_actually_used"
    usecols = list(range(4))
    names = ["station", "network", "x", "z"]
    metadata = pd.read_csv(metadata_path, delim_whitespace=True,
                           usecols=usecols, names=names)

    return metadata


metadata = get_metadata()


def get_fname(row, channel):
    return f'{row["network"]}.{row["station"]}.{channel}.semv'


def get_dataset(metadata, get_fname):
    dataset = xr.Dataset.from_dataframe(metadata)

    channels = ["BXX", "BXZ"]
    for channel in channels:
        data = []
        for _, row in metadata.iterrows():
            path = Path("OUTPUT_FILES") / get_fname(row, channel)
            time, waveform = np.loadtxt(path).T
            data.append(waveform)
        data = np.stack(data, axis=-1)
        dataset.coords["time"] = time
        dataset[channel] = (("time", "index"), data)
    return dataset


def to_dataarray(dataset, channel):
    xarr = dataset[channel]
    xarr.coords["index"] = dataset["x"]
    xarr = xarr.rename({"index": "offset"})
    return xarr


def symmetry(xarr):
    time = xarr.coords["time"].values
    offset = xarr.coords["offset"].values
    offset = np.concatenate((-offset[1:][::-1], offset))
    coords = {"time": time, "offset": offset}
    values = np.concatenate(
        (-xarr.values[:, 1:][:, ::-1], xarr.values), axis=1)
    return xr.DataArray(values, coords, xarr.dims)


if __name__ == "__main__":
    dataset = get_dataset(metadata, get_fname)
    xarr = to_dataarray(dataset, "BXX")
    # xarr = symmetry(xarr)  # <<<<<<<<<< uncomment for axisymetric simulation
    xarr.to_netcdf("simulation.nc")
