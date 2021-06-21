#!/users/kshahin/miniconda2/bin/python

import numpy as np
import matplotlib.pyplot as plt
import hera_cal.abscal as abscal
import uvtools.dspec as dspec
import hera_pspec
from hera_pspec import pspecbeam
from hera_pspec import pspecdata
import itertools
import scipy
from scipy import signal
import pickle
import copy
from hera_cal.utils import polnum2str, polstr2num, jnum2str, jstr2num
from hera_cal.io import HERAData, HERACal
from hera_cal.io import DataContainer
from hera_cal import apply_cal
from hera_cal import io
from hera_cal import smooth_cal
from hera_cal import vis_clean
from hera_cal import redcal
from pyuvdata import UVFlag
import glob
import tqdm
import os
import shutil
from hera_cal import frf
import imp
from hera_pspec import utils as uvp_utils
from hera_pspec import plot as pspec_plot
from hera_pspec import grouping
import argparse

ap = argparse.ArgumentParser(description="Calibration")
ap.add_argument("flag_file", type=str, help="Flag File")
ap.add_argument("data_file", type=str, help="Data File")
ap.add_argument("model_file", type=str, help="Model File")
args = ap.parse_args()
flag_files = [args.flag_file]
data_files = [args.data_file]
model_files = [args.model_file]

ntimes = []
for flags in flag_files:
    uvf = UVFlag(flags)
    ntimes.append(uvf.Ntimes)
nflagstimes = np.min(ntimes)
nfiles = []
for day,flagfile in enumerate(flag_files):
    if not os.path.exists(f"Day_{day}"):
        os.mkdir(f"Day_{day}")
    hdlst = io.HERAData(data_files)
    lsts = np.unique(np.hstack(list(hdlst.lsts)))
    dlst = np.mean(np.diff(lsts))
    flags = UVFlag(flagfile)
    flags.select(frequencies = flags.freq_array[(flags.freq_array>=115*1e+6) & (flags.freq_array<175*1e+6)])
    flags.select(times=flags.time_array[(flags.lst_array>=lsts.min() - dlst/2) & (flags.lst_array<=lsts.max()+dlst/2)])
    flags_chunks = [flags.flag_array.squeeze()[i*60:(i+1)*60] for i in range(nflagstimes//60)]
    file_number=0

    for datafile, modelfile, flags_chunk in tqdm.tqdm(zip(data_files, model_files, flags_chunks)):
        hd_data = HERAData(datafile)
        freqs = hd_data.freqs[(hd_data.freqs>=115*1e+6) & (hd_data.freqs<175*1e+6)]
        data, flag, nsample = hd_data.read(polarizations=["nn"], frequencies=freqs)
        for bl in data:
            flag[bl] = flags_chunk
        hd_data.update(flags=flag)
        hd_data.write_uvh5(f"Day_{day}/data_{day}_{file_number}.uvh5", clobber=True)
        del data, flag, nsample, hd_data
        redcal.redcal_run(input_data=f"Day_{day}/data_{day}_{file_number}.uvh5", clobber=True, solar_horizon=90, verbose=True)
        abscal.post_redcal_abscal_run(data_file=f"Day_{day}/data_{day}_{file_number}.uvh5", redcal_file=f"Day_{day}/data_{day}_{file_number}.omni.calfits", model_files=[f"Day_{day}/data_{day}_{file_number}.uvh5"], clobber=True, data_solar_horizon=90, model_solar_horizon=90)
        file_number += 1
    nfiles.append(file_number)

    cs=smooth_cal.CalibrationSmoother(calfits_list=sorted(glob.glob(f"Day_{day}/data_{day}_*.abs.calfits")))
    cs.time_freq_2D_filter(time_scale=21600)
    cs.write_smoothed_cal(clobber=True, output_replace=(".abs.",".smooth_abs."))

    for file_number in tqdm.tqdm(range(nfiles[-1])):
        apply_cal.apply_cal(data_infilename=f"Day_{day}/data_{day}_{file_number}.uvh5", data_outfilename=f"Day_{day}/data_{day}_{file_number}_smoothcal.uvh5", new_calibration=f"Day_{day}/data_{day}_{file_number}.smooth_abs.calfits", clobber=True)

        vc = vis_clean.VisClean(f"Day_{day}/data_{day}_{file_number}_smoothcal.uvh5")
        vc.read()
        vc.vis_clean(standoff=100, min_dly=600, mode="dpss_leastsq", flag_model_rms_outliers=True, max_contiguous_edge_flags=1)
        vc.write_filtered_data(filled_outfilename=f"Day_{day}/data_{day}_{file_number}_filtered.uvh5", clobber=True)
        os.remove(f"Day_{day}/data_{day}_{file_number}.uvh5")
