#!/lustre/aoc/projects/hera/kshahin/miniconda3/envs/hera3/bin/python
import pdb
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
ap.add_argument("job_ID", type=int, help="Job ID")
ap.add_argument('model_directory', type=str, help="Model Directory")
ap.add_argument('data_directory', type=str, help="Data Directory")
args = ap.parse_args()
job_ID = args.job_ID
day = job_ID//32
chunk = np.mod(job_ID,32)
model_directory = args.model_directory
data_directory = args.data_directory


pdb.set_trace()
data_file = sorted(glob.glob(f'{data_directory}/*.uvh5'))[chunk]
model_file = sorted(glob.glob(f'{model_directory}/*.uvh5'))[chunk]
flag_files = [f"/users/kshahin/HERA_Calibration/H1C_Flags/{jd}.flags.h5" for jd in [2458098,2458099,2458101,2458102,2458103,2458104,2458105,2458106,
                                                                                                        2458107,2458108,2458109,2458110,2458111,2458112,2458113,2458114,
                                                                                                        2458115,2458116]]
flag_file = flag_files[day]
if not os.path.exists(f"Day_{day}"):
    os.mkdir(f"Day_{day}")
hdlst = io.HERAData(data_file)
lsts = np.unique(np.hstack(list(hdlst.lsts)))
dlst = np.mean(np.diff(lsts))
flags = UVFlag(flag_file)
flags.select(frequencies = flags.freq_array[(flags.freq_array>=115*1e+6) & (flags.freq_array<175*1e+6)])
flags.select(times=flags.time_array[(flags.lst_array>=lsts.min() - dlst/2) & (flags.lst_array<=lsts.max()+dlst/2)])

hd_data = HERAData(data_file)
freqs = hd_data.freqs[(hd_data.freqs>=115*1e+6) & (hd_data.freqs<175*1e+6)]
data, flag, nsample = hd_data.read(polarizations=["nn"], frequencies=freqs)
for bl in data:
    flag[bl] = flags.flag_array.squeeze()

hd_data.update(flags=flag)
hd_data.write_uvh5(f"Day_{day}/data_{day}_{chunk}.uvh5", clobber=True)
del data, flag, nsample, hd_data
redcal.redcal_run(input_data=f"Day_{day}/data_{day}_{chunk}.uvh5", clobber=True, solar_horizon=90, verbose=True)
abscal.post_redcal_abscal_run(data_file=f"Day_{day}/data_{day}_{chunk}.uvh5", redcal_file=f"Day_{day}/data_{day}_{chunk}.omni.calfits", model_files=[f"Day_{day}/data_{day}_{chunk}.uvh5"], clobber=True, data_solar_horizon=90, model_solar_horizon=90)

cs=smooth_cal.CalibrationSmoother(calfits_list=sorted(glob.glob(f"Day_{day}/data_{day}_*.abs.calfits")))
cs.time_freq_2D_filter(time_scale=21600)
cs.write_smoothed_cal(clobber=True, output_replace=(".abs.",".smooth_abs."))


apply_cal.apply_cal(data_infilename=f"Day_{day}/data_{day}_{chunk}.uvh5", data_outfilename=f"Day_{day}/data_{day}_{chunk}_smoothcal.uvh5", new_calibration=f"Day_{day}/data_{day}_{chunk}.smooth_abs.calfits", clobber=True)

vc = vis_clean.VisClean(f"Day_{day}/data_{day}_{chunk}_smoothcal.uvh5")
vc.read()
vc.vis_clean(standoff=100, min_dly=600, mode="dpss_leastsq", skip_if_flag_within_edge_distance=1, flag_model_rms_outliers=True, max_contiguous_edge_flags=1)
vc.write_filtered_data(filled_outfilename=f"Day_{day}/data_{day}_{chunk}_filtered.uvh5", clobber=True)
os.remove(f"Day_{day}/data_{day}_{chunk}.uvh5")
