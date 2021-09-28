#!/lustre/aoc/projects/hera/kshahin/miniconda3/envs/hera3/bin/python
import pdb
import numpy as np
import matplotlib.pyplot as plt
import hera_cal.abscal as abscal
import itertools
import scipy
from scipy import signal
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
import argparse

ap = argparse.ArgumentParser(description="Calibration")
ap.add_argument("chunks", type=int, help="Chunks per Day")
ap.add_argument('model_directory', type=str, help="Model Directory")
ap.add_argument('data_directory', type=str, help="Data Directory")
args = ap.parse_args()

chunks = args.chunks
chunk = np.mod(chunks,1)
day = chunks//1
model_directory = args.model_directory
data_directory = args.data_directory


data_file = data_directory
model_file = model_directory
flag_files = [f"/lustre/aoc/projects/hera/aewallwi/H1C_flags/{jd}.flags.h5" for jd in [2458098,2458099,2458101,2458102,2458103,2458104,2458105,2458106,
                                                                                                        2458107,2458108,2458109,2458110,2458111,2458112,2458113,2458114,
                                                                                                        2458115,2458116]]
bad_ants = [np.loadtxt(f"/users/kshahin/kshahin/HERA_Calibration/hera_pipelines/pipelines/h1c/idr2/v2/bad_ants/{ba}.txt") for ba in [2458098,2458099,2458101,2458102,2458103,2458104,2458105,2458106,2458107,2458108,
                                                                                                                          2458109,2458110,2458111,2458112,2458113,2458114,2458115,2458116]]
flag_file = flag_files[day]
bad_ant = bad_ants[day]
if not os.path.exists(f"/users/kshahin/kshahin/HERA_Calibration/DayfxP_{day}"):
    os.mkdir(f"/users/kshahin/kshahin/HERA_Calibration/DayfxP_{day}")
flags = UVFlag(flag_file)
flags.select(frequencies = flags.freq_array[(flags.freq_array>=115*1e+6) & (flags.freq_array<175*1e+6)])
flags.select(times=flags.time_array[2600:2660])
hd_data = HERAData(data_file)
freqs = hd_data.freqs[(hd_data.freqs>=115*1e+6) & (hd_data.freqs<175*1e+6)]
data, flag, nsample = hd_data.read(polarizations=["ee"], frequencies=freqs)
for bl in data:
    if (bl[0] == bad_ant).any() or (bl[1] == bad_ant).any():
                flag[bl] = np.ones_like(flag[bl])
    flag[bl] = flags.flag_array.squeeze()
hd_data.update(flags=flag)
hd_data.write_uvh5(f"/users/kshahin/kshahin/HERA_Calibration/DayfxP_{day}/data_{day}_{chunk}.uvh5", clobber=True)
del data, flag, nsample, hd_data
redcal.redcal_run(input_data=f"/users/kshahin/kshahin/HERA_Calibration/DayfxP_{day}/data_{day}_{chunk}.uvh5", clobber=True, solar_horizon=90, verbose=True)
abscal.post_redcal_abscal_run(data_file=f"/users/kshahin/kshahin/HERA_Calibration/DayfxP_{day}/data_{day}_{chunk}.uvh5", redcal_file=f"/users/kshahin/kshahin/HERA_Calibration/DayfxP_{day}/data_{day}_{chunk}.omni.calfits", model_files=[data_file], clobber=True, data_solar_horizon=90, model_solar_horizon=90)

cs=smooth_cal.CalibrationSmoother(calfits_list=sorted(glob.glob(f"/users/kshahin/kshahin/HERA_Calibration/DayfxP_{day}/data_{day}_*.abs.calfits")))
cs.time_freq_2D_filter(time_scale=21600)
cs.write_smoothed_cal(clobber=True, output_replace=(".abs.",".smooth_abs."))


apply_cal.apply_cal(data_infilename=f"/users/kshahin/kshahin/HERA_Calibration/DayfxP_{day}/data_{day}_{chunk}.uvh5", data_outfilename=f"/users/kshahin/kshahin/HERA_Calibration/DayfxP_{day}/data_{day}_{chunk}_smoothcal.uvh5", new_calibration=f"/users/kshahin/kshahin/HERA_Calibration/DayfxP_{day}/data_{day}_{chunk}.smooth_abs.calfits", clobber=True)

vc = vis_clean.VisClean(f"/users/kshahin/kshahin/HERA_Calibration/DayfxP_{day}/data_{day}_{chunk}_smoothcal.uvh5")
vc.read()
vc.vis_clean(standoff=100, min_dly=600, mode="dpss_leastsq", skip_if_flag_within_edge_distance=1, flag_model_rms_outliers=True, max_contiguous_edge_flags=1)
vc.write_filtered_data(filled_outfilename=f"/users/kshahin/kshahin/HERA_Calibration/DayfxP_{day}/data_{day}_{chunk}_filtered.uvh5", clobber=True)
os.remove(f"/users/kshahin/kshahin/HERA_Calibration/DayfxP_{day}/data_{day}_{chunk}.uvh5")
