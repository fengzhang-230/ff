# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 14:12:56 2019

@author: fzhang
if aa.empty:
    print('NO')
else:
    print('OK')
"""

import pandas as pd
import numpy as np
from os import path
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.dates import DateFormatter, MonthLocator
from One_Dimensional import *

# ==== SETTING-UP =============================================================
dir_data = r'E:\Model\python\Ebo'
dir_src = r'E:\Model\python\Ebo\EB_Py'
dir_pic = r'E:\Model\python\Ebo\Pic'

# execfile(path.join(dir_src, 'One_Dimensional.py'))

n_snow = 5
config_file = r'E:\Model\python\Ebo\EB_Py\config_file.cfg'
sim = 'BND2012-2019.csv'
depth = 'Depth.csv'
SoilT = 'relationship2012-2019.csv'
grad = 'rtt_12.CSV'
bore_t = 'Bore_tem.csv'
bore_d = 'bore_depth.csv'
bcl = 'ncl_default.rgb'

simfile = path.join(dir_data, sim)
depthfile = path.join(dir_data, depth)
SoilTfile = path.join(dir_data, SoilT)
gradfile = path.join(dir_data, grad)
boretfile = path.join(dir_data, bore_t)
boredfile = path.join(dir_data, bore_d)
colorfile = path.join(dir_data, bcl)

simdata = pd.read_csv(simfile, skiprows=1, delimiter=',',
                      names=['Date', 'Ta', 'Ts', 'Snd', 'Sden'])
simdata.index = pd.to_datetime(simdata['Date'])
years = simdata.index.year.unique()

SoilT = pd.read_csv(SoilTfile)
SoilT.Date = pd.to_datetime(SoilT['Date'])
SoilT.index = pd.to_datetime(SoilT['Date'])

gradient = np.loadtxt(gradfile, encoding='UTF-8-SIG')
xyn = pd.read_csv(depthfile, header=None)
node = np.size(xyn) + n_snow

date0 = '10/7/2012'
date1 = '4/9/2019'
# ========================== Run  =====================================
one_d = Permafrost_one_dimesional(n_snow, config_file, simdata)
one_d.update(node)
one_d.load_ini_tmp(xyn, gradient)


Comb = Combine_obs_sim(SoilT)
calib_depths = one_d.load_Init_para()
# =====================Plot compare of sim and obs soil temperatres==================================
matplotlib.rcdefaults()
resu_pic = path.join(dir_pic, 'comp_temp.png')
fig = plt.figure(figsize=(7.2, 6.2))
for i in np.arange(len(calib_depths)):
    # ax depend on the number of calib_depths.
    ax = fig.add_subplot(3, 2, i+1)
    Comb.plot_obs_sim_data(ax, calib_depth=depth[i], ylim=[-15, 15],
                           calib_period_onset=date0, calib_period_end=date1)
    months = MonthLocator(range(1, 13), bymonthday=1, interval=12)

    monthsFmt = DateFormatter("%b")
    ax.xaxis.set_major_locator(months)

    ax.xaxis.set_major_formatter(monthsFmt)
    if i % 2 == 0:
        plt.ylabel(r'Soil Temperature ($^\circ$C)', fontsize=10)
    if i != 4 and i != 5:
        plt.xticks(color='w')
    if i == 4 or i == 5:
        #        plt.ylim([-15,15])
        #        self.plot_temperature_field(calib_period_onset, calib_period_end, ylim)
        plt.xticks(fontsize=9)
        plt.yticks(fontsize=9)

        months = MonthLocator(range(1, 13), bymonthday=1, interval=12)

        monthsFmt = DateFormatter("%Y-%m")

        ax.xaxis.set_major_locator(months)

        ax.xaxis.set_major_formatter(monthsFmt)
        for ticketlabel in ax.xaxis.get_ticklabels():
            ticketlabel.set_rotation(30)
        ax.autoscale_view()

    if i == 5:
        ax.legend(fontsize=7, loc=0, shadow=True, frameon=True)

plt.subplots_adjust(hspace=0.25, left=0.125, right=0.8, wspace=0.20)
plt.savefig(resu_pic, dpi=500, bbox_inches='tight')

# ===================plot profile =============================================
colors = np.loadtxt(colorfile, skiprows=3)
Pro_pic = path.join(dir_pic, 'Profile_soilTem.png')
Comb.plot_temperature_field(onset=date0, ending=date1, depth_bot=2.5, colors)

# =======Plot soil temperatures according to sim and bore hole's  =============
bore_tem = pd.read_csv(boretfile)
bore_dep = pd.read_csv(boredfile)

fig1 = plt.figure(figsize=(7, 4.1))
ax1 = fig1.add_subplot(1, 2, 1)

Comb.plot_boretemp_profile(
     ax1, bore_tem, bore_dep,onset=date0, ending=date1, ylim=[0, 2.5])

ax2 = fig1.add_subplot(1, 2, 2)

Comb.plot_obstemp_profile(
    ax2, calib_depths,onset=date0, ending=date1, ylim=[0.05, 0.77])

plt.subplots_adjust(hspace=0.25, left=0.125, right=0.8, wspace=0.25)
plt.savefig(Pro_pic, dpi=500, bbox_inches='tight')
plt.show()
