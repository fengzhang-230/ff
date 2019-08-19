# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 15:48:07 2019

@author: fzhang
"""
import pandas as pd
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter, MonthLocator
from matplotlib.colors import LinearSegmentedColormap



class Permafrost_one_dimesional(object):
    """
    One-dimensional heat transfer model
    """

    def __init__(self, n_snow_layer, Init_para, simdata):
        """
        time step is 1 day, error accepted is 0.001, nowlayer is set to 5,        
        and groung surface started from 5 layer.
        """

        self.nstep = 1
        self.eet = 0.001
        self.ees = 0.001
        self.rtf = 273.15
        self.rhowater = 1000.
        self.ngrnd = int(n_snow_layer)

        self.Init_para = Init_para

        self.load_Init_para()

        self.simdata = simdata
        self.years = sorted(list(set(self.simdata["Year"])))

        x0 = self.simdata.loc[self.date_onset_sim]
        if x0.empty:
            print('\033[3;30;41m' +
                  "can not find the 'Running start date'" + '\033[0m')
        else:
            self.n_start=x0.Date

        x1 = self.simdata.loc[self.date_final_sim]
        if x1.empty:
            print('\033[3;30;41m'+"can not find the 'Forcing by ground surface date'" + '\033[0m')
            print('\033[3;30;41m'+"It was set as ZERO" + '\033[0m')
            self.n_final_calibration = 0
        else:
            self.n_final_calibration =x1.Date
        x2 = self.simdata.loc[self.date_final_running]
        if x2.empty:
            print('\033[3;30;41m' +
                  "can not find the 'Running ending date'" + '\033[0m')
        else:
            self.n_final = x2.Date
        
        self.snow = self.simdata[self.date_onset_sim:self.date_final_sim]['Snd']
        self.airt = self.simdata[self.date_onset_sim:self.date_final_sim]['Ta']
        self.groundt = self.simdata[self.date_onset_sim:self.date_final_sim]['Ts']
        self.sden = self.simdata[self.date_onset_sim:self.date_final_sim]['Sden']
        print("====  Loading initial permafrost temperature profile ...")
        self.load_ini_tmp()
        print("====  Loading soil thermal parameters ...")

        f0 = open('geo_par')
        str0 = f0.readline()
        str0 = str0.split(" ")
        self.i_grids = int(str0[0])
        self.n_levels = int(str0[1])
        self.geo_par = np.loadtxt(self.geopar_filename, skiprows=1)

        self.theta0 = self.geo_par[:, 0]
        self.uw_a1 = self.geo_par[:, 3]
        self.uw_b1 = self.geo_par[:, 4]
        self.uw_c1 = self.geo_par[:, 5]
        self.uw_d1 = self.geo_par[:, 6]
        self.uw_a2 = self.geo_par[:, 7]
        self.uw_b2 = self.geo_par[:, 8]
        self.uw_c2 = self.geo_par[:, 9]
        self.uw_d2 = self.geo_par[:, 10]
        self.k_para = self.geo_par[:, 13]
        self.high = self.geo_par[:, 14]

    def update(self, node):

        print("====  Running ...")
        print("====  Since " + self.date_onset_sim)
        self.node = node
        self.tn = len(self.simdata)
        self.result = np.zeros((self.tn, self.node))
        for ii in np.arange(self.tn):
            print(f'{float(ii) /self.tn:.3f}')

            self.mntime = ii
            self.month = self.groundt.index[self.mntime].month
            self.deltat = 1.0 * self.nstep

            self.snowh = self.snow[self.mntime] + 0.0

            self.snowden = self.sden[self.mntime] + 0.0

            self.snow_k_c()

            self.kcom[0:self.ngrnd] = self.ksnow
            self.cap[0:self.ngrnd] = self.csnow

            self.snowh = round(self.snowh, 2)
            # Following is original diversion in fortran, 5 layers:

            if self.snowh <= 0.003:
                self.lmn = self.ngrnd - 0
            elif self.snowh <= 0.07:
                self.lmn = self.ngrnd - 1
            elif self.snowh <= 0.15:
                self.lmn = self.ngrnd - 2
            elif self.snowh <= 0.24:
                self.lmn = self.ngrnd - 3
            elif self.snowh <= 0.35:
                self.lmn = self.ngrnd - 4
            else:
                self.lmn = self.ngrnd - 5

            for iz in np.arange(self.lmn, self.ngrnd):
                self.xyn[iz] = -1.0*self.snowh * \
                    (self.ngrnd - iz) / (self.ngrnd - self.lmn)
                # self.rtt[iz] = self.rairt[self.mntime]
                self.dx[iz] = self.snowh/(self.ngrnd - self.lmn)
            if self.lmn == self.ngrnd:
                self.rtt[self.lmn] = self.groundt[self.mntime]
            else:
                self.rtt[self.lmn] = self.airt[self.mntime]
            if ii == 0:

                self.snowh == 0.0

                self.lmn = self.ngrnd + 0

                self.rtt[self.lmn] = self.ini_TTOP + 0.

            while 1:
                self.unfrozenwater()

                self.heat_capacity()

                self.heat_conduc()

                self.abcd()

                self.tram()

                if ii != 0:
                    self.rtt[self.lmn +
                             1:self.node] = self.dtt[self.lmn+1:self.node]
                    break
                else:
                    xdf = np.abs(
                        self.dtt[self.lmn+1:self.node] - self.rtt[self.lmn+1:self.node])
                    reee = np.max(xdf)
                    if reee <= self.eet:

                        self.rtt[self.lmn +
                                 1:self.node] = self.dtt[self.lmn+1:self.node]

                        break
                    else:

                        self.rtt[self.lmn +
                                 1:self.node] = self.dtt[self.lmn+1:self.node]

                self.result[self.mntime, ] = self.rtt

        print('\x1b[6;30;42m' + "====  Finish" + '\x1b[0m')

        print("====  ")

    def unfrozenwater(self):
        a = np.zeros(self.node)
        b = np.zeros(self.node)
        c = np.zeros(self.node)
        d = np.zeros(self.node)
        theta = np.zeros(self.node)
        self.thetau = np.zeros(self.node)
        self.thetai = np.zeros(self.node)
        self.rll = np.zeros(self.node)
        self.dxthetau = np.zeros(self.node)

        def fomula_uw(a, b, c, d, t):
            uw = d-(d/(a+b*np.exp(c*t)))
            dxu = (b*c*d*np.exp(c*t))/((a+b*np.exp(c*t))**2)
            return uw, dxu
        vfomula_uw = np.vectorize(fomula_uw)

        # 平衡过程
        # 融化过程

        if self.month in [2, 3, 4, 5, 6, 7, 8, 9, 10]:
            for i in np.arange(self.n_levels):
                if i == 0:

                    idx0 = np.where(self.xyn[self.ngrnd:] <= self.high[i])[
                        0] + self.ngrnd
                else:

                    idx0 = np.where((self.xyn[self.ngrnd:] <= self.high[i]) & (
                        self.xyn[self.ngrnd:] > self.high[i-1]))[0] + self.ngrnd
                a[idx0] = self.uw_a2[i]
                b[idx0] = self.uw_b2[i]
                c[idx0] = self.uw_c2[i]
                d[idx0] = self.uw_d2[i]
                theta[idx0] = 1.0 - self.theta0[i]
            self.thetau[self.ngrnd:], self.dxthetau[self.ngrnd:] = vfomula_uw(a[self.ngrnd:],
                                                                              b[self.ngrnd:], c[self.ngrnd:], d[self.ngrnd:], self.rtt[self.ngrnd:])
            self.thetai[self.ngrnd:] = 1.0 - \
                theta[self.ngrnd:] - self.thetau[self.ngrnd:]
            self.rll[self.ngrnd:] = self.rhowater * (333.20 + 4.955 * self.rtt[self.ngrnd:]
                                                     + 0.02987 * self.rtt[self.ngrnd:]**2.) / 1000.
            idx2 = np.where(self.thetai[self.ngrnd:] < 0.)[0] + self.ngrnd
            self.thetai[idx2] = 0
            self.rkk = np.array([0.25, 0.25, 0.35, 0.1, 1.1, 1.1, 1.1, 2.])

        else:
            for i in np.arange(self.n_levels):

                if i == 0:
                    idx0 = np.where(self.xyn[self.ngrnd:] <= self.high[i])[
                        0] + self.ngrnd
                else:
                    idx0 = np.where((self.xyn[self.ngrnd:] <= self.high[i]) &
                                    (self.xyn[self.ngrnd:] > self.high[i-1]))[0] + self.ngrnd

                a[idx0] = self.uw_a1[i]
                b[idx0] = self.uw_b1[i]
                c[idx0] = self.uw_c1[i]
                d[idx0] = self.uw_d1[i]
                theta[idx0] = 1.0 - self.theta0[i]
              # np.savetxt('para_a.csv',a,delimiter=',')
            self.thetau[self.ngrnd:], self.dxthetau[self.ngrnd:] = vfomula_uw(
                a[self.ngrnd:], b[self.ngrnd:], c[self.ngrnd:], d[self.ngrnd:], self.rtt[self.ngrnd:])
            self.thetai[self.ngrnd:] = 1.0 - \
                theta[self.ngrnd:] - self.thetau[self.ngrnd:]
            self.rll[self.ngrnd:] = self.rhowater * \
                (333.20 + 4.955 * self.rtt[self.ngrnd:] +
                 0.02987 * self.rtt[self.ngrnd:]**2.) / 1000.
            idx2 = np.where(self.thetai[self.ngrnd:] < 0.)[0] + self.ngrnd
            self.thetai[idx2] = 0
            self.rkk = np.array([0.25, 0.25, 0.35, 0.1, 1.1, 1.1, 1.1, 2.])

    def cal_rmse(self, x, y):

        x = np.reshape(x, np.size(x))
        y = np.reshape(y, np.size(y))

        idx_nan = np.where(~np.isnan(x-y))
        return np.sqrt(((x[idx_nan] - y[idx_nan]) ** 2).mean())

    def heat_conduc(self):
        kice = np.zeros(self.node)
        kwater = np.zeros(self.node)
        ksoil = np.zeros(self.node)
        self.kcom = np.zeros(self.node)

        for i in np.arange(self.n_levels):
            if i == 0:
                idx0 = np.where(self.xyn[self.ngrnd:] <= self.high[i])[
                    0] + self.ngrnd
            else:
                idx0 = np.where((self.xyn[self.ngrnd:] <= self.high[i]) & (
                    self.xyn[self.ngrnd:] > self.high[i-1]))[0] + self.ngrnd

            ksoil[idx0] = self.k_para[i]

        kice[self.ngrnd:] = 0.4685 + 488.19 / \
            (self.rtf + self.rtt[self.ngrnd:])

        kwater[self.ngrnd:] = 0.11455 + 1.6318E-3 * \
            (self.rtf + self.rtt[self.ngrnd:])

        self.kcom[self.ngrnd:] = kice[self.ngrnd:] ** self.thetai[self.ngrnd:] *\
            ksoil[self.ngrnd:]**(1.0-self.thetau[self.ngrnd:]-self.thetai[self.ngrnd:]) *\
            kwater[self.ngrnd:] ** self.thetau[self.ngrnd:]

        self.kcom[self.ngrnd:] = self.kcom[self.ngrnd:] * 24.0 * 3600.0

    def heat_capacity(self):
        capice = np.zeros(self.node)
        capwater = np.zeros(self.node)
        cap_v = np.zeros(self.node)
        self.cap = np.zeros(self.node)
        csoil = self.geo_par[:, 1]
        # 青藏高原不同半径热融湖下融区发展差异的非线性分析    令锋、吴青柏
        for i in np.arange(self.n_levels):
            if i == 0:
                idx0 = np.where(self.xyn[self.ngrnd:] <= self.high[i])[
                    0] + self.ngrnd
            else:
                idx0 = np.where((self.xyn[self.ngrnd:] <= self.high[i]) & (
                    self.xyn[self.ngrnd:] > self.high[i-1]))[0] + self.ngrnd
            csoil[idx0] = csoil[i]
        capice[self.ngrnd:] = 1.94 + 7.14E-3 * self.rtt[self.ngrnd:]
        capwater[self.ngrnd:] = 4.20843 + 1.11362E-1 * self.rtt[self.ngrnd:] + \
            5.12142E-3 * self.rtt[self.ngrnd:]**2 + \
            9.3482E-5 * self.rtt[self.ngrnd:]**3
        cap_v[self.ngrnd:] = self.thetau[self.ngrnd:] * capwater[self.ngrnd:] + \
            self.thetai[self.ngrnd:] * capice[self.ngrnd:] + \
            (1.0 - self.thetau[self.ngrnd:] -
             self.thetai[self.ngrnd:]) * csoil[self.ngrnd:]
        self.cap[self.ngrnd:] = (cap_v[self.ngrnd:] +
                                 self.rll[self.ngrnd:]*self.dxthetau[self.ngrnd:]) * 1.0E6

    def snow_k_c(self):
        self.ksnow = self.snowden ** 2 * 2.9E-6
        self.csnow = self.snowden * 2.09E3
        # self.rksnow =0.0688*np.exp(0.0088*self.airt+4.6682*self.rosnow)
        self.ksnow = self.ksnow * 24.0 * 3600

    def load_Init_para(self):
        '''
       Init_para,including Bnds_filename：Ta, Tg, Depsnow, Densnow
       Bnds_filename: soil depth; geopar_filename: soil thermal parameters and 
       soil VWC parameters at different depths; date_onset_sim：star year of 
       this simulation, date_final_sim: end year of simulation; 
       ini_TTOP : first day's temperature(soil or air); T_Gradient: Tem Gradient
       from site; hflux: Low boundary, heat flux ; 
       calib_dep: Number of Calibration soil depth; calib_depths: Calibration depths
       Calibration_Data_File: Soil tempeartures of every day.

        '''

        f1 = open(self.Init_para)
        print('==============Init_para==============================')
        print("====  Init_para" + self.Init_para + ' .....')
        print('=====================================================')
        # line-1 :
        # thermal parameter of Soil
        ss = self.extract_cfg(f1)
        self.geopar_filename = ss

        # line-2 :
        ss = self.extract_cfg(f1)
        self.date_onset_sim = ss

        # line-3 :
        ss = self.extract_cfg(f1)
        self.date_final_sim = ss

        # line-4 :
        ss = self.extract_cfg(f1)
        self.date_final_running = ss

        ss = self.extract_cfg(f1)
        self.ini_TTOP = float(ss)

        # line-5 :
        ss = self.extract_cfg(f1)
        self.hflux = float(ss)

        # line-6 :
        ss = self.extract_cfg(f1)
        self.calib_dep = int(ss)

        # line-7 :
        calib_depths = np.zeros(self.calib_dep)
        ss = self.extract_cfg(f1)
        ss = ss.split(',')
        cont = -1

        if np.size(ss) != self.calib_dep:
            print("!!! Error:: please check calibration input ")
        else:
            for ii in ss:
                cont = cont + 1
                calib_depths[cont] = float(ii)
        return calib_depths
    def extract_cfg(self, f1):

        ss = f1.readline()
        ss = ss.split("|")
        ss = ss[1]

        return ss.strip()

    def abcd(self):
        self.A = np.zeros(self.node)
        self.C = np.zeros(self.node)

        self.B = np.zeros(self.node)
        self.D = np.zeros(self.node)

        self.rkw = np.zeros(self.node)
        self.rke = np.zeros(self.node)

        self.rcw = np.zeros(self.node)
        self.rce = np.zeros(self.node)
        self.rcccp = np.zeros(self.node)
        self.rdx = np.zeros(self.node)
        self.apzero = np.zeros(self.node)

        self.rkw[self.lmn+1:self.ngrnd+1] = self.kcom[self.lmn:self.ngrnd]
        self.rkw[self.ngrnd + 1:(self.node - 1)] = 2.0 * self.kcom[self.ngrnd:(self.node - 2)] * \
            self.kcom[(self.ngrnd+1):(self.node-1)] / (self.kcom[self.ngrnd:(self.node-2)] +
                                                       self.kcom[(self.ngrnd+1):(self.node-1)])

        self.rke[self.lmn+1:self.ngrnd] = self.kcom[self.lmn+1:self.ngrnd]
        self.rke[self.ngrnd:(self.node - 1)] = 2.0 * self.kcom[self.ngrnd:(self.node - 1)] *
        self.kcom[(self.ngrnd+1):(self.node)] / (self.kcom[self.ngrnd:(self.node-1)] +
                                                 self.kcom[(self.ngrnd+1):(self.node)])

        self.A[self.lmn+1:self.node-1] = self.rkw[self.lmn +
                                                  1:self.node-1] / self.dx[self.lmn:self.node-2]
        self.C[self.lmn+1:self.node-1] = self.rke[self.lmn +
                                                  1:self.node-1] / self.dx[self.lmn+1:self.node-1]

        self.rcw[self.lmn+1:self.node-1] = self.rccc[self.lmn:self.node -
                                                     2] * self.dx[self.lmn:self.node-2]
        self.rce[self.lmn+1:self.node-1] = self.rccc[self.lmn +
                                                     2:self.node] * self.dx[self.lmn+1:self.node-1]

        self.rcccp[self.lmn + 1:self.node - 1] = (self.rcw[self.lmn + 1:self.node - 1] +
                                                  self.rce[self.lmn+1:self.node-1]) / (self.dx[self.lmn+1:self.node-1] +
                                                                                       self.dx[self.lmn:self.node-2])

        self.rdx[self.lmn + 1:self.node - 1] = (self.dx[self.lmn + 1:self.node - 1] +
                                                self.dx[self.lmn:self.node-2])/2.0

        self.apzero[self.lmn + 1:self.node - 1] = self.rcccp[self.lmn + 1:self.node - 1] * \
            self.rdx[self.lmn+1:self.node-1] / self.rdeltat

        self.B[self.lmn+1:self.node-1] = -1.0*(self.A[self.lmn+1:self.node-1] +
                                               self.C[self.lmn+1:self.node-1] + self.apzero[self.lmn+1:self.node-1])

        self.D[self.lmn+1] = -1.0*self.apzero[self.lmn+1] * \
            self.rtt[self.lmn+1] - self.A[self.lmn+1] * self.rtt[self.lmn]
        self.D[self.lmn+2:self.node-1] = -1.0*self.apzero[self.lmn +
                                                          2:self.node-1] * self.rtt[self.lmn+2:self.node-1]

        self.A[self.node-1] = self.kcom[self.node-2]/self.dx[self.node-2]
        self.C[self.node-1] = 0.0
        self.B[self.node-1] = -1.0*self.A[self.node-1]
        self.D[self.node-1] = -1.0*self.qq

    def tram(self):

        rp = np.zeros(self.node)
        rq = np.zeros(self.node)

        rp[self.lmn+1] = -1.0*self.C[self.lmn+1]/self.B[self.lmn+1]
        rq[self.lmn+1] = self.D[self.lmn+1]/self.B[self.lmn+1]

        for i in np.arange(self.lmn + 2, self.node):
            pp = self.A[i] * rp[i-1] + self.B[i]
            rp[i] = -1.0*self.C[i]/pp
            rqq = self.D[i] - self.A[i] * rq[i-1]
            rq[i] = rqq/pp

        self.dtt = np.zeros(self.node)
        self.dtt[self.node-1] = rq[self.node-1]
        for i in np.arange(self.node-2, self.lmn+0, -1):
            self.dtt[i] = rp[i] * self.dtt[i + 1] + rq[i]

    def load_ini_tmp(self, xyn, gradient):
        # Number of soil grids, # of soil + 5 layers of snow
        self.rtt = np.zeros(self.node)
        self.dx = np.zeros(self.node)
        self.xyn = np.zeros(self.node)

        self.xyn[self.ngrnd:] = xyn

        self.dx[self.ngrnd:self.node - 1] = self.xyn[(self.ngrnd + 1):self.node] \
            - self.xyn[self.ngrnd:self.node-1]
        self.rtt[self.ngrnd:self.node] = gradient

    def finalize(self):
        print "====  Reorgnizing model results ..."
        # Extract soil temperatures:
        self.model_date = self.simdata[self.n_start:self.n_final+1]['Date']
        self.model_depth = self.xyn[self.ngrnd:]

        self.model_input_air = self.rairt[self.n_start:self.n_final+1]
        self.model_input_gsf = self.rmt[self.n_start:self.n_final+1]
        self.model_input_snd = self.rsnow[self.n_start:self.n_final+1]

        self.model_t_soil = self.result[self.n_start:self.n_final+1,self.ngrnd:]
        self.model_date_matrix = self.datevec(self.model_date)

        print("====  ")
        print('\x1b[5;30;42m' + "====  The end ..." + '\x1b[0m')
        return self.model_date, self.model_t_soil, self.model_depth
    def datevec(self, date_ordinal):

        n_date = np.size(date_ordinal)
        date = np.zeros([n_date, 3], dtype='i4')

        cont = -1
        for xx in date_ordinal:
            cont = cont + 1
            x = pd.to_datetime(xx)
            date_0[cont, 0] = x.year
            date_0[cont, 1] = x.month
            date_0[cont, 2] = x.day
        return date

class Combine_obs_sim(object):
    '''
    Combine the observed and simulated results.
    '''

    def __init__(self, obs_soilt):

        self.obs_soilt = obs_soilt #  observed soil temperatures at depths

    def combine_obsim(self, calib_period_onset, calib_period_end, calib_depth):

        one_d = Permafrost_one_dimesional()
        self.model_date, self.model_t_soil, self.model_depth = one_d.finalize()
        calib_depths = one_d.load_Init_para()
        x0 = self.model_date.loc[calib_period_onset]

        if x0.empty:
            print('\x1b[6;30;42m' + \
                "can not find the 'Calibration start date' in Simulated results" + '\x1b[0m')
        else:
            print('OK')
            x0 = x0
            
        x1 = self.model_date.loc[calib_period_end]
        if x1.empty:
            print('\x1b[6;30;42m' +
                "can not find the 'Calibration rnd date' in Simulated results" + '\x1b[0m') 
        else:
            x1=x1
    
        self.obs_date = self.obs_soilt.index[0]
        y0 = self.obs_date.loc[calib_period_onset]
        if y0.empty:
            print("can not find the 'Calibration start date' in obs Soil Temp")
        else:
            y0 = y0.Date
        y1 = self.obs_date.loc[calib_period_end]
        if y1.empty:
            print("can not find the 'Calibration end date' in obs Soil Temp")
        else:
            y1 = y1.Date
        idx1 = np.where(abs(calib_depths - calib_depth) <= 0.0001)[0] + 1
        self.col = self.obs_soilt.columns[idx]
        obs = self.obs_soilt[y0:y1][col]

        self.model_t_soil = pd.DataFrame(self.model_t_soil)
        self.model_t_soil.index = pd.to_datetime(x0, x1)
        
        self.idx0 = np.where(abs(self.model_depth - calib_depth) <= 0.0001)[0]
        
        sim = self.model_t_soil[x0:x1][idx0]
        
        comp_data = pd.concat((obs1, sim), axis=1)

        return comp_data

    def plot_obs_sim_data(self, ax,calib_period_onset, calib_period_end, calib_depth, ylim):

        calib_depth = float(calib_depth)
        comp_data = self.combine_obsim(calib_period_onset, calib_period_end, calib_depth)
        calib_depth = float(calib_depth)
        comp_data = self.confrontations_obs_sim_data(
            calib_period_onset, calib_period_end, calib_depth)

        fig, ax = plt.subplots(figsize=(6, 4))

        obs = comp_data[self.col]
        sim = comp_data[self.idx0]

        plt.plot(obs, 'r.', lw=0.8,label='Obs')
        plt.plot(sim, 'k-', lw=1,label='Sim')
        ax.legend(loc=0, shadow=True)
        print(f'{calib_depth:.2f}' + ' cm')
        plt.ylim(ylim)

    def plot_temperature_field(self, onset, ending, depth_bot, colors):

        year0 = onset[-4:]
        year1 = ending[-4:]

        if int(year0) != int(year1):
            date_label = year0+'-'+year1
        else:
            date_label = year0


        idx0 = np.where(self.model_date == date0)[0][0]
        idx1 = np.where(self.model_date == date1)[0][0]

        date_selected = self.model_date[onset:ending]

        depth_index = np.where(self.model_depth <= depth_bot)[0]

        depth_selected = self.model_depth.iloc[depth_index]

        t_soil_selected = self.model_t_soil[onset:ending].iloc[:, depth_index]

        t_max = np.nanmax(t_soil_selected)
        t_min = np.nanmin(t_soil_selected)

        max_abs = np.round(np.max(np.abs([t_max,t_min])))
        max_abs = max_abs - np.mod(max_abs, 5)

        levels = np.linspace(-1*max_abs, max_abs, num=11)

        colorm = LinearSegmentedColormap.from_list('NCL_Default', colors)

        fig2, ax2 = plt.subplots(figsize=(6, 3))
#        levels = [-15, -10, -5, -2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1, 1.5, 2.0, 2.5, 5, 10, 15]
        cs3 = plt.contourf(date_selected, depth_selected,
                           np.transpose(t_soil_selected),
                           levels,
                           extend='both', cmap=self.cm)
        cs0 = plt.contour(date_selected, depth_selected,
                          np.transpose(t_soil_selected), [-1, 0, 1], colors='k')

        plt.clabel(cs0, fmt='%d')

        months = MonthLocator(range(1, 13), bymonthday=1, interval=10)
        monthsFmt = DateFormatter("%Y-%m")
        ax2.xaxis.set_major_locator(months)
        ax2.xaxis.set_major_formatter(monthsFmt)
        for ticketlabel in ax2.xaxis.get_ticklabels():
            ticketlabel.set_rotation(30)
        ax2.autoscale_view()
        ax2.minorticks_on()
#        plt.xlabel('date',fontsize=16)
#        plt.xlabel(date_label,fontsize=11)

        plt.ylabel('Depth (m)', fontsize=11)
        plt.grid(axis='both')
        plt.xticks(fontsize=9)
        plt.yticks(fontsize=9)
        fig2.gca().invert_yaxis()
        font = {'family': 'serif',
                'color': 'black',
                'weight': 'normal',
                'size': 11,
                }

        cbar = plt.colorbar(
            cs3, label='Soil Temperature ($^\circ$C)', ticks=levels)
        cbar.set_label('Soil Temperature ($^\circ$C)', fontdict=font)

    def get_temp_profile(self, onset, ending):
        # =====Temperatures profile from onset-end at depths     
        profile_data = self.model_t_soil[onset:ending]

        profile_data = np.transpose(profile_data)
        profile      = np.zeros((np.size(profile_data, axis=0), 4))
        profile[:,0] = self.model_depth
        profile[:,1] = np.nanmin(profile_data, axis = 1)
        profile[:,2] = np.nanmin(profile_data, axis = 1)
        profile[:,3] = np.nanmin(profile_data, axis = 1) 
                    
        return profile

    def get_obstemp_profile(self, calib_depths, onset, ending):

        profile_data = self.obs_soilt[onset, ending].iloc[:,1:]
        profile_data = np.transpose(profile_data)
        profile      = np.zeros((np.size(profile_data, axis=0), 7))
        profile[:, 0] = calib_depths
        profile[:,1] = np.nanmin(profile_data, axis = 1)
        profile[:,2] = np.nanmean(profile_data, axis = 1)
        profile[:, 3] = np.nanmax(profile_data, axis=1)

        return profile

    def plot_boretemp_profile(self, ax, bore_tem, bore_dep, onset, ending, ylim):

        profile = self.get_temp_profile(onset, ending)

        bore = pd.DataFrame()
        bore['min'] = np.nanmin(bore_t, axis=1)
        bore['mean'] = np.nanmean(bore_t, axis=1)
        bore['max'] = np.nanmax(bore_t, axis=1)
        bore['dep'] = dep

        title = onset[-4:]+'-'+ending[-4:]

        plt.plot(profile[:, 1], profile[:, 0], color='g',
                 linestyle='solid', label='Min_sim')
        plt.plot(profile[:, 2], profile[:, 0], color='darkslateblue',
                 linestyle='solid',  label='Avg_sim')
        plt.plot(profile[:, 3], profile[:, 0],
                 color='darkorange', linestyle='solid', label='Max_sim')
        plt.gca().invert_yaxis()
        
        plt.plot(bore['min'], bore['dep'],  'b--', label='Min_obs')
        plt.plot(bore['mean'], bore['dep'], 'k--', label='Avg_obs')
        plt.plot(bore['max'], bore['dep'], 'r--', label='Max_obs')
        plt.plot([0, 0], [0, np.max(model_depth)], 'k-')

        plt.ylim(0, 16)
        plt.legend(shadow=False, loc=3, bbox_to_anchor=(
            0.46, 0.06), frameon=False)
        plt.xlabel('Soil Temperature ($^\circ$C)', fontsize=11)
        plt.xlim(-14, 14)
        plt.ylabel('Depth (m)', fontsize=11)
        plt.xticks(fontsize=9)
        plt.yticks(fontsize=9)
        plt.title(title, fontsize=11)
        ax.minorticks_on()
        #  Invert y axis
        plt.gca().invert_yaxis()

    def plot_obstemp_profile(self, ax, calib_depths, onset, ending, ylim):
        profile0 = self.get_temp_profile(onset, ending)
        profile1 = self.get_obstemp_profile(calib_depths, onset, ending)

        plt.plot(profile0[:,1], profile0[:,0], color='g',linestyle='solid', label='Min_sim')
        plt.plot(profile0[:,2], profile0[:,0],  color='darkslateblue',linestyle='solid', label='Avg_sim')
        plt.plot(profile0[:,3], profile0[:,0],color='darkorange',linestyle='solid', label='Max_sim')
        plt.gca().invert_yaxis()
        
        plt.plot(profile1[:,1], profile1[:,0], 'b--',label='Min_obs')
        plt.plot(profile1[:,2], profile1[:,0], 'k--',label='Avg_obs')
        plt.plot(profile1[:,3], profile1[:,0], 'r--',label='Max_obs')
        plt.plot([0,0], [0,np.max(self.model_depth)],'k-')               
        
        plt.ylim(ylim)
        plt.xlabel('Soil Temperature ($^\circ$C)',fontsize=11)
        plt.xlim(-14,14)
        plt.xticks(fontsize=9)
        plt.yticks(fontsize=9)
        plt.title(title,fontsize=11)
        ax.minorticks_on()
        plt.gca().invert_yaxis()      
        
    





sim_result = pd.read_csv(r'E:\Model\EB\self.model_t_soil.csv',header=None)
        
        
