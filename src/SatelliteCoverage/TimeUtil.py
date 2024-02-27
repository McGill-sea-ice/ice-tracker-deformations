
import os
import sys
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime, date, time, timedelta
import subprocess
import numpy.ma as ma
import shutil

import pyproj

from scipy.interpolate import griddata
import matplotlib.path as mpath


class TimeUtil:
    def __init__(self,
                 TimeMeta=None, config = None):

        if TimeMeta != None:
            self.StartYear = TimeMeta.StartYear
            self.StartMonth =TimeMeta.StartMonth
            self.StartDay = TimeMeta.StartDay
            self.StartHour = TimeMeta.StartHour
            self.EndYear = TimeMeta.EndYear
            self.EndMonth = TimeMeta.EndMonth
            self.EndDay = TimeMeta.EndDay
            self.EndHour = TimeMeta.EndHour
            self.tstep = TimeMeta.tstep

        elif config != None:
            self.StartYear = int(config['start_year'])
            self.StartMonth = int(config['start_month'])
            self.StartDay = int(config['start_day'])
            self.StartHour = int(config['start_hour'])
            self.EndYear = int(config['end_year'])
            self.EndMonth = int(config['end_month'])
            self.EndDay = int(config['end_day'])
            self.EndHour = int(config['end_hour'])
            self.tstep = int(config['tstep'])


        self.StartDate = datetime(self.StartYear,self.StartMonth,self.StartDay,hour=self.StartHour)

        self.ThisTime = self.StartDate

        if self.EndYear == self.StartYear:
            self.EndDate = datetime(self.EndYear,self.EndMonth,self.EndDay,hour=self.EndHour)
        else:
            self.EndDate = datetime(self.StartYear,self.EndMonth,self.EndDay,hour=self.EndHour)

        self.NextTime = self.ThisTime + timedelta(seconds=self.tstep*60*60)

        self.LastTime = self.ThisTime - timedelta(seconds=self.tstep*60*60)

        self.ndays = int((self.EndDate - self.StartDate).days+1)

        self.nstep = int(int((self.EndDate - self.StartDate).days)*24.0/self.tstep + int((self.EndDate - self.StartDate).seconds/(self.tstep*60*60) +1))

        if self.EndYear != self.StartYear:
            self.nstep = int(self.nstep * (self.EndYear-self.StartYear))

        self.year = self.StartDate.strftime("%Y")
        self.dailyclock = 0
        self.StartDate.strftime("RCMS1SUD_%Y%m%d_000")
        ThisTime_str = self.ThisTime.strftime("%Y%m%d")
        NextTime_str = self.NextTime.strftime("%Y%m%d")
        self.ThisTimeFile = "RCMS1SID_%s_%s_dt72_tol72_dx.nc" % (ThisTime_str, NextTime_str)

        self.MLDP_init_date = self.StartDate - timedelta(days=1)
        self.NextTimeFile = self.NextTime.strftime("%Y%m%d%H_000")
        self.LastTimeFile = self.LastTime.strftime("%Y%m%d%H_000")

    def daterange(self, nmax = None):
#------------------------------------------------
#  Creates the list of dates in the given time period
#------------------------------------------------
        if nmax != None:
            for n in range(0, nmax):
                yield self.StartDate + timedelta(hours = n*self.tstep)
        else:
            for n in range(self.dailyclock, int(self.nstep)):
                yield self.StartDate + timedelta(hours = n*self.tstep)

    def step(self):
        self.dailyclock = self.dailyclock + 1

