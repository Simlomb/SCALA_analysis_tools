#! /usr/bin/env python
# -*- coding: utf-8 -*-
from scipy import interpolate
import numpy             as N
import random
import os
DATAPATH = os.path.realpath(__file__).split('calib_library.py')[0]+'refdata/'

p0     = N.loadtxt(DATAPATH+"CLAP0_cal_new.txt")
p1     = N.loadtxt(DATAPATH+"CLAP1_cal_new.txt")
random.seed()

class Simulate():

    
    def __init__(self,clap_number):

        self.clap_number = clap_number
        random.seed()
        if clap_number == None:
            self.interpolate_mirror = self.A18_calib()
        else:
            self.interpolate_array = self.Iterate_interp()
           
    def _randomize_me_(self,val,err,Npoints=1):
        """
        """
        return N.random.normal(val,err,Npoints)
    
    def Interp_clap0(self):
        
        p0_new  = N.zeros((len(p0),2))
        p0_new  = p0[:,1:]/10.**9.
        self.y0 = N.zeros((len(p0)))

        for i in range(len(p0)):
            self.y0[i] = self._randomize_me_(p0_new[i,0],p0_new[i,1])

        return interpolate.interp1d(p0[:,0], self.y0*10.**9., kind='linear')

         
    def Interp_clap1(self):

        p1_new  = N.zeros((len(p1),2))
        p1_new  = p1[:,1:]/10.**10.
        self.y1 = N.zeros((len(p1)))
        
        for i in range(len(p1)):
            self.y1[i] = self._randomize_me_(p1_new[i,0],p1_new[i,1])
    
        return interpolate.interp1d(p1[:,0], self.y1*10.**10., kind='linear')

        

    def Iterate_interp(self):

        self.interp_func0 = []
        self.interp_func1 = []
        if self.clap_number == 1:
            for i in range(100):
                self.interp_func1.append(self.Interp_clap1())

            return self.interp_func1
        else:
            for i in range(100):
                self.interp_func0.append(self.Interp_clap0())
            return self.interp_func0



    def A18_calib(self):

        MIRROR_CALIB = N.loadtxt("Calibration_mirror_ratio_new.txt")
        
        
        self.interp_func1 = []
        for j in range(100):
            self.m = N.zeros((len(MIRROR_CALIB)))
            for i in range(len(MIRROR_CALIB)):
                self.m[i] = self._randomize_me_(MIRROR_CALIB[i,1],MIRROR_CALIB[i,2])
            self.interp_func1.append(interpolate.interp1d(MIRROR_CALIB[:,0], MIRROR_CALIB[:,1], kind='linear',bounds_error=False, fill_value=0.))
        return self.interp_func1

