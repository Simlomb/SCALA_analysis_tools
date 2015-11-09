#! /usr/bin/env python
# - Analysis of SCALA data by Simona Lombardo
# - 2015
# -*- coding: utf-8 -*-

import numpy              as N
import pyfits             as F
import calib_library as I
from scipy import integrate

from os import listdir
from os.path import isfile,isdir, join
import scipy.constants as SC
import CLAP_data_process as Cdp
from scipy import interpolate

class SnifsData:
    """
    Open SNIFS datacubes and saves data and header
    """
    def __init__(self,cubename):
        """
        Open SNIFS datacubes and saves data and header
        """
        self.data       = F.open(cubename)[0].data
        self.variance   = F.open(cubename)[1].data
        self.header     = F.open(cubename)[0].header

        self.channel    = self.header.get('CHANNEL')
        self.W_spix     = self.header.get('NAXIS1')
        self.N_spix     = self.header.get('NAXIS2')
        self.L_npix     = self.header.get('NAXIS3')
        
        self.W_step     = self.header.get('CDELT1')
        self.N_step     = self.header.get('CDELT2')
        self.L_step     = self.header.get('CDELT3')
        
        self.W_start    = self.header.get('CRVAL1')
        self.N_start    = self.header.get('CRVAL2')
        self.L_start    = self.header.get('CRVAL3')

        self.exptime    = self.header.get('EFFTIME')
        
        self.lbda       = N.arange(self.L_npix) * self.L_step + self.L_start
            


class SCALA_Calib:
    """
    """
    def __init__(self,list_SCALA, list_CLAP, throughput=True, clap_number=1,spaxel_sub=0):
        """
        """
        self.list_SCALA   = sorted(list_SCALA)
        self.list_CLAP    = sorted(list_CLAP)
        self.throughput   = throughput
	        
        self._clap_number = clap_number
        self._spaxel_sub  = spaxel_sub
        
        self.clap_files0  = [i for i in self.list_CLAP if i[2] == '0']
        self.clap_files1  = [i for i in self.list_CLAP if i[2] == '1']
        if self._clap_number == 0:
            self._clap_files = self.clap_files0
        else:
            self._clap_files = self.clap_files1       
	       
        self.scala_B_channel, self.scala_R_channel =[],[]
        for s in self.list_SCALA :
            if s.split(".")[0][-1] == 'B':
                self.scala_B_channel = N.append(self.scala_B_channel, s)
            else:
                self.scala_R_channel = N.append(self.scala_R_channel, s)
        self.Inter1      = self.interp_profile(1)
        self.snifs_data_B, self.snifs_data_R, self.clap_data, self.integrated_clap = [],[],[],[]
        print "First I load clap and snifs data and I fit all the CLap data"
        for i in range(len(self._clap_files)):
            self.snifs_data_B.append(SnifsData(self.scala_B_channel[i]))
            self.snifs_data_R.append(SnifsData(self.scala_R_channel[i]))
            if clap_number == 0:
                print "WARNING the next step will not work because you are not using the class Clap0_Data for clap 0"
                break
            else:
                self.clap_data.append(Cdp.Clap1_Data(self._clap_files[i]))
                self.clap_data[i].sig_back()
                self.integrated_clap.append(self.clap_data[i].light)
        print "Clap data fitting performed. Now I compute the throughput for each spaxel"
        self.A18         = I.Simulate(None).interpolate_mirror
        self.clap1_simul = I.Simulate(1).interpolate_array



    def interp_profile(self, arm_number):
        """
        perform the interpolation of the theoretical bandpass of SCALA
        for one of the arm the arm_number indicate wich one
        """
        if arm_number == 1:
            arm = N.loadtxt('arm4.txt')
        return interpolate.interp1d(arm[:,0],arm[:,1], kind='cubic')


    def Clean_array(self,Arr):
        """
        from useful_tools. 
        Look for N.NaN is Arr and remove them
        ---
        return array 
        """
        Arr_C = N.asarray(Arr).copy()
        FlagNaN = N.isnan(Arr_C) | N.isinf(Arr_C)
        return Arr_C[-FlagNaN]


    def _background_(self, delta, lbdacent, spxx, spxy,line,linewidth):
        """
        Compute the background taking the mean value of an interval between two lines
        observed in SNIFS.
        INPUT :
               delta   : half separation in index between two central wavelength
                         in SNIFS (usually is 250 A)
               lbdacent: the central wavelength from SNIFS
               spxx    : x position of the spaxel examined
               spxy    : y position of the spaxel examined
               line    : wavelength number observed (usually they are in tot
                         4 for B and 10 for R)
        """
        first_bkg      = lbdacent-delta[0]
        second_bkg     = lbdacent+delta[1]
        bkg_step1 = delta[0]-linewidth-5
        bkg_step2 = delta[1]-linewidth-5
        background = self.Clean_array(self.snifs_data.data[first_bkg-bkg_step1:first_bkg+bkg_step1, spxx,spxy])
        background = N.append(background,self.Clean_array(self.snifs_data.data[second_bkg-bkg_step2:second_bkg+bkg_step2,spxx,spxy]))
        mean_bkg = N.mean(background)
        return mean_bkg


    def Clap_info_and_light_level(self, file_number):
        """
        use the info contained in the previously run CLAP
        class which produce some objects with info
        about the data and run the method to compute the
        light level for every wavelength observed in the
        file considered

        INPUT: file_number: the number of the file used within the list
                            self.list_CLAP or self.list_SCALA

        OUTPUT: every object created from the class CLAP and the method
                sig_back
        """
        self.Clap       = self.clap_data[file_number]
        #self.Clap.sig_back()



    def weight_CLAP_Data(self, lbda_cent):
        """
        Compute the weighting function for the Clap data
        to assign a weight due to the particular shape for
        the interested fiber arm (which is always the same
        for CLAP 1). With this is possible to convert the
        CLAP light level from ADUs to W.
        
        INPUT:
                lbda_cent: the central wavelength as observed by SNIFS
                           in Angstrom
        OUTPUT:
                weight_funct: the weighting function to  apply to CLAP light
                              level to convert the ADU in W 
        """
        lbd_weight   = N.linspace(-46.,53.,1000)
        clap1_int    = N.mean([d(lbd_weight+lbda_cent) for d in self.clap1_simul], axis=0)
        arm1_int     = self.Inter1(lbd_weight)
        integ_weight = integrate.simps(arm1_int*clap1_int, lbd_weight+lbda_cent)
        weight_funct = integrate.simps(arm1_int, lbd_weight+lbda_cent)/integ_weight
        # We have to compute the error on that!!!
        return (weight_funct/43.5)


    def select_SNIFS_data_single_line(self, start, end, spxx, spxy):
        """
        Properly cut the SNIFS data around the line to examine
        INPUT:
              start: starting index from which the line starts
              end  : ending index where the line ends
              spxx : x position of the spaxel examined
              spxy : y position of the spaxel examined
        OUTPUT:
              lbda_line   : array containing the lbda belonging to the
                            wavelength line observed by SNIFS
              line_profile: array containing the data belonging to the
                            line observed by SNIFS
              variance    : array containing the variance of the
                            line_profile data
                  
        """
        lbda_line     = self.snifs_data.lbda[start:end]
        line_profile  = self.snifs_data.data[start:end, spxx, spxy]
        variance      = self.snifs_data.variance[start:end, spxx, spxy]
        return lbda_line, line_profile, variance
        
    def find_delta(self,line,len_Cl_data,spxx, spxy):
        """
        """
       
        if line == 0:
            delta1 = (self.Clap.lbda[line]-self.snifs_data.lbda[2])/(2*self.snifs_data.L_step)
            delta2 = (self.Clap.lbda[line+1]-self.Clap.lbda[line])/(2*self.snifs_data.L_step)
        elif line == len_Cl_data-1: 
            delta1 = (self.Clap.lbda[line]-self.Clap.lbda[line-1])/(2*self.snifs_data.L_step)
            delta2 = (self.snifs_data.lbda[-3]-self.Clap.lbda[line])/(2*self.snifs_data.L_step)
        else:
            delta1 = (self.Clap.lbda[line]-self.Clap.lbda[line-1])/(2*self.snifs_data.L_step)
            delta2 = (self.Clap.lbda[line+1]-self.Clap.lbda[line])/(2*self.snifs_data.L_step)
        return int(delta1), int(delta2)


        
    def select_condition_for_line_analysis(self,channel,line,len_Cl_data,index,delta,linewidth,spxx,spxy):
        """
        Select the first and last wavelength to cut the SNIFS data to examine a
        specific line. The selection depends on the channel and the line position
        respect to the wavelength range observed by SNIFS
        INPUT :
               channel    : SNIFS channel
               line       : wavelength number observed (usually they are in tot
                            4 for B and 10 for R)
               len_Cl_data: length of the subRun (4 for B and 10 for R)
               index      : position of the clap wavelength in SNIFS array
               delta      : half separation in index between two central wavelength
                            in SNIFS (usually is 250 A)
               linewidth  : the width of the SNIFS line observed in index unit
               spxx       :  x position of the spaxel examined
               spxy       : y position of the spaxel examined
        OUTPUT:
              index of the element which correspond to the start of the line observed
              by SNIFS, and index of the element which correspond to the end (start,end)
        """
        if channel == 'R' and len_Cl_data == 4:
            index_snifs = N.argmax(self.snifs_data.data[:index+delta[1],spxx,spxy])
            return 0, index_snifs+linewidth, index_snifs
        elif line == 0:
                index_snifs = N.argmax(self.snifs_data.data[:index+delta[1],spxx,spxy])
                if index_snifs < linewidth:
                    return 0, index_snifs+linewidth, index_snifs
                else:
                    return index_snifs-linewidth,index_snifs+linewidth, index_snifs
        else:
            #print self.snifs_data.data[index-delta[0],spxx,spxy], self.snifs_data.data[index+delta[1],spxx,spxy ],index, delta[0], delta[1]
            index_snifs = N.argmax(self.snifs_data.data[index-delta[0]:index+delta[1],spxx,spxy])+index-delta[0]
            return index_snifs-linewidth,index_snifs+linewidth, index_snifs


    def integ_snifs_line(self,obs_line,lbda,variance):
        """
        Compute the integral over the line observed by SNIFS using
        the simpson method (IS IT OK TO USE IT??) and the related
        error
        INPUT:
               obs_line: array data of the line observed by SNIFS
                         after they have been background subtsracted
                         the data must be in unit of energy (erg)
               lbda    : array wavelengths belonging to the line
               variance: variance of the data
        OUTPUT:
               int_profile: line profile integrated
               snifs_err  : err due to the integration
        
        """
        if self.throughput:
            line_profile = N.asarray(obs_line)*self.photon_energy(lbda)
        else:
            line_profile =  N.asarray(obs_line)
        int_profile = integrate.simps(line_profile,lbda)
        snifs_err   = (self.snifs_data.L_step*self.snifs_data.L_step)*(variance[0]+4*N.sum(variance[2:-2:2])+16*N.sum(variance[1::2])+(variance[-1]))/(9*int_profile*int_profile)
        return int_profile, snifs_err
    
    def photon_energy(self, lbda):
        """
        Compute the enrgy of photons
        INPUT :
               lbda : array of wavelengths for the photons
        OUTPUT:
               E : array of energy of photons
        """
        E =  (SC.h*SC.c)*10**(10)/lbda #photon energy
        return N.array(E)
        
    
    def convert_factor(self,lbda_cent):
        """
        INPUT:
                lbda_cent: the central wavelength as observed by SNIFS
                           in Angstrom
        OUTPUT:
                the correcting factor to consider the geometry and the shape
                of CLAP, SNIFS, and the SCALA beam:
                I HAVE TO WRITE IT DOWN HERE!!!!!
        """
        A18_value = N.mean([f(lbda_cent) for f in self.A18]) #Check the A18 value scaling factor!!!!
        return 93549.95727539062/A18_value

    
    def analysis_result(self,snifs_light,clap_light,geom_factor):
        """
        Compute the result of the analysis combining the light measured by CLAP
        with that measured by SNIFS and the related error
        INPUT:
               snifs_light: light level for one wavelength measured by SNIFS #[erg/A/s/cm**2/spaxel]
               clap_light : light level measured by CLAP #[W]
               geom_factor: converting factor for different dimension between
                            CLAP, SNIFS and the SCALA beam
        OUTPUT:
               Calib_value: value computed as result of the calibration analysis: if
                            the parameter throughput is set on True then it compute
                            the throughput estimation for a wavelength and a spaxel
                            if it is set on False, then it computes the comparison
                            between the SNIFS pipeline calib and the SCALA one for a
                            given wavelength and spaxel
                Tot_error : is the Calib_value error 
        """
        if self.throughput:
            Calib_value = snifs_light*geom_factor/clap_light
            #print snifs_light, clap_light
            #Tot_error   = N.sqrt(Clap_error +integ_clap_err+A18_error+snifs_err)*Calib_value
            ###Compute the error!!!!!
        else:
            Calib_value = snifs_light*geom_factor*((0.58**2)*(10**-7)/clap_light)*(self.snifs_data.exptime)
            #Tot_error   = N.sqrt(Clap_error +integ_clap_err+A18_error+snifs_err)*Calib_value

        return Calib_value#, Tot_error


    def channel_analysis(self, file_number, channel, line, spxx,spxy):
        """
        INPUT :
               channel    : SNIFS channel
               file_number: index of file used for analysis
               
        """
        if channel == 'B':
            self.snifs_data = self.snifs_data_B[file_number]
        else:
            self.snifs_data = self.snifs_data_R[file_number]
        #print channel, file_number,line
        lbdacent_clap = self.Clap.lbda[line]
        index         = int((lbdacent_clap-self.snifs_data.L_start)/self.snifs_data.L_step)
        delta         = self.find_delta(line,len(self.Clap.lbda), spxx,spxy)
        linewidth     = int(70./self.snifs_data.L_step)
        #print file_number, channel,line, spxx, spxy,lbdacent_clap
        #considering that we have datacubes for all the data in B and R
        #we need to distinguish between the ones with data or without
        if lbdacent_clap > self.snifs_data.lbda[-1-5] or lbdacent_clap < (self.snifs_data.L_start-4.):
            self.empty_line     = True
            self.lbdacent_snifs = lbdacent_clap
            self.lbda_line      = lbdacent_clap
        else:
            self.empty_line = False
            start, end, index_snifs = self.select_condition_for_line_analysis(channel,line,len(self.Clap.lbda),index,delta,linewidth,spxx,spxy)
            snifs_line              = self.select_SNIFS_data_single_line(start,end,spxx,spxy)
            self.lbdacent_snifs     = self.snifs_data.lbda[index_snifs]
            mean_bkg                = self._background_(delta,index_snifs,spxx,spxy,line,linewidth)
            self.line_bkg_sub       = snifs_line[1] - mean_bkg
            self.line_bkg_sub[N.isnan(self.line_bkg_sub)] = 0.
            variance_Snifs          = snifs_line[2]
            variance_Snifs[N.isnan(variance_Snifs)] = 0.
            
        if self.empty_line:
            self.line_int     = 0.
            self.line_profile = 0.
            Calibration = 0.
        else:
            weight_func = self.weight_CLAP_Data(self.lbdacent_snifs)
            clap_light  = self.integrated_clap[file_number][line][0]*weight_func
            if clap_light == 0.:
                Calibration = 0.
            else:
                geom_factor = self.convert_factor(self.lbdacent_snifs)
                snifs_light,err_snifs = self.integ_snifs_line(self.line_bkg_sub,snifs_line[0], variance_Snifs)
                Calibration = self.analysis_result(snifs_light,clap_light,geom_factor)

        return Calibration
