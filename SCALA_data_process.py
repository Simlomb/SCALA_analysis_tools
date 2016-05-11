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
    def __init__(self,list_SCALA, list_CLAP, analys_type, clap_number=1,spaxel_sub=0):
        """
        
        """
    
        self.list_SCALA = sorted(list_SCALA)
        self.list_CLAP = sorted(list_CLAP)
        self.analys_type = analys_type
        self._clap_number = clap_number
        self._spaxel_sub = spaxel_sub
        
        self.clap_files0  = [i for i in self.list_CLAP if i[2] == '0']
        self.clap_files1  = [i for i in self.list_CLAP if i[2] == '1']
        if self._clap_number == 0:
            self._clap_files = self.clap_files0
        else:
            self._clap_files = self.clap_files1       

         ## files to be used for the clap data cross talk correction
        cross_talk1 = N.loadtxt('crosstalk_clap1.txt')
        self.cross_talk0 = interpolate.interp1d(cross_talk1[:,0], cross_talk1[:,1], kind='linear') #light measurement from clap0
        self.cross_talk1 = interpolate.interp1d(cross_talk1[:,0], cross_talk1[:,3], kind='linear') # cross talk measurement from clap1
        
        self.scala_B_channel, self.scala_R_channel =[],[]
        for s in self.list_SCALA :
            if s.split(".")[0][-1] == 'B':
                self.scala_B_channel = N.append(self.scala_B_channel, s)
            else:
                self.scala_R_channel = N.append(self.scala_R_channel, s)
        self.Inter1, self.var_Inter1 = self.interpolate_arm_curve(4)
        self.snifs_data_B, self.snifs_data_R, self.clap_data, self.integrated_clap, self.integrated_clap2 = [],[],[],[],[]
        print "First I load clap and snifs data and I fit all the CLap data"
        for i in range(len(self._clap_files)):
            self.snifs_data_B.append(SnifsData(self.scala_B_channel[i]))
            self.snifs_data_R.append(SnifsData(self.scala_R_channel[i]))
            if clap_number == 0:
                print "WARNING the next step will not work because you are not using the class Clap0_Data for clap 0"
                break
            else:
                self.clap_data.append(Cdp.Clap1_Data(self._clap_files[i]))
                self.clap_data[i].sig_back(N.zeros((len(self.clap_data[i].lbda))),no_out=True)
                self.clap_data[i].light[N.isnan(self.clap_data[i].light)] = 0.
                self.clap_data[i].other_light[N.isnan(self.clap_data[i].other_light)] = 0.
                if self.clap_data[i].light[0,0] == 0.:
                    print " This file: %s will be skipped because of shutter failure" %self._clap_files[i]
                    self.integrated_clap.append(self.clap_data[i].light)
                    self.integrated_clap2.append(self.clap_data[i].other_light)
                else:
                    print "Applying cross talk correction to file %s" %self._clap_files[i]
                    clap0_class = Cdp.Clap0_Data(self._clap_files[i][0:2]+'0'+self._clap_files[i][3:])
                    clap0_light_level = clap0_class.sig_back(self.clap_data[i].mask_list, no_snifs=True)
                    cross_talk_correction = self.cross_talk1(self.clap_data[i].lbda)*(-clap0_light_level[0])/self.cross_talk0(self.clap_data[i].lbda)
                    self.clap_data[i].sig_back(cross_talk_correction,no_out=True)
                    self.integrated_clap.append(self.clap_data[i].light)
                    self.integrated_clap2.append(self.clap_data[i].other_light)
        #print N.shape(self.integrated_clap),  N.shape(self.integrated_clap2)
        
        print "Clap data fitting performed. Now I compute the throughput for each spaxel"
        self.A18         = I.Simulate(None).interpolate_mirror
        self.clap1_simul = I.Simulate(1).interpolate_array


    
    def interpolate_arm_curve(self,arm_number):
        """
        perform the interpolation of the theoretical bandpass of SCALA
        for one of the arm the arm_number and return the averaged intrerpolation
        curve computed on the lbd_weight datapoints
        """
        arm_data = N.loadtxt('arm%s_profile.txt' %arm_number)
        arm_data[:,0] = arm_data[:,0]*10.
        interpolation = []
        for j in range(100):
            m = N.zeros((len(arm_data[:,0])))
            for i in range(len(arm_data[:,0])):
                m[i] = N.random.normal(arm_data[i,1],arm_data[i,2],size=1)
            interpolation.append(interpolate.interp1d(arm_data[:,0], arm_data[:,1], kind='cubic'))

        lbd_weight   = N.linspace(0.,99.940, 1001)
        apply_interp = []
        apply_interp = [N.append(apply_interp, interpolation[i](lbd_weight)) for i in range(len(interpolation[:]))]
        mean_interp = N.mean(apply_interp, axis=0)
        var_interp = N.var(apply_interp, axis=0)
        return mean_interp, var_interp

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

    def clipped_array(self,array):
        """
        Apply sigma clipping to the data
        INPUT:
              array: array on which compute the sigma clipping
              sigma: sigma value of the clipping. Default:5
        OUTPUT:
              flagout: array with Boolean values to be applied
                       on the array to have the clipped data
        """
        array = N.asarray(array)
        maximum  = N.argmax(array)
        mask = int(35/self.snifs_data.L_step)
        if len(array) > 2*mask:
            new_array = N.append(array[:maximum-mask],array[maximum+mask+1:])
        else:
            new_array = array
        #print maximum,mask
        return new_array

    def find_nearest(self,array,value):
        """
        """
        idx = (N.abs(array-value)).argmin()
        return idx


    def _backgroundold_(self,snifs_data,channel,cube):
        """
        Compute the background taking the mean value of an interval between two lines
        observed in SNIFS.
        INPUT :
               delta   : half separation in index between two central wavelength
                         in SNIFS (usually is 250 A)
        """
        mask,mean_bkg,var_bkg = [],[],[]
        mask_back = N.array(([False]*len(self.snifs_data.lbda)))
        if cube == True :
            # integrate snifs data to have a better measurement of the maximum of light
            integrated_data = self.datacube_to_spectrum(self.snifs_data.data)
        else:
            integrated_data = N.copy(snifs_data)
        for i in range(len(self.Clap.lbda)):
            lbdacent_clap,index,delta,linewidth = self.set_up_analysis(i,channel)
            #print lbdacent_clap
            if lbdacent_clap > self.snifs_data.lbda[-1-5] or lbdacent_clap < (self.snifs_data.L_start-4.):
                self.full_line = N.append(self.full_line, False)
                self.lbdacent_snifs = N.append(self.lbdacent_snifs,lbdacent_clap)
                self.max_snifs = N.append(self.max_snifs, 0.0001)
                mask.append(mask_back)
                mean_bkg = N.append(mean_bkg,0.)
                var_bkg = N.append(var_bkg,0.)
                #print 'False'
            else:
                self.full_line = N.append(self.full_line, True)
                start,end,index_snifs = self.select_condition_for_line_analysis(channel,i,len(self.Clap.lbda),index,delta,linewidth,integrated_data)
                lbda_cent_mean = self.lbda_mean_snifs(index_snifs, integrated_data)
                self.lbdacent_snifs = N.append(self.lbdacent_snifs, lbda_cent_mean)
                self.max_snifs = N.append(self.max_snifs, integrated_data[index_snifs]/225.)
                mask_tmp = (self.snifs_data.lbda>=self.snifs_data.lbda[start])*(self.snifs_data.lbda<self.snifs_data.lbda[end+1])
                mask.append(mask_tmp)
                extra_mask = int(100/self.snifs_data.L_step)
                first_bkg = index_snifs-delta[0]
                second_bkg = index_snifs+delta[1]
                bkg_step1 = delta[0]-(index_snifs-start)-extra_mask
                bkg_step2 = delta[1]-(end-index_snifs)-extra_mask
                background1 = self.Clean_array(snifs_data[first_bkg-bkg_step1:first_bkg+bkg_step1])
                background2 = self.Clean_array(snifs_data[second_bkg-bkg_step2:second_bkg+bkg_step2])
                if len(background1) == 0 and len(background2) != 0:
                    background = self.clipped_array(background2)
                elif len(background2) == 0 and len(background1) != 0:
                    background = self.clipped_array(background1)
                elif len(background1) == 0 and len(background2) == 0:
                    bkg_step1 = delta[0]-(index_snifs-start)-extra_mask/3
                    bkg_step2 = delta[1]-(end-index_snifs)-extra_mask/3
                    background1 = self.Clean_array(snifs_data[first_bkg-bkg_step1:first_bkg+bkg_step1])
                    background2 = self.Clean_array(snifs_data[second_bkg-bkg_step2:second_bkg+bkg_step2])
                    if len(background1) == 0 and len(background2) != 0:
                        background = self.clipped_array(background2)
                    elif len(background2) == 0 and len(background1) != 0:
                        background = self.clipped_array(background1)
                else:
                    background = N.append(self.clipped_array(background1),self.clipped_array(background2))
                #we want to remove part of the high cross talk contribution
                #so we can reduce the bias on the background
                mean_bkg = N.append(mean_bkg,N.mean(background))
                var_bkg = N.append(var_bkg,N.var(background))
                #print i,first_bkg, bkg_step1,second_bkg,bkg_step2

        if N.any(self.full_line):
            return  mean_bkg,N.array(mask),var_bkg
        else:
            """ return nothing """
            return 0.,0.,0.


    def lbda_mean_snifs(self,index_snifs, data):
        """
        """
        
        lbd_array = self.snifs_data.lbda[index_snifs-5:index_snifs+6]
        clean_data = self.Clean_array(data[index_snifs-5:index_snifs+6])
        lbda_cent_mean = N.sum(lbd_array*clean_data)/N.sum(clean_data) 
        return lbda_cent_mean

    def datacube_to_spectrum(self,datacube):
        """
        Collapse a datacube into a spectrum
        """
        datacube[N.isnan(datacube)] = 0.
        spectrum = N.sum(datacube, axis=1)
        spectrum = N.sum(spectrum, axis=1)
        return spectrum
    
    def _background_(self,snifs_data, file_number, channel, cube=True):
        """
        Compute the background taking the mean value of the masked array of data
        where the mask excludes the lines observed in SNIFS.
        INPUT :
               snifs_data : SNIFS data (is a spectrum)
               channel    : SNIFS channel
               file_number: index of file used for analysis
               cube    : if the imput data is from a cube (True)
                         or only a spectrum (False)
        OUPUT :
               mean_bkg : mean background
               mask : mask for the lines it has the shape
                      (number of lines, snifs_data)    
        """
        #print channel,file_number,spxx,spxy
        mask = []
        mask_back = N.array(([False]*len(self.snifs_data.lbda)))
        if cube == True :
            # integrate snifs data to have a better measurement of the maximum of light
            self.snifs_data.data[N.isnan(self.snifs_data.data)] = 0.
            integrated_data = N.sum(self.snifs_data.data, axis=1)
            integrated_data = N.sum(integrated_data, axis=1)
        else:
            integrated_data = N.copy(snifs_data)
        for i in range(len(self.Clap.lbda)):
            #print i
            lbdacent_clap,index,delta,linewidth = self.set_up_analysis(i,channel)
            if lbdacent_clap > self.snifs_data.lbda[-1-5] or lbdacent_clap < (self.snifs_data.L_start-4.):
                self.full_line = N.append(self.full_line, False)
                self.lbdacent_snifs = N.append(self.lbdacent_snifs,lbdacent_clap)
                self.max_snifs = N.append(self.max_snifs, 0.0001)
                mask.append(mask_back)
            else:
                self.full_line = N.append(self.full_line, True)
                start,end,index_snifs = self.select_condition_for_line_analysis(channel,i,len(self.Clap.lbda),index,delta,linewidth,integrated_data)
                lbda_cent_mean = self.lbda_mean_snifs(index_snifs, integrated_data)
                self.lbdacent_snifs = N.append(self.lbdacent_snifs, lbda_cent_mean)
                self.max_snifs = N.append(self.max_snifs, integrated_data[index_snifs]/225.)
                #self.lbdacent_snifs = N.append(self.lbdacent_snifs, self.snifs_data.lbda[index_snifs])
                #if channel == 'R' and i == (len(self.Clap.lbda)-3):
                #    mask_tmp = (self.snifs_data.lbda>=self.snifs_data.lbda[start])
                #    mask.append((self.snifs_data.lbda>=self.snifs_data.lbda[start])*(self.snifs_data.lbda<self.snifs_data.lbda[end+1]))
                #    mask_back = N.ma.mask_or(mask_back,mask_tmp)
                #else:
                mask_tmp = (self.snifs_data.lbda>=self.snifs_data.lbda[start])*(self.snifs_data.lbda<self.snifs_data.lbda[end+1])
                mask.append(mask_tmp)
                extra_mask = int(100/self.snifs_data.L_step)
                if (self.snifs_data.lbda[start]-61.) <= self.snifs_data.lbda[0]:
                    mask_tmp1 = (self.snifs_data.lbda<self.snifs_data.lbda[end+1+extra_mask])
                    mask_back = N.ma.mask_or(mask_back,mask_tmp1)
                elif (self.snifs_data.lbda[end+1]+61.) >= self.snifs_data.lbda[-1]:
                    mask_tmp1 = (self.snifs_data.lbda>=self.snifs_data.lbda[start-extra_mask])
                    mask_back = N.ma.mask_or(mask_back,mask_tmp1)
                else:
                    mask_tmp1 = (self.snifs_data.lbda>=self.snifs_data.lbda[start-extra_mask])*(self.snifs_data.lbda<self.snifs_data.lbda[end+1+extra_mask])
                    mask_back = N.ma.mask_or(mask_back,mask_tmp1)
                    
        back_ground_mask = snifs_data[-mask_back]
        lbd_masked = self.snifs_data.lbda[-mask_back]
       
        # we have to check in case there are ghost lines and avoid to consider them in the background
        if self.analys_type == 't':
            if N.max(back_ground_mask) >= 200. and channel == 'B':
                index_ghost = N.argmax(back_ground_mask)
                if lbd_masked[index_ghost] > self.snifs_data.lbda[-1-40]:
                    mask_tmp = (lbd_masked>=lbd_masked[index_ghost-40])
                else:
                    mask_tmp = (lbd_masked>=lbd_masked[index_ghost-40])*(lbd_masked<=lbd_masked[index_ghost+35])
                back_ground_mask = back_ground_mask[-mask_tmp]
                lbd_masked = lbd_masked[-mask_tmp]
            if channel == 'R' and N.max(back_ground_mask) >= 2*N.mean(back_ground_mask) and N.any(self.full_line):
                index_second_order = N.argmax(back_ground_mask)
                if lbd_masked[index_second_order] > 9000.:
                    if lbd_masked[index_second_order] > self.lbdacent_snifs[-1]:
                        mask_tmp = (lbd_masked>=lbd_masked[index_second_order-40])
                    else:
                        mask_tmp = (lbd_masked>=lbd_masked[index_second_order-40])*(lbd_masked<=lbd_masked[index_second_order+40])
                    back_ground_mask = back_ground_mask[-mask_tmp]
                    lbd_masked = lbd_masked[-mask_tmp]
        else:
            if N.max(back_ground_mask) >= N.abs(N.mean(snifs_data[-20:])*200.) and channel == 'B':
                index_ghost = N.argmax(back_ground_mask)
                if lbd_masked[index_ghost] > self.snifs_data.lbda[-1-40]:
                    mask_tmp = (lbd_masked>=lbd_masked[index_ghost-40])
                else:
                    mask_tmp = (lbd_masked>=lbd_masked[index_ghost-40])*(lbd_masked<=lbd_masked[index_ghost+35])
                back_ground_mask = back_ground_mask[-mask_tmp]
                lbd_masked = lbd_masked[-mask_tmp]
        if N.any(self.full_line):
            clean_background = self.Clean_array(back_ground_mask)
            mean_bkg = N.median(clean_background*1)
            #print len(clean_background), mean_bkg
            return mean_bkg,N.array(mask),N.var(clean_background), back_ground_mask,lbd_masked
        else:
            """ return nothing """
            return 0.,0.,0.
                


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
        lbd_weight = N.linspace(0.,99.940, 1001)
        lbd_line = N.linspace(lbd_weight[0]-lbd_weight[len(lbd_weight)/2],lbd_weight[-1]-lbd_weight[len(lbd_weight)/2],len(lbd_weight))+lbda_cent
        clap1_int = N.mean([d(lbd_line) for d in self.clap1_simul], axis=0)
        err_clap1_int = N.var([d(lbd_line) for d in self.clap1_simul], axis=0)
        integ_weight = integrate.simps(self.Inter1*clap1_int, lbd_line)
        sqr_Inter1 = self.Inter1**4
        sqr_clap1 = clap1_int**2
        Int_func = integrate.simps(self.Inter1, lbd_line)
        Int_sqr = integ_weight**2
        
        weight_funct = Int_func/integ_weight
        weight_funct_err = ((self.var_Inter1[0]*(Int_sqr+sqr_clap1[0]-2*integ_weight*clap1_int[0])+(sqr_Inter1[0]*err_clap1_int[0]))\
                            +4*N.sum((self.var_Inter1[2:-2:2]*(Int_sqr+sqr_clap1[2:-2:2]-2*integ_weight*clap1_int[2:-2:2])\
                            +(sqr_Inter1[2:-2:2]*err_clap1_int[2:-2:2])))+16*N.sum((self.var_Inter1[1:-1:2]*(Int_sqr\
                            +sqr_clap1[1:-1:2]-2*integ_weight*clap1_int[1:-1:2])+(sqr_Inter1[1:-1:2]*err_clap1_int[1:-1:2])))\
                            +((self.var_Inter1[-1]*(Int_sqr+sqr_clap1[-1]-2*integ_weight*clap1_int[-1])\
                            +(sqr_Inter1[-1]*err_clap1_int[-1]))))*(99.940/1001.)**2/(1892.25*9.*Int_sqr**2)
        # We have to compute the error on that!!!
        return (weight_funct/43.5),weight_funct_err


    def select_SNIFS_data_single_line(self, start, end):
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
        line_profile  = self.snifs_data.data[start:end]
        variance      = self.snifs_data.variance[start:end]
        return lbda_line, line_profile, variance
        
    def find_delta(self,line,len_Cl_data,channel):
        """
        """
       
        if line == 0:
            if channel == 'B' and len_Cl_data > 5 :
                delta1 = (self.Clap.lbda[line+1]-self.Clap.lbda[line])/(2*self.snifs_data.L_step)
                delta2 = (self.snifs_data.lbda[-3]-self.Clap.lbda[line])/(2*self.snifs_data.L_step)
            else:
                delta1 = (self.Clap.lbda[line]-self.snifs_data.lbda[2])/(2*self.snifs_data.L_step)
                delta2 = (self.Clap.lbda[line+1]-self.Clap.lbda[line])/(2*self.snifs_data.L_step)
        elif line == len_Cl_data-1: 
            delta1 = (self.Clap.lbda[line]-self.Clap.lbda[line-1])/(2*self.snifs_data.L_step)
            delta2 = (self.snifs_data.lbda[-3]-self.Clap.lbda[line])/(2*self.snifs_data.L_step)
        else:
            delta1 = (self.Clap.lbda[line]-self.Clap.lbda[line-1])/(2*self.snifs_data.L_step)
            delta2 = (self.Clap.lbda[line+1]-self.Clap.lbda[line])/(2*self.snifs_data.L_step)
        return int(delta1), int(delta2)


        
    def select_condition_for_line_analysis(self,channel,line,len_Cl_data,index,delta,linewidth,integrated_data):
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
        OUTPUT:
              index of the element which correspond to the start of the line observed
              by SNIFS, and index of the element which correspond to the end (start,end)
        """
        
        if channel == 'R' and len_Cl_data == 4:
            index_snifs = N.argmax(integrated_data[:index+delta[1]])
            return 0, index_snifs+linewidth, index_snifs
        elif line == 0:
                index_snifs = N.argmax(integrated_data[:index+delta[1]])
                if index_snifs < linewidth:
                    return 0, index_snifs+linewidth, index_snifs
                elif channel == 'B' and len_Cl_data > 5:
                    index_snifs = N.argmax(integrated_data[index-delta[0]:index+delta[1]])+index-delta[0]
                    return index_snifs-linewidth,N.min((index_snifs+linewidth,len(self.snifs_data.lbda)-2)), index_snifs
                else:
                    return index_snifs-linewidth,index_snifs+linewidth, index_snifs
        else:
            #print self.snifs_data.data[index-delta[0],spxx,spxy], self.snifs_data.data[index+delta[1],spxx,spxy ],index, delta[0], delta[1]
            index_snifs = N.argmax(integrated_data[index-delta[0]:index+delta[1]])+index-delta[0]
            if index_snifs+linewidth > len(self.snifs_data.lbda)-2:
                return index_snifs-linewidth,N.min((index_snifs+linewidth,len(self.snifs_data.lbda)-2)), index_snifs
            else:
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
        if self.analys_type == 't':
            energy = self.photon_energy(lbda)
            line_profile = N.asarray(obs_line)*energy
        else:
            #line_profile =  N.asarray(obs_line)
            line_profile = N.asarray(obs_line)
        int_profile = integrate.simps(line_profile,lbda)
        if self.analys_type == 't':
            sqrt_energy = energy**2
            snifs_err = (self.snifs_data.L_step**2)*((variance[0]*sqrt_energy[0])+4*N.sum(variance[2:-2:2]*sqrt_energy[2:-2:2])+16*N.sum(variance[1:-1:2]*sqrt_energy[1:-1:2])+(variance[-1]*sqrt_energy[-1]))/9.
        else:
            snifs_err   = (self.snifs_data.L_step*self.snifs_data.L_step)*(variance[0]+4*N.sum(variance[2:-2:2])+16*N.sum(variance[1:-1:2])+(variance[-1]))/9.
        return int_profile, snifs_err
    
    def photon_energy(self, lbda):
        """
        Compute the enrgy of photons
        INPUT :
               lbda : array of wavelengths for the photons in A
        OUTPUT:
               E : array of energy of photons in J
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
        Parameters:

        D_p=2.22   		# Diameter of primary mirror in m
        D_s=0.625		# Diameter of secondary mirror in m
        D_o=0.86		# Central obstruction in m
        l_pd=0.0058		# Length of the photodiode active area in m
        d_pl=1			# diameter of the artificila planet produced by SCALA in deg
        D_fov=6.4		# filed of view of SNIFS in arcsec
        N_spaxel=15		# 15x15 spaxels fill field of view
        D_m_SCALA=0.2	# diameter of the scala mirrors in m
        D_mask=0.16		# diameter of the holes in SCALA mask in m
        #collecting area of the telescope
        A_tel=(D_p**2-D_o**2)*cons.pi/4.
        #collecting area of the photodiode
        A_pd=l_pd**2
        #Area illuminated by SCALA
        A_SCALA=D_m_SCALA**2*cons.pi/4.
        #Area illuminated with SCALA mask
        A_SCALA=D_mask**2*cons.pi/4.
        # solid angle illuminated by SCALA
        sa_SCALA=(d_pl*cons.pi/180)**2*cons.pi/4.
        # solid angle of one spaxel (15x15 Spaxels)
        sa_spaxel=((D_fov/N_spaxel)*cons.pi/(180.*3600.))**2
        #scaling Factor S
        S=A_pd*sa_SCALA/(A_SCALA*sa_spaxel) # 93549.95727539062
        diameters=(221.8, 66.0)
        """
        A18_value = N.mean([f(lbda_cent) for f in self.A18]) #Check the A18 value scaling factor!!!!
        var_A18 = N.var([f(lbda_cent) for f in self.A18])
        if self.analys_type == 't':
            geom_fact_var = (93549.95727539062**2/A18_value**4)*var_A18
            return 93549.95727539062/A18_value,geom_fact_var,A18_value,var_A18
        else:
            # factor*(222**2 - 66**2)*pi/4*10**(-7)= 330.103828289671
            geom_fact_var = (330.103828289671**2/A18_value**4)*var_A18
            return 330.103828289671/A18_value,geom_fact_var,A18_value,var_A18

    
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
        
        snifs, snifs_err = snifs_light
        clap, clap_var = clap_light
        factor, factor_var = geom_factor
        #print snifs, clap, factor,self.snifs_data.exptime
        if self.analys_type == 't':
            Calib_value = snifs*factor/clap
            #print snifs_light, clap_light
            #Tot_error   = N.sqrt(Clap_error +integ_clap_err+A18_error+snifs_err)*Calib_value
            ###Compute the error!!!!!
            Tot_error = (snifs_err/(snifs)**2+factor_var/(factor)**2+clap_var/(clap)**2)*Calib_value**2
            #print snifs_err/(snifs)**2,factor_var/(factor)**2,clap_var/(clap)**2, N.sqrt(Tot_error)/Calib_value
        else:
            Calib_value = snifs*factor*self.snifs_data.exptime/clap
            Tot_error = (snifs_err/(snifs)**2+factor_var/(factor)**2+clap_var/(clap)**2+(0.1/self.snifs_data.exptime)**2)*Calib_value**2
        #print Calib_value
        return Calib_value, Tot_error



    def set_up_analysis(self, line,channel):
        """
        Compute some fundamental quantities that are used in the later step
        """
        lbdacent_clap = self.Clap.lbda[line]
        index         = int((lbdacent_clap-self.snifs_data.L_start)/self.snifs_data.L_step)
        delta         = self.find_delta(line,len(self.Clap.lbda),channel)
        linewidth     = int(((110./self.snifs_data.L_step)+1)/2.)
        return lbdacent_clap,index,delta,linewidth


    def select_snifs_FWHM(self,obs_line,lbda,variance):
        """
        """
        snifs_line_func = interpolate.interp1d(lbda,obs_line, kind='cubic', bounds_error=False, fill_value=0.)
        new_lbda = N.linspace(lbda[0], lbda[-1],1000)
        snifs_line_max = N.max(snifs_line_func(new_lbda))
        FWHM = snifs_line_max/2.
        lbda_start = new_lbda[self.find_nearest(snifs_line_func(new_lbda[:len(new_lbda)/2]),FWHM)]
        lbda_stop = new_lbda[len(new_lbda)/2:][self.find_nearest(snifs_line_func(new_lbda[len(new_lbda)/2:]),FWHM)]
        lbd_width = lbda_stop-lbda_start
        Full_line = lbd_width/2.
        #lbda_start_full = new_lbda[self.find_nearest(snifs_line_func(new_lbda[:len(new_lbda)/2]),Full_line)]
        #lbda_stop_full = new_lbda[len(new_lbda)/2:][self.find_nearest(snifs_line_func(new_lbda[len(new_lbda)/2:]),Full_line)]
        #lbd_fwhm = N.linspace(lbda_start-Full_line, lbda_stop+Full_line,1000)
        lbd_fwhm = N.linspace(lbda_start, lbda_stop,500)
        snifs_line_fwhm = snifs_line_func(lbd_fwhm)
        var_line_func = interpolate.interp1d(lbda,variance, kind='linear', bounds_error=False, fill_value=0.)
        var_line_fwhm = var_line_func(lbd_fwhm)
        return lbd_fwhm,snifs_line_fwhm,var_line_fwhm

class Cube_analysis(SCALA_Calib):
    """
    This class is a child of SCALA_Calib
    and it is meant to run the analysis on
    datacubes of SCALA
    """
    
    def channel_analysis(self, file_number, channel, spxx,spxy,snifs_fwhm=True,back_old=True):
        """
        INPUT :
               channel    : SNIFS channel
               file_number: index of file used for analysis
               
        """
        if channel == 'B':
            self.snifs_data = self.snifs_data_B[file_number]
        else:
            self.snifs_data = self.snifs_data_R[file_number]
        self.lbdacent_snifs, Calibration, Calibration2, self.full_line,error_Cal,error_Cal2 = [],[],[],[],[],[]
        snifs_int,err_snifs_int,tot_clap_light,err_tot_clap_light,int_clap_try,err_int_clap_try = [],[],[],[],[],[]
        if back_old:
            mean_bkg,mask,var_bkg = self._backgroundold_(self.snifs_data.data[:,spxx,spxy],channel,cube=True)
        else:
            mean_bkg,mask,var_bkg = self._background_(self.snifs_data.data[:,spxx,spxy],file_number, channel,cube=True)
        snifs_time = N.ones((len(self.Clap.lbda)))*self.snifs_data.exptime
        for line in range(len(self.Clap.lbda)):
            #considering that we have datacubes for all the data in B and R
            #we need to distinguish between the ones with data or without
            if self.full_line[line] == False:
                Calibration = N.append(Calibration,0.)
                error_Cal = N.append(error_Cal,0.)
                Calibration2 = N.append(Calibration2,0.)
                error_Cal2 = N.append(error_Cal2,0.)
                # for test purpose
                snifs_int = N.append(snifs_int,0.)
                err_snifs_int = N.append(err_snifs_int,0.)
                tot_clap_light = N.append(tot_clap_light, 0.)
                err_tot_clap_light = N.append(err_tot_clap_light,0.)
                int_clap_try = N.append(int_clap_try,0.)
                err_int_clap_try = N.append(err_int_clap_try,0.)
            else:
                #return mean_bkg, mask[1]
                variance_Snifs = self.snifs_data.variance[mask[line], spxx, spxy] # variance corresponding to the snifs line
                variance_Snifs[N.isnan(variance_Snifs)] = 0.
                if back_old:
                    line_bkg_sub = self.snifs_data.data[mask[line], spxx, spxy] - mean_bkg[line]
                    var_snifs = variance_Snifs+var_bkg[line]
                else:
                    line_bkg_sub = (self.snifs_data.data[mask[line], spxx, spxy]) - mean_bkg
                    var_snifs = variance_Snifs+var_bkg
                line_bkg_sub[N.isnan(line_bkg_sub)] = 0.
                weight_func,weight_func_err = self.weight_CLAP_Data(self.lbdacent_snifs[line])
                # clap_light is the light measure by clap computed integrating the fit function to the data,
                # clap_light2 is the light measured by clap computed the data directly with simpsons method
                clap_light  = self.integrated_clap[file_number][line][0]*weight_func
                var_Clap_light = (self.integrated_clap[file_number][line][1]*weight_func**2)+(weight_func_err*self.integrated_clap[file_number][line][0]**2)
                clap_light2 = self.integrated_clap2[file_number][line][0]*weight_func
                var_Clap_light2 = (self.integrated_clap2[file_number][line][1]*weight_func**2)+(weight_func_err*self.integrated_clap2[file_number][line][0]**2)
                int_clap_try = N.append(int_clap_try,self.integrated_clap2[file_number][line][0])
                err_int_clap_try = N.append(err_int_clap_try,self.integrated_clap2[file_number][line][1])
                if clap_light == 0.:
                    Calibration = N.append(Calibration,0.)
                    error_Cal = N.append(error_Cal,0.)
                    Calibration2 = N.append(Calibration2,0.)
                    error_Cal2 = N.append(error_Cal2,0.)
                    #for test purpose
                    snifs_int = N.append(snifs_int,0.)
                    err_snifs_int = N.append(err_snifs_int,0.)
                    tot_clap_light = N.append(tot_clap_light, 0.)
                    err_tot_clap_light = N.append(err_tot_clap_light,0.)
                else:
                    geom_factor,geom_factor_var,A18,var_A18 = self.convert_factor(self.lbdacent_snifs[line])
                    tot_clap_light = N.append(tot_clap_light, clap_light2*A18)
                    err_tot_clap_light = N.append(err_tot_clap_light,var_Clap_light2*A18 + var_A18*clap_light2)
                    if snifs_fwhm:
                        lbd_fwhm,snifs_line_fwhm,var_line_fwhm = self.select_snifs_FWHM(line_bkg_sub,self.snifs_data.lbda[mask[line]], var_snifs)
                        snifs_light,err_snifs = self.integ_snifs_line(snifs_line_fwhm,lbd_fwhm, var_line_fwhm)
                    else:
                        snifs_light,err_snifs = self.integ_snifs_line(line_bkg_sub,self.snifs_data.lbda[mask[line]], var_snifs)
                    cal, err_cal = self.analysis_result((snifs_light,err_snifs),(clap_light,var_Clap_light),(geom_factor,geom_factor_var))
                    cal2, err_cal2 = self.analysis_result((snifs_light,err_snifs),(clap_light2,var_Clap_light2),(geom_factor,geom_factor_var))
                    Calibration = N.append(Calibration,cal)
                    error_Cal = N.append(error_Cal,err_cal)
                    Calibration2 = N.append(Calibration2,cal2)
                    error_Cal2 = N.append(error_Cal2,err_cal2)
                    snifs_int = N.append(snifs_int,snifs_light)
                    err_snifs_int = N.append(err_snifs_int,err_snifs)
                    #if spxx == 7 and spxy == 7: 
                    #    print N.sqrt(err_snifs)/snifs_light,N.sqrt(var_Clap_light2)/clap_light2,N.sqrt(geom_factor_var)/geom_factor
                    #    print N.sqrt(err_cal2)/cal2,cal2
        
        return Calibration,error_Cal,Calibration2,error_Cal2,snifs_int,err_snifs_int,int_clap_try,err_int_clap_try,tot_clap_light,err_tot_clap_light,snifs_time

class Spectrum_analysis(SCALA_Calib):
    """
    This class is a child of SCALA_Calib
    and it is meant to run the analysis on
    spectra from collapsed datacubes of SCALA
    """

    def channel_analysis(self, file_number, channel):
        """
        your stuff
        """
        #print channel
        if channel == 'B':
            self.snifs_data = self.snifs_data_B[file_number]
        else:
            self.snifs_data = self.snifs_data_R[file_number]
        integrated_data = self.datacube_to_spectrum(self.snifs_data.data)  
        variance_Snifs = self.datacube_to_spectrum(self.snifs_data.variance) # variance corresponding to the snifs 
        self.lbdacent_snifs,self.max_snifs, Calibration, Calibration2, self.full_line,error_Cal,error_Cal2 = [],[],[],[],[],[],[]
        snifs_int,err_snifs_int,tot_clap_light,err_tot_clap_light,int_clap_try,err_int_clap_try = [],[],[],[],[],[]
        mean_bkg,mask,var_bkg= self._backgroundold_(integrated_data,channel,cube=False)
        snifs_time = N.ones((len(self.Clap.lbda)))*self.snifs_data.exptime
        max_snifs = self.max_snifs*snifs_time/self.clap_data[file_number].sub_exp
        for line in range(len(self.Clap.lbda)):
            #considering that we have datacubes for all the data in B and R
            #we need to distinguish between the ones with data or without
            #print line
            if self.full_line[line] == False:
                Calibration = N.append(Calibration,0.)
                error_Cal = N.append(error_Cal,0.)
                Calibration2 = N.append(Calibration2,0.)
                error_Cal2 = N.append(error_Cal2,0.)
                # for test purpose
                snifs_int = N.append(snifs_int,0.)
                err_snifs_int = N.append(err_snifs_int,0.)
                tot_clap_light = N.append(tot_clap_light, 0.)
                err_tot_clap_light = N.append(err_tot_clap_light,0.)
                int_clap_try = N.append(int_clap_try,0.)
                err_int_clap_try = N.append(err_int_clap_try,0.)
            else:
                #return mean_bkg, mask[1]
                line_bkg_sub = integrated_data[mask[line]] - mean_bkg[line]
                line_bkg_sub[N.isnan(line_bkg_sub)] = 0.
                #print N.shape(line_bkg_sub)
                variance_Snifs_line = variance_Snifs[mask[line]] # variance corresponding to the snifs line
                variance_Snifs_line[N.isnan(variance_Snifs_line)] = 0.
                var_snifs = variance_Snifs_line+var_bkg[line]
                #print N.shape(var_snifs)
                weight_func,weight_func_err = self.weight_CLAP_Data(self.lbdacent_snifs[line])
                # clap_light is the light measure by clap computed integrating the fit function to the data,
                # clap_light2 is the light measured by clap computed the data directly with simpsons method
                clap_light  = self.integrated_clap[file_number][line][0]*weight_func
                var_Clap_light = (self.integrated_clap[file_number][line][1]*weight_func**2)+(weight_func_err*self.integrated_clap[file_number][line][0]**2)
                clap_light2 = self.integrated_clap2[file_number][line][0]*weight_func
                var_Clap_light2 = (self.integrated_clap2[file_number][line][1]*weight_func**2)+(weight_func_err*self.integrated_clap2[file_number][line][0]**2)
                int_clap_try = N.append(int_clap_try,self.integrated_clap2[file_number][line][0])
                err_int_clap_try = N.append(err_int_clap_try,self.integrated_clap2[file_number][line][1])
                if clap_light == 0.:
                    Calibration = N.append(Calibration,0.)
                    error_Cal = N.append(error_Cal,0.)
                    Calibration2 = N.append(Calibration2,0.)
                    error_Cal2 = N.append(error_Cal2,0.)
                    #for test purpose
                    snifs_int = N.append(snifs_int,0.)
                    err_snifs_int = N.append(err_snifs_int,0.)
                    tot_clap_light = N.append(tot_clap_light, 0.)
                    err_tot_clap_light = N.append(err_tot_clap_light,0.)
                else:
                    geom_factor,geom_factor_var,A18,var_A18 = self.convert_factor(self.lbdacent_snifs[line])
                    tot_clap_light = N.append(tot_clap_light, clap_light2*A18)
                    err_tot_clap_light = N.append(err_tot_clap_light,var_Clap_light2*A18 + var_A18*clap_light2)
                    #print N.shape(self.snifs_data.lbda[mask[line]])
                    lbd_fwhm,snifs_line_fwhm,var_line_fwhm = self.select_snifs_FWHM(line_bkg_sub,self.snifs_data.lbda[mask[line]], var_snifs)
                    snifs_light,err_snifs = self.integ_snifs_line(snifs_line_fwhm,lbd_fwhm, var_line_fwhm)
                    cal, err_cal = self.analysis_result((snifs_light,err_snifs),(clap_light,var_Clap_light),(geom_factor,geom_factor_var))
                    cal2, err_cal2 = self.analysis_result((snifs_light,err_snifs),(clap_light2,var_Clap_light2),(geom_factor,geom_factor_var))
                    Calibration = N.append(Calibration,cal)
                    error_Cal = N.append(error_Cal,err_cal)
                    Calibration2 = N.append(Calibration2,cal2)
                    error_Cal2 = N.append(error_Cal2,err_cal2)
                    snifs_int = N.append(snifs_int,snifs_light)
                    err_snifs_int = N.append(err_snifs_int,err_snifs)
        
        return Calibration,error_Cal,Calibration2,error_Cal2,snifs_int,err_snifs_int,int_clap_try,err_int_clap_try,tot_clap_light,err_tot_clap_light,snifs_time,max_snifs
            
