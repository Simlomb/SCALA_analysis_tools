#! /usr/bin/env python
# - Analysis of CLAPs data by Simona Lombardo
# - 2015
# -*- coding: utf-8 -*-

import numpy                    as N
import pyfits                   as F
import scipy.special.orthogonal as OP
from scipy import integrate
from scipy.special import erf
from scipy.optimize import leastsq


        
class Clap:
    """
    """
    
    def __init__(self,clap_file):
        """
        Open CLAPs file and save data info
        INPUT :
               clap_file: file name of clap data
        """
        
        self._Fits      = F.open(clap_file)
        self._Primary   = self._Fits[0]
        self._Tables    = self._Fits[1]
        
        self.lbda       = self._Tables.data["LBDA"][:]
        self.expotime   = self._Tables.data["EXPTIME"][:]
        self.time_start = self._Tables.data["TIMESTAR"][:]
        self.time_on    = self._Tables.data["TIMEON"][:]
        self.time_off   = self._Tables.data["TIMEOFF"][:]
        self.time_stop  = self._Tables.data["TIMESTOP"][:]
        self.frequency  = self._Primary.header["CLAPFREQ"]*1000.
        self.data       = self._Primary.data



    def clipped_array(self,array,std_=None,mean_=None,sigma=5):
        """
        Apply sigma clipping to the data
        INPUT:
              array: array on which compute the sigma clipping
              std_ : std_ of the array data. Defoult: None
              mean_: mean of the array data. Default: None
              sigma: sigma value of the clipping. Default:5
        OUTPUT:
              flagout: array with Boolean values to be applied
                       on the array to have the clipped data
        """
        array = N.asarray(array)
        if std_ == None:
            std_  = array.std()
        if mean_ == None:
            mean_ = array.mean()
            
        flagout = N.abs(array-mean_)>sigma*std_
        #array_clipped = array[-flagout]
        return flagout


    def residuals(self,p, y, x):
        """
        the function that describes CLAPs data has to be minimized using this new function
        the parameters are:
        1) a is the amplitude of the signal
        2) b is the steepness of the curve from the background to the signal
        3) w is the half of the width of the top hat
        4) xe is the center of the top hat
        5) back is the background level
        INPUT:
               p : array containing the values for the fitting function
                   (b,a,w,xe,back)
               y : array of data to fit
               x : array of corresponding x values to the y
        OUTPUT:
               err: residuals of the fit
        """
        
        b,a,w,xe,back = p
        err = y - ((a/4.)*(erf(b*(w+x-xe))+1)*(erf(b*(w-x+xe))+1)+back)
        return err

    def data_shape(self,x,p):
        """
        Function describing the shape of the CLAPs data
        INPUT :
               x : array of the x axis of the data
               p : array containing the parameters
                   describing the function. p is the same as
                   described in residuals
        OUTPUT:
               the y array computed according to the function
        """

        return ((p[1]/4.)*(erf(p[0]*(p[2]+x-p[3]))+1)*(erf(p[0]*(p[2]-x+p[3]))+1))

    def data_integral(self, x,y,data_clipped, sub_exp, x0, back):
        """
        It fits the data using the residuals function and integrate
        over the area covered by that function
        INPUT :
               x           : array of x
               y           : array of data
               data_clipped: light level (signal amplitude)
               sub_expo    : exposure time for one wavelength
               x0          : time of beginning of the exposure
               back        : background level of the exposure
        OUTPUT:
               light       : integrated value of the light
                             exposure for one wavelength
               err         : error on the light parameter
        """
        # first guess of fit parameters
        slope = data_clipped/(223.591-223.584) # the b parameter
        w     = (sub_exp/2.)+0.0025            # the w parameter
        xe    = w - 0.0025 + x0     # the xe parameter
        p     = [slope,data_clipped,w,xe,0.]
        # the fit of the data
        fit  = leastsq(self.residuals, p, args=(-(y-back), x))
        # the integral of the total amount of light produced by SCALA
        light, err = integrate.quad(self.data_shape, x[0], x[-1], args=(fit[0]))
        return light, err



class Clap1_Data(Clap):
    """
    This Class handle the data from Clap1 and compute the light level
    """
            
    
    def sig_back(self, no_out=True):
        """
        Compute the signal and background level of CLAP1 data
        and also compute the overall amount of light produced by SCALA in
        the different wavelegths
        OUTPUT:
               data_clipped : array of amplitudes in which are stored the
                              exposures values for every wavelength observed
                              in the file considered. This is used only for test.
               x_data_clipped:corresponding time for the data_clipped
               light         : array containing all the integrated exposures for the
                               wavelengths observed in a file. This is what we use
                               for light level estimation from CLAP
               mask_list     : lists of masks for each wavelength examined in a file.
                               This is the mask applied to distinguish between the
                               signal and background level for CLAP1
               clip_back     : array containing only the background data cleaned by
                               cosmics for all the wavelengths observed
               x_back        : corresponding times to clip_back
        """
        
        self.sub_exp,self.data_clipped,self.x_data_clipped, self.clip_back,self.x_back = [],[],[],[],[]
        # we have to consider the case we need also the data from CLAP0
        # since the S/N for this CLAP is much small it is easier to process
        # the data if we use the mask from the CLAP1 analysis and therefore
        # we have to save this information 
        self.mask_list = [1.] # mask to separate the signal data from the background
        self.light     = [0.,0.]
        #-- Data clipped --#
        # Loop over the different wavelengths in a file
        # analysing every single block of back-signal-back data 
        for i in range(len(self.lbda)):
            #print i
            # first we remove the nan from the data (clap1)
            data_new      = N.asarray([l for l in self.data[i,:] if N.isnan(l)==False])
            x_array       = N.linspace(self.time_start[i],self.time_stop[i],len(data_new))
            # we cut the data in a small part where we know there is only signal 
            data_sub      = N.asarray(data_new[int((self.time_on[i]-self.time_start[i]+0.5)*self.frequency):int((self.time_off[i]-self.time_start[i]-0.5)*self.frequency)])
            # we make a first rough mask to this smaller dataset to clean it from cosmics
            flagout       = self.clipped_array(data_sub,sigma =5)
            # we apply this mask
            masked_data   = data_sub[-flagout]
            # we can now obtain the mean and std of the signal to create a refined mask for the whole dataset 
            mask          = self.clipped_array(data_new,masked_data.std(), masked_data.mean(),sigma=5)
            # we now have to remove the possible cosmics which were in the background
            # and end up in the signal by mistake, the time separation must be smaller then 0.010 s
            x_array_sig   = x_array[-mask]
            mask_rem_cosm = x_array_sig<=0.010
            
            # we apply this mask to get the entire signal dataset
            data_new_sig  = data_new[-mask]
            # this method actually cut also the rising and decreasing part of the signal which must be considered as well
            # we use the information got untill now only to compute the initial guess of the parameters of the fit performed later
            
            # we compute a first estimate of the exposure time of the single wavelength counting how many data are in the signal
            # we use this value to compute the half width of the top hat (w) and its center (xe) 
            self.sub_exp  = N.append(self.sub_exp,  x_array_sig[-mask_rem_cosm][-1] -  x_array_sig[-mask_rem_cosm][0])
            # we can now compute the mean value of the signal
            raw_mean_clip = N.mean(data_new_sig[-mask_rem_cosm])
            # applying again the mask we get everything that was excluded before
            clip_back      = data_new[mask]
            x_b            =  x_array[mask]
            # to obtain the background we have to clean it from cosmics and signal related data 
            flagout_b      = self.clipped_array(clip_back,clip_back[:1000].std(), clip_back[:1000].mean(),sigma=4)
            reclip_back    = clip_back[-flagout_b]
    
            # The next two lines will be useful for daytime light leaks removal
            self.clip_back = N.append(self.clip_back,reclip_back)
            self.x_back    = N.append(self.x_back, x_b[-flagout_b])
            # we can compute the average background level
            clip_back_mean = N.mean(reclip_back)
            self.x_data_clipped    = N.append(self.x_data_clipped,(-self.time_start[i]+self.time_stop[i])/2. + self.time_start[i])
            
            if -(raw_mean_clip-N.mean(reclip_back[-1000:])) <= 5.:
                print "shutter issue at wave %d, Please check file '%s'" %(self.lbda[i], clap_file)
            else:
                # the amplitude of the signal (parameter a in the fit)
                self.data_clipped  = N.append(self.data_clipped,-(raw_mean_clip - clip_back_mean))
                #print i, self.data_clipped
                #-- Error clipped --#
                #data_error         = N.sqrt(raw_error**2 + ((N.std(clip_back1)/N.sqrt(len(clip_back1)-1.))**2+(N.std(clip_back2)/N.sqrt(len(clip_back2)-1.))**2)/4.)
                #self.data_error    = N.append(self.data_error, data_error)

            integral_result = self.data_integral(x_array,data_new,self.data_clipped[i],self.sub_exp[i],x_array_sig[-mask_rem_cosm][0],clip_back_mean)
            # the integral of the total amount of light produced by SCALA
            self.light = N.vstack((self.light, integral_result))
            self.mask_list.append(mask)

        self.light = N.delete(self.light,0,0)
        self.mask_list.remove(1.)
        if no_out:
           """return nothing"""
        else:
            return self.data_clipped, self.x_data_clipped,self.light, self.mask_list, self.clip_back,self.x_back
     

class Clap0_Data(Clap):
    """
    This Class handle the data from Clap0 and compute the light level
    """
            
    
    def sig_back(self, mask):
        """
        Compute the signal and background level of CLAP0 data
        and also compute the overall amount of light produced by SCALA in
        the different wavelegths using the mask computed for Clap1 to
        distinguish the signal data from the background data
        INPUT :
               mask : this is a list of masks from the CLAP1 analysis
                      which must be applied on CLAP0 data to distinguish
                      between signal and background
        OUTPUT:
              data_clipped: array of amplitudes in which are stored the
                            exposures values for every wavelength observed
                            in the file considered. This is used only for test.
               light      : array containing all the integrated exposures for the
                            wavelengths observed in a file. This is what we use
                            for light level estimation from CLAP0
               clip_back  : array containing only the background data cleaned by
                            cosmics for all the wavelengths observed
               x_back     : corresponding times to clip_back

        """

        self.light = [0.,0.]
        self.sub_exp,self.data_clipped, self.clip_back,self.x_back = [],[],[],[]
        #-- Data clipped --#
        # Loop over the different wavelengths in a file
        # analysing every single block of back-signal-back data
        mask       = N.array(mask)
        for i in range(len(self.lbda)):
            #print i
            # first we remove the nan from the data (clap0)
            data_new           = N.asarray([l for l in self.data[i,:] if N.isnan(l)==False])
            x_array            = N.linspace(self.time_start[i],self.time_stop[i],len(data_new))
            # we compute a first estimate of the exposure time of the single wavelength counting how many data are in the signal
            # we use this value to compute the half width of the top hat (w) and its center (xe)
            # we now have to remove the possible cosmics which were in the background
            # and end up in the signal by mistake, the time separation must be smaller then 0.010 s
            x_array_sig   = x_array[-mask[i]]
            mask_rem_cosm = x_array_sig<=0.010
            data_new_sig  = data_new[-mask]
            self.sub_exp       = N.append(self.sub_exp,  x_array_sig[-mask_rem_cosm][-1] -  x_array[-mask_rem_cosm][0])
            # first we have to compute the background level as before
            clip_back          = data_new[mask[i]]
            flagout_b          = self.clipped_array(clip_back,clip_back[:1000].std(), clip_back[:1000].mean(),sigma=4)
            reclip_back        = clip_back[-flagout_b]
            # The next two lines will be useful for daytime light leaks removal
            self.clip_back = N.append(self.clip_back,reclip_back)
            self.x_back    = N.append(self.x_back, x_b[-flagout_b])
            # we can compute the average background level
            clip_back_mean     = N.mean(reclip_back)
            # now we compute the signal level
            self.data_clipped  = N.append(self.data_clipped, -(N.mean(data_new_sig[-mask_rem_cosm])-clip_back_mean))
            integral_result    = self.data_integral(x_array,data_new,self.data_clipped[i],self.sub_exp[i],x_array_sig[-mask_rem_cosm][0],clip_back_mean)
            # the integral of the total amount of light produced by SCALA
            self.light = N.vstack((self.light, integral_result))
            
        self.light = N.delete(self.light,0,0)

        return self.data_clipped, self.light, self.clip_back, self.x_back

'''
class Mirror_Ratio():
    """
    This class take care of computing the integral of the
    Clap data for both Claps using the previous classes
    """

    def __init__(self,clap_file):
        """
        Open CLAPs file and save data
        """
        self.file_name = clap_file
        self.MIRROR    = "/Users/simonalombardo/new_backup/SCALA_CLAP/DATA/mirror_calib_2015/"
        self.clap1     = Clap1_Data(self.MIRROR+self.file_name)
        self.clap0     = Clap0_Data(self.MIRROR+self.file_name[0:2]+'0'+self.file_name[3:])

    def save_ratio_files(self,save_file):
        """
        This method must be used only when the data don't present
        daytime light leaks pollution or have been cleaned from that
        """

        data1  = self.clap1.sig_back()
        data0  = self.clap0.sig_back(data1[3])
        N.savetxt('%s.txt' % save_file, N.hstack((data0[1],data1[2])))


    def leaks_removal(self):
        """
        This method is meant to remove the light leaks effect from the data
        fitting Legendre polynoms (order 100) to the background level
        To make everything faster the background has been averaged in block of
        300 datapoints
        """
        
        data_1      = self.clap1.sig_back()[4:]
        # Before fitting the Legendre polynoms we have to reduce the dataset to make it fatser
        # we take the eaveraged values every 300 datapoints of the background dataset
        resamp_back = [N.mean(data_1[0][i*300:300*i+300]) for i in range(len(data_1[0])/300)]
        # Now we have to create the matrix of the Polynoms order 130
        Model       = N.array(self.Legendre_poly(data_1[1][149::300], degree=130))
        fit         = N.linalg.lstsq(Model.T, resamp_back)
        leaks_mean  = [0]*len(self.clap1.data[0,:])
        self.fit_param = fit[0]
        # after computing the fit and its parameters we iterate over the wavelengths
        # and we apply the fit parameters to every sub block of data and then we
        # subtract it so we remove the light leaks effect from the data
        
        for i in range(len(self.clap1.lbda)):
            x_array    = N.linspace(self.clap1.time_start[i],self.clap1.time_stop[i],len(self.clap1.data[0,:]))
            fit_Model  = N.array(self.Legendre_poly(x_array,x_min=self.clap1.time_start[0],x_max=self.clap1.time_stop[-1], degree=130))
            y_fit      = [N.sum(fit_Model.T[d,:]*fit[0]) for d in range(len(x_array))]
            leaks_mean = N.vstack((leaks_mean,y_fit))

        leaks_mean   = N.delete(leaks_mean,0,0)
        correct_data = [self.clap1.data[i,:] - leaks_mean[i,:] for i in range(len(self.clap1.lbda))]
        # now we replace the data with the corrected ones
        self.clap1.data = correct_data



    def Legendre_poly(self,x,x_min=None, x_max=None, degree=5):
        """
        From Mickael. This create a M*N array (M = Size of the given x-axis and N the degree of the polysome you want)
        """
        x = N.asarray(x)
        if x_min == None:
            x_min = x.min()
        if x_max == None:
            x_max = x.max()
        X = ((x - x_min ) / x_max )* 2 - 1
        continuum_degree = degree
        return [OP.legendre(i)(X) for i in range(continuum_degree)]
''' 
    

        
        
    
