#! /usr/bin/env python
# - Analysis of CLAPs data by Simona Lombardo
# - 2015
# -*- coding: utf-8 -*-

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
from scipy import interpolate


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
        self.clap_file  = clap_file
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
        err = (y - ((a/4.)*(erf(b*(w+x-xe))+1)*(erf(b*(w-x+xe))+1)+back))
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
               sigma       : uncertainty on the y data
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

    def line_two_point(self,x0, y0,x,err_y0,full_line=True):
        """
        """
        x1,x2 = x0
        y1,y2 = y0
        y1_err,y2_err = err_y0
        m = (y2 - y1)/(x2 - x1)
        err_m = (y2_err+y1_err)/(x2 - x1)
        q = -x1*((y2 - y1)/(x2 -x1)) + y1
        err_q = (x1/(x2 -x1))*(y2_err+y1_err)+y1_err
        if full_line:
            #lets consider only the error in the middle of the slope as background error
            err_middle = err_m*(x1+(x2-x1)/2.)+err_q
            return (m*x + q),err_middle
        else:
            return m,q
    
    def line_intersection(self,m,q):
        """
        Compute the intersection point between two lines
        """
        m1,m2 = m
        q1,q2 = q
        return (q2-q1)/(m1-m2)

    def compute_background(self,clip_back,x_b,start,end):
        """
        compute the mask to apply to background left and right
        """
        flagout_b_l = self.clipped_array(clip_back,clip_back[:1000].std(), clip_back[:1000].mean(),sigma=4)
        mask_l = (-flagout_b_l)*(x_b<start) #mask background left
        
        flagout_b_r = self.clipped_array(clip_back,clip_back[-1000:].std(), clip_back[-1000:].mean(),sigma=4)
        mask_r = (-flagout_b_r)*(x_b>end) #mask background right
    
        return mask_l, mask_r

    def remove_back(self,clip_back,x_b,mask_l,mask_r,all_data,x_all):
        """
        This is the background removal in the case of linear changes
        No strong light leaks
        """
        
        reclip_back_l = clip_back[mask_l][:-9] #clipped background left
        x_back_l = x_b[mask_l][:-9]
        reclip_back_r = clip_back[mask_r][10:] #clipped background right
        x_back_r = x_b[mask_r][10:]

        # Now we compute the two avereged point for the linear fit
        y0_l = N.array((N.mean(reclip_back_l[:len(reclip_back_l)/2]),N.mean(reclip_back_l[len(reclip_back_l)/2:]))) # for the left side
        y0_lerr = N.array((N.var(reclip_back_l[:len(reclip_back_l)/2])/(len(reclip_back_l)/2.),N.var(reclip_back_l[len(reclip_back_l)/2:])/(len(reclip_back_l)/2.)))
        y0_r = N.array((N.mean(reclip_back_r[:len(reclip_back_r)/2]),N.mean(reclip_back_r[len(reclip_back_r)/2:]))) # for the right side
        y0_rerr = N.array((N.var(reclip_back_r[:len(reclip_back_r)/2])/(len(reclip_back_r)/2.),N.var(reclip_back_r[len(reclip_back_r)/2:])/(len(reclip_back_r)/2.)))
        x0_l = N.array((N.mean(x_back_l[:len(x_back_l)/2]),N.mean(x_back_l[len(x_back_l)/2:])))
        x0_r = N.array((N.mean(x_back_r[:len(x_back_r)/2]),N.mean(x_back_r[len(x_back_r)/2:])))
        # compute the line for those two points in the back left and right
        m_l,q_l = self.line_two_point(x0_l,y0_l,x_back_l,y0_lerr,full_line=False) # parameters for left line
        m_r,q_r = self.line_two_point(x0_r,y0_r,x_back_r,y0_rerr,full_line=False) # parameters for left line
        # lets find the interception between the two
        x_intercept = self.line_intersection(N.array((m_l,m_r)), N.array((q_l,q_r)))
        #if x_intercept >= x_back_l[-1] and x_intercept <= x_back_r[0]:
        #    # Now we separete the all dataset and apply the linear fit consequently
        #    mask_slope = x_all<=x_intercept 
        #    data_corrected = N.append((all_data[mask_slope]- self.line_two_point(x0_l,y0_l,x_all[mask_slope])), (all_data[-mask_slope]- self.line_two_point(x0_r,y0_r,x_all[-mask_slope])))
        #else:
        y0 = N.array((N.mean(reclip_back_l),N.mean(reclip_back_r)))
        #if y0[0]>y0[1]+300.:
        #    mean_slope,err_slope =  self.line_two_point(x0_l,y0_l,x_all,y0_lerr)
        #else:
        y0_err = N.array((N.var(reclip_back_l)/len(reclip_back_l),N.var(reclip_back_r)/len(reclip_back_r)))
        x0 = N.array((N.mean(x_back_l),N.mean(x_back_r)))
        mean_slope,err_slope =  self.line_two_point(x0,y0,x_all,y0_err)
        data_corrected = (all_data-mean_slope)
        return data_corrected,err_slope

    def separate_sig_back(self, all_data,x_array,line,error_data=None,mask_data=None):
        """
        Separate the sig from the background using sigma clipping
        """
        if mask_data == None:
            # we cut the data in a small part where we know there is only signal 
            data_sub = N.asarray(all_data[int((self.time_on[line]-self.time_start[line]+0.5)*self.frequency):int((self.time_off[line]-self.time_start[line]-0.5)*self.frequency)])
            # we make a first rough mask to this smaller dataset to clean it from cosmics
            flagout = self.clipped_array(data_sub,sigma =5)
            # we apply this mask
            masked_data = data_sub[-flagout]
            # we can now obtain the mean and std of the signal to create a refined mask for the whole dataset 
            mask = self.clipped_array(all_data,masked_data.std(), masked_data.mean(),sigma=5)
            # we now have to remove the possible cosmics which were in the background
            # and end up in the signal by mistake, the time separation must be smaller then 0.010 s
        else:
            mask = mask_data
        x_array_sig = x_array[-mask]
        index_rem_cosm = x_array_sig[1:]-x_array_sig[:-1]
        mask_rem_cosm = index_rem_cosm>0.010
        mask_rem_cosm = N.insert(mask_rem_cosm,len(mask_rem_cosm)/2,False)
        dim = len(x_array_sig)-len(x_array_sig[-mask_rem_cosm])
        # we apply this mask to get the entire signal dataset
        # this method actually cut also the rising and decreasing part of the signal which must be considered as well
        data_new_sig = all_data[-mask]
        x_corrected = x_array_sig[-mask_rem_cosm]
        data_corrected = data_new_sig[-mask_rem_cosm]
        # the following loop helps eliminating residual cosmics
        while dim>0:
            index_rem_cosm = x_corrected[1:]-x_corrected[:-1]
            mask_rem_cosm = index_rem_cosm>0.010
            mask_rem_cosm = N.insert(mask_rem_cosm,len(mask_rem_cosm)/2,False)
            dim = len(x_corrected)-len(x_corrected[-mask_rem_cosm])
            x_corrected = x_corrected[-mask_rem_cosm]
            data_corrected = data_corrected[-mask_rem_cosm]
        # we use the information got untill now only to compute the initial guess of the parameters of the fit performed later
        # we compute a first estimate of the exposure time of the single wavelength counting how many data are in the signal
        # we use this value to compute the half width of the top hat (w) and its center (xe) 
        #self.sub_exp  = N.append(self.sub_exp,  x_corrected[-1] -  x_corrected[0])
        # applying again the mask we get everything that was excluded before
        clip_back = all_data[mask]
        x_b =  x_array[mask]
        mask_clean = (x_b>=x_corrected[0]) * (x_b<= x_corrected[-1])
        x_b = x_b[-mask_clean]
        clip_back = clip_back[-mask_clean]
        if error_data == None and mask_data == None:
            return x_corrected, data_corrected, x_b, clip_back
        elif error_data != None and mask_data == None:
            #err_sig = error_data[-mask]
            #err_sig_corr = err_sig[-mask_rem_cosm]
            #err_back = error_data[mask]
            #err_clip_back = err_back[-mask_clean]
            #return x_corrected, data_corrected,err_sig_corr, x_b, clip_back,err_clip_back, mask
            return x_corrected, data_corrected, x_b, clip_back, mask
        else:
            return x_corrected, data_corrected, x_b, clip_back


    def weighted_avg_and_stderr(self,values, weights):
        """
        Return the weighted average and variance.

        values, weights -- Numpy ndarrays with the same shape.
        """
        average = N.average(values, weights=weights)
        variance = N.average((values-average)**2, weights=weights)/(len(values)-1)
        return (average, variance)
        

    def snifs_wave_from_clap(self,clap_wave):
        """
        function compiuted with linear fit in fit_clap_snifs_wave.py
        add new data there when ready. results of that fit y=x*m+c:
        for wavelength <=5000.
        array([  0.9985409 ,  19.55074642]),
        array([[  5.42283511e-08,  -2.24787767e-04],
        [ -2.24787767e-04,   9.46129317e-01]])
       
        for wavelength >5000.:
        (array([  0.99826321,  23.70693881]),
        array([[  2.28328894e-09,  -1.70968873e-05],
        [ -1.70968873e-05,   1.33349329e-01]])
       
        This function is needed to compute the central wavelength for the
        a18 study considering that the wavelength from monochromator
        is not reproducible and sliglthly different from snifs's
        """
        if clap_wave <= 5000.:
            return (clap_wave* 0.9985409 + 19.55074642)
        else:
            return (clap_wave* 0.99826321 + 23.70693881)




class Clap1_Data(Clap):
    """
    This Class handle the data from Clap1 and compute the light level
    """
            
    
    def sig_back(self,cross_talk_correction, no_out=False, no_snifs=False,white_light=None):
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
        
        self.sub_exp,self.data_clipped,self.x_data_clipped, self.clip_back,self.x_back,data_error_tot = [],[],[],[],[],[]
        # we have to consider the case we need also the data from CLAP0
        # since the S/N for this CLAP is much small it is easier to process
        # the data if we use the mask from the CLAP1 analysis and therefore
        # we have to save this information 
        self.mask_list = [1.] # mask to separate the signal data from the background
        self.light, self.other_light = [0.,0.],[0.,0.]
        #-- Data clipped --#
        # Loop over the different wavelengths in a file
        # analysing every single block of back-signal-back data 
        for i in range(len(self.lbda)):
            #print i
            # first we remove the nan from the data (clap1)
            data_new = N.asarray([l for l in self.data[i,:] if N.isnan(l)==False])
            x_array = N.linspace(self.time_start[i],self.time_stop[i],len(data_new))
            x_corrected, sig_corrected, x_b, clip_back = self.separate_sig_back(data_new,x_array,i)
            
            # we can now compute the mean value of the signal
            raw_mean_clip = N.mean(sig_corrected)
            # lets first compute the two mask to separete the clean background left and right
            mask_l, mask_r = self.compute_background(clip_back,x_b,x_corrected[0],x_corrected[-1])
            self.clip_back = N.append(self.clip_back,N.append(clip_back[mask_l],clip_back[mask_r]))
            self.x_back    = N.append(self.x_back, N.append(x_b[mask_l],x_b[mask_r]))
            if -(raw_mean_clip-N.mean(clip_back[mask_r])) <= 5.:
                print "shutter issue at wave %d, Please check file '%s'" %(self.lbda[i],self.clap_file)
                integral_result = N.array((0.,0.))
                self.light = N.vstack((self.light, integral_result))
                self.other_light = N.vstack((self.other_light,integral_result))
            else:
                # lets apply the linear fit to remove the background
                data_new_corrected, err_back_corrected = self.remove_back(clip_back,x_b,mask_l,mask_r,data_new,x_array)
                x_new_correct, sig_new_correct,x_b_correct, clip_back_correct,mask = self.separate_sig_back(data_new_corrected,x_array,i,error_data=err_back_corrected)
                mask_l_new, mask_r_new = self.compute_background(clip_back_correct,x_b_correct,x_new_correct[0],x_new_correct[-1])
                back_corrected = N.append(clip_back_correct[mask_l_new][:-9],clip_back_correct[mask_r_new][10:])
                #err_back_corrected = N.append(err_back_corr[mask_l_new][:-9],err_back_corr[mask_r_new][10:])
                #mean_back,err_mean_back =  self.weighted_avg_and_stderr(back_corrected,1/err_back_corrected)
                mean_back = N.mean(back_corrected)
                #print i, mean_back,N.sqrt(err_back_corrected)
                # lets compose the cleaned data set so that we can integrate all data with simpsons without any cosmics
                data_cleaned, x_cleaned = [],[]
                data_cleaned = N.append(data_cleaned, -(clip_back_correct[mask_l_new]))
                if white_light == None:
                    data_cleaned = N.append(data_cleaned, -(sig_new_correct)+cross_talk_correction[i])
                else:
                    data_cleaned = N.append(data_cleaned, -(sig_new_correct)+cross_talk_correction[i]-white_light[i])
                    print 1-(N.mean(-(sig_new_correct)+cross_talk_correction[i]-white_light[i])/ N.mean(-(sig_new_correct)+cross_talk_correction[i]))
                data_cleaned = N.append(data_cleaned, -(clip_back_correct[mask_r_new]))
                x_cleaned = N.append(x_cleaned, x_b_correct[mask_l_new])
                x_cleaned = N.append(x_cleaned, x_new_correct)
                x_cleaned = N.append(x_cleaned, x_b_correct[mask_r_new])
                # The next two lines will be useful for daytime light leaks removal
                #self.clip_back = N.append(self.clip_back,reclip_back)
                #self.x_back    = N.append(self.x_back, x_b[-flagout_b])
                # we can compute the average background level
                #clip_back_mean = N.mean(reclip_back)
                #self.x_data_clipped    = N.append(self.x_data_clipped,(-self.time_start[i]+self.time_stop[i])/2. + self.time_start[i])
                # we compute a first estimate of the exposure time of the single wavelength counting how many data are in the signal
                # we use this value to compute the half width of the top hat (w) and its center (xe) 
                self.sub_exp = N.append(self.sub_exp,  x_new_correct[-1] -  x_new_correct[0])
                #raw_mean_clip_correct,err_mean_clip_corr = self.weighted_avg_and_stderr(sig_new_correct[3:-2],1/err_sig_corr[3:-2])
                raw_mean_clip_correct = N.mean(sig_new_correct[3:-2])
                err_mean_clip_corr = N.var(sig_new_correct[3:-2])/len(sig_new_correct[3:-2])
                # the amplitude of the signal (parameter a in the fit)
                self.data_clipped = N.append(self.data_clipped,-raw_mean_clip_correct+cross_talk_correction[i])
                self.x_data_clipped = N.append(self.x_data_clipped,(x_new_correct[-1] -  x_new_correct[0])/2. + x_new_correct[0]+self.time_start[i])
                #print i, self.data_clipped
                #-- Error clipped --#(
                data_error         = N.sqrt(err_mean_clip_corr+err_back_corrected)
                data_error_tot     = N.append(data_error_tot, data_error)
                #print self.data_clipped[i], data_error, data_error/self.data_clipped[i]
                # We now compute the integral of the data in two ways: integral will be computed by integrating the
                # fit function of the data, and other_integral will be the simple simpsons integration of the clap
                # data directly with corresponding error propagation
                integral_result = self.data_integral(x_array,data_new_corrected,self.data_clipped[i],self.sub_exp[i],x_new_correct[0],mean_back)
                integral_result_other = integrate.simps(data_cleaned,x_cleaned,even='first')
                error_other_integral = (((x_cleaned[-1]-x_cleaned[0])/len(x_cleaned))**2)*(2+4*len(x_cleaned[2:-2:2])+16*len(x_cleaned[1:-1:2]))*err_back_corrected/9.
                # the integral of the total amount of light produced by SCALA
                other_integral = N.array((integral_result_other,error_other_integral))
                self.light = N.vstack((self.light, integral_result))
                self.other_light = N.vstack((self.other_light,other_integral))
                self.mask_list.append(mask)
                #print self.light[i+1][0],self.other_light[i],error_other_integral,len(x_cleaned)
                
        #print N.shape(self.other_light)

        self.light = N.delete(self.light,0,0)
        self.other_light = N.delete(self.other_light,0,0)
        self.mask_list.remove(1.)
        if no_out:
           """return nothing"""
        elif no_snifs:
            
            return self.data_clipped, data_error_tot
        else:
            return self.data_clipped, self.x_data_clipped,self.light, self.mask_list, self.clip_back,self.x_back,self.other_light
     

class Clap0_Data(Clap):
    """
    This Class handle the data from Clap0 and compute the light level
    """
            
    
    def sig_back(self, mask,no_snifs=False):
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
        self.sub_exp,self.data_clipped, self.clip_back,self.x_back,data_error_tot = [],[],[],[],[]
        #-- Data clipped --#
        # Loop over the different wavelengths in a file
        # analysing every single block of back-signal-back data
        mask       = N.array(mask)
        for i in range(len(self.lbda)):
            #print i
            # first we remove the nan from the data (clap0)
            data_new           = N.asarray([l for l in self.data[i,:] if N.isnan(l)==False])
            x_array            = N.linspace(self.time_start[i],self.time_stop[i],len(data_new))
            
            x_corrected, sig_corrected, x_b, clip_back = self.separate_sig_back(data_new,x_array,i,mask_data=mask[i])
            # lets first compute the two mask to separete the clean background left and right
            mask_l, mask_r = self.compute_background(clip_back,x_b,x_corrected[0],x_corrected[-1])
            self.clip_back = N.append(self.clip_back,N.append(clip_back[mask_l],clip_back[mask_r]))
            self.x_back    = N.append(self.x_back, N.append(x_b[mask_l],x_b[mask_r]))
            # lets apply the linear fit to remove the background
            data_new_corrected, err_back_corrected = self.remove_back(clip_back,x_b,mask_l,mask_r,data_new,x_array)
            x_new_correct, sig_new_correct, x_b_correct, clip_back_correct = self.separate_sig_back(data_new_corrected,x_array,i,mask_data=mask[i])
            mask_l_new, mask_r_new = self.compute_background(clip_back_correct,x_b_correct,x_new_correct[0],x_new_correct[-1])
            back_corrected = N.append(clip_back_correct[mask_l_new][:-9],clip_back_correct[mask_r_new][10:])
            mean_back = N.mean(back_corrected)
            #print i, mean_back
            raw_mean_clip_correct = N.mean(sig_new_correct[3:-2])
            err_mean_clip_corr = N.var(sig_new_correct[3:-2])/len(sig_new_correct[3:-2])
            # the amplitude of the signal (parameter a in the fit)
            self.data_clipped = N.append(self.data_clipped,-(raw_mean_clip_correct - mean_back))
            self.sub_exp = N.append(self.sub_exp,  x_new_correct[-1] -  x_new_correct[0])
            #self.x_data_clipped = N.append(self.x_data_clipped,(x_new_correct[-1] -  x_new_correct[0])/2. + x_new_correct[0]+self.time_start[i])
            #print i, self.data_clipped
            #-- Error clipped --#(
            data_error         = N.sqrt(err_mean_clip_corr+err_back_corrected)
            data_error_tot     = N.append(data_error_tot, data_error)

            integral_result = self.data_integral(x_array,data_new_corrected,self.data_clipped[i],self.sub_exp[i],x_new_correct[0],mean_back)
            # the integral of the total amount of light produced by SCALA
            self.light = N.vstack((self.light, integral_result))
            
        self.light = N.delete(self.light,0,0)
        if no_snifs:
            return self.data_clipped, data_error_tot
        else: 
            return self.data_clipped, self.light, self.clip_back, self.x_back
