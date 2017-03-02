#! /usr/bin/env python
# - Analysis of CLAP data by Simona Lombardo
# - 2015
# -*- coding: utf-8 -*-

import numpy              as N
import optparse
import SCALA_data_process as SC
from os import listdir
from os.path import isfile,isdir, join
from scipy import interpolate
##########################################################################
###############################  MAIN  ###################################
##########################################################################
if __name__ == '__main__':
    ## This function only analyses the CLAP data:
    # Example: python CLAP_main.py -path '/user_forlder/' -c SC1_filename.fits,SC0_filename.fits -n 'new_file_name'
    # the options are :
    #                  -p path for the data folder
    #                  -c list of file names to useseparated by ','
    #                  -n name of the txt files to save, without .txt
    #                  if only path is specified then all the CLAP's files
    #                  in that directory are used. If also -c is added then
    #                  only the selected files are used.
    #The files produced are:
    #                  name.txt a file with the lambda in the first column
    #                  and the C1*Es*D*G in the other column, where G is the
    #                  geometric factor for CLAP only and C1 is integrated over time.
    #                  name_a18.txt is the file with lambda and Es
    #                  name_integ.txt is the file with lambda and the integrated
    #                  clap measurement over time.
    #                  name_no_integ.txt is the file with lambda and the clap data
    #                  not integrated over time,  C1*Es*D*G
    #  WARNINGS:
    #           This implementation works for CLAP1 data only, needs to be checked
    #           and possibly modified for CLAP0 data!!
    #           This analysis does not include errors!!! Do it!
    #           The wavelengths are not sorted 
    ##########################################################################
    #######################  OPTION ##########################################
    ##########################################################################
    parser = optparse.OptionParser()
    parser.add_option("-p", "--path",
                      help="path where to find the data [%default]", default='./')
    parser.add_option("-c", "--clap_list",
                      help="List of CLAP files [%default]", default='empty')
    parser.add_option("-n", "--name",
                      help="Name for the file saved", default='clap_data_process')
    parser.add_option("-l", "--clap_cal",
                      help="clap calibration curve to use", default='1')
    
    
    opts,args = parser.parse_args()
    c_list = sorted([str(item) for item in opts.clap_list.split(',')])
    if str(opts.clap_list) != 'empty':
        c_list = sorted([str(item) for item in opts.clap_list.split(',')])
    else:
        c_list = sorted([ f for f in listdir(str(opts.path)) if isfile(join(str(opts.path),f)) and f.split(".")[0][:2]=="SC"])
    
   
    scala = SC.SCALA_Calib([], c_list, 't', path=str(opts.path))
    geometric_factor =314.16/(0.336*2.39e-4)
    clap_lbd,integ_clap1,a18,only_int,no_integ_clap1 = [],[],[],[],[]
    # integ_clap1 is the complete value of clap integration
    # a18 is an array with all the a18 values from the wavelength scanned
    # only int is an array with the clap values integrated over the time
    alternative_clap_C = N.loadtxt('/Users/simonalombardo/Greg_data/response_test.txt')
    intrp = interpolate.interp1d(alternative_clap_C[:,0], alternative_clap_C[:,1], kind='linear',bounds_error=False, fill_value=N.nan)
    for n in range(len(scala.clap_data)):
        print n
        lbd_new = scala.clap_data[n].snifs_wave_from_clap(scala.clap_data[n].lbda)
        clap_lbd = N.append(clap_lbd,lbd_new)
        A18_value = scala.A18[0](lbd_new)#N.mean([f(clap_lbd[n]) for f in scala.A18[0]])
        #print N.shape(A18_value)
        #print A18_value
        if N.shape(A18_value) != (1,):
            for j in range(len(scala.clap_data[n].lbda)):
                print j
                a18.append(A18_value[j])
                only_int.append(scala.integrated_clap[n][j][0])
                if str(opts.clap_cal) == '1':
                    clap1_int = scala.clap1_simul[1](lbd_new[j])
                else:
                    #print 'Heiiii'
                    clap1_int = intrp(lbd_new[j])
                weight_func = 1/(clap1_int*43.5)
                #weight_func,weight_func_err = scala.weight_CLAP_Data(lbd_new[j])
                integ_clap1 = N.append(integ_clap1,scala.integrated_clap[n][j][0]*A18_value[j]*geometric_factor*weight_func)
                no_integ_clap1 = N.append(no_integ_clap1,scala.data_clap_clipped[n][j]*A18_value[j]*geometric_factor*weight_func)
        else:
            if str(opts.clap_cal) == '1':
                clap1_int = scala.clap1_simul[1](lbd_new)
            else:
                #print 'Heiiii2'
                clap1_int = intrp(lbd_new)
            #clap1_int = scala.clap1_simul[1](lbd_new)
            weight_func = 1/(clap1_int*43.5)
            #weight_func,weight_func_err = scala.weight_CLAP_Data(lbd_new)
            integ_clap1 = N.append(integ_clap1,scala.integrated_clap[n][0][0]*A18_value*geometric_factor*weight_func)
            no_integ_clap1 = N.append(no_integ_clap1,scala.data_clap_clipped[n][0]*A18_value*geometric_factor*weight_func)
            a18.append(A18_value)
            only_int.append(scala.integrated_clap[n][0][0])
            
    cal_tot = N.array((clap_lbd,integ_clap1))
    a18_tot = N.array((clap_lbd,a18))
    int_tot = N.array((clap_lbd,only_int))
    cal_tot1 = N.array((clap_lbd,no_integ_clap1))
    N.savetxt('%s.txt' %str(opts.name), cal_tot.T)
    N.savetxt('%s_a18.txt' %str(opts.name), a18_tot.T)
    N.savetxt('%s_integ.txt' %str(opts.name), int_tot.T)
    N.savetxt('%s_no_integ.txt' %str(opts.name), cal_tot1.T)
    
    
