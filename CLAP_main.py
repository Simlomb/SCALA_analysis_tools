#! /usr/bin/env python
# - Analysis of CLAP data by Simona Lombardo
# - 2015
# -*- coding: utf-8 -*-

import numpy              as N
import optparse
import SCALA_data_process as SC
from os import listdir
from os.path import isfile,isdir, join
##########################################################################
###############################  MAIN  ###################################
##########################################################################
if __name__ == '__main__':
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
    
    
    opts,args = parser.parse_args()
    c_list = sorted([str(item) for item in opts.clap_list.split(',')])
    if str(opts.clap_list) != 'empty':
        c_list = sorted([str(item) for item in opts.clap_list.split(',')])
    else:
        c_list = sorted([ f for f in listdir(str(opts.path)) if isfile(join(str(opts.path),f)) and f.split(".")[0][:2]=="SC"])
    
   
    scala = SC.SCALA_Calib([], c_list, 't', path=str(opts.path))
    geometric_factor = 201.06/(0.33*2.39e-4)
    clap_lbd,integ_clap1,a18,only_int = [],[],[],[]
    # integ_clap1 is the complete value of clap integration
    # a18 is an array with all the a18 values from the wavelength scanned
    # only int is an array with the clap values integrated over the time

    for n in range(len(scala.clap_data)):
        clap_lbd = N.append(clap_lbd,scala.clap_data[n].lbda)
        A18_value = N.mean([f(clap_lbd[n]) for f in scala.A18])
        if N.shape(A18_value) != ():
            for j in range(len(scala.clap_data[n].lbda)):
                a18.append(A18_value[j])
                only_int.append(scala.integrated_clap2[n][j][0])
                integ_clap1 = N.append(integ_clap1,scala.integrated_clap2[n][j][0]*A18_value[j]*geometric_factor)
        else:
            integ_clap1 = N.append(integ_clap1,scala.integrated_clap2[n][0][0]*A18_value*geometric_factor)
            a18.append(A18_value)
            only_int.append(scala.integrated_clap2[n][0][0])
            
    cal_tot = N.array((clap_lbd,integ_clap1))
    a18_tot = N.array((clap_lbd,a18))
    int_tot = N.array((clap_lbd,only_int))
    N.savetxt('%s.txt' %str(opts.name), cal_tot.T)
    N.savetxt('%s_a18.txt' %str(opts.name), a18_tot.T)
    N.savetxt('%s_integ.txt' %str(opts.name), int_tot.T)
    
    
