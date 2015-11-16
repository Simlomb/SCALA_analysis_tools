#! /usr/bin/env python
# - Analysis of SCALA data by Simona Lombardo
# - 2015
# -*- coding: utf-8 -*-

import numpy              as N
import optparse
import SCALA_data_process as SC
import pickle             as pk

##########################################################################
###############################  MAIN  ###################################
##########################################################################
if __name__ == '__main__':
    ##########################################################################
    #######################  OPTION ##########################################
    ##########################################################################
    parser = optparse.OptionParser()
    parser.add_option("-s", "--scala_list",
                      help="List of SNIFS files from a given night")
    parser.add_option("-c", "--clap_list",
                      help="List of CLAP files from the same night as scala's")
    parser.add_option("-k", "--kind",
                      help="Kind of analysis to perform. 't' for throughput or "
                      "'f' for flux calibration [%default]", default='t')

    opts, args = parser.parse_args()

    if opts.kind not in ['t', 'f']:
        raise IOError("Analysis kind (option --kind) must be throughput ('t') "
                      "or flux calibration ('f')")
    elif opts.scala_list is None or opts.clap_list is None:
        raise IOError("Please use opts.scala_list and opts.clap_list properly")

    s_list = sorted([str(item) for item in opts.scala_list.split(',')])
    c_list = sorted([str(item) for item in opts.clap_list.split(',')])
    if len(s_list) != len(c_list):
        for i in range(len(s_list)/2):
            while (s_list[i*2][-13:-10]) != (c_list[i][-8:-5]):
                print ("WARNING one SCALA file is missing for blue and red!!"
                       "The corresponding CLAP data would be "
                       "%s and it will be removed by the analysis!" % c_list[i])
                del c_list[i]
                del c_list[i+(len(c_list)/2)]

    scala = SC.SCALA_Calib(s_list, c_list, opts.kind == 't')
    all_wave = N.sum([len(scala.clap_data[k].lbda)
                      for k in range(len(scala.clap_data))])
    matrix_cal_B = N.zeros((2, all_wave, 15, 15))
    matrix_cal_R = N.zeros((2, all_wave, 15, 15))
    for i in range(15):
        for j in range(15):
            print "Computing the calibration of spaxel %ix%i" %(i, j)
            calib_sp_B, calib_sp_R, lbda_B, lbda_R = [], [], [], []
            for k in range(len(scala.clap_data)):
                dimension = len(scala.clap_data[k].lbda)
                calib_file_B, calib_file_R = [], []
                scala.Clap_info_and_light_level(k)
                for d in range(dimension):
                    calib_file_B = N.append(calib_file_B,
                                            scala.channel_analysis(k, 'B',
                                                                   d, i, j))
                    calib_file_R = N.append(calib_file_R,
                                            scala.channel_analysis(k, 'R',
                                                                   d, i, j))
                    lbda_B = N.append(lbda_B, scala.lbdacent_snifs)
                    lbda_R = N.append(lbda_R, scala.lbdacent_snifs)
                calib_sp_B = N.append(calib_sp_B, calib_file_B)
                calib_sp_R = N.append(calib_sp_R, calib_file_R)

            index_sortedB = N.argsort(lbda_B)
            index_sortedR = N.argsort(lbda_R)
            matrix_cal_B[:, :, i, j] = N.array((lbda_B[index_sortedB],
                                                calib_sp_B[index_sortedB]))
            matrix_cal_R[:, :, i, j] = N.array((lbda_R[index_sortedR],
                                                calib_sp_R[index_sortedR]))
	#print "DONE!!!"
name = c_list[0][4:14]
with open('%s_B.pickle' %name, 'wb') as transmiss_B:
    pk.dump(matrix_cal_B, transmiss_B)

with open('%s_R.pickle' %name, 'wb') as transmiss_R:
    pk.dump(matrix_cal_R, transmiss_R)
