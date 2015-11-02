#! /usr/bin/env python
# - Analysis of SCALA data by Simona Lombardo
# - 2015
# -*- coding: utf-8 -*-


import optparse
import new_analysis_SCALA_6 as SC

##########################################################################
###############################  MAIN  ###################################
##########################################################################
if __name__ == '__main__':
    ##########################################################################
    #######################  OPTION ##########################################
    ##########################################################################
    parser = optparse.OptionParser()
    parser.add_option("-s", "--scala_list",
                      help="list of files from SNIFS from a same night",
                      )
    parser.add_option("-c", "--clap_list",
                      help="list of  files from clap from the same night as scala's",
                      )
    
    parser.add_option("-t", "--throughput",
                      help="kind of analysis to perform, is True if is the throughput or False otherwise",
                      default='True')
                      
    
    opts,args = parser.parse_args()

    
    scala = SC.SCALA_Calib(opts.scala_list, opts.clap_list, opts.throughput)
    all_wave     = N.sum([len(scala.clap_data[k].lbda) for k in range(len(scala.clap_data))])
    matrix_cal_B = N.zeros((all_wave,all_wave,15,15))
    matrix_cal_R = N.zeros((all_wave,all_wave,15,15)) 
    for i in range(15):
            
        for j in range(15):
            print ("Computing the calibration of spaxel %ix%i" %(i,j))
            calib_sp_B,calib_sp_R,lbda_B,lbda_R = [],[],[],[]
            for k in range(len(scala.clap_data)):
                dimension = len(scala.clap_data[k].lbda)
                calib_file_B,calib_file_R = [],[]
                scala.Clap_info_and_light_level(k)
                for d in range(dimension):
                    calib_file_B = N.append(calib_B,scala.channel_analysis(k, 'B', d, i,j))
                    calib_file_R = N.append(calib_R, scala.channel_analysis(k, 'R', d, i,j))
                    lbda_B       = N.append(lbda_B, scala.lbdacent_snifs)
                    lbda_R       = N.append(lbda_R, scala.lbdacent_snifs)
                calib_sp_B = N.append(calib_sp_B, calib_file_B)
                calib_sp_R = N.append(calib_sp_R, calib_file_R)
                
            index_sortedB         = N.argsort(lbda_B)
            index_sortedR         = N.argsort(lmbd_R)
            matrix_cal_B[:,:,i,j] = N.array((lbda_B[index_sortedB],calib_sp_B[index_sortedB]))
            matrix_cal_R[:,:,i,j] = N.array((lbda_R[index_sortedR],calib_sp_R[index_sortedR]))
