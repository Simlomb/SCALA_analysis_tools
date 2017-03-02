#! /usr/bin/env python
# - Analysis of SCALA data by Simona Lombardo
# - 2015
# -*- coding: utf-8 -*-

import numpy              as N
import optparse
import SCALA_data_process as SC
import pickle             as pk



def SCALA_cube(s_list, c_list,analys_type):
    """
    perform the cube analysis computing the analysis on the entire
    datacube, hence spaxel by spaxel
    """
    scala = SC.Cube_analysis(s_list, c_list,analys_type)
    name = c_list[0][4:14]
    clap_lbd=[]
    for n in range(len(scala.clap_data)):
        clap_lbd = N.append(clap_lbd,scala.clap_data[n].lbda)

    N.savetxt('%s_clap_lbd.txt' %name,clap_lbd)

    all_wave     = N.sum([len(scala.clap_data[k].lbda) for k in range(len(scala.clap_data))])
    matrix_cal_B = N.zeros((2,all_wave,15,15))
    matrix_cal_R = N.zeros((2,all_wave,15,15))
    
    #matrix_var_cal_R2 = N.zeros((2,all_wave,15,15))
    # for test
    matrix_snifs_int_B = N.zeros((2,all_wave,15,15))
    matrix_snifs_int_R = N.zeros((2,all_wave,15,15))
    matrix_clap_int_B = N.zeros((2,all_wave,15,15))
    matrix_clap_int_R = N.zeros((2,all_wave,15,15))
    matrix_clap_tot_B = N.zeros((2,all_wave,15,15))
    matrix_clap_tot_R = N.zeros((2,all_wave,15,15))
    lbd_tot_b = N.zeros((all_wave))
    lbd_tot_r = N.zeros((all_wave))
    calib_tot_b,var_tot_b = [],[]
    calib_tot_r,var_tot_r = [],[]
    for i in range(15):
            
        for j in range(15):
            print ("Computing the calibration of spaxel %ix%i" %(i,j))
            calib_sp_B,calib_sp_R,var_calib_sp_B,var_calib_sp_R,lbda_B,lbda_R = [],[],[],[],[],[]
            #for test purpose
            snifs_int_sp_B,err_snifs_int_sp_B,clap_int_sp_B,err_clap_int_sp_B,tot_clap_sp_B,err_tot_clap_sp_B,snifs_exp_time_B = [],[],[],[],[],[],[]
            snifs_int_sp_R,err_snifs_int_sp_R,clap_int_sp_R,err_clap_int_sp_R,tot_clap_sp_R,err_tot_clap_sp_R,snifs_exp_time_R = [],[],[],[],[],[],[]
            for k in range(len(scala.clap_data)):
                #print c_list[k]
                scala.Clap_info_and_light_level(k)
                result_calib_B = scala.channel_analysis(k, 'B',i,j)
                calib_sp_B = N.append(calib_sp_B,result_calib_B[0])
                var_calib_sp_B = N.append(var_calib_sp_B,result_calib_B[1])
                lbda_B = N.append(lbda_B, scala.lbdacent_snifs)
                ### for test purpose in B
                snifs_int_sp_B = N.append(snifs_int_sp_B,result_calib_B[2])
                err_snifs_int_sp_B = N.append(err_snifs_int_sp_B,result_calib_B[3])
                clap_int_sp_B = N.append(clap_int_sp_B,result_calib_B[4])
                err_clap_int_sp_B = N.append(err_clap_int_sp_B,result_calib_B[5])
                tot_clap_sp_B = N.append(tot_clap_sp_B,result_calib_B[6])
                err_tot_clap_sp_B = N.append(err_tot_clap_sp_B,result_calib_B[7])
                snifs_exp_time_B = N.append(snifs_exp_time_B,result_calib_B[8])
                
                result_calib_R= scala.channel_analysis(k, 'R',i,j)
                calib_sp_R = N.append(calib_sp_R,result_calib_R[0])
                var_calib_sp_R = N.append(var_calib_sp_R,result_calib_R[1])
                lbda_R       = N.append(lbda_R, scala.lbdacent_snifs)
                ### for test purpose in R
                snifs_int_sp_R = N.append(snifs_int_sp_R,result_calib_R[2])
                err_snifs_int_sp_R = N.append(err_snifs_int_sp_R,result_calib_R[3])
                clap_int_sp_R = N.append(clap_int_sp_R,result_calib_R[4])
                err_clap_int_sp_R = N.append(err_clap_int_sp_R,result_calib_R[5])
                tot_clap_sp_R = N.append(tot_clap_sp_R,result_calib_R[6])
                err_tot_clap_sp_R = N.append(err_tot_clap_sp_R,result_calib_R[7])
                snifs_exp_time_R = N.append(snifs_exp_time_R,result_calib_R[8])
            if i == 7 and j == 7:
                index_sortedB         = N.argsort(lbda_B)
                index_sortedR         = N.argsort(lbda_R)
                N.savetxt('%s_clap_lbd_B_sorted.txt' %name,clap_lbd[index_sortedB])
                N.savetxt('%s_clap_lbd_R_sorted.txt' %name,clap_lbd[index_sortedR])
                N.savetxt('%s_snifs_lbd_B_sorted.txt' %name,lbda_B[index_sortedB])
                N.savetxt('%s_snifs_lbd_R_sorted.txt' %name,lbda_R[index_sortedR])
                N.savetxt('exp_time_snifs_B.txt', snifs_exp_time_B)
                N.savetxt('exp_time_snifs_R.txt', snifs_exp_time_R)
            matrix_cal_B[:,:,i,j] = N.array((lbda_B,calib_sp_B))
            matrix_cal_R[:,:,i,j] = N.array((lbda_R,calib_sp_R))
            # for test purpose
            matrix_snifs_int_B[:,:,i,j] = N.array((snifs_int_sp_B,err_snifs_int_sp_B))
            matrix_snifs_int_R[:,:,i,j] = N.array((snifs_int_sp_R,err_snifs_int_sp_R))
            matrix_clap_int_B[:,:,i,j] = N.array((clap_int_sp_B,err_clap_int_sp_B))
            matrix_clap_int_R[:,:,i,j] = N.array((clap_int_sp_R,err_clap_int_sp_R))
            matrix_clap_tot_B[:,:,i,j] = N.array((tot_clap_sp_B,err_tot_clap_sp_B))
            matrix_clap_tot_R[:,:,i,j] = N.array((tot_clap_sp_R,err_tot_clap_sp_R))
            lbd_tot_b = lbd_tot_b + lbda_B
            lbd_tot_r = lbd_tot_r + lbda_R
            calib_tot_b.append(calib_sp_B)
            calib_tot_r.append(calib_sp_R)
            var_tot_b.append(var_calib_sp_B)
            var_tot_r.append(var_calib_sp_R)
        
    # save the txt file of averaged calibration over all spaxels
    cal_B_tot=[calib_tot_b[i][index_sortedB] for i in range(len(calib_tot_b[:]))]
    cal_B_tot = N.array((cal_B_tot))
    cal_R_tot = [calib_tot_r[i][index_sortedR] for i in range(len(calib_tot_r[:]))]
    cal_R_tot = N.array((cal_R_tot))
    var_B_tot=[var_tot_b[i][index_sortedB] for i in range(len(var_tot_b[:]))]
    var_B_tot = N.array((var_B_tot))
    var_R_tot=[var_tot_r[i][index_sortedR] for i in range(len(var_tot_r[:]))]
    var_R_tot = N.array((var_R_tot))

    N.savetxt('%s_calibration_B_new.txt' %name, cal_B_tot)
    N.savetxt('%s_calibration_R_new.txt' %name, cal_R_tot)
    N.savetxt('%s_var_calibration_B_new.txt' %name, var_B_tot)
    N.savetxt('%s_var_calibration_R_new.txt' %name, var_R_tot)
    
    #save all the other useful info in dacube format

    with open('%s_B_new.pickle' %name, 'wb') as transmiss_B:
        pk.dump(matrix_cal_B, transmiss_B)

    with open('%s_R_new.pickle' %name, 'wb') as transmiss_R:
        pk.dump(matrix_cal_R, transmiss_R)
    
    '''
    with open('var_%s_B2_new.pickle' %name, 'wb') as var_transmiss_B2:
        pk.dump(matrix_var_cal_B2, var_transmiss_B2)

    with open('var_%s_R2_new.pickle' %name, 'wb') as var_transmiss_R2:
        pk.dump(matrix_var_cal_R2, var_transmiss_R2)
    '''
    #### for test
    with open('%s_B_snifs_int.pickle' %name, 'wb') as snifs_B:
        pk.dump(matrix_snifs_int_B, snifs_B)

    with open('%s_R_snifs_int.pickle' %name, 'wb') as snifs_R:
        pk.dump(matrix_snifs_int_R, snifs_R)

    with open('%s_B_clap_int.pickle' %name, 'wb') as clap_B:
        pk.dump(matrix_clap_int_B, clap_B)

    with open('%s_R_clap_int.pickle' %name, 'wb') as clap_R:
        pk.dump(matrix_clap_int_R, clap_R)

    with open('%s_B_clap_tot.pickle' %name, 'wb') as clap_tot_B:
        pk.dump(matrix_clap_tot_B, clap_tot_B)

    with open('%s_R_clap_tot.pickle' %name, 'wb') as clap_tot_R:
        pk.dump(matrix_clap_tot_R, clap_tot_R)

    
def SCALA_spec(s_list, c_list,analys_type):
    """
    perform the spectrum analysis collapsing all the spaxels
    together and computing the analysis on the resulting spectrum
    """
      
    scala = SC.Spectrum_analysis(s_list, c_list,analys_type)
    name = c_list[0][4:14]
    clap_lbd=[]
    for n in range(len(scala.clap_data)):
        clap_lbd = N.append(clap_lbd,scala.clap_data[n].lbda)

    N.savetxt('%s_clap_lbd.txt' %name,clap_lbd)

    all_wave     = N.sum([len(scala.clap_data[k].lbda) for k in range(len(scala.clap_data))])
    
    # for test
    matrix_snifs_int_B = N.zeros((2,all_wave))
    matrix_snifs_int_R = N.zeros((2,all_wave))
    matrix_clap_int_B = N.zeros((2,all_wave))
    matrix_clap_int_R = N.zeros((2,all_wave))
    matrix_clap_tot_B = N.zeros((2,all_wave))
    matrix_clap_tot_R = N.zeros((2,all_wave))
    matrix_max_snifs_B = N.zeros((2,all_wave))
    matrix_max_snifs_R = N.zeros((2,all_wave))
    
   
            
    calib_sp_B,calib_sp_R,var_calib_sp_B,var_calib_sp_R,lbda_B,lbda_R = [],[],[],[],[],[]
    #for test purpose
    snifs_int_sp_B,err_snifs_int_sp_B,clap_int_sp_B,err_clap_int_sp_B,tot_clap_sp_B,err_tot_clap_sp_B,snifs_exp_time_B,max_snifs_B = [],[],[],[],[],[],[],[]
    snifs_int_sp_R,err_snifs_int_sp_R,clap_int_sp_R,err_clap_int_sp_R,tot_clap_sp_R,err_tot_clap_sp_R,snifs_exp_time_R,max_snifs_R = [],[],[],[],[],[],[],[]
    for k in range(len(scala.clap_data)):
        print ("Computing the calibration of file %i" %k)
        #print c_list[k]
        scala.Clap_info_and_light_level(k)
        result_calib_B = scala.channel_analysis(k, 'B')
        calib_sp_B = N.append(calib_sp_B,result_calib_B[0])
        var_calib_sp_B = N.append(var_calib_sp_B,result_calib_B[1])
        lbda_B = N.append(lbda_B, scala.lbdacent_snifs)
        ### for test purpose in B
        snifs_int_sp_B = N.append(snifs_int_sp_B,result_calib_B[2])
        err_snifs_int_sp_B = N.append(err_snifs_int_sp_B,result_calib_B[3])
        clap_int_sp_B = N.append(clap_int_sp_B,result_calib_B[4])
        err_clap_int_sp_B = N.append(err_clap_int_sp_B,result_calib_B[5])
        tot_clap_sp_B = N.append(tot_clap_sp_B,result_calib_B[6])
        err_tot_clap_sp_B = N.append(err_tot_clap_sp_B,result_calib_B[7])
        snifs_exp_time_B = N.append(snifs_exp_time_B,result_calib_B[8])
        max_snifs_B = N.append(max_snifs_B,result_calib_B[9])
        
        result_calib_R= scala.channel_analysis(k, 'R')
        calib_sp_R = N.append(calib_sp_R,result_calib_R[0])
        var_calib_sp_R = N.append(var_calib_sp_R,result_calib_R[1])
        lbda_R       = N.append(lbda_R, scala.lbdacent_snifs)
        ### for test purpose in R
        snifs_int_sp_R = N.append(snifs_int_sp_R,result_calib_R[2])
        err_snifs_int_sp_R = N.append(err_snifs_int_sp_R,result_calib_R[3])
        clap_int_sp_R = N.append(clap_int_sp_R,result_calib_R[4])
        err_clap_int_sp_R = N.append(err_clap_int_sp_R,result_calib_R[5])
        tot_clap_sp_R = N.append(tot_clap_sp_R,result_calib_R[6])
        err_tot_clap_sp_R = N.append(err_tot_clap_sp_R,result_calib_R[7])
        snifs_exp_time_R = N.append(snifs_exp_time_R,result_calib_R[8])
        max_snifs_R = N.append(max_snifs_R,result_calib_R[9])
        
    index_sortedB         = N.argsort(lbda_B)
    index_sortedR         = N.argsort(lbda_R)
        
    N.savetxt('%s_clap_lbd_B_sorted.txt' %name,clap_lbd[index_sortedB])
    N.savetxt('%s_clap_lbd_R_sorted.txt' %name,clap_lbd[index_sortedR])
    N.savetxt('%s_snifs_lbd_B_sorted.txt' %name,lbda_B[index_sortedB])
    N.savetxt('%s_snifs_lbd_R_sorted.txt' %name,lbda_R[index_sortedR])
    N.savetxt('exp_time_snifs_B.txt', snifs_exp_time_B)
    N.savetxt('exp_time_snifs_R.txt', snifs_exp_time_R)
        
    matrix_cal_B = N.array((lbda_B[index_sortedB],calib_sp_B[index_sortedB],var_calib_sp_B[index_sortedB]))
    matrix_cal_R = N.array((lbda_R[index_sortedR],calib_sp_R[index_sortedR],var_calib_sp_R[index_sortedR]))
    
    # for test purpose
    matrix_snifs_int_B[:,:] = N.array((snifs_int_sp_B,err_snifs_int_sp_B))
    matrix_snifs_int_R[:,:] = N.array((snifs_int_sp_R,err_snifs_int_sp_R))
    matrix_clap_int_B[:,:] = N.array((clap_int_sp_B,err_clap_int_sp_B))
    matrix_clap_int_R[:,:] = N.array((clap_int_sp_R,err_clap_int_sp_R))
    matrix_clap_tot_B[:,:] = N.array((tot_clap_sp_B,err_tot_clap_sp_B))
    matrix_clap_tot_R[:,:] = N.array((tot_clap_sp_R,err_tot_clap_sp_R))
    matrix_max_snifs_R[:,:] = N.array((lbda_R[index_sortedR],max_snifs_R[index_sortedR]))
    matrix_max_snifs_B[:,:] = N.array((lbda_B[index_sortedB],max_snifs_B[index_sortedB]))
	#print "DONE!!!"

    N.savetxt('%s_calibration_B.txt' %name, matrix_cal_B.T)
    N.savetxt('%s_calibration_R.txt' %name, matrix_cal_R.T)
    '''
    with open('%s_B_new.pickle' %name, 'wb') as transmiss_B:
        pk.dump(matrix_cal_B, transmiss_B)

    with open('%s_R_new.pickle' %name, 'wb') as transmiss_R:
        pk.dump(matrix_cal_R, transmiss_R)


    with open('var_%s_B2_new.pickle' %name, 'wb') as var_transmiss_B2:
        pk.dump(matrix_var_cal_B2, var_transmiss_B)

    with open('var_%s_R2_new.pickle' %name, 'wb') as var_transmiss_R2:
        pk.dump(matrix_var_cal_R2, var_transmiss_R)
    ''' 
    #### for test
    with open('%s_B_snifs_int.pickle' %name, 'wb') as snifs_B:
        pk.dump(matrix_snifs_int_B, snifs_B)

    with open('%s_R_snifs_int.pickle' %name, 'wb') as snifs_R:
        pk.dump(matrix_snifs_int_R, snifs_R)

    with open('%s_B_clap_int.pickle' %name, 'wb') as clap_B:
        pk.dump(matrix_clap_int_B, clap_B)

    with open('%s_R_clap_int.pickle' %name, 'wb') as clap_R:
        pk.dump(matrix_clap_int_R, clap_R)

    with open('%s_B_clap_tot.pickle' %name, 'wb') as clap_tot_B:
        pk.dump(matrix_clap_tot_B, clap_tot_B)

    with open('%s_R_clap_tot.pickle' %name, 'wb') as clap_tot_R:
        pk.dump(matrix_clap_tot_R, clap_tot_R)

    with open('%s_B_max_snifs.pickle' %name, 'wb') as max_B:
        pk.dump(matrix_max_snifs_B, max_B)

    with open('%s_R_max_snifs.pickle' %name, 'wb') as max_R:
        pk.dump(matrix_max_snifs_R, max_R)



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
    parser.add_option("-d", "--data_type",
                      help="Kind of format of data. 'cube' for cube analysis or "
                      "'spec' for spectrum analysis [%default]", default='cube')

    opts, args = parser.parse_args()

    if opts.kind not in ['t', 'f']:
        raise IOError("Analysis kind (option --kind) must be throughput ('t') "
                      "or flux calibration ('f')")
    elif opts.scala_list is None or opts.clap_list is None:
        raise IOError("Please use opts.scala_list and opts.clap_list properly")
                      
    
    opts,args = parser.parse_args()
    s_list = sorted([str(item) for item in opts.scala_list.split(',')])
    c_list = sorted([str(item) for item in opts.clap_list.split(',')])
    if len(s_list) != len(c_list):
        for i in range(len(s_list)/2):
            while (s_list[i*2][-13:-10]) != (c_list[i][-8:-5]):
                print ("WARNING one SCALA file is missing for blue and red!! The corresponding CLAP data would be %s and it will be removed by the analysis!" %c_list[i])
                del c_list[i]
                del c_list[i+(len(c_list)/2)]
                
    # with open('/Users/simonalombardo/new_SCALA_analysis/test_SCALA/159_flux_Calib/psf_sum_norm_B6.pickle', 'rb') as psf_info:
    #            psf_B = pk.load(psf_info)

    #with open('/Users/simonalombardo/new_SCALA_analysis/test_SCALA/159_flux_Calib/psf_sum_norm_R6.pickle', 'rb') as psf_info:
    #            psf_R = pk.load(psf_info)
                
    analys_type = str(opts.kind)
    data_format = str(opts.data_type)

    if data_format == 'cube':
        SCALA_cube(s_list, c_list,analys_type)
        print "Cube analysis performed"
    else:
        SCALA_spec(s_list, c_list,analys_type)
        print "Spectrum analysis performed"
