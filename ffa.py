#!/usr/bin/env python
import numpy as np
import math
import time
import sys
import os
import ast
import re
from matplotlib import rc
import matplotlib.pyplot as plt
import scipy.signal
import glob
import subprocess
import optparse
import argparse
import pickle

import ffa_tools as ft
import ffa_stages as fs
import FFA_cy as FFA
import ffa_help 

time_total = time.time()

def func(List,IDs):
    return [List[i] for i in IDs]
    
    
def main():
    """
    Runs the Fast Folding Algorithm on a de-dispersed time series
    The parameters, set in Config.py, are : 
                 the mininum duty-cycle to look for,
                 the Signal-to-Noise treshold to select candidates 
                 the ranges of trial periods
                 the minimum sampling intervals corresponding to period ranges
    
    """

    parser = optparse.OptionParser(prog="ffa.py", \
                        version="Emilie Parent (Spring, 2016)", \
                        description="Applies the Fast Folding Algorithm on "\
                " a de-dispersed time series")

    parser.add_option('--h',dest='show_help',action='store_true', \
	    default = False)

    parser.add_option('--mindc',dest='mindc',type = 'float', \
        help="Minimum duty-cycle to look for. Default is 0.5%. "\
            "Options are: 0.5, 1., 1.5. It multiplies the list "\
            "of minimum sampling intervals that has to be tested in each "\
            "subranges of periods by 2 (mindc=1) or 3 (mindc =1.5)",\
            default = None)

    parser.add_option('--SN_tresh',dest='SN_tresh', type='float', \
            help="Signal-to-noise treshold for picking candidates", \
            default = None)

    parser.add_option('--win_detrend',dest='win_detrend', type='float', \
            help="Time window for detrending (in sec)", \
            default = None)

    parser.add_option('--plot',dest='periodograms', \
             action='store_true',  help="Produces periodograms "\
             "for few different duty cycles, default is False", \
	     default=False)

    parser.add_option('--save',dest='save', \
             action='store_true',  help="Dump SNR values into pickle",
	     default = False)

    parser.add_option('--noxwin',dest='noxwin', \
             action='store_true',  help="Don't show the periodogram",\
	     default = False)

    
    options,remains = parser.parse_args()
    if len(sys.argv[1:]) == 0:
    	print "No argument were given"+'\n'+'You must at least give the ffa a time series'+'\n'+'\n'
	ffa_help.print_usage()
	sys.exit()

    show_help = options.show_help
    periodograms = options.periodograms
    save = options.save 
    noxwin = options.noxwin
    SN_tresh = options.SN_tresh
    win_detrend = options.win_detrend
    mindc = options.mindc
    ffa_time=time.time()
    if show_help:
	ffa_help.print_usage()
    
    # -------------     Read configuration file        ------------------

    import Config_ffa
    p_ranges = Config_ffa.p_ranges
    dt_list  = Config_ffa.dt_list

    if options.SN_tresh == None: 
        SN_tresh = Config_ffa.SN_tresh
    if options.win_detrend == None: 
        win_detrend = Config_ffa.win_detrend
    if options.mindc == None : 
        mindc = Config_ffa.mindc
    # -------------     Select beam    ------------------
    args = sys.argv[1:]
    beam = [arg for arg in sys.argv[1:] if '.dat' in arg][0]
    try:
    	ts,name = ft.get_timeseries(beam)
    except: 
	print "Couldn't read the time series. "+'\n',\
	      "Make sure that the time series is the last element in your command line,"+'\n',\
	      "and that it is in the right format"
	sys.exit()

    T,dt,DM = ft.get_info_beam(name+'.inf')
    N = int(T/dt)
    
    # -------------     Downsampling    ------------------
    # mindc : minimum duty-cycle tested, set in cfg file 
    if mindc > 0.5:
        if mindc ==1	: dt_list = 2*dt_list 
        if mindc ==1.5	: dt_list = 3*dt_list
    dwn_ideal = int(dt_list[0]/dt)
    
    # min,max :give a range of accepted downsampling amount - for initial downsampling
    minimum_dwn,maximum_dwn  = dwn_ideal-int(dwn_ideal*0.10),dwn_ideal+int(dwn_ideal*0.20)
    ts,dwn = ft.select_factor(ts,minimum_dwn,maximum_dwn)
    
    print 'Downsampling the time series ..' 
    ts = ft.downsample(ts, dwn)
    
    # -------------     Detrending    ------------------
    
    print "Detrending .."
    #window size : over which statistics are computed 
    window_size = win_detrend*int(len(ts)/T)        
    break_points = np.arange(0,len(ts),window_size) 
    ts = scipy.signal.detrend(ts,bp=break_points)
    
    # Normalize w.r.t. maximum
    ts = ts/max(ts)
    sigma_total= np.std(ts)
    
    #===================================    FFA    ===================================
    print "Entering FFA"
    dt= T/len(ts)
    from ffa_var import *

    
    # count_lim: used in stage 2 and 3; how many consecutive downsamplings 
    count_lim = 2        
    
    # Going through subsequent sub-ranges of periods (set in config_ffa.cfg)
    # Each range of periods has it own initial sampling interval
    for num in range(len(p_ranges)):
        if num > 0:
            dwn_ideal = int(dt_list[num]/dt)
            if (dwn_ideal ==1) or (dwn_ideal == 0) :
                dwn_ideal =2 

            minimum_dwn,maximum_dwn = dwn_ideal-int(dwn_ideal*0.05),dwn_ideal+int(dwn_ideal*0.15)
            ts,dwn = ft.select_factor(ts,minimum_dwn,maximum_dwn)
            ts = ft.downsample(ts, dwn)
            sigma_total=sigma_total*np.sqrt(dwn)
            dt = T/len(ts)

        print "  Folding, period range of ", p_ranges[num], " ..."
        all_SNs_x1, all_Ps_x1, dts_x1 =  fs.ffa_code_stage1(ts, dt,T, sigma_total,p_ranges[num][0],\
                			 p_ranges[num][1], count_lim,name)
        all_SNs_x2, all_Ps_x2 , dts_x2 = fs.ffa_code_stage2(ts, dt,T, sigma_total,p_ranges[num][0],\
                            		  p_ranges[num][1], count_lim,name)
        all_SNs_x3, all_Ps_x3 , dts_x3 = fs.ffa_code_stage3(ts, dt,T, sigma_total,p_ranges[num][0],\
                			  p_ranges[num][1], count_lim,name)

        all_SNs1.append(all_SNs_x1), all_Ps1.append(all_Ps_x1), dts_1.append(dts_x1)
        all_SNs_2.append(all_SNs_x2), all_Ps_2.append(all_Ps_x2), dts_2.append(dts_x2)
        all_SNs_3.append(all_SNs_x3), all_Ps_3.append(all_Ps_x3), dts_3.append(dts_x3)

    
    
    # -------------         end of FFA        ------------------
    
    # Format the lists of S/N, periods, sampling intervals
    for i in range(len(all_Ps1)):
	SNs1.extend(all_SNs1[i])
        SNs2_phase1.extend(all_SNs_2[i][0]), SNs2_phase2.extend(all_SNs_2[i][1])
        SNs4_phase1.extend(all_SNs_2[i][2]), SNs4_phase2.extend(all_SNs_2[i][3])
        SNs4_phase3.extend(all_SNs_2[i][4]), SNs4_phase4.extend(all_SNs_2[i][5])
        SNs3_phase1.extend(all_SNs_3[i][0]), SNs3_phase2.extend(all_SNs_3[i][1])
        SNs3_phase3.extend(all_SNs_3[i][2])
        SNs9_phase1.extend(all_SNs_3[i][3]), SNs9_phase2.extend(all_SNs_3[i][4])
        SNs9_phase3.extend(all_SNs_3[i][5]), SNs9_phase4.extend(all_SNs_3[i][6])
        SNs9_phase5.extend(all_SNs_3[i][7]), SNs9_phase6.extend(all_SNs_3[i][8])
        SNs9_phase7.extend(all_SNs_3[i][9]), SNs9_phase8.extend(all_SNs_3[i][10])
        SNs9_phase9.extend(all_SNs_3[i][11])

        Ps1.extend(all_Ps1[i])    , dt1s.extend(dts_1[i])
        Ps2.extend(all_Ps_2[i][0]), dt2s.extend(dts_2[i][0])
        Ps4.extend(all_Ps_2[i][1]), dt4s.extend(dts_2[i][1])
        Ps3.extend(all_Ps_3[i][0]), dt3s.extend(dts_3[i][0])
        Ps9.extend(all_Ps_3[i][1]), dt9s.extend(dts_3[i][1])

        if count_lim ==2:
            SNs8_phase1.extend(all_SNs_2[i][6]), SNs8_phase2.extend(all_SNs_2[i][7])
            SNs8_phase3.extend(all_SNs_2[i][8]), SNs8_phase4.extend(all_SNs_2[i][9])
            SNs8_phase5.extend(all_SNs_2[i][10]), SNs8_phase6.extend(all_SNs_2[i][11])
            SNs8_phase7.extend(all_SNs_2[i][12]), SNs8_phase8.extend(all_SNs_2[i][13])
            SNs27_phase1.extend(all_SNs_3[i][12]), SNs27_phase2.extend(all_SNs_3[i][13])
            SNs27_phase3.extend(all_SNs_3[i][14]), SNs27_phase4.extend(all_SNs_3[i][15])
            SNs27_phase5.extend(all_SNs_3[i][16]), SNs27_phase6.extend(all_SNs_3[i][17])
            SNs27_phase7.extend(all_SNs_3[i][18]), SNs27_phase8.extend(all_SNs_3[i][19])
            SNs27_phase9.extend(all_SNs_3[i][20]), SNs27_phase10.extend(all_SNs_3[i][21])
            SNs27_phase11.extend(all_SNs_3[i][22]), SNs27_phase12.extend(all_SNs_3[i][23])
            SNs27_phase13.extend(all_SNs_3[i][24]), SNs27_phase14.extend(all_SNs_3[i][25])
            SNs27_phase15.extend(all_SNs_3[i][26]), SNs27_phase16.extend(all_SNs_3[i][27])
            SNs27_phase17.extend(all_SNs_3[i][28]), SNs27_phase18.extend(all_SNs_3[i][29])
            SNs27_phase19.extend(all_SNs_3[i][30]), SNs27_phase20.extend(all_SNs_3[i][31])
            SNs27_phase21.extend(all_SNs_3[i][32]), SNs27_phase22.extend(all_SNs_3[i][33])
            SNs27_phase23.extend(all_SNs_3[i][34]), SNs27_phase24.extend(all_SNs_3[i][35])
            SNs27_phase25.extend(all_SNs_3[i][36]), SNs27_phase26.extend(all_SNs_3[i][37])
            SNs27_phase27.extend(all_SNs_3[i][38])

            Ps8.extend(all_Ps_2[i][2]), dt8s.extend(dts_2[i][2])
            Ps27.extend(all_Ps_3[i][2]), dt27s.extend(dts_3[i][2])

    Ps1, SNs1 = np.concatenate(Ps1), np.concatenate(SNs1)

    print "Folding done "
    time_tot =  (time.time() - ffa_time)
    print ( " --- %.7s seconds is the FFA time ---" % time_tot),'\n'


    # ==============================     Post - FFA    =========================================
    # Calculate the mode and the MAD of array of S/N, for 7 (or 5) diff. duty cycles.
    # Assumes that the statistics are the same when duty cycle is the same
    
    loc1, scale1 = ft.param_sn_uniform(SNs1)
    loc2, scale2 = ft.param_sn_uniform(SNs2_phase1)
    loc4, scale4 = ft.param_sn_uniform(SNs4_phase1)
    loc3, scale3 = ft.param_sn_uniform(SNs3_phase1)
    loc9, scale9 = ft.param_sn_uniform(SNs9_phase1)

    if count_lim==2:
        loc8,  scale8  = ft.param_sn_uniform(SNs8_phase1)
        loc27, scale27 = ft.param_sn_uniform(SNs27_phase1)
        
    # Making lists of arrays ; all S/Ns for each duty cycle in one object
    # Don't know why if I remove this I get runtime error dictionary changed size during iteration
    x2,x4,x3,x9 = [var for var in vars() if 'SNs2_phase' in var],\
		  [var for var in vars() if 'SNs4_phase' in var],\
		  [var for var in vars() if 'SNs3_phase' in var],\
		  [var for var in vars() if 'SNs9_phase' in var]

    list_SNS = [SNs1]+ 	[locals()[var] for var in vars() if 'SNs2_phase' in var] + \
			[locals()[var] for var in vars() if 'SNs4_phase' in var] + \
			[locals()[var] for var in vars() if 'SNs3_phase' in var] + \
			[locals()[var] for var in vars() if 'SNs9_phase' in var] 
    
    list_PS = [Ps1]+ 2*[Ps2]+ 4*[Ps4]+ 3*[Ps3]+ 9*[Ps9]

    # -------------    Save the data into a pickle file    ------------------
    
    if save :
	SNR_to_dump = {}
	list_snr_str = ['SNs1']+  [str(var) for var in vars() if 'SNs2_phase' in var] + \
				  [str(var) for var in vars() if 'SNs4_phase' in var] + \
				  [str(var) for var in vars() if 'SNs3_phase' in var] + \
				  [str(var) for var in vars() if 'SNs9_phase' in var] 

	for i in range(len(list_snr_str)):
		SNR_to_dump[list_snr_str[i]]=np.array((list_PS[i],list_SNS[i]))
	g = open(name+'_periodogram_data.txt','w')
	pickle.dump(SNR_to_dump,g)
	print "Saved periodogram data to file : ",str(name+'_periodogram_data.txt')

     # -------------    Merging all the lists for candidate selection    ------------------

   

    list_DTS =   [dt1s]+   2*[dt2s]+   4*[dt4s]+   3*[dt3s]+   9*[dt9s]  
    list_locs =   [loc1]+   2*[loc2]+   4*[loc4]+   3*[loc3]+   9*[loc9]
    list_scales = [scale1]+ 2*[scale2]+ 4*[scale4]+ 3*[scale3]+ 9*[scale9]
    
    if count_lim ==2:
	list_SNS = list_SNS + [locals()[var] for var in vars() if 'SNs8_phase' in var] + \
		              [locals()[var] for var in vars() if 'SNs27_phase' in var]	

	list_PS =     list_PS 	  + 8*[Ps8]    + 27*[Ps27]
	list_DTS =    list_DTS 	  + 8*[dt8s]   + 27*[dt27s]
        list_locs =   list_locs   + 8*[loc8]   + 27*[loc27] 
        list_scales = list_scales + 8*[scale8] + 27*[scale27] 
    
    # ======================   Picking Candidates    ==================================
    print "Picking candidates ..."

    # write cands: only True when you get to the end. 
    write_cands, write_cands[-1] = [False]*len(list_SNS), True
    
     # -------------   Making all SNs baseline uniform + Candidates selection  ----- 
    list_SNS = [(np.array(list_SNS[i]).astype(float) - np.array(list_locs[i]).astype(float))/np.array(list_scales[i]).astype(float) for i in range(len(list_SNS))]
    
    j = [list_SNS[i] >= SN_tresh for i in range(len(list_SNS))]
    good_id = [np.nonzero(j[i]) for i in range(len(j))]
    good_id = np.concatenate(np.array(good_id))

    good_p = [func(list_PS[i],good_id[i]) for i in range(len(good_id))]
    good_sn = [func(list_SNS[i],good_id[i]) for i in range(len(good_id))]
    good_dt = [func(list_DTS[i],good_id[i]) for i in range(len(good_id))]

    good_p  = np.concatenate(good_p)
    good_sn = np.concatenate(good_sn)
    good_dt = np.concatenate(good_dt)

    good_p = [round(good_p[i],4) for i in range(len(good_p))]
    good_sn = [round(good_sn[i],4) for i in range(len(good_sn))]
    good_dt = [round(good_dt[i],5) for i in range(len(good_dt))]

    cands = ft.ffa_cands()
    cands.add_cand(good_p,good_sn,good_dt)
    ft.cands_to_file(cands,name,'_precands.ffa')

    # Check if candidates were selected
    if len(cands.periods)==0:
        print "No cands with S/N > ",SN_tresh," were detected"
        sys.exit()


    # ======================   Sifting Candidates    ==================================
    print 'Sifting the candidates ...'
    #making precands list (must be sifted)
    candsfile_str = name+'.txt'
    candsfile_list = open(candsfile_str ,'w')
    candsfile_list.write(name+'_precands.ffa')
    candsfile_list.close()
    if os.stat(name+'_precands.ffa').st_size >0:
        ft.apply_sifting(candsfile_str,name+'_cands.ffa')    

    
    print "Completed ", name
    print "Total time for FFA: ",time.time()-time_total        
    subprocess.call(["rm",candsfile_str])
    subprocess.call(["rm",name+'_precands.ffa'])
    
    # ======================   Generating Periodograms    ==================================
    # this makes 3 subplots : for delta_t, 2 x delta_t and 9 x delta_t (downsampled data, at different phases)
    if periodograms :
        rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
        plt.rc('font', family='serif')
        plt.figure(figsize=(16,14),facecolor='w')
        xmax = p_ranges[-1][1]
        xmin = p_ranges[0][0]
        ymin=0
        ymax = np.array([list_SNS[0].max(), list_SNS[1].max(), list_SNS[2].max(), list_SNS[10].max(),\
            list_SNS[12].max(), list_SNS[15].max(), list_SNS[18].max()]).max()
        if (round(xmax)-round(xmin))<12:
            steps = 1 
        else:
            steps = int((round(xmax)-round(xmin))/12)
            
        plt.subplot(311)
	# Initial sampling interval (shortest)
        plt.suptitle("Periodograms"+'\n'+name, fontsize=20)
        plt.plot(Ps1,list_SNS[0],color='mediumpurple',linewidth=1.3,label=r'$\Delta t \sim $ '+str(round(list_DTS[0][0],4))+'s')
        plt.ylabel(' S/N ' ,fontsize=20)
        plt.fill_between(Ps9, 5, 0,alpha = '0.2',color='grey')
        plt.xlim(xmin=xmin,xmax=xmax)
        plt.ylim(ymin,ymax)
        plt.xticks(np.arange(round(xmin), round(xmax), steps),fontsize=20)
        plt.yticks(fontsize=20)
	plt.minorticks_on()
        plt.legend(frameon=False,prop={'size':14})
    
        plt.subplot(312)
        # sampling interval : 2X initial
	plt.plot(Ps2,list_SNS[1],color='slateblue',linewidth=1.0,label=r'$\Delta t \sim $ '+str(round(list_DTS[1][0],4))+'s')
        plt.plot(Ps2,list_SNS[2],color='steelblue',linewidth=1.0)
        plt.ylabel(' S/N ' ,fontsize=20)
        plt.fill_between(Ps9, 5, 0,alpha = '0.2',color='grey')
        plt.xlim(xmin=xmin,xmax=xmax)
        plt.ylim(ymin,ymax)
        plt.xticks(np.arange(round(xmin), round(xmax), steps),fontsize=20)
        plt.yticks(fontsize=20)
	plt.minorticks_on()
        plt.legend(frameon=False,prop={'size':14})
    
	plt.subplot(313)
        # sampling interval : 9X initial, only  4 phases plotted out of 9
        plt.plot(Ps9,list_SNS[10],color='mediumseagreen',linewidth=1.0,label=r'$\Delta t \sim $ '+str(round(list_DTS[10][0],4))+'s')
        plt.plot(Ps9,list_SNS[12],color='lime',linewidth=1.0)
        plt.plot(Ps9,list_SNS[15],color='forestgreen',linewidth=1.0)
        plt.plot(Ps9,list_SNS[18],color='darkolivegreen',linewidth=1.0)
        plt.ylabel(' S/N ' ,fontsize=20)
        plt.xlabel('Period (s)',fontsize=20)
        plt.fill_between(Ps9, 5, 0,alpha = '0.2',color='grey')
        plt.xlim(xmin=xmin,xmax=xmax)
        plt.ylim(ymin,ymax)
        plt.xticks(np.arange(round(xmin), round(xmax), steps),fontsize=20)
        plt.yticks(fontsize=20)
	plt.minorticks_on()
        plt.legend(frameon=False,prop={'size':14})
	if not noxwin:
		plt.show()

        plt.savefig(name+'.png')
        print "Periodogram available : ", name+'.png'


if __name__=='__main__':
    main()
