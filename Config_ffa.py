import numpy as np
#configuration file FFA

#[FFA_settings]
p_ranges = np.array([[0.5,1.],[1.,2.],[2.,5.],[5.,10.],[10.,15.],[15.,30.]])
dt_list = np.array([0.002,0.005,0.01,0.02,0.050,0.075])
win_detrend = p_ranges[-1][-1]*2. # twice the largest trial period
SN_tresh    = 6.
mindc 	    = 0.5     	# min_dc = 0.5 , 1 or 1.5 (this is the min duty-cycle in %)
numdms      = 2         # relevant for sifting with ffa_final.py 
metric      = 'A'       # 'A' : SNR = max-med/std, or 'B': same but exluding a 20% window around the peak
