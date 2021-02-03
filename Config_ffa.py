import numpy as np
#configuration file FFA

#[FFA_settings]
p_ranges = np.array([[1.,2.],[2.,5.],[5.,10.],[10.,15.],[15.,30.]]) #[0.5,1]
#p_ranges = np.array([[3.0,3.5]])
dt_list = np.array([0.005,0.01,0.02,0.050,0.075]) #[0.002]
#dt_list = np.array([[0.01]])
win_detrend = p_ranges[-1][-1]*2. # twice the largest trial period
SN_tresh    = 6. #6.
mindc 	    = 0.5     	# min_dc = 0.5 , 1 or 1.5 (this is the min duty-cycle in %)
numdms      = 2         # relevant for sifting with ffa_final.py 
metric      = 'A'       # 'A', 'B', or 'C', see Parent et al. 2018 for details on the metrics.
