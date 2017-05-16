import sys 
def print_usage():
    help = "        Runs the Fast Folding Algorithm on a de-dispersed time series"+'\n'+\
    	"	    	 The parameters, set in Config.py, are : "+'\n'+\
    	"                 the mininum duty-cycle to look at "+'\n'+\
    	"                 the signal-to-noise treshold to select candidates "+'\n'+\
    	"                 the ranges of trial periods "+'\n'+\
    	"                 the minimum sampling intervals corresponding to period ranges"+'\n'+\
    	'		 \n'+\
    	"      	 version = Emilie Parent, McGill University (Spring, 2016) "+'\n'+'\n'+\
	"    		--h    	 :	Displays this help message"+'\n'+\
	"				"+'\n'+\
	"    		--mindc  :	Minimum duty-cycle to look for. Default is 0.5%. "+'\n'+\
	"            			Options are: 0.5, 1., 1.5. It multiplies the list "+'\n'+\
	"            			of minimum sampling intervals that has to be tested in each "+'\n'+\
	"            			subranges of periods by 2 (mindc=1) or 3 (mindc =1.5)"+'\n'+\
	"            			default is set in Config.py "+'\n'+\
	"				"+'\n'+\
	"   		--SN_tresh: 	Signal-to-noise treshold for picking candidates"+'\n'+\
	"			        default is set in Config.py"+'\n'+\
	"            			"+'\n'+\
	"   		--win_detrend: 	Time window for detrending (in sec)"+'\n'+\
	"			        default is set in Config.py"+'\n'+\
	"            			"+'\n'+\
	"   		--plot    :	Produces periodograms for a few different duty cycles,"+'\n'+\
	"				default is False "+'\n'+\
	"				"+'\n'+\
	"    		--save    :	Dump SNR values into pickle"+'\n'+\
	"			        default is False"+'\n'+\
	"				"+'\n'+\
	"    		--noxwin  :	Don't show the periodogram"+'\n'\
	"			        default is False"

    print help
    sys.exit()

