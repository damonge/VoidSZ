import os
import sys
import numpy as np

fname_y         =sys.argv[1]
fname_msk       =sys.argv[2]
fname_data_voids=sys.argv[3]
prefix_mocks    =sys.argv[4]
suffix_mocks    =sys.argv[5]
prefix_out      =sys.argv[6]

for n in [15,20,30] :
    print n
    dirname=prefix_out+"_%d"%n
    os.system("mkdir -p "+dirname)
    os.system("./FieldXCorr "+fname_y+" "+fname_msk+" "+fname_data_voids+" "+
              dirname+"/wth_voids.txt 3 %d 0 0. NO_WEIGHT NO_CUT theta_eff 1 > "%n+dirname+"/log_voids")

    for i in np.arange(1000)+1 :
        a="%04d"%i
        fname_cat=prefix_mocks+"_"+a+"_"+suffix_mocks
        fname_wth=dirname+"/wth_mock_"+a+".txt"
        fname_log=dirname+"/wth_mock_"+a+".log"
        print fname_cat
        os.system("./FieldXCorr "+fname_y+" "+fname_msk+" "+fname_cat+" "+fname_wth+
                  " 3 %d 0 0. NO_WEIGHT NO_CUT theta_eff 1 > "%n+fname_log)
